# Infrastructure Update Design

**Date:** 2026-05-20  
**Branch:** `parker/update-infra`  
**Scope:** Python 3.12 upgrade, uv migration (hatchling), bug fixes, test improvements

---

## 1. Goals

1. Upgrade Python requirement from 3.9 to 3.12
2. Replace `setup.py` + `requirements.txt` with `pyproject.toml` (hatchling backend) + `uv.lock`
3. Update all three CI workflows (Linux, macOS, Windows) to use `uv`
4. Fix 11 confirmed bugs across `guide.py`, `pool.py`, `batch_retest.py`, and `check_seq.py`
5. Improve test suite to catch the above bugs going forward

---

## 2. Branch & File Changes

| Action | File |
|--------|------|
| Remove | `setup.py` |
| Remove | `requirements.txt` (contained circular self-reference `pyguide~=0.1.0`) |
| Add | `pyproject.toml` |
| Add | `uv.lock` (generated via `uv lock`) |
| Update | `.github/workflows/python_package-linux.yml` |
| Update | `.github/workflows/python_package-macos.yml` |
| Update | `.github/workflows/python_package-windows.yml` |
| Add | `docs/superpowers/specs/2026-05-20-infra-update-design.md` (this file) |

---

## 3. `pyproject.toml`

Build backend: **hatchling** (uv default; no configuration needed for this package layout).

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "pyguide"
version = "0.1.1"
requires-python = ">=3.12"
dependencies = [
    "pandas>=2.0.0",
    "numpy>=1.26.0",
    "biothings-client>=0.2.6",
    "mygene>=3.2.2",
]

[project.optional-dependencies]
dev = ["pytest>=7.0"]

[project.scripts]
pyguide-order         = "pyguide.guide:main"
pyguide-collate       = "pyguide.pool:main"
pyguide-batch-retest  = "pyguide.batch_retest:main"
pyguide-check-seq     = "pyguide.check_seq:main"

[project.urls]
Homepage = "https://github.com/pgrosjean/pyguide"

[tool.hatch.build.targets.wheel]
packages = ["pyguide"]
```

Dependency version rationale:
- `pandas>=2.0.0` — pandas 1.x does not support Python 3.12
- `numpy>=1.26.0` — numpy 1.23 does not support Python 3.12; 1.26 is the first 3.12-compatible release
- `biothings-client>=0.2.6` — minimum version retained from original; no 3.12 incompatibility known
- `mygene>=3.2.2` — unchanged from original
- `streamlit` — kept out of declared deps (README instructs separate install)
- `setuptools` — removed from runtime deps (it is a build tool, not a runtime requirement)

---

## 4. CI Workflows

All three workflows (`python_package-linux.yml`, `python_package-macos.yml`, `python_package-windows.yml`) get the same treatment:

- Bump `actions/checkout` from `v3` → `v4`
- Replace `actions/setup-python@v2` + `pip install -e ".[dev]"` with `astral-sh/setup-uv@v5` + `uv sync --extra dev`
- Replace `pytest -v` with `uv run pytest -v`
- Set Python version to `3.12` in the uv setup step

Resulting steps block (identical across all three):

```yaml
steps:
  - uses: actions/checkout@v4

  - name: Set up uv and Python 3.12
    uses: astral-sh/setup-uv@v5
    with:
      python-version: "3.12"

  - name: Install dependencies
    run: uv sync --extra dev

  - name: Test with pytest
    run: uv run pytest -v
```

Each workflow retains its existing `runs-on` value (`ubuntu-latest`, `macos-latest`, `windows-latest`) and trigger configuration.

---

## 5. Bug Fixes

### 5.1 `guide.py` — `filter_database` (lines 151–161)

**Bug:** If both `guide_ids` and `gene_names` are empty, `filtered_df` is never assigned, causing `NameError` on return.

**Fix:** Add an early return for the empty case:
```python
if len(guide_ids) == 0 and len(gene_names) == 0:
    return df.iloc[0:0]
```

---

### 5.2 `guide.py` — `make_query_map` (line 1013)

**Bug:** `alias = list(alias)` on a string produces `['A','P','O','E']` (list of characters) instead of `['APOE']`.

**Fix:**
```python
if pd.notna(alias):
    if not isinstance(alias, list):
        alias = [alias]
```

---

### 5.3 `guide.py` — `filter_cloned_guides` (line 463)

**Bug:** After `sort_values`, `gene_filter.head(num).index.values` returns pandas label indices (e.g. `[3, 7, 12]`), which are then used as positional indices into `bool_ind = np.ones(N)`. If any label ≥ N, this raises `IndexError`.

**Fix:** Reset the index immediately after the two `sort_values` calls at the top of the function:
```python
filtered_df = filtered_df.sort_values('score', ascending=False).sort_values('gene').reset_index(drop=True)
```

---

### 5.4 `guide.py` — `order_guides` reset_index not captured (lines 1177, 1183, 1189)

**Bug:** `collated_df.reset_index(drop=True)` returns a new DataFrame; the result is discarded three times, leaving the index unreset.

**Fix:** Assign result on each line:
```python
collated_df = collated_df.reset_index(drop=True)
```

---

### 5.5 `guide.py` — `order_guides` `gene_list` possibly undefined (line 1195)

**Bug:** `gene_list` is only assigned inside `if len(gene_names) > 0` but is referenced unconditionally at line 1195 (`get_local_db(gene_list)`). If `len(gene_names) == 0`, this raises `NameError`.

**Fix:** Initialize before the conditional block:
```python
gene_list = []
```

---

### 5.6 `guide.py` — `split_dataframe` off-by-one (line 550)

**Bug:** `num_chunks = len(df) // chunk_size + 1` always produces one extra empty chunk when `len(df)` is exactly divisible by `chunk_size` (e.g. 96 guides → 2 chunks, second empty).

**Fix:** Use ceiling division:
```python
import math
num_chunks = math.ceil(len(df) / chunk_size)
```

---

### 5.7 `guide.py` — `main` dead code (lines 1366–1369)

**Bug:** The inner `if/else` has the same condition as the outer `if`, making the inner `else` branch unreachable.

**Fix:** Collapse to a single assert:
```python
if args.order_format in ("pooled", "batch-retest"):
    assert primer_df is not None, "You must run pyguide-collate (or pyguide-batch-retest) to generate primers before pooled ordering."
```

---

### 5.8 `pool.py` — `generate_pooled_list` hardcoded path separator (line 46)

**Bug:** `f"{base_dir}/{file_name}"` hardcodes `/` after `base_dir` is computed with `\\` on Windows, producing malformed paths.

**Fix:**
```python
file = os.path.join(base_dir, file_name)
```

---

### 5.9 `batch_retest.py` — `lib_num` never incremented (lines 52–54)

**Bug:** `lib_num` is set to `0` before the loop and never incremented, so every guide is written with `lib_num=0`.

**Fix:** Use `enumerate`:
```python
for idx, (guide_id, (left_primer, right_primer)) in enumerate(zip(guide_ids, primer_list)):
    outfile.write(f"{guide_id}\t{left_primer}\t{right_primer}\t{idx}\n")
```

---

### 5.10 `check_seq.py` — `pd.read_csv` missing separator (line 190)

**Bug:** `pd.read_csv(file)` on a tab-delimited `.txt` file parses each full line as one column.

**Fix:**
```python
df = pd.read_csv(file, sep="\t")
```

---

### 5.11 `check_seq.py` — incomplete `update_db` branch (lines 186–194)

**Bug:** When `--update_db` is passed, `df` is loaded but nothing is done with it. The log file write is also skipped entirely in this branch.

**Fix:** Move the log file write outside the `if/else` so it always runs, and add a `# TODO` comment marking the incomplete database-update logic clearly:
```python
# Always write the log file
with open(log_file_name, "w") as f:
    for file, guide_id in sorted(zip(file_list, list(guide_id_arr)), key=lambda x: x[1]):
        f.write(f"{file} \t {guide_id}\n")

if args.update_db:
    # TODO: implement database update logic
    pass
```

---

## 6. Test Improvements

The following improvements bring the test suite from smoke-testing to behaviorally verified:

### 6.1 Fix fragile path construction (all test files)

Replace:
```python
file_path_1 = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
```
With:
```python
file_path_1 = os.path.dirname(os.path.abspath(__file__))
```

### 6.2 Add output file existence assertions (test_single.py, test_arrayed.py, test_pooled.py)

After each `order_guides()` call, verify the expected output file was created using `glob.glob` on the output directory. Clean up generated files after each test to avoid cross-test pollution.

### 6.3 Add `batch-retest` end-to-end test (test_pooled.py)

Add `test_batch_retest_order_i()` covering the full flow: `batch_retest.main()` → `guide.order_guides(..., order_format="batch-retest", ...)`.

### 6.4 Add mouse + CRISPRa coverage to check_seq tests (test_check_seq.py)

Add `test_check_seq_crispra()` and a test that explicitly verifies that files containing only `NNNNN` sequences return `None` guide IDs.

### 6.5 Add input validation tests (new test_validation.py)

Test that `order_guides` raises `AssertionError` for:
- `guides_per_gene > 10`
- Both `guide_ids` and `gene_names` empty
- Invalid `order_format`

---

## 7. Out of Scope

- `check_seq.py` `update_db` full implementation (marked TODO in bug fix 5.11)
- Adding `streamlit` to declared dependencies
- Any functional changes to the guide ordering logic beyond the bug fixes above
- Increasing the mouse gene list in test fixtures

---

## 8. Success Criteria

1. `uv sync --extra dev && uv run pytest -v` passes on Python 3.12 locally
2. All three CI workflows pass on Linux, macOS, Windows
3. All 11 bugs listed in Section 5 are fixed with corresponding test coverage
4. `pyguide-order`, `pyguide-collate`, `pyguide-batch-retest`, `pyguide-check-seq` entry points all work via `uv run`
