"""
Streamlit UI for pyguide Tools
- Creates a NEW temp folder each time "Run" is pressed
- Provides a cleanup button for that folder
"""

import os
import io
import shutil
import zipfile
import subprocess
import streamlit as st
from typing import List
from datetime import datetime

def main():
    st.title("pyguide Streamlit App (New Temp Folder Per Run)")

    # -------------------------
    # 1) Sidebar Configuration
    # -------------------------
    st.sidebar.header("Configuration")
    order_format = st.sidebar.selectbox(
        "Order Format",
        ["single", "arrayed", "pooled", "batch-retest"]
    )
    user_name = st.sidebar.text_input("User Name", value="JohnDoe")
    ai_status = st.sidebar.selectbox("CRISPR Type", ["i", "a"])  # CRISPRi or CRISPRa
    guides_per_gene = st.sidebar.number_input("Guides per gene", value=5, min_value=1, max_value=10)
    organism = st.sidebar.selectbox("Organism", ["human", "mouse"])
    check_db = st.sidebar.checkbox("Check local db for cloned guides?", value=False)
    ntc_frac = st.sidebar.slider("NTC fraction (pooled/batch-retest)", 0.0, 1.0, 0.2)
    
    st.write(
        f"**Selected Options**:\n\n"
        f"- Order Format: `{order_format}`\n"
        f"- CRISPR: `{ai_status}`\n"
        f"- Guides per gene: `{guides_per_gene}`\n"
        f"- Organism: `{organism}`\n"
        f"- Check DB: `{check_db}`\n"
        f"- NTC fraction: `{ntc_frac}`\n"
    )

    # ----------------------
    # 2) File Uploads
    # ----------------------
    st.subheader("File Uploads")
    wishlist_files = []
    primer_file = None

    if order_format in ["single", "arrayed"]:
        wf = st.file_uploader("Upload .txt wishlist file", type=["txt"])
        if wf:
            wishlist_files.append(wf)
    elif order_format == "pooled":
        wfs = st.file_uploader(
            "Upload one or more .txt wishlist files (for collation)",
            type=["txt"],
            accept_multiple_files=True
        )
        if wfs:
            wishlist_files = list(wfs)

        primer_file = st.file_uploader(
            "Optional primer file (txt/csv/tsv). If omitted, a default is used.",
            type=["txt", "csv", "tsv"]
        )
    else:  # batch-retest
        wf = st.file_uploader("Upload a .txt wishlist of guide IDs", type=["txt"])
        if wf:
            wishlist_files.append(wf)

        primer_file = st.file_uploader(
            "Optional primer file (txt/csv/tsv). If omitted, a default is used.",
            type=["txt", "csv", "tsv"]
        )

    # ----------------------
    # 3) Run Pipeline
    # ----------------------
    st.subheader("Run Pipeline")
    run_button = st.button("Run")

    # We'll store the path to the newly created folder in session_state
    if "latest_run_folder" not in st.session_state:
        st.session_state["latest_run_folder"] = None

    if run_button:
        # 3a) Create a unique temp folder for THIS run
        tmpdir = create_temp_output_directory()
        st.session_state["latest_run_folder"] = tmpdir

        # 3b) Pipeline logic
        if order_format in ["single", "arrayed"]:
            if not wishlist_files:
                st.error("Please upload a wishlist file.")
                return

            wishlist_path = save_uploaded_file(wishlist_files[0], tmpdir)

            cmd = [
                "pyguide-order",
                "--wishlist_file", wishlist_path,
                "--name", user_name,
                "--ai", ai_status,
                "--guides_per_gene", str(guides_per_gene),
                "--order_format", order_format,
                "--organism", organism
            ]
            if check_db:
                cmd.append("--check_db")

            st.info(f"Running: `{' '.join(cmd)}`")
            run_cli_command(cmd)
            st.success("Done ordering (single/arrayed).")

        elif order_format == "pooled":
            if not wishlist_files:
                st.error("Please upload one or more wishlist files.")
                return

            wishlist_paths = [save_uploaded_file(f, tmpdir) for f in wishlist_files]

            primer_path = ""
            if primer_file:
                primer_path = save_uploaded_file(primer_file, tmpdir)

            cmd_collate = ["pyguide-collate", "--wishlist_files"] + wishlist_paths
            if primer_path:
                cmd_collate += ["--primer_file", primer_path]

            st.info(f"Running Collate: `{' '.join(cmd_collate)}`")
            run_cli_command(cmd_collate)
            st.success("Collation completed.")

            collated_file_path = find_collated_file(tmpdir)
            if not collated_file_path:
                st.error("Could not find the collated wishlist file.")
                return

            cmd_order = [
                "pyguide-order",
                "--wishlist_file", collated_file_path,
                "--name", user_name,
                "--ai", ai_status,
                "--guides_per_gene", str(guides_per_gene),
                "--order_format", order_format,
                "--organism", organism,
                "--ntc_frac", str(ntc_frac)
            ]
            if check_db:
                cmd_order.append("--check_db")

            st.info(f"Running Order: `{' '.join(cmd_order)}`")
            run_cli_command(cmd_order)
            st.success("Done ordering (pooled).")

        else:  # batch-retest
            if not wishlist_files:
                st.error("Please upload a wishlist file containing guide IDs.")
                return

            wishlist_path = save_uploaded_file(wishlist_files[0], tmpdir)

            primer_path = ""
            if primer_file:
                primer_path = save_uploaded_file(primer_file, tmpdir)

            cmd_batch = ["pyguide-batch-retest", "--wishlist_file", wishlist_path]
            if primer_path:
                cmd_batch += ["--primer_file", primer_path]

            st.info(f"Running batch-retest: `{' '.join(cmd_batch)}`")
            run_cli_command(cmd_batch)
            st.success("Batch retest wishlist generated.")

            batch_file_path = find_batch_retest_file(tmpdir)
            if not batch_file_path:
                st.error("Could not find the batch_retest wishlist file.")
                return

            cmd_order = [
                "pyguide-order",
                "--wishlist_file", batch_file_path,
                "--name", user_name,
                "--ai", ai_status,
                "--guides_per_gene", str(guides_per_gene),
                "--order_format", order_format,
                "--organism", organism,
                "--ntc_frac", str(ntc_frac)
            ]
            if check_db:
                cmd_order.append("--check_db")

            st.info(f"Running Order: `{' '.join(cmd_order)}`")
            run_cli_command(cmd_order)
            st.success("Done ordering (batch-retest).")

    # ----------------------
    # 4) Show/Download results for the LATEST run folder
    # ----------------------
    latest_dir = st.session_state.get("latest_run_folder", None)
    if latest_dir and os.path.exists(latest_dir):
        out_files = os.listdir(latest_dir)
        if out_files:
            st.subheader("Output Files from the Latest Run")

            # Individual file downloads
            st.write("Download individual files:")
            for f in out_files:
                path = os.path.join(latest_dir, f)
                if os.path.isfile(path):
                    with open(path, "rb") as gf:
                        st.download_button(
                            label=f"Download {f}",
                            data=gf.read(),
                            file_name=f
                        )

            # "Download All" ZIP button
            st.write("Or download everything as a single ZIP:")
            if st.button("Download All Files as ZIP"):
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
                    for f in os.listdir(latest_dir):
                        fpath = os.path.join(latest_dir, f)
                        if os.path.isfile(fpath):
                            with open(fpath, "rb") as gf:
                                zipf.writestr(f, gf.read())
                zip_buffer.seek(0)

                st.download_button(
                    label="Click to Download ZIP",
                    data=zip_buffer,
                    file_name="pyguide_outputs.zip"
                )

            # Cleanup button
            if st.button("Cleanup This Temp Folder"):
                cleanup_tempdir(latest_dir)
                st.session_state["latest_run_folder"] = None
                st.success("Latest run folder removed.")
        else:
            st.info("The latest run folder is empty.")
    else:
        st.info("No run folder yet or it has been cleaned up.")


# ----------------------------------------------------------------------
# Helper Functions
# ----------------------------------------------------------------------

def create_temp_output_directory() -> str:
    """
    Creates a unique output folder, e.g. 'temp_pyguide_outputs_YYYYMMDD-HHMMSS'
    and returns its path.
    """
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    folder_name = f"temp_pyguide_outputs_{ts}"
    tmpdir = os.path.join(os.getcwd(), folder_name)
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir


def run_cli_command(cmd: List[str]):
    """
    Runs a console command using subprocess.run and prints stdout/stderr in Streamlit.
    """
    process = subprocess.run(cmd, capture_output=True, text=True)
    if process.returncode != 0:
        st.error(f"Command failed with return code {process.returncode}")
        st.error(process.stderr)
    else:
        st.text("Command succeeded.")
        if process.stdout:
            st.text(process.stdout)
        if process.stderr:
            st.text(process.stderr)


def save_uploaded_file(uploaded_file, directory: str) -> str:
    """
    Saves 'uploaded_file' into 'directory' and returns the local path.
    """
    fname = uploaded_file.name
    local_path = os.path.join(directory, fname)
    with open(local_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return local_path


def find_collated_file(directory: str) -> str:
    """
    Looks for 'collated_pooled_wish_list_YY_MM_DD.txt' in 'directory'.
    """
    for f in os.listdir(directory):
        if "collated_pooled_wish_list" in f and f.endswith(".txt"):
            return os.path.join(directory, f)
    return ""


def find_batch_retest_file(directory: str) -> str:
    """
    Looks for 'batch_retest_wishlist_YY_MM_DD.txt' in 'directory'.
    """
    for f in os.listdir(directory):
        if "batch_retest_wishlist" in f and f.endswith(".txt"):
            return os.path.join(directory, f)
    return ""


def cleanup_tempdir(dir_path: str):
    """Removes the specified temp directory entirely."""
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path, ignore_errors=True)


if __name__ == "__main__":
    main()
