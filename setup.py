from setuptools import setup, find_packages

setup(
    name='pyguide',
    version='0.1.1',
    packages=find_packages(),
    python_requires=">=3.9.0",
    install_requires=["pandas>=1.4.1",
                      "numpy>=1.22.2",
                      "biothings-client==0.2.6",
                      "setuptools>=57.4.0",
                      "mygene>=3.2.2"],
    extras_require={
        "dev": [
            "pytest>=3.7",
        ]
    },
    entry_points={
        'console_scripts': ['pyguide-order=pyguide.guide:main',
                            'pyguide-collate=pyguide.pool:main',
                            'pyguide-batch-retest=pyguide.batch_retest:main',
                            'pyguide-check-seq=pyguide.check_seq:main'],
    },
    url='https://github.com/pgrosjean/pyguide',
    license='MIT',
    author='Parker Grosjean',
    author_email='parker.grosjean@gmail.com',
    description='Tools for CRISPRi/a gRNA work'
)
