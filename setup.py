from setuptools import setup, find_packages

setup(
    name='pyguide',
    version='0.1.0',
    packages=find_packages(),
    install_requires=["pandas>=1.4.1",
                      "numpy>=1.22.2",
                      "setuptools>=57.4.0"],
    extras_require={
        "dev": [
            "pytest>=3.7",
        ]
    },
    entry_points={
        'console_scripts': ['pyguide-order=pyguide.guide:main'],
    },
    url='https://github.com/pgrosjean/pyguide',
    license='MIT',
    author='Parker Grosjean',
    author_email='parker.grosjean@gmail.com',
    description='Tools for CRISPRi/a gRNA work'
)
