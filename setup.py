# Note: Because this project benchmarks PySparse under Python 2.7, pyproject.toml doesn't work

from setuptools import find_packages, setup

setup(
    name='fipy_solver_benchmarking',
    version='0.0.1',
    author="Jonathan Guyer",
    author_email="guyer@nist.gov",
    license="NIST Public Domain",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2",
        "License :: Public Domain",
        "Operating System :: OS Independent",
    ],
)
