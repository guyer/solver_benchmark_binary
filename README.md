Timing and solver convergence info developed while working on #701 (Add 
support for PETSc solvers)

## INSTALLATION

- Create conda environment:
    - To benchmark PySparse (and other Python 2.7 solvers)

      `mamba env create --name <ENVNAME> --file environment_27.yml`
      
    - To benchmark Python 3 solvers

      `mamba env create --name <ENVNAME> --file environment_3x.yml`
      
- Activate the environment:

  `mamba activate <ENVNAME>`

- Install a working copy of FiPy:

  `python -m pip install --editable <PATH_TO_FIPY>`
  
- Install the supporting package for this project:

    - for Python 2.7:

      `python setup.py develop`

    - for Python 3.x:

      `python -m pip install --editable .`
      

## USAGE

## DATA & FILE OVERVIEW

```
./
├── README.md
│    { Description of `solvers_and_timings/`.
├── TODO.md
│    { Tasks to work on
├── codes/
│    { Runable programs
│   ├── config/
│   │    { Python logging configuration files.
│   ├── notebooks/
│   │    { [Jupyter](https://jupyter.org/) notebook files.
│   └── scripts/
│        { Python and bash programs.
└── results/
    ├── diffusion/
    │    { Simulation results from `{fipy}/examples/benchmarking/solvers/diffusion.py`.
    ├── nucleation/
    │    { Simulation results using Phase Field Benchmark 8(b)
    └── old/
         { old stuff
```

(generated with `tree -FL 2 --info --filesfirst .`)

### DATA-SPECIFIC INFORMATION FOR: *.csv

[Comma-separated value](https://en.wikipedia.org/wiki/Comma-separated_values)-formatted
delimited text.

### DATA-SPECIFIC INFORMATION FOR: .info

[`tree` command file comment file](https://en.wikipedia.org/wiki/Tree_(command))-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.ipynb

[Jupyter notebook](https://nbformat.readthedocs.io/en/latest/format_description.html#notebook-file-format)-formatted
text.

### DATA-SPECIFIC INFORMATION FOR: *.json

[JavaScript Object Notation](https://www.json.org/)-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.log

[Python event log](https://docs.python.org/3/library/logging.html).

### DATA-SPECIFIC INFORMATION FOR: *.md

[Markdown](https://daringfireball.net/projects/markdown/)-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.mtx

[Matrix Market](https://math.nist.gov/MatrixMarket/formats.html)-formatted
text.

### DATA-SPECIFIC INFORMATION FOR: *.npy

[NumPy array](https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format)-formatted
binary.

### DATA-SPECIFIC INFORMATION FOR: *.npz

A [NumPy zipped archive of files](https://numpy.org/doc/stable/reference/generated/numpy.savez.html)
named after the variables they contain.  The archive is not compressed and
each file in the archive contains one variable in
[`.npy`](https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html#module-numpy.lib.format)
format.

### DATA-SPECIFIC INFORMATION FOR: *.py

[Python code](https://python.org)-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.sh

[Bash script](https://en.wikipedia.org/wiki/Bash_(Unix_shell))-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.txt

Should be `.tsv` or `.csv`?
