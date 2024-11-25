Timing and solver convergence info developed while working on
[#701 (Add support for PETSc solvers)](https://github.com/usnistgov/fipy/pull/701).

## INSTALLATION

```bash
mamba env create --file workflow/envs/snakemake.yml
```

## USAGE

### Configuration

Edit `config/config.yml` to change the revisions to be benchmarked, the 
benchmarks to run, the solver suites to test, and the simulation 
parameters.

### Run

Notebooks for simulating different problems can be found in
`workflow/notebooks/*.py.ipynb`.  These notebooks are intended to be
batch-processed by [Snakemake](https://snakemake.readthedocs.io/).

#### Bootstrap

The revisions to be benchmarked must be cloned and appropriate conda
environments created for them before simulation variants can be run.  From
the base directory of this repository:

```bash
snakemake --use-conda all_envs
```

Progress can be tracked by launching
[panoptes](https://github.com/panoptes-organization/panoptes) and invoking
`snakemake` with the `--wms-monitor` option.

#### Simulation and Analysis

Once the FiPy revision environments have been established by the 
`all_envs` rule, all permutations can be run an plotted with:

```bash
snakemake --use-conda all
```

### Reporting

```bash
snakemake --report report.html
```

### Retrieving results

Bring results from CTCMS to laptop with

```bash
rsync -avzh  mr-french.nist.gov:/data/guyer/solver_benchmark_binary/results .
```

## DATA & FILE OVERVIEW

```
./
├── LICENSE.md
│    { Terms of use.
├── README.md
│    { Description of this workflow.
├── config/
│    { Snakemake workflow configuration.
│   └── config.yml
├── images/
│    { Workflow graphs.
├── logs/
│    { Workflow progress reports.
├── resources/
│    { Data used by workflows.
│   └── t=300.0.npz
├── results/
│    { Workflow results.
└── workflow/
     { Snakemake workflow definition.
    ├── Snakefile
    │    { Snakemake workflow entrypoint.
    ├── envs/
    │    { Conda environment files.
    ├── notebooks/
    │    { [Jupyter](https://jupyter.org/) notebook files.
    ├── report/
    │    ⎧ [Snakemake report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html)
    │    ⎩ configuration files.
    ├── rules/
    │    { Snakemake modules.
    └── scripts/
         { Python and bash programs.
```

(generated with `tree -FL 2 --info --filesfirst --gitignore .`)

### DATA-SPECIFIC INFORMATION FOR: *.csv

[Comma-separated value](https://en.wikipedia.org/wiki/Comma-separated_values)-formatted
delimited text.

### DATA-SPECIFIC INFORMATION FOR: .info

[`tree` command file comment file](https://en.wikipedia.org/wiki/Tree_(command))-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.py.ipynb

[Jupyter notebook](https://nbformat.readthedocs.io/en/latest/format_description.html#notebook-file-format)-formatted
text containing [Python code](https://python.org) code.

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

### DATA-SPECIFIC INFORMATION FOR: *.rst

[ReStructuredText](https://docutils.sourceforge.io/docs/user/rst/quickstart.html)-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.sh

[Bash script](https://en.wikipedia.org/wiki/Bash_(Unix_shell))-formatted text.

### DATA-SPECIFIC INFORMATION FOR: *.txt

Should be `.tsv` or `.csv`?
