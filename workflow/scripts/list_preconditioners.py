"""Introspect FiPy for names of preconditioners in active suite
"""

import os

if "FIPY_SOLVERS" not in os.environ:
    os.environ["FIPY_SOLVERS"] = snakemake.wildcards.suite

import fipy as fp

if __name__ == "__main__":
    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('list_preconditioners')
    fh = logging.FileHandler(str(snakemake.log[0]))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        preconditioners = []

        for k, v in fp.solvers.__dict__.items():
            if isinstance(v, type) and issubclass(v, fp.preconditioner.Preconditioner):
                preconditioners.append(k)

        with open(snakemake.output[0], 'w') as f:
            f.write(" ".join(sorted(preconditioners)))
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e
