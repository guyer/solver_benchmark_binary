import os
import pandas as pd

def concat_csv(input, output, log):
    try:
        li = [pd.read_csv(fname, index_col=False) for fname in input]
        if li:
            df = pd.concat(li, ignore_index=True)
        else:
            df = pd.DataFrame()
        df.to_csv(output, index=False)
    except Exception as e:
        with open(log, 'w') as f:
            f.write(repr(e))
        raise e

def concat_json(input, output, log):
    try:
        li = [pd.read_json(fname) for fname in input]
        if li:
            df = pd.concat(li, ignore_index=True)
        else:
            df = pd.DataFrame()
        df.to_json(output)
    except Exception as e:
        with open(log, 'w') as f:
            f.write(repr(e))
        raise e

def get_permutation_ids(wildcards):
#     path = checkpoints.all_permutations.get(**wildcards).output[0]
    path = "config/all_permutations.csv"
    df = pd.read_csv(path)
    return df.index.map("{:07d}".format)

def extract_config_by_id(wildcards, input, output, log):
    import logging

    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('make_config')
    fh = logging.FileHandler(str(log))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        permutations = pd.read_csv(input)
        permutations.loc[int(wildcards.id)].to_json(output)
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e

def read_config(path):
    import json

    with open(path, 'r') as f:
        return json.load(f)
