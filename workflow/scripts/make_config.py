import logging

# https://stackoverflow.com/a/55849527/2019542
logger = logging.getLogger('make_config')
fh = logging.FileHandler(str(snakemake.log[0]))
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

try:
    permutations = get_all_permutations(snakemake.wildcards)
    permutations.loc[snakemake.wildcards.id].to_json(snakemake.output[0])
except Exception as e:
    logger.error(e, exc_info=True)
    pass
