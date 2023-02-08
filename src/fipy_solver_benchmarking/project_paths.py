# Inspired by https://medium.com/analytics-vidhya/data-directory-in-jupyter-notebooks-dc46cd79eb2f

from pathlib import Path

module_path = Path(__file__)
project_path = (module_path / "../../").resolve()
data_path = (project_path / "data").resolve()
results_path = (project_path / "results").resolve()
