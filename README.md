# ADNMS
Analysis of de novo molecular structures


Automated tool for analysis of results of the CReM de novo generation tool: https://github.com/ci-lab-cz/crem-pharm.

## Dependencies
RDKit [https://www.rdkit.org/docs/Install.html]

Meeko: `pip install meeko`

Python bindings for vina: `pip install vina`

## Usage
The tool is implemented as a pbs file that can be run in terminal using:
```
qsub eval.pbs /path/to/input/database/file /path/to/output/directory/ /path/to/reference/molecule/pdb/file
```

This needs to be run from the directory containing all the helper python scripts in order to work.
If you're using anaconda, initialize your environment before running the tool.

## Explanation
The tool extracts all terminal molecules from the database (visited_ids_count = max(visited_ids_count)) into an sdf file. (Extraction condition can be changed in get_sdf.sql file.) Next, the sdf file is split according to number of pharmacophore features the molecules matched. For each category, physicochemical properties, Tanimoto similarity scores to the reference molecule and synthetic accessibility score are calculated and gathered in one csv file.
