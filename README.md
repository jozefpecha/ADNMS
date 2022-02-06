# ADNMS
Analysis of de novo molecular structures


Automated tool for analysis of results of the CReM de novo generation tool: https://github.com/ci-lab-cz/crem-pharm.

## Dependencies
RDKit [https://www.rdkit.org/docs/Install.html]

Meeko: `pip install meeko`

Python bindings for vina: `pip install vina`

Helper rdkit scripts accesible from https://github.com/DrrDom/rdkit-scripts:
physchem_calc.py
vina_dock.py
rmsd_rdkit.py
read_input.py

## Usage
The tool is implemented as a pbs file that can be run in terminal using:
```
qsub eval.pbs
```


This needs to be run from the directory containing all the helper python scripts in order to work.
If you're using anaconda, initialize your environment before running the tool.


The tool takes several arguments:
`input_db` - path to database created by CReM tool
`output_dir` - path to ouput directory
`ref_mols` - path to a smi file with tab separated smiles and names of reference molecules for similarity calculations
`prot_pdbqt` - pdbqt file of protein for docking
`compare_dir` - optional argument for comparison of generated molecules with already existing ligands for the given protein

## Explanation
The tool extracts all terminal molecules from the database (visited_ids_count = max(visited_ids_count)) into an sdf file. (Extraction condition can be changed in get_sdf.sql file.) Next, the sdf file is split according to number of pharmacophore features the molecules matched. For each category, physicochemical properties, Tanimoto similarity scores to the reference molecule, synthetic accessibility score, docking score and rmsd between generated molecules before and after docking are calculated and gathered in one csv file.
