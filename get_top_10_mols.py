import pandas as pd
import argparse
import os

from rdkit import Chem
from rdkit.Chem import AllChem

def get_mol_names(results_file):
	df = pd.read_csv(results_file, sep="	", index_col=False)

	df_sorted = df.sort_values(by='docking_score')
	df_sorted = df_sorted['docking_score'] > -20
	df_top = df_sorted.head(10)

	mol_names = df_top['Name'].tolist()

	return mol_names

def get_top_mols(sdf_file, mol_names):
	suppl = Chem.SDMolSupplier(sdf_file)
	mols = [mol for mol in suppl]

	top_mols = []
	[top_mols.append(x) for x in mols if x.GetProp('_Name') in mol_names]

	return top_mols

def main():
	parser = argparse.ArgumentParser(description='Extract sdf file with top 1O molecules')
	parser.add_argument('-i', '--input', metavar='DIRECTORY', required=True, help='path to input directory')
	args = parser.parse_args()

	mol_names = get_mol_names(os.path.join(args.input, 'results.csv'))
	top_mols = get_top_mols(os.path.join(args.input, 'docking_results.sdf'), mol_names)

	with Chem.SDWriter(os.path.join(args.input, "top_10.sdf")) as w:
		for m in top_mols:
			w.write(m)

if __name__ == '__main__':
	main()