from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
import os
import pandas as pd

from similarity_search import calc_similarity


def load_mols(in_file, out_file):
	"""

	param in_file: sdf file with mols extracted from db
	return: lists of RDKit Mols, smiles, matched ids and matched ids count
	"""
	suppl = Chem.SDMolSupplier(str(in_file))
	mols = [mol for mol in suppl]
	names, smiles, matched_ids, matched_ids_count = [], [], [], []
	for mol in mols:
		names.append(new_name)
		smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True))		
		matched_ids.append(mol.GetProp("matched_ids"))
		matched_ids_count.append(mol.GetProp("matched_ids_count"))
	

	return mols, names, smiles, matched_ids, matched_ids_count

def load_refs(ref_file):
	"""

	param ref_file: smi file containing smiles of reference molecules
	return: list of reference molecules
	"""
	ref_mols = []
	ref_mol_names = []
	with open(ref_file, "r") as f:
		for line in f:
			ref_mol_name = line.strip().split('\t')[1]
			ref_mol = Chem.MolFromSmiles(line.strip().split('\t')[0])
			ref_mol.SetProp('_Name', ref_mol_name)

			ref_mols.append(ref_mol)
			ref_mol_names.append(ref_mol_name)
			

	return ref_mols, ref_mol_names


def get_fp(mols):
	"""

	param mols: list of RDKit Mols
	return: list of Morgan FPs5
	"""
	mfp_2 = [AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]

	return mfp_2
		

def get_pdbs(mols):
	"""

	param mols: list of RDKit Mols
	param mol_names: list of names of molecules
	return: pdb files of generated molecules
	"""
	for mol in mols:
		file_name = f"docking_results/pdb_files/{mol.GetProp('_Name')}.pdb"
		file = os.path.join(args.out, file_name)
		Chem.MolToPDBFile(mol, file)
	
	return None


def main():
	parser = argparse.ArgumentParser(description='Extract fragment smiles')
	parser.add_argument('-i', '--input', metavar='STRING', required=True, help='path to input sdf file')
	parser.add_argument('-r', '--ref', metavar='FILE', required=True, help='path to file with reference molecules smiles and their tab separated names') # input as a file with reference smiles
	parser.add_argument('-o', '--out', metavar='STRING', required=True, help='path to output directory')
	parser.add_argument('-p', '--pdb', action='store_true', default=False, help='creates pdb files of ligands if specified')
	#parser.add_argument('-g', '--graph', action='store_true', default=False, help='creates a histogram from similarity values and saves it as svg')
	args = parser.parse_args()

	mols, names, smiles, matched_ids, matched_ids_count = load_mols(args.input, args.out)

	ref_mols, ref_names = load_refs(args.ref)

	ref_fps = get_fp(ref_mols)
	ligands = get_fp(mols)

	sim_scores = {}
	sim_scores['Name'] = names
	sim_scores['smiles'] = smiles
	sim_scores['matched_ids_count'] = matched_ids_count
	sim_scores['matched_ids'] = matched_ids

	for ref_fp in list(zip(ref_names, ref_fps)):
		if not ref_fp[1]:
			return None
		sim_scores[f'{ref_fp[0]}_similarity_scores'] = calc_similarity(ref_fp[1], ligands)

	df_results = pd.DataFrame.from_dict(sim_scores)
	df_results.to_csv(os.path.join(args.out, 'similarity_scores.csv'), index=False, sep='\t', na_rep='NaN')


	# write out a file with all smiles
	with open(os.path.join(args.out, 'smiles.smi'), "w") as f1:
		for i in list(zip(smiles, names)):
			f1.write(i[0] + "\t" + i[1] + "\n")
	
	if args.pdb:
		get_pdbs(mols)

if __name__ == '__main__':
	main()