from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import argparse
import os

from similarity_search import calc_similarity


def load_mols(in_file):
	"""

	param in_file: sdf file with mols extracted from db
	return: lists of RDKit Mols, smiles, matched ids and matched ids count
	"""
	suppl = Chem.SDMolSupplier(str(in_file))
	mols = [mol for mol in suppl]
	smiles, matched_ids, matched_ids_count = [], [], []
	for mol in mols:
		smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True))
		matched_ids.append(mol.GetProp("matched_ids"))
		matched_ids_count.append(mol.GetProp("matched_ids_count"))

	return mols, smiles, matched_ids, matched_ids_count


def get_fp(mols):
	"""

	param mols: list of RDKit Mols
	return: list of Morgan FPs
	"""
	mfp_2 = [AllChem.GetMorganFingerprintAsBitVect(x, 2) for x in mols]

	return mfp_2



def names(mols):
	"""

	param mols: list of RDKit Mols
	return: sequential names of generated mols
	"""
	id_names = []
	mol_names = [mol.GetProp("_Name") for mol in mols]
	for i in range(1, len(mols)+1):
		id_names.append(f"F{i:06d}_{mol_names[i-1]}")

	return id_names


def get_pdbs(mols, mol_names):
	"""

	param mols: list of RDKit Mols
	param mol_names: list of names of molecules
	return: pdb files of generated molecules
	"""
	for mol, name in list(zip(mols, mol_names)):
		mol.SetProp('_Name', name)
		file_name = f"docking_results/pdb_files/{name}.pdb"
		file = os.path.join(args.out, file_name)
		Chem.MolToPDBFile(mol, file)
	
	return None


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Extract fragment smiles')
	parser.add_argument('-i', '--input', metavar='STRING', required=True, help='path to input sdf file')
	parser.add_argument('-r', '--ref', metavar='STRING', required=True, help='path to reference ligand pdb file')
	parser.add_argument('-o', '--out', metavar='STRING', required=True, help='path to output directory')
	args = parser.parse_args()

	pp_output = args.out
	os.makedirs(pp_output, exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "pdb_files"), exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "pdbqt_files"), exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "pdbqt_models"), exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "docking_logs"), exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "rsmd", "smiles_files"), exist_ok=True)
	os.makedirs(os.path.join(pp_output, "docking_results", "rsmd", "rsmd_results"), exist_ok=True)

	mols, smiles, matched_ids, matched_ids_count = load_mols(args.input)
	ref = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromPDBFile(args.ref), 2)
	ligands = get_fp(mols)
	sim_scores = calc_similarity(ref, ligands)
	mol_names = names(mols)
	zipped = list(zip(smiles, mol_names))

	# write out a file with all smiles
	with open(os.path.join(args.out, 'smiles.smi'), "w") as f1:
		for n, i in enumerate(zipped):
			f1.write(zipped[n][0] + "\t" + zipped[n][1] + "\n")

	# write out smiles files for each molecule (needed for rmsd calculation)
	for n, i in enumerate(zipped):
		with open(os.path.join(args.out, f"docking_results/rsmd/smiles_files/{i[1]}.smi"), "w") as f:
			f.write(i[0] + "\t" + i[1])

	# write out a temporary results file with combined data
	with open(os.path.join(args.out, 'tmp_res.csv'), "w") as f2:
		f2.write("Name" + "\t" + "Smiles" + "\t" + "matched_ids_count" + "\t" + "matched_ids" + "\t" + "sim_score" + "\n")
		for n in range(len(mols)):
			f2.write(mol_names[n] + '\t' + smiles[n] + '\t' + str(matched_ids_count[n]) + '\t' + matched_ids[n] + '\t' + str(round(sim_scores[n], 4)) + '\n')

	get_pdbs(mols, mol_names)

