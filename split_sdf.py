from rdkit import Chem
import os
import argparse

def get_keys(in_file):
	suppl = Chem.SDMolSupplier(in_file)
	mols = [mol for mol in suppl]
	matched_ids_count = []
	[(matched_ids_count.append(mol.GetProp("matched_ids_count"))) for mol in mols  if mol.GetProp("matched_ids_count") not in matched_ids_count]

	return mols, {key:[] for key in matched_ids_count}

def split_sdf(mols, keys):
	for n, mol in enumerate(mols):
		new_name = f"F{(n+1):06d}_{mol.GetProp('_Name')}"
		mol.SetProp("_Name", new_name)
		keys[mol.GetProp("matched_ids_count")].append(mol)
    
	return keys

def write_sdfs(mol_dicts, output_dir):
	for key, val in mol_dicts.items():
		os.makedirs(os.path.join(args.out, key), exist_ok=True)
		path = f'{os.path.join(args.out, key)}/partial.sdf'
		with Chem.SDWriter(path) as w:
			for n, m in enumerate(val):
				w.write(m)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Divide sdf based on matching_ids_count')
	parser.add_argument('-i', '--input', metavar='STRING', required=True, help='path to input sdf file') 
	parser.add_argument('-o', '--out', metavar='STRING', required=True, help='path to output directory')
	args = parser.parse_args()

	os.makedirs(os.path.join(args.out), exist_ok=True)
	mols, dicts = get_keys(args.input)
	splitted_sdf = split_sdf(mols, dicts)
	write_sdfs(splitted_sdf, args.out)

	print("Finished spliting sdf files!")