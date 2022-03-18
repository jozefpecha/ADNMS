import os
import argparse
import pandas as pd
import json
import sqlite3

def get_data(out_folder, res_dir):
	out_file=os.path.join(out_folder, 'combined_results.csv')
	feature_folder_names = list(os.listdir(res_dir))
	collected_results = []

	for feature_folder in feature_folder_names:
		feature_folder_path = os.path.join(res_dir, feature_folder)
		if os.path.isdir(feature_folder_path):
			results_file_path=os.path.join(res_dir, feature_folder, 'results.csv')
			if feature_folder == min(feature_folder_names):
				results = pd.read_csv(results_file_path, sep="	", index_col=False)
				collected_results.append(results)
			else:
				results = pd.read_csv(results_file_path, sep="	", index_col=False, skiprows=0)
				collected_results.append(results)
		else:
			continue

	combined_results = pd.concat(collected_results, ignore_index=True)
	combined_results.fillna('None')

	combined_results.to_csv(out_file, sep='\t', index=False)

	return combined_results

def extract_enumerated_mols(db_file):
	conn = sqlite3.connect(db_file)
	cur = conn.cursor()
	
	query = cur.execute("SELECT sum(nmols) FROM mols WHERE rowid in (SELECT min(rowid) FROM mols WHERE parent_mol_id is null group by id)")
	record = query.fetchall()
	
	cur.close()
	
	return record

def analyze_data(combined_results, out_folder, db_dir=None):
	json_list = []

	if db_dir != None:
		if os.path.isdir(db_dir):
			sim_time_file_path=os.path.join(db_dir, 'execution_time.txt')
			if os.path.exists(sim_time_file_path):
				with open(sim_time_file_path, 'r') as f:
					time = 0
					for line in f.readlines():
						time += int(line.strip().split('\t')[1])
			else:
				time=None
				print(f'File {sim_time_file_path} not found! Skipping run_time evaluation.')

			db_file_path=os.path.join(db_dir, 'res.db')
			if os.path.exists(db_file_path):
				enumerated_mols=extract_enumerated_mols(db_file_path)[0][0]
			else:
				enumerated_mols=None
				print(f'{sim_time_dir} is not a directory! Skipping enumerated mols evaluation.')
		else:
			raise Exception(f'{db_dir} is not a directory!')
	else:
		time=None
		enumerated_mols=None

	max_features=combined_results.matched_ids_count.max()
	n_features=max_features-combined_results.matched_ids_count.min()

	generated_molecules=0
	
	list_of_dicts = []
	for n in range(n_features+1):
		dict_name = f'feature_stats_N-{n}'
		tmp_dict = {}
		tmp_dict[f'N-{n}_mols']=int((combined_results.matched_ids_count==(max_features-n)).sum())
		tmp_dict[f'N-{n}_MW<450']=int(((combined_results.matched_ids_count==(max_features-n) & (combined_results.MW<450)).sum()))
		tmp_dict[f'N-{n}_MW<450_&_RMSD<=2']=int(((combined_results.matched_ids_count==(max_features-n) & (combined_results.MW<450) & (combined_results.rmsd_score<=2)).sum()))
		tmp_dict[f'N-{n}_MW<450_&_logP<=4']=int(((combined_results.matched_ids_count==(max_features-n) & (combined_results.MW<450) & (combined_results.logP<=4)).sum()))
		tmp_dict[f'N-{n}_MW<450_&_RMSD<=2_&_logP<=4']=int(((combined_results.matched_ids_count==(max_features-n) & (combined_results.MW<450) & (combined_results.logP<=4) & (combined_results.rmsd_score<=2)).sum()))
		
		generated_molecules += tmp_dict[f'N-{n}_mols']
		
		dict_list = {dict_name: tmp_dict}
		list_of_dicts.append(dict_list)
	
	json_list.append({'max_features': max_features, 'simulation_time': time, 'enumerated_mols': enumerated_mols, 'generated_molecules': generated_molecules})
	json_list.append(list_of_dicts)

	out_file_summary=os.path.join(out_folder, 'summary.json')
	with open(out_file_summary, 'w') as f:
		json.dump(json_list, f, indent=2)

def main():
	parser = argparse.ArgumentParser(description='Extract sdf file with top 1O molecules')
	parser.add_argument('-i', '--res_dir', metavar='DIRECTORY', required=True, help='path to input directory')
	parser.add_argument('-o', '--out_folder', metavar='DIRECTORY', required=True, help='path to output directory')
	parser.add_argument('-db', '--db_dir', metavar='DIRECTORY', required=False, help='path to directory with a txt file with simulation run time')
	args = parser.parse_args()

	combined_results = get_data(args.out_folder, args.res_dir)
	analyze_data(combined_results, args.out_folder, args.db_dir)

if __name__ == '__main__':
	main()
