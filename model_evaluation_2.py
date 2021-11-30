import os
import pandas as pd
import argparse

def get_dock_and_rmsd(names, logs_dir, rmsd_dir):
	"""
	
	param logs_dir: directory with docking log files
	param rmsd_dir: directory with rmsd scores files
	return: dataframes with docking and rmsd scores
	"""
	dock_scores = []
	rmsd_scores = []
	for name in names:
		with open(os.path.join(logs_dir, f'{name}.txt'), "r") as f1:
			line = f1.readlines()[26]
			dock_scores.append(float(line.split()[1]))
		
		with open(os.path.join(rmsd_dir, f'{name}.txt'), "r") as f2:
			rmsd_scores.append(float(f2.readline().strip().split()[-1]))
	
	df_dock = pd.DataFrame(dock_scores, columns=['Docking_score'])
	df_rmsd = pd.DataFrame(rmsd_scores, columns=['RMSD_score'])

	return df_dock, df_rmsd


def load_data(logs_dir, rmsd_dir, physchem_file, sim_scores_file, sa_scores_file, output_file):
	"""

	param logs_dir: directory with docking log files
	param rmsd_dir: directory with rmsd scores files
	param physchem_file: path to file containing physchem values of mols
	param sim_scores_file: path to file containing similarity scores
	param sa_scores_file: path to file containing synthtic accessibility scores
	param outpzt_file: path to file to output combined results
	return: csv file containing combined results
	"""
	physchem = pd.read_csv(physchem_file, sep="\t", index_col=False)
	sim_scores = pd.read_csv(sim_scores_file, sep="\t", index_col=False)
	sa_scores = pd.read_csv(sa_scores_file, sep="\t", index_col=False, names = ['Name', 'SA_score'], usecols=[1])

	dock_scores, rmsd_scores = get_dock_and_rmsd(sim_scores['Name'].tolist(), logs_dir, rmsd_dir)
	
	data_frames = [physchem, rmsd_scores, sa_scores]
	df = pd.concat(data_frames, join='outer', axis=1).fillna('NaN')

	df.to_csv(output_file, index=False, sep="\t")

	return None

def main():
	parser = argparse.ArgumentParser(description='Combine data into one file')
	parser.add_argument('-p', '--physchem', metavar='STRING', required=True, help='path to physchem props input file')
	parser.add_argument('-s', '--sim', metavar='STRING', required=True, help='path to similarity scores input file')
	parser.add_argument('-d', '--dock', metavar='DIRECTORY', required=True, help='path to docking scores input directory')
	parser.add_argument('-r', '--rmsd', metavar='DIRECTORY', required=True, help='path to rmsd scores input directory')
	parser.add_argument('-sa', '--sa', metavar='STRING', required=True, help='path to sa scores input file')
	parser.add_argument('-o', '--out', metavar='STRING', required=True, help='path to output file')
	args = parser.parse_args()

	load_data(args.dock, args.rmsd, args.physchem, args.sim, args.sa, args.out)

if __name__ == '__main__':
	main()
