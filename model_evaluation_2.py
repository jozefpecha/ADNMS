import os
import pandas as pd
import argparse

def get_dock_score(logs_dir):
	"""

	param logs_dir: directory with docking log files
	return: list of best docking scores of molecules
	"""
	dock_scores = []
	for log_file in os.listdir(logs_dir):
		with open(os.path.join(logs_dir, log_file), "r") as f:
			for n, i in enumerate(f.readlines()):
				# line 26 contains best docking score in all log files
				if n == 26:
					dock_scores.append(float(i.split()[1]))

	return dock_scores

def get_rmsd_score(rmsd_dir):
	"""

	param rmsd_dir: directory with rmsd scores files
	return: list of rmsd scores corresponding to the best bidning poses
	"""
	rmsd_scores = []
	for rmsd_log_file in os.listdir(rmsd_dir):
		with open(os.path.join(rmsd_dir, rmsd_log_file), "r") as f:
			rmsd_scores.append(float(f.readline().strip().split()[-1]))

	return rmsd_scores

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
	sa_scores = pd.read_csv(sa_scores_file, sep="\t", index_col=False, header=None, usecols=[1])

	dock_scores = get_dock_score(logs_dir)
	rmsd_scores = get_rmsd_score(rmsd_dir)

	df = physchem.merge(sim_scores, how="outer")
	df["docking_scores"] = dock_scores
	df["rmsd_scores"] = rmsd_scores
	df["sa_scores"] = sa_scores

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
