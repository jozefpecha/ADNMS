import os
import pandas as pd
import argparse
import sqlite3

def extract_docking_scores(db_file):
	conn = sqlite3.connect(db_file)
	cur = conn.cursor()
    
	query = cur.execute("SELECT id, docking_score FROM mols")
	cols = ['Name', 'docking_score']
	record = pd.DataFrame.from_records(data = query.fetchall(), columns=cols)
    
	cur.close()
    
	return record


def load_data(logs_db, rmsd_file, physchem_file, sim_scores_file, sa_scores_file, output_file):
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
	sa_scores = pd.read_csv(sa_scores_file, sep="\t", index_col=False, names = ['Name', 'SA_score'])
	rmsd_scores = pd.read_csv(rmsd_file, sep="\t", index_col=False, names = ['Name', 'rmsd_score'], usecols=[0, 2])
	
	dock_scores = extract_docking_scores(logs_db)

	df1 = physchem.merge(sim_scores, how='outer', on='Name')
	df2 = df1.merge(sa_scores, how='outer', on='Name')
	df3 = df2.merge(rmsd_scores, how='outer', on='Name')
	df4 = df3.merge(dock_scores, how='outer', on='Name')

	df4 = df4[['Name', 'Smiles', 'sim_score', 'docking_score', 'rmsd_score', 'SA_score', 'matched_ids', 'MW', 'MR', 'HBA', 'HBD',
        'complexity', 'NumRings', 'RTB', 'TPSA', 'logP', 'Csp3', 'fmf', 'QED', 'HAC', 'NumRingsFused', 'N_unique_hba_hbd_atoms',
        'max_ring_size', 'ChiralCenters']]

	df4.to_csv(output_file, index=False, sep="\t")

	return None

def main():
	parser = argparse.ArgumentParser(description='Combine data into one file')
	parser.add_argument('-p', '--physchem', metavar='STRING', required=True, help='path to physchem props input file')
	parser.add_argument('-s', '--sim', metavar='STRING', required=True, help='path to similarity scores input file')
	parser.add_argument('-d', '--dock', metavar='FILE', required=True, help='path to docking scores database')
	parser.add_argument('-r', '--rmsd', metavar='FILE', required=True, help='path to rmsd scores input file')
	parser.add_argument('-sa', '--sa', metavar='STRING', required=True, help='path to sa scores input file')
	parser.add_argument('-o', '--out', metavar='STRING', required=True, help='path to output file')
	args = parser.parse_args()

	load_data(args.dock, args.rmsd, args.physchem, args.sim, args.sa, args.out)

if __name__ == '__main__':
	main()
