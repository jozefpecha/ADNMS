import pandas as pd
import argparse
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns


def extract_docking_scores(db_file):
	conn = sqlite3.connect(db_file)
	cur = conn.cursor()
    
	query = cur.execute("SELECT id, docking_score FROM mols")
	cols = ['Name', 'docking_score']
	record = pd.DataFrame.from_records(data = query.fetchall(), columns=cols)
    
	cur.close()
    
	return record


def load_data(logs_db, rmsd_file, physchem_file, sim_scores_file, sa_scores_file, output_dir):
	"""

	param logs_dir: directory with docking log files
	param rmsd_dir: directory with rmsd scores files
	param physchem_file: path to file containing physchem values of mols
	param sim_scores_file: path to file containing similarity scores
	param sa_scores_file: path to file containing synthtic accessibility scores
	param output_dir: outpur directory for combined results
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

	df4.to_csv(os.path.join(output_dir, 'results.csv'), index=False, sep="\t")

	return df4

def plot_docking_scores(crem_mols_data, pp_output, comaprison_data_dir = None):
	"""
	plots docking scores as a histogram
	param crem_mols_data: a pandas dataframe with docking scores
	param pp_output: path to output directory
	param comparison_data_dir: path to directory with docking results to compare with generated molecules' dock scores
	return: None
	"""

	def get_dockscore(pp_dockscores):
    """
    extracts the best docking score
    param pp_dockscores: path to folder. The script expect that there are folders with docking results, `log` folder
    return: list of tuples where there is molecule name on the first position and docking score on the second one
    """
    df_dock = []
    for ll in os.listdir(pp_dockscores):
        with open(os.path.join(pp_dockscores, ll)) as f:
            ds = f.readlines()[26]
            df_dock.append((os.path.splitext(ll)[0], float(ds.split()[1])))
    return df_dock

    def plotting(pp_dock_res, pp_output):
    """
	plots docking scores if comparison data directory if specified
	param pp_dock_res: directory with comparison molecules' docking scores in log files
	param pp_out: path to output directory
    """
    df = pd.DataFrame(columns=['name', 'dockscore', 'status'])
    for ll in os.listdir(pp_dock_res):
        dft = pd.DataFrame(get_dockscore(os.path.join(pp_dock_res, ll, 'log')), columns=['name', 'dockscore'])
        dft['status'] = [ll] * dft.shape[0]
        df = df.append(dft, sort=False)

    sns.histplot(df, x="dockscore", hue="status", element="step")
    plt.title("Distribution of docking score for active, inactive and generated molecules")
    plt.xlabel("Docking score")
    plt.savefig(os.path.join(pp_output, 'dock_scores_histogram.svg'))

    return None

    # if no comparison data, plots a simple histogram of dock scores of generated molecules
	if comaprison_data_dir = None:
		sns.histplot(data=crem_mols_data['docking_score'], kde=True, element='step')
		plt.title('Distribution of docking score for generated molecules')
		plt.xlabel('Docking score')
		plt.savefig(os.path.join(pp_output, 'dock_scores_histogram.svg'))
		return None

	elif comaprison_data_dir is not None:
		if os.path.isdir(comaprison_data_dir):
			if not os.listdir(comaprison_data_dir):
				print('Specified directory is empty.')
				return None
			else:
				pp_dock_res = get_dockscore(comaprison_data_dir)
				plotting(pp_dock_res, pp_out)
		else:
			print('Provided path is not a directory.')
			return None

def main():
	parser = argparse.ArgumentParser(description='Combine data into one file')
	parser.add_argument('-p', '--physchem', metavar='STRING', required=True, help='path to physchem props input file')
	parser.add_argument('-s', '--sim', metavar='STRING', required=True, help='path to similarity scores input file')
	parser.add_argument('-d', '--dock', metavar='FILE', required=True, help='path to docking scores database')
	parser.add_argument('-r', '--rmsd', metavar='FILE', required=True, help='path to rmsd scores input file')
	parser.add_argument('-sa', '--sa', metavar='STRING', required=True, help='path to sa scores input file')
	parser.add_argument('-o', '--out', metavar='DIRECTORY', required=True, help='path to output directory')
	parser.add_argument('-c', '--compare', metavar='DIRECTORY', required=False, default=None, help='path to directory with docking logs for comparison')
	args = parser.parse_args()

	data = load_data(args.dock, args.rmsd, args.physchem, args.sim, args.sa, args.out)
	plot_docking_scores(data, args.out, comaprison_data_dir=args.compare)

if __name__ == '__main__':
	main()
