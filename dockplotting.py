import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_dockscore(pp_dockscores):
    """
    extract the best docking score
    :param pp_dockscores: path to folder. The script expect that there are folders with docking results, `log` folder
    :return: list of tuples where there is molecule name on the first position and docking score on the second one
    """
    df_dock = []
    for ll in os.listdir(pp_dockscores):
        with open(os.path.join(pp_dockscores, ll)) as f:
            ds = f.readlines()[26]
            df_dock.append((os.path.splitext(ll)[0], float(ds.split()[1])))
    return df_dock


def plotting(pp_dock_res, pp_output):
    df = pd.DataFrame(columns=['name', 'dockscore', 'status'])
    for ll in os.listdir(pp_dock_res):
        dft = pd.DataFrame(get_dockscore(os.path.join(pp_dock_res, ll, 'log')), columns=['name', 'dockscore'])
        dft['status'] = [ll] * dft.shape[0]
        df = df.append(dft, sort=False)

    sns.histplot(df, x="dockscore", hue="status", element="step")
    plt.title("Distribution of docking score for active, inactive and generated molecules")
    plt.xlabel("Docking score")
    plt.savefig(pp_output)