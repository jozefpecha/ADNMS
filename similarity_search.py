#!/usr/bin/env python3

__author__ = 'Alina Kutlushina'

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import argparse


def get_fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(mol), 2)


def calc_similarity(ref_fp, fps):
    if not ref_fp:
        return None
    sim_scores = DataStructs.BulkTanimotoSimilarity(ref_fp, fps)
    return sim_scores


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Tanimoto index between two vectors of molecules')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=True, type=str,
                        help='path to input smi file. No header. Space- or tab-separated.')
    parser.add_argument('-r', '--ref', metavar='FILENAME', required=True, type=str,
                        help='path to smi file of reference molecules. No header. Space- or tab-separated.')
    parser.add_argument('-o', '--out', metavar='FILENAME', required=True, type=str, help='path to output file.')
    args = parser.parse_args()

    df_res = pd.DataFrame(columns=[i.strip().split()[1] for i in open(args.ref).readlines()],
                          index=[i.strip().split()[1] for i in open(args.input).readlines()])

    rfps = [get_fp(i.strip().split()[0]) for i in open(args.ref).readlines()]
    if len(rfps) != df_res.shape[1]:
        list_rfps = []
        list_rnames = []
        with open(args.ref) as fr:
            for line in fr.readlines():
                rsml, rname = line.strip().split()
                rfp = get_fp(rsml)
                if not rfp:
                    continue
                list_rfps.append(rfp)
                list_rnames.append(rname)
        df_res.columns = list_rnames

    with open(args.input) as fi:
        for line in fi.readlines():
            isml, iname = line.strip().split()
            bulk_tanimoto = calc_similarity(get_fp(isml), rfps)
            if bulk_tanimoto:
                df_res.loc[iname] = bulk_tanimoto
    df_res.to_csv(args.out)