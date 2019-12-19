import glob
import pandas as pd
import re
import numpy as np
import cooler
import bioframe


def get_files(cooldir):
    cools = glob.glob(cooldir + '/**/*mcool', recursive=True)
    return cools


def get_df(cools, resolution = 100000):
    df = pd.DataFrame([np.array(re.split("[\_\./\-]+", cool))[[-8, -5]].tolist() + [cool]
                       for cool in cools],
                      columns=['cell_line', 'assembly', 'path'])
    c_list = list()
    for i in range(len(df)):
        c_list.append(cooler.Cooler(df.iloc[i]['path'] + '::/resolutions/' + str(resolution)))

    df['cooler'] = pd.Series(c_list)

    return df

def get_genecov(df):
    genecov_dict = dict()
    for assembly in df.assembly.unique():
        bins = df[df.assembly == assembly]['cooler'].iloc[0].bins()[:]
        genecov = bioframe.tools.frac_gene_coverage(bins, assembly)
        genecov_dict[assembly] = genecov['gene_coverage']

    return genecov_dict
