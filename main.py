import sys
import cooler
import fileprocessing
import eigendecomp


def main(argv):
    if (argv == None):
        input_data = sys.argv
    else:
        input_data = argv.split()

    if input_data[1].find('.mcool') > -1:
        return get_track(input_data[1], input_data[2])
    else:
        return multi_track(input_data[1], input_data[2])

def get_track(coolfile, n=3, resolution = 100000):
    df = fileprocessing.get_df([coolfile], resolution)
    genecov_dict = fileprocessing.get_genecov(df)
    bins = df['cooler'].iloc[0].bins()[:]
    bins['gene_count'] = genecov_dict[df.assembly.iloc[0]]
    regions = bins.chrom.unique()[:-3]
    track = eigendecomp.cooler_eigendecomp(df['cooler'].iloc[0], bins, regions, operations = ['clip', 'log','cliplog', 'tanh', 'tanhlog'], trackname = 'gene_count', n = n)
    return track

def multi_track(cooldir, n=3, resolution = 100000):
    """
    This analysis will do the analysis of all the coolfiles in a folder.
    """
    cools = fileprocessing.get_files(cooldir)
    df = fileprocessing.get_df(cools, resolution)
    genecov_dict = fileprocessing.get_genecov(df)
    tracks = list()
    for i in range(len(df)):
        print('Analyzing', df['cell_line'].iloc[i])
        bins = df['cooler'].iloc[i].bins()[:]
        bins['gene_count'] = genecov_dict[df.assembly.iloc[i]]
        regions = bins.chrom.unique()[:-3]
        track = eigendecomp.cooler_eigendecomp(df['cooler'].iloc[i], bins, regions, operations = ['clip', 'log','cliplog', 'tanh', 'tanhlog'], trackname = 'gene_count', n = n)
        tracks.append(track)

    return df, tracks


if __name__ == "__main__":
    sys.exit(main(None))
