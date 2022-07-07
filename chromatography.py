# script to plot chromatogram 
# author: Lam Nguyen


def load_fsa_file(file):
    from Bio import SeqIO
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    filename = os.path.basename(file)
    record = SeqIO.read(file, 'abi')
    
    return(filename, record)


def build_linear_model(record):
    from sklearn.linear_model import LinearRegression
    import numpy as np
    import pandas as pd
    
    peaks_as_base = record.annotations['abif_raw']['Peak12']
    matched_peaks = [True if x == '1' else False for x in 
                     record.annotations['abif_raw']['Peak20'].decode().split(",")]
    matched_points = [x for x in record.annotations['abif_raw']['Peak2']]
    sizes = np.array(peaks_as_base)[matched_peaks].reshape(-1, 1)
    points = np.array(matched_points)[matched_peaks].reshape(-1, 1)
    df = pd.DataFrame(np.array([np.sqrt(points), sizes]).reshape(2,-1)).T
    df.columns = ['points', 'sizes']
    df['points2'] = df['points']**2
    df['points3'] = df['points']**3

    # print(df)
    
    # build linear model 
    # model = LinearRegression().fit(df[['points', 'points2']], df['sizes'])
    model = LinearRegression().fit(df[['points', 'points2', 'points3']], df['sizes'])
#     print(model.score(points, sizes))
    return model, np.array([np.min(sizes), np.max(sizes)])

def load_snapshot_file(snapshot_file):
    import pandas as pd
    import numpy as np
    df = pd.read_csv(snapshot_file, sep='\t')
    df['Allele 2'] = df['Allele 2'].fillna('')
    df['Allele 2'] = np.where(df['Allele 2'] == "", df['Allele 1'], df['Allele 2'])
    df['Genotype'] = df['Allele 1'] + df['Allele 2']
    df['Gene'] = [x[0] for x in df['Marker'].str.split("_")]
#     df['Customer'] = [x[-1].strip() for x in df['Sample Name'].str.split("-")]
    df['Customer'] = [x[:-3] for x in df['Sample Name']]
    
    df = df[['Sample File', 'Sample Name', 'Customer', 'Panel', 
             'Gene', 'Marker', 'Genotype', 'Allele 1', 'Allele 2', 
             'Size 1', 'Size 2', 'Height 1', 'Height 2']]
    df['x'] = np.where(np.isnan(df['Size 2']), df['Size 1'], 
                                (df['Size 1'] + df['Size 2'])/2)

    df['y'] = df[['Height 1', 'Height 2']].apply(np.nanmax, axis=1)
    df['text'] = df['Marker'] + ':' + df['Genotype']
    return df


def plot_chromatogram(file, bin_data, plot_folder=None, c_limit=[20, 100], iscustom=False, ref=False):
    from collections import defaultdict
    import matplotlib.pyplot as plt
    import numpy as np
    # load fsa file
    filename, record = load_fsa_file(file)
    # build model from record
    # model, d_limit = build_linear_model(record)

    channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
    colors = ['blue', 'green', 'black', 'red']
    labels = ['G', 'A', 'C', 'T']

    trace = defaultdict(list)

    for chanel in channels:
        trace[chanel] = record.annotations['abif_raw'][chanel]

    # predict size
    fontsize = 15
    # total_points = record.annotations['abif_raw']['Scan1']
    # data_points = np.arange(0, total_points).reshape(-1, 1)
    # data_points = np.sqrt(data_points)
    # data = pd.DataFrame(np.array([data_points]).reshape(1, -1)).T
    # data.columns = ['points']
    # data['points2'] = data['points'] ** 2
    # data['points3'] = data['points'] ** 3
    # size_pred = model.predict(data)
    size_pred = record.annotations['abif_raw']['SMap2']

    # get size and height
    annotation = bin_data[bin_data['Sample File'] == filename]
    # print(annotation)

    # check custom range
    if iscustom:
        limit = c_limit.copy()
    else:
        limit = [0, 120]

    # plot chromagraphy
    fig, ax = plt.subplots()
    ax.set_xlim(limit)
    ax.grid(linestyle='--', alpha=0.5, color='gray')
    for data, color, label in zip(channels, colors, labels):
        ax.plot(size_pred, trace[data], color=color, label=label)
    for x, y, s in zip(annotation['x'], annotation['y'], annotation['Marker']):
        ax.text(x, y + 10, s, ha='center', size=8, bbox=dict(fc='white', alpha=0.5))

    if ref:
        ref_channel = record.annotations['abif_raw']['DATA105']
        ax.plot(size_pred, ref_channel, color='orange', alpha=0.3, label='RefSize')

    ax.set_title(f'{filename}\n', size=fontsize, fontdict={'fontweight': 'bold'})

    plot_file = f'{filename}.png'
    ax.legend(prop={'size': fontsize}, loc='upper right', ncol=5)
    ax.margins(x=0, y=0)
    fig.set_size_inches([16, 4])
    fig.savefig(f'{plot_folder}/{plot_file}', dpi=150)

    return plot_file


def main(file, fsa_folder = 'fsas'):
    import os
    bin_data =  load_snapshot_file(file)
    fsa_files = bin_data['Sample File'].unique()
    for fsa in fsa_files:
        full_fsa = f'{fsa_folder}/{fsa}'

        if os.path.exists(full_fsa):
            plot_chromatogram(full_fsa, bin_data, iscustom=True, ref=True, plot_folder='results/chromatogram_plot')
        else:
            print(f'{full_fsa} not existed!')

if __name__ == '__main__':
    import sys
    # bin_file = sys.argv[1]
    # out_folder = sys.argv[2]
    bin_file = '/Volumes/GoogleDrive/My Drive/PhD/Works/SPMED/Genotyping/SNAPshot/fsa/CYP2C19_BIN SET/CYP2C19_BIN SET.txt'
    out_folder = '/Volumes/GoogleDrive/My Drive/PhD/Works/SPMED/Genotyping/SNAPshot/fsa/CYP2C19_BIN SET/CYP2C19_fsa_file/'
    main(bin_file, out_folder)




