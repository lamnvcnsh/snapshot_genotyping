import pandas as pd
import numpy as np
import re
import itertools
import pickle

def load_snapshot_file(snapshot_file):
    df = pd.read_csv(snapshot_file, sep='\t').fillna("")
    df['Allele 2'] = np.where(df['Allele 2'] == "", df['Allele 1'], df['Allele 2'])
    df['Genotype'] = df['Allele 1'] + df['Allele 2']
    df['Gene'] = [x[0] for x in df['Marker'].str.split("_")]
    df['ExtractedSample'] = [extract_sample_name(x) for x in df['Sample Name']]
    df['SampleStatus'] = [is_valid_sample_name(x) for x in df['Sample Name']]
    
    df = df[['Sample File', 'Sample Name', 'ExtractedSample', 'SampleStatus', 'Panel', 
             'Gene', 'Marker', 'Genotype', 'Allele 1', 'Allele 2']]
    return df

def is_valid_sample_name(sample_name, p = '_S[0-9]+$'):
    pattern = re.compile(p)
    if pattern.search(sample_name):
        return(True)
    else:
        return(False)
    
def extract_sample_name(sample_name, p = '_S[0-9]+$'):
    pattern = re.compile(p)
    if pattern.search(sample_name):
        return(pattern.sub('', sample_name))
    else:
        print(f'{sample_name} is not valid! Keep as origin')
        return(sample_name)
    
def is_valid_alleles(a,b,x,y):
    # a and b are allele 1 and allele 2
    # x and y are wildtype and mutant alleles
    return True if (a in x or a in y) and (b in x or b in y) else False

def format_genotype(a,b,x,y):
    if (a in y) and (b in x):
        return f'{b}/{a}'
    else:
        return f'{a}/{b}'
    
def merge_bin_and_definition(bin_text_file, tables):
    bin_data = load_snapshot_file(bin_text_file)
    bin_sample = bin_data.merge(tables['marker_table'], how='left', 
                                    left_on='Marker', right_on='marker')
    bin_sample['AlleleStatus'] = [is_valid_alleles(a,b,x,y)
                                  for a,b,x,y in zip(bin_sample['Allele 1'], 
                                                     bin_sample['Allele 2'],
                                                     bin_sample['wildtype'],
                                                     bin_sample['mutant'])]
    bin_sample['FormatedGenotype'] = [format_genotype(a,b,x,y)
                                      for a,b,x,y in zip(bin_sample['Allele 1'], 
                                                         bin_sample['Allele 2'],
                                                         bin_sample['wildtype'],
                                                         bin_sample['mutant'])]
    bin_sample.loc[bin_sample.Genotype == '', ('FormatedGenotype')] = 'error'
    bin_sample['IsValidGenotype'] = np.where(bin_sample['Genotype'] != "", True, False)
    bin_sample['IsValidMarker'] = [x in tables['marker_table']['marker'].tolist() for x in bin_sample['Marker']]
    bin_sample['IsCalled'] = np.where(bin_sample['Genotype'] != "", True, False)
    return(bin_sample)

def is_heterozygous(allele):
    if '/' in allele:
        allele = '{}{}'.format(allele[-3], allele[-1])
#     print(allele)
    homozygous = ['AA', 'CC', 'GG', 'TT']
    return(allele not in homozygous)

def generate_pattern(allele):
    if is_heterozygous(allele):
        pattern = f"({allele}|{allele[1]}{allele[0]})"
    else:
        pattern = allele
    return pattern

def generate_all_pattern(alleles):
    return '_'.join([generate_pattern(allele) for allele in alleles])

def get_sample_marker_by_gene(tables, sample_data, gene):
    reference_marker = get_reference_marker_by_gene(tables, gene)
    sample_marker = sample_data[sample_data['Gene'] == gene][['Marker', 'Genotype']]
    sample_marker_final = reference_marker.merge(sample_marker, how='left', 
                                                 left_on="marker", right_on='Marker').fillna("")
    sample_marker_final['Genotype'] = np.where(sample_marker_final['Genotype'] == "", 
                                               sample_marker_final['wildtype'] + sample_marker_final['wildtype'],
                                               sample_marker_final['Genotype'])
    return(sample_marker_final[['marker', 'Genotype']])

def get_reference_marker_by_gene(tables, gene):
    gene_table = tables['gene_table']
    marker_table = tables['marker_table']
    gene_id = gene_table[gene_table['gene'] == gene]['gene_id'].values[0]
    return(marker_table[marker_table['gene_id'] == gene_id][['marker', 'wildtype']])

def generate_marker_target(tables, gene):
    genotype_marker_table = tables['genotype_marker_table']
    df = genotype_marker_table[genotype_marker_table['gene'] == gene]
    genotype_id = df['genotype_id'].tolist()[0]
    maker_names = genotype_marker_table[genotype_marker_table['genotype_id'] == genotype_id]['marker']
    
    df = df[['gene', 'genotype', 'marker', 'value']].pivot(index='genotype', columns='marker', values='value')
    df = df[maker_names]
    df['target'] = df[maker_names].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    return(df)

def is_valid_marker(gene_data, tables, gene):
    sample_marker = sorted(gene_data['Marker'].tolist())
    defined_marker = sorted(tables['genotype_marker_table']
                            [tables['genotype_marker_table']['gene'] == gene]['marker'].unique().tolist())
    return sample_marker == defined_marker
        
def is_valid_genotype(gene_data):
    sample_genotype = gene_data['Genotype'].tolist()
    return '' not in sample_genotype

def is_valid_gene_data(gene_data, tables, gene):
    return is_valid_marker(gene_data, tables, gene) and is_valid_genotype(gene_data)


def non_star_gene_list():
    # return genes are not required to call star alleles
    gene_list = ['ABCB1', 'ANKK1', 'COMT', 'DRD2', 'DRD3', 'EPHX1', 
                 'FKBP5', 'HTR1A','HTR2A', 'HTR2C', 'MC4R', 'OPRM1', 'SCN1A']
    return(gene_list)

def rename_columns(df):
    df.columns = ['sample_name', 'gene', 'rsid', 'genotype', 
                  'location', 'file_name', 'allele_status', 
                  'is_valid_marker', 'is_valid_genotype']
    return(df)

def genotype_variant(gene_data):
    df = gene_data[['ExtractedSample', 'Gene', 'rsid', 'FormatedGenotype', 
                    'nt_change', 'Sample File', 'AlleleStatus', 
                    'IsValidMarker', 'IsValidGenotype', 'IsCalled']]
    
    df.columns = ['sample_name', 'gene', 'rsid', 'genotype', 
                  'location', 'file_name', 'allele_status', 
                  'is_valid_marker', 'is_valid_genotype', 'is_called']
    return(df)

def genotype_star_allele(sample_data, gene_data, tables, sample, gene):
    file_name = '|'.join(list(dict.fromkeys(gene_data['Sample File'])))
    
    allele_status = all(gene_data['AlleleStatus'])
    is_val_marker = False
    is_val_genotype = False
    is_called = False
    genotype = 'error'

    if is_valid_genotype(gene_data):
        is_val_genotype = True
        if is_valid_marker(gene_data, tables, gene):
            is_val_marker = True
            marker_target = generate_marker_target(tables, gene)
            pattern = generate_all_pattern(
                get_sample_marker_by_gene(tables=tables,
                                          sample_data=sample_data, gene=gene)['Genotype'])
            matched = [bool(re.fullmatch(pattern=pattern, string=target)) 
                       for target in marker_target['target']]
            genotype = '|'.join(marker_target.iloc[np.where(matched)].index.values.tolist())
            
    
            if sum(matched) >= 1:
                is_called = True
        
            if (genotype == ''):
                genotype = 'other/other'
        
    df = pd.DataFrame({'sample_name': sample, 
                       'gene': gene, 
                       'rsid': [''], 
                       'genotype': genotype, 
                       'location': [''], 
                       'file_name': file_name, 
                       'allele_status': allele_status, 
                       'is_valid_marker': is_val_marker,
                       'is_valid_genotype': is_val_genotype ,
                       'is_called': is_called})
    
    return(df)

def genotype_calling(bin_sample, tables):
    samples = bin_sample['ExtractedSample'].unique().tolist()
    genotype_df = pd.DataFrame()

    for sample in samples:
        # extract sample data
        sample_data = bin_sample[bin_sample['ExtractedSample'] == sample]

        genes = sample_data['Gene'].unique().tolist()

        for gene in genes:
            # extract gene data
            gene_data = sample_data[sample_data['Gene'] == gene]

            if gene in non_star_gene_list():
                tmp_df = genotype_variant(gene_data)
            else:
                tmp_df = genotype_star_allele(sample_data, gene_data, tables, sample, gene)

            genotype_df = pd.concat([genotype_df, tmp_df])

    genotype_df.sort_values(['sample_name','gene'], inplace=True)
    genotype_df.reset_index(drop=True, inplace=True)
    genotype_df['QC'] = genotype_df[['allele_status', 'is_valid_genotype', 'is_valid_marker', 'is_called']].sum(axis=1)
    genotype_df['QC'] = np.where(genotype_df['QC'] == 4, 'Pass', 'Warning')
    
    return genotype_df

def main(bin_text_file, definition_pickle = 'resource/tables.pdata', excel =False):
    import os
    result = 'result/genotype_calling'
    if not (os.path.exists(result) and os.path.isdir(result)):
        os.mkdir(result)

    tables = pickle.load(open(definition_pickle, 'rb'))
    bin_sample = merge_bin_and_definition(bin_text_file, tables)
    genotype_df = genotype_calling(bin_sample, tables)

    if excel:
        base = os.path.basename(bin_text_file)
        outfile = os.path.splitext(base)[0]
        genotype_df.to_excel(f'{result}/{outfile}_result.xlsx', index=False)

    return genotype_df


if __name__ == '__main__':
    import sys
    df = main(bin_text_file=sys.argv[1], excel=True)
