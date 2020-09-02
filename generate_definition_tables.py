import pandas as pd
import numpy as np
import re
import itertools
import pickle

def get_combinations(star_allele):
    homozygous = [f'{x}/{x}' for x in star_allele]
    heterozygous = ["/".join(map(str, comb)) for comb in itertools.combinations(star_allele, 2)]
    combinations = homozygous + heterozygous
    return(pd.DataFrame({"genotype_id": np.arange(1,len(combinations)+1), "genotype": combinations}))

def get_variant_combines(data, genotype):
    
    allele_1 = genotype.split("/")[0]
    allele_2 = genotype.split("/")[1]
    allele_1_idx = np.where(data.index == allele_1)[0][0]
    allele_2_idx = np.where(data.index == allele_2)[0][0]
    
    variant_combines = [f'{x}{y}' for x,y in zip(data.iloc[allele_1_idx], data.iloc[allele_2_idx])]
    variant_name = [var for var in data.columns]
    
    return pd.DataFrame({'genotype': [genotype]*len(variant_name),
                         'marker': variant_name,
                         'value': variant_combines})

def get_marker_info(excel_file, sheetname, start_col=3, nrows=6):
    tmp_df = pd.read_excel(excel_file, sheet_name=sheetname, nrows=nrows)
    # remove first three columns and transform data
    tmp_df = tmp_df.iloc[:,np.arange(start_col, tmp_df.shape[1])].T
    
    # set column names and reset index
    columns = ["star_allele","aa_change", "nt_change", 'rsid', "wildtype", "mutant"]
    tmp_df.columns = columns
    tmp_df['marker'] = tmp_df.index
    tmp_df.reset_index(drop=True, inplace=True)
    marker_df = tmp_df[["marker"]+ columns]
    marker_df.replace({"\n":"_"}, regex=True, inplace=True)
    
    return(marker_df)

def get_definition_table(excel_file, gene):
    marker_names = get_marker_info(excel_file, gene)['marker']
    tmp = pd.read_excel(excel_file, sheet_name=gene, skiprows=7, header=None)
    tmp.index = tmp.iloc[:,1]
    tmp.drop([0,1,2], axis=1, inplace=True)
    tmp.columns = [x for x in marker_names]
    tmp.fillna("", inplace=True)
    
    return(tmp)


def generate_gene_table(excel_file):
    sheet_names = pd.ExcelFile(excel_file).sheet_names
    sheet_names = [x for x in sheet_names if x not in ['etc.']]
    
    gene_table = pd.DataFrame({"gene_id": np.arange(1,len(sheet_names)+1), 
                             "gene": sheet_names})
    return(gene_table)

def generate_marker_table(excel_file):
    
    gene_table = generate_gene_table(excel_file)
    marker_table = pd.DataFrame()

    for idx, gene in zip(gene_table['gene_id'], gene_table['gene']):
        marker_tmp = get_marker_info(excel_file, gene)
        marker_tmp['gene_id'] = [idx]*marker_tmp.shape[0]
        marker_table = pd.concat([marker_table, marker_tmp], sort=False)
        marker_table.reset_index(drop=True, inplace=True)
        marker_table.fillna("-", inplace=True)
        marker_table['marker_id'] = marker_table.index + 1
    
    return(marker_table)

def generate_genotype_table(excel_file):
    gene_table = generate_gene_table(excel_file)
    
    genotype_table = pd.DataFrame()
    
    for idx, gene in zip(gene_table['gene_id'], gene_table['gene']):
        definition_table = get_definition_table(excel_file, gene)
        genotype_tmp = get_combinations(definition_table.index)
        genotype_tmp['gene_id'] = [idx]*genotype_tmp.shape[0]
        genotype_table = pd.concat([genotype_table, genotype_tmp])
    genotype_table['genotype_id'] = np.arange(1, genotype_table.shape[0]+1)
    
    return(genotype_table)

def generate_genotype_by_marker(data, gene):
    
    combines = get_combinations(data.index)
    final_df = pd.DataFrame()
    for i in combines.index:
        genotype_id = combines.iloc[i][0]
        genotype = combines.iloc[i][1]
        
        tmp_df = get_variant_combines(data, genotype)
        tmp_df['gene'] = gene
        
        final_df = pd.concat([final_df, tmp_df])
    final_df['id'] = np.arange(1, final_df.shape[0]+1)   
    return(final_df)

def generate_genotype_marker_table(excel_file, gene_table):
    genotype_marker_table = pd.DataFrame()

    for idx, gene in zip(gene_table['gene_id'], gene_table['gene']):
        defination_table = get_definition_table(excel_file, gene)
        genotype_marker_tmp = generate_genotype_by_marker(defination_table, gene)
        
        genotype_marker_table = pd.concat([genotype_marker_table, genotype_marker_tmp])
    genotype_marker_table['id'] = np.arange(1, genotype_marker_table.shape[0]+1)
    
    return(genotype_marker_table)

def load_all_reference_tables(excel_file):
    print('loading gene tables...')
    gene_table = generate_gene_table(excel_file)
    print('loading marker tables...')
    marker_table = generate_marker_table(excel_file)
    print('loading genotype tables...')
    genotype_table = generate_genotype_table(excel_file)
    print('loading genotype marker tables...')
    genotype_marker_table = generate_genotype_marker_table(excel_file, gene_table)
    
    genotype_marker_table = genotype_marker_table.merge(gene_table, how='left', on='gene')
    genotype_marker_table = genotype_marker_table.merge(marker_table[['marker', 'marker_id']], 
                                                        how='left', on='marker')
    genotype_marker_table = genotype_marker_table.merge(genotype_table, 
                                                        how='left', on=['gene_id', 'genotype'])
#     genotype_marker_table = genotype_marker_table[['id', 'gene_id', 'genotype_id', 'marker_id', 'value']]
    print('------DONE-----')
    
    return {'gene_table':gene_table, 
            'marker_table':marker_table, 
            'genotype_table':genotype_table, 
            'genotype_marker_table':genotype_marker_table}

def main(excel_file, output='tables.pdata'):

	import os
	outfolder='resource'

	if not (os.path.exists(outfolder) and os.path.isdir(outfolder)):
		os.mkdir(outfolder)

	tables = load_all_reference_tables(excel_file)

	# It will overwrite the previous tables
	file = open(f'{outfolder}/{output}', 'wb')
	pickle.dump(tables, file)

if __name__ == '__main__':
	import sys
	main(sys.argv[1])
