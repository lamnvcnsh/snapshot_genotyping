import pickle
import numpy as np
import re
import pandas as pd

tables = pickle.load(open('resource/tables.pdata', 'rb'))


def get_alleles(tables, gene):
    g_tables = tables['genotype_marker_table']
    genotypes = g_tables[g_tables['gene'] == gene]
    alleles = [x.split("/")[0] for x in genotypes['genotype']]
    alleles = [re.sub('_[0-9]', '', x) for x in alleles]
    alleles = np.unique(np.array(alleles))
    alleles = [x for x in alleles if x != 'other']

    return alleles


def format_hla_genotype(genotype, gene, tables_file='resource/tables.pdata'):
    tables = pickle.load(open(tables_file, 'rb'))
    star_alleles = get_alleles(tables, gene)
    genotypes = ['']*len(star_alleles)
    genes = [gene]*len(star_alleles)

    for i, allele in enumerate(star_alleles):
        counts = genotype.count(allele)
        if counts == 2:
            genotypes[i] = genotype
        elif counts == 1:
            genotypes[i] = f'other/{allele}'
        else:
            genotypes[i] = 'other/other'

    print(genotypes)
    return pd.DataFrame.from_dict({'gene': genes,
                                   'rsid': star_alleles,
                                   'genotype': genotypes,
                                   'location': ''})

