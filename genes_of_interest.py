#!/usr/bin/env python3
import pickle

with open('data/gene_id.pkl', 'rb') as f:
    id_to_gene = pickle.load(f)

WANTED = ['KRAS', 'NRAS', 'DIS3', 'TP53', 'FAM46C', 'BRAF',
          'TRAF3', 'PRDM1', 'IRF4', 'INTS12', 'IDH1',
          'FGFR3', 'PTPN11', 'EZR']

for G in WANTED:
    for key, value in id_to_gene.items():
        if value == G:
            print("%s\t%s" % (value, key))
