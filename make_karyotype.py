#!/usr/bin/env python3

import Bio
from Bio import SeqIO
import pandas as pd

gbFile = '/work/data/NCBI_plasmids/plasmid/ncbi_plasmid.gbff'

colors = {
    "Staphylococcus" : "84,48,5",
    "Enterococcus" : "140,81,10",
    "Lactobacillus" : "191,129,45",
    "Bacillus" : "223,194,125",
    "Enterobacter" : "199,234,229",
    "Escherichia" : "128,205,193",
    "Klebsiella" : "53,151,143",
    "Salmonella" : "1,102,94",
    "Yersinia" : "0,60,48",
    "Borrelia" : "246,232,195",
    "Borreliella" : "246,232,195"
}

data = []

print("reading genbank file")
for record in SeqIO.parse(gbFile, "genbank"):
    id = record.id.split('.')[0]
    color = "255,255,255"
    genus = record.annotations["organism"].split(" ")[0]
    genus = genus.replace("[", "").replace("]", "")
    taxonomy = '_'.join(record.annotations["taxonomy"])
    if genus in colors.keys():
        color = colors[genus]
    data.append(["chr", "-", id, id, 1, len(record), color, taxonomy])
    #print("chr - {} {} 1 {} {} {}".format(id, id, len(record), color, taxonomy))

print("making dataframe")
df = pd.DataFrame(data, 
        columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb", "taxstring"])
print("saving txt file")
df.to_csv("karyotype_rand.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)
print("sorting by taxonomy and saving")
df.sort_values(by="taxstring", inplace=True)
df.to_csv("karyotype_tax.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)
print("sorting by length and saving")
df.sort_values(by="seqlen", inplace=True)
df.to_csv("karyotype_len.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)

