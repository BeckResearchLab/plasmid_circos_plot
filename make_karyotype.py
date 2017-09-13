#!/usr/bin/env python3

import Bio
from Bio import SeqIO
import pandas as pd

gbFile = '/work/data/NCBI_plasmids/plasmid/ncbi_plasmid.gbff'

colors = {
    "Yersinia" : "50,0,15",
    "Escherichia" : "103,0,31",
    "Klebsiella" : "178,24,43",
    "Staphylococcus" : "214,96,77",
    "Salmonella" : "244,165,130",
    "Acinetobacter" : "253,219,199",
    "Enterococcus" : "209,229,240",
    "Enterobacter" : "146,197,222",
    "Shigella" : "67,147,195",
    "Citrobacter" : "33,102,172",
    "Borrelia" : "5,48,97",
    "Lactobacillus" : "5,20,45"
}

data = []

print("reading genbank file")
for record in SeqIO.parse(gbFile, "genbank"):
    id = record.id.split('.')[0]
    color = "255,255,255"
    genus = record.annotations["organism"].split(" ")[0]
    taxonomy = '_'.join(record.annotations["taxonomy"])
    if genus in colors.keys():
        color = colors[genus]
    data.append(["chr", "-", id, id, 1, len(record), color, taxonomy])
    #print("chr - {} {} 1 {} {} {}".format(id, id, len(record), color, taxonomy))

print("making dataframe")
df = pd.DataFrame(data, 
        columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb", "taxstring"])
print("saving txt file")
df.to_csv("karyotype.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)
print("sorting by taxonomy and saving")
df.sort_values(by="taxstring", inplace=True)
df.to_csv("karyotype_tax.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)
print("sorting by length and saving")
df.sort_values(by="seqlen", inplace=True)
df.to_csv("karyotype_len.txt", columns=["chr", "-", "id", "id2", "1", "seqlen", "rgb"],
        sep=' ', header=False, index=False)

