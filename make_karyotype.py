#!/usr/bin/env python3

import Bio
from Bio import SeqIO

gbFile='/work/data/NCBI_plasmids/plasmid/ncbi_plasmid.gbff'

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

for record in SeqIO.parse(gbFile, "genbank"):
    id=record.id.split('.')[0]
    color="255,255,255"
    genus=record.annotations["organism"].split(" ")[0]
    if genus in colors.keys():
        color=colors[genus]
    print("chr - {} {} 1 {} {}".format(id, id, len(record), color))
