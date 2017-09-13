#!/usr/bin/env python3

import Bio
from Bio import SeqIO

gbFile='/work/data/NCBI_plasmids/plasmid/ncbi_plasmid.gbff'

for record in SeqIO.parse(gbFile, "genbank"):
    id=record.id.split('.')[0]
    print("{}\t{}".format(id, record.annotations["organism"]), end='')
    #for level in record.annotations["taxonomy"]:
    #    print("\t{}".format(level), end='')
    print("")
