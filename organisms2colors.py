#!/usr/bin/env python3

import csv

import pandas as pd

print("reading in organism names")
odf = pd.read_csv("organisms.txt", sep='\t', header=None, 
        names=["name", "organism"])

print("reading in taxonomy lineages")
mdf = pd.read_csv("lineages-2017-05-04.csv", dtype=str)

print("merging data frames")
df = pd.merge(odf, mdf, left_on="organism", right_on="species")

df.to_csv('/tmp/o.tsv', sep='\t', quoting=csv.QUOTE_NONE, index=False)
