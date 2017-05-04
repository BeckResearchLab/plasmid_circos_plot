#!/usr/bin/env python3

import csv

import pandas as pd

print("reading in contig names and lengths")
# read the file that contains lengths and names
#chr - NC_003894 NC_003894 1 9253 145,255,145
kdf = pd.read_csv("karyotype.txt", sep=' ', header=None, 
        names=["chr", "ignore", "contig", "name",
                "start", "end", "color"])
kdf["length"] = (kdf.end - kdf.start) + 1
kdf.drop(["chr", "ignore", "contig", "color", "start", "end"], 1, inplace=True)

print("read filtered matches")
# read the matches and fix up the names
mdf = pd.read_csv("filtered_matches.txt", sep='\t')
fixed = mdf["qseqid"].apply(lambda x : x.split("|")[1].split(".")[0])
mdf["qseqid"] = fixed
fixed = mdf["sseqid"].apply(lambda x : x.split("|")[1].split(".")[0])
mdf["sseqid"] = fixed
mdf.rename(columns = { "length" : "alength" }, inplace=True)

print("merging data frames")
# join mdf twice, once with lengths on qseqid and then on sseqid
df = pd.merge(mdf, kdf, left_on="qseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "length" : "qlength"}, inplace=True)
df = pd.merge(df, kdf, left_on="sseqid", right_on="name")
df.drop("name", 1, inplace=True);
df.rename(columns = { "length" : "slength"}, inplace=True)

print("filtering out redundant plasmids")
# remove any matches that are > 50% length of either
qscrit = df["alength"] < df["slength"] * .5
qqcrit = df["alength"] < df["qlength"] * .5
scrit = df["alength"] > 5000
print(df.size)
df = df[(qscrit) & (qqcrit) & (scrit)]
print(df.size)

print("writing to 'links.txt'")
links = df[["qseqid", "qstart", "qend", "sseqid", "sstart", "send"]].copy()
links.to_csv("links.txt", sep='\t', quoting=csv.QUOTE_NONE, header=False, index=False)
