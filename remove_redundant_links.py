import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Cut out non-unique plasmid connections in links file")

parser.add_argument("Infile", help="Path to existing links file")
parser.add_argument("Outfile", help="Desired name of output links file")

args = parser.parse_args()

columnNames = ['qlocus', 'qstart', 'qend', 'slocus', 'sstart', 'send', 'color']
df = pd.read_table(args.Infile, sep="\t", names=columnNames, dtype={'qstart': int, 'qend': int})

outDF = df.drop_duplicates(['qlocus', 'slocus'])

outDF = outDF[(outDF['qend'] - outDF['qstart']) > 5000].sample(frac=1)

outDF.to_csv(args.Outfile, sep="\t", index=False)