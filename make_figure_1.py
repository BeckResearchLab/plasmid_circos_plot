import pandas as pd
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

def get_lengths(genbank):
    lenDict = {}
    for record in SeqIO.parse(genbank, "genbank"):
        recordName = record.id.split(".")[0]
        lenDict[recordName] = len(record.seq)
    return lenDict

def expand_covered(old, new, direct):
    if (direct == 'start') and (new < old):
        return new
    elif (direct == 'end') and (new > old):
        return new
    else:
        return old


def coverage(locus, links, length):
    try:
        covered = np.ones(int(length))
    except KeyError:
        print("Locus:{}\tlength:{}\n".format(locus, length))
        exit()
    t = open("test.txt", 'w')
    t.write(str(len(str(covered))))
    regions = links[links['plasmid1'] == locus][['start1', 'end1']]
    t.write(regions.to_string())
    t.write("\n")
    for i, row in regions.iterrows():
        covered[row['start1']-1:row['end1']-1] = np.zeros(row['end1']-row['start1'])
        unique, counts = np.unique(covered, return_counts=True)
        x = dict(zip(unique, counts))
        t.write(str(x[0]))
        t.write("\n")
    exit()
    return 1-(covered.sum()/length)

def scoverage(locus, links):
    return(links[links['plasmid1'] == locus]['start1'].min().item())

def ecoverage(locus, links):
    return(links[links['plasmid1'] == locus]['end1'].max().item())

def numCov(locus, links):
    return(links[links['plasmid1'] == locus]['plasmid2'].count().item())

def import_links(linkFile):
    linkCols = ['plasmid1', 'start1', 'end1', 'plasmid2', 'start2', 'end2', 'relationship']
    linkDF = pd.read_table(linkFile, names=linkCols, sep="\t")
#    linkDF = linkDF.replace(to_replace="color=red", value="interspecies")
#    linkDF = linkDF.replace(to_replace="color=black", value="intraspecies")
    linkDF['length'] = linkDF['end1'] - linkDF['start1']
    plasDF = pd.DataFrame(linkDF.plasmid1.unique(), columns=['locus'])
    return linkDF, plasDF


def plotter(plasDF, linkDF, outfile):
#    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(6, 6))
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(223)
    ax3 = plt.subplot(222)
    ax4 = plt.subplot(224)

    ax1.hist(plasDF['cprop'], bins=20, color='c')
    n2, bins2, patches2 = ax2.hist(plasDF['links'], bins=20, color='m')
    linkMax = plasDF['links'].max()
    ax3data = linkDF.groupby('relationship').count()['plasmid1']
    deduplicated = linkDF.drop_duplicates(subset=['plasmid1', 'plasmid2'])
    connectionCounts = deduplicated.groupby('plasmid1').count()
    bar1 = ax3.bar([0.3, 0.3], height=[ax3data['color=black'], ax3data['color=black']], color='c', width=[0.2, 0.2],
            align='center')
    bar2 = ax3.bar([0.3, 0.3], height=[ax3data['color=red'], ax3data['color=red']],
            bottom=[ax3data['color=black'], ax3data['color=black']], color='m', width=[0.2, 0.2], align='center')
    ax3.legend((bar1, bar2), ("Intraspecies", "Interspecies"))
    n4, bins4, patches4 = ax4.hist(connectionCounts['plasmid2'], bins=20, color='c')
    print(map(lambda x: int(x), bins4))

#    bxDta = [linkDF[linkDF['relationship'] == 'intraspecies']['length'], linkDF[linkDF['relationship'] == 'interspecies']['length']]
#    ax3.boxplot(bxDta, labels=['intraspecies', 'interspecies'], boxprops=dict(color='c'),
#                flierprops=dict(marker='.', markerfacecolor='white', markeredgecolor='black'),
#                whiskerprops=dict(linestyle='-', color='c'), medianprops=dict(color='black'))

    ax1.set_xlabel("Mosaic proportion of plasmid")
    ax1.set_ylabel("Number of plasmids")

    ax2.set_yscale("symlog")
    ax2.set_xlabel("Links per plasmid")
    ax2.set_xticks(bins2)
    counter = 0
    for label in ax2.xaxis.get_ticklabels():
        if counter % 5 != 0:
            label.set_visible(False)
        counter += 1
    ax2.set_ylabel("Number of plasmids")
    ax2.set_xlim([0, linkMax])

#    ax3.set_yscale("log")
    ax3.set_ylabel("Number of Links")
    ax3.set_xlim(left=0, right=2)
#    ax3.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax3.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2e'))
    ax3.set_xticklabels(["Plasmid links"])
    ax3.set_xticks([0.3])

    ax4.set_xlabel("Unique plasmids linked")
    ax4.set_xticks(map(lambda x: int(x), bins4))
    for label in ax4.xaxis.get_ticklabels():
        if counter % 5 != 0:
            label.set_visible(False)
        counter += 1
    ax4.set_ylabel("Number of plasmids")

    plt.tight_layout()
    plt.savefig(outfile)

    return connectionCounts['plasmid2'].median()

parser = argparse.ArgumentParser(description="Make figure 1")

parser.add_argument("LinkFile", help="Links file created for circos")
parser.add_argument("GenBank", help="Genbank format file of original plasmids")
parser.add_argument("-o", "--out", default="Figure_1.png", help="Name for plot")

args = parser.parse_args()

lenDict = get_lengths(args.GenBank)
links, plasmids = import_links(args.LinkFile)


plasmids['length'] = plasmids['locus'].map(lenDict)
#plasmids['cstart'] = plasmids['locus'].map(lambda x: scoverage(x, links))
#plasmids['cend'] = plasmids['locus'].map(lambda x: ecoverage(x, links))
#plasmids['cprop'] = (plasmids['cend'] - plasmids['cstart'])/plasmids['length']
plasmids['cprop'] = plasmids['locus'].map(lambda x: coverage(x, links, plasmids[plasmids['locus'] == x]['length'].item()))
plasmids['links'] = plasmids['locus'].map(lambda x: numCov(x, links))
print(plasmids[plasmids['links'] > 10000])

print("Proportion median:{}\nLinks median:{}\n".format(plasmids['cprop'].median(), plasmids['links'].median()))

medianUnique = plotter(plasmids, links, args.out)

print("Unique linked plasmids median:{}".format(medianUnique))
