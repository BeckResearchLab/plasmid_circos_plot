# plasmid_circos_plot
Circos plot and summary charts describing possible recombination sites from NCBI all plasmid database

### You will need:
* circos
* biopython
* pandas


### To run
* ./make_organisms.py
* ./make_karyotype.py
* ./filtered2links.py
* ./make_chromolist
* edit circos.conf and update the chromosomes list
* select one of the karyotype_X.txt files and copy to karyotype.txt
* ./run
* python make_figure_1.py -o Figure1.png links.txt ncbi_plasmids.genbank
