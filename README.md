# RinsilicoPCR - Test Primers and Extract Amplicons using BLAST and Biostrings (R package)



## Installation
###Requirements:
* Install NCBI [BLAST+](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/) suite
* Install [Bioconductor](http://www.bioconductor.org/install/) and the Bioconductor package `Biostrings`
* Install packages `plyr` and `data.table`

###Installation from github through `devtools`
```
library(devtools)
install_github("zachth/RinsilicoPCR")
```


## Example
```
##Run local primer blast

#load query sequence
query <- readDNAStringSet("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/010/525/GCA_000010525.1_ASM1052v1/GCA_000010525.1_ASM1052v1_genomic.fna.gz")

#define primers
primers <- c("GGCTGGATCACCTCCTT","GCCWAGGCATCCDCC")

#make blast database from primers
db <- makeprimerdb(primers)

#Define blast arguments: these are the same that can be found in `RinsilicoPCR::blast.params$primers`
args <- list(
             outfmt="'6 qacc sacc pident length qlen slen qstart sstart qend send evalue'",
             max_target_seqs=20,
             word_size=8
             )

#run blast
blastres <- blast(query,db,args)

#determine primer pairs
pairs <- primer.insilico(blastres)

#extract amplicons from query
amplicons <- extract(pairs)

>amplicons
  A DNAStringSet instance of length 3
    width seq                                                                                                                               names
[1]   829 TCTAAGGATGATCCTTCAGTCTTCGGGCCTTTCGGGCTTCGAGCTATCGGATCTCTTGGAAAC...TGAGAACCTCAGCCGAGGAGTGGGCATGGACGATGAGAACGATCAAGTGTCTTAAGGGCATTC AP009384.1.681988...
[2]   829 TCTAAGGATGATCCTTCAGTCTTCGGGCCTTTCGGGCTTCGAGCTATCGGATCTCTTGGAAAC...TGAGAACCTCAGCCGAGGAGTGGGCATGGACGATGAGAACGATCAAGTGTCTTAAGGGCATTC AP009384.1.477848...
[3]   829 TCTAAGGATGATCCTTCAGTCTTCGGGCCTTTCGGGCTTCGAGCTATCGGATCTCTTGGAAAC...TGAGAACCTCAGCCGAGGAGTGGGCATGGACGATGAGAACGATCAAGTGTCTTAAGGGCATTC AP009384.1.504494...
```

