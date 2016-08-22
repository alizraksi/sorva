SORVA (Significance Of Rare VAriants)
=====================================
SORVA is a stand-alone program for querying the SORVA database and calculating the significance of your NGS sequencing findings. 

The SORVA dataset contains calculations on the number of individuals who have a rare variant in a given gene for numerous filtering threshold scenarios, which may be used for calculating the significance of an observed rare variant being causal for disease. Run SORVA to answer the question: How often do individuals carry a mutation in a given gene? You can specify the type of variant, e.g. only count an individual if they are homozygous for a loss-of-function variant, after filtering out common variants with a minor allele frequency >= 5%, by specifying options while running the script.

Requirements
============
* Python. Python version 3.5.1 has been tested.
* SciPy. Version 0.13.3 has been tested.

Reference
=========
Aliz R Rao et al. Calculating the statistical significance of rare variants causal for Mendelian disorders

Usage
=====
    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD

    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD --protdomains

    #The output of the first command becomes the -f parameter of the following command to calculate p-value:
    python statistics.py -f 0.00239616613419 --n1 10 --s1 3

Installation
============
    tar xvzf sorva-1.x.tar.gz;
    cd sorva-1.x;

    #for options run:
    python sorva.py --help
    python statistics.py --help


SORVA
=====
Options:
  --consequence {nonsyn,lof}, -c {nonsyn,lof}
                        Variant consequence filtering threshold. nonsyn =
                        missense or more severe, also includes LOF. lof =
                        potential loss-of-function (LOF), which includes stop
                        loss gain, splice site mutation, and frameshift indel.
  --population {ALL,EAS,EUR,AFR,AMR,SAS}, -p {ALL,EAS,EUR,AFR,AMR,SAS}
                        Superpopulation as specified by 1000 Genomes Project.
                        EAS=East Asian. EUR=European. AFR=African. AMR=ad-
                        mixed American. SAS=South Asian. ALL=all populations.
  --maf {0.05,0.01,0.005,0.001,0.0005}, -m {0.05,0.01,0.005,0.001,0.0005}
                        Minor allele frequency (MAF) threshold.
  --binarity {binary,countvariants}, -b {binary,countvariants}
                        Count the number of individuals (binary), or the total
                        number of variants seen in individuals
                        (countvariants). The latter option will count an
                        individual multiple times if they have multiple
                        variants in the gene.
  --zygosity {het,hom,both}, -z {het,hom,both}
                        Zygosity of variants to include. het=heterozygous
                        variants only. hom=homozygous variants only. both=both
                        types of variants.
  --gene GENE, -g GENE  Gene name. Use standard gene name from HGNC or specify
                        Ensembl gene ID starting with ENSG.
  --protdomains         Output protein domain info.
  --genelist {ensembl75}
                        Gene list and version.
  --genomebuild {hg19}  Human genome build.
  --verbose, -v         Gene mutational burden output should be in the form of
                        a full sentence.


STATISTICS
==========
Options:
  -p P           background proportion of individuals in population who have a
                 mutation in the gene of interest (output of sorva.py)
  --n1 N1        number of singletons (unrelated individuals) sequenced
  --s1 S1        number of singletons who have a mutation in the gene of
                 interest
  --n1_2 N1_2    number of families where individuals share 1/2 of their
                 genome (e.g. parent-child pair or siblings)
  --s1_2 S1_2    number of families sharing 1/2 of their genome who share
                 mutations in the gene of interest
  --n1_4 N1_4    number of families where individuals share 1/4 of their
                 genome (e.g. grandparent-grandchild pair, aunt/uncle-
                 niece/nephew pair, half siblings, two first cousins or three
                 siblings)
  --s1_4 S1_4    number of families sharing 1/4 of their genome who share
                 mutations in the gene of interest
  --n1_8 N1_8    number of families where individuals share 1/8 of their
                 genome
  --s1_8 S1_8    number of families sharing 1/8 of their genome who share
                 mutations in the gene of interest
  --n1_16 N1_16  number of families where individuals share 1/16 of their
                 genome
  --s1_16 S1_16  number of families sharing 1/16 of their genome who share
                 mutations in the gene of interest
  --n1_32 N1_32  number of families where individuals share 1/32 of their
                 genome
  --s1_32 S1_32  number of families sharing 1/32 of their genome who share
                 mutations in the gene of interest
  --verbose, -v  output additional information on how P-value was calculated


FILES
=====
[distribution]
	data		[Directory contains 840 mutational burden data files]
	LICENSE.md	[GNU General Public License v3]
	README.md		[This file]
	sorva.py	[Python script to query data]
	statistics.py	[Python script to calculate significance P-value]


