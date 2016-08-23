SORVA (Significance Of Rare VAriants)
=====================================
SORVA is a stand-alone program for querying the SORVA database and calculating the significance of your NGS sequencing findings. 

The SORVA dataset contains calculations on the number of individuals who have a rare variant in a given gene for numerous filtering threshold scenarios, which may be used for calculating the significance of an observed rare variant being causal for disease. Run SORVA to answer the question: How often do individuals carry a mutation in a given gene? You can specify the type of variant, e.g. only count an individual if they are homozygous for a loss-of-function variant, after filtering out common variants with a minor allele frequency >= 5%, by specifying options while running the script.

## Requirements

* Python. Python version 3.5.1 has been tested.
* SciPy. Version 0.13.3 has been tested.

## Reference

Aliz R Rao et al. Calculating the statistical significance of rare variants causal for Mendelian disorders

## Installation

After extracting the downloaded zip file, navigate into the directory and for options, run:

    python sorva.py --help
    python statistics.py --help


## Usage

### Quick start

Let's go through the following example:

What fraction of the general ('ALL') population has a heterozygous LOF variant anywhere in the gene DMD, after filtering out all variants with MAF > 0.05 ?

    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD

How about anywhere within described protein domains in DMD?

    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD --protdomains

We had sequenced 10 unrelated individuals and 3 of them had heterozygous variants in DMD. How likely is this to occur by random chance? The output of the first command becomes the -f parameter of the following command to calculate our p-value:

    python statistics.py -f 0.0023961661341853034 --n1 10 --s1 3
    
The resulting value is the nominal P-value, which we have to correct for multiple testing. A simple but conservative method is the Bonferroni correction, where we multiply the P-value by the number of genes we had sequenced:

    P-value = 1.63029317901e-06 * 24000 = 0.0391

### sorva.py

```Options:
  -h, --help            show this help message and exit
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
  --zygosity {het,hom,both,hom_or_compoundhet}, -z {het,hom,both,hom_or_compoundhet}
                        Zygosity of variants to include. het=heterozygous
                        variants only. hom=homozygous variants only. both=both
                        types of variants. hom_or_compoundhet=homozygous
                        variants or potential compound heterozygous variants
                        (two heterozygous variants at different loci within
                        the gene)
  --gene GENE, -g GENE  Gene name. Use standard gene name from HGNC or specify
                        Ensembl gene ID starting with ENSG.
  --protdomains         Output protein domain info.
  --genelist {ensembl75}
                        Gene list and version.
  --genomebuild {hg19}  Human genome build.
  --verbose, -v         Gene mutational burden output should be in the form of
                        a full sentence.
```

### statistics.py

```Options:
  -h, --help           show this help message and exit
  -f F                 background frequency / fraction of individuals in
                       population who have a mutation in the gene of interest
                       (output of sorva.py)
  --n1 N1              number of singletons (unrelated individuals) sequenced
  --s1 S1              number of singletons who have a mutation in the gene of
                       interest
  --n1_2 N1_2          number of families where individuals share 1/2 of their
                       genome (e.g. parent-child pair or full siblings)
  --s1_2 S1_2          number of families sharing 1/2 of their genome who
                       share mutations in the gene of interest
  --n1_4 N1_4          number of families where individuals share 1/4 of their
                       genome (e.g. grandparent-grandchild pair, aunt/uncle-
                       niece/nephew pair, half siblings, two first cousins or
                       three siblings)
  --s1_4 S1_4          number of families sharing 1/4 of their genome who
                       share mutations in the gene of interest
  --n1_8 N1_8          number of families where individuals share 1/8 of their
                       genome (e.g. second cousins)
  --s1_8 S1_8          number of families sharing 1/8 of their genome who
                       share mutations in the gene of interest
  --n1_16 N1_16        number of families where individuals share 1/16 of
                       their genome (e.g. third cousins)
  --s1_16 S1_16        number of families sharing 1/16 of their genome who
                       share mutations in the gene of interest
  --n1_32 N1_32        number of families where individuals share 1/32 of
                       their genome
  --s1_32 S1_32        number of families sharing 1/32 of their genome who
                       share mutations in the gene of interest
  --n_custom N_CUSTOM  number of families where individuals share a fraction
                       of their genome that is a custom value, specified by
                       --custom
  --s_custom S_CUSTOM  number of families sharing user-specified fraction of
                       their genome who share mutations in the gene of
                       interest
  --custom CUSTOM      fraction of genome shared by family members sequenced
  --verbose, -v        output additional information on how P-value was
                       calculated
```

## Files

```[distribution]
	data		[Directory contains 840 mutational burden data files]
	COPYING.txt	[Copy of GNU General Public License v3]
	LICENSE.md	[Copyright and license information]
	README.md	[This file]
	sorva.py	[Python script to query data]
	statistics.py	[Python script to calculate significance P-value]
```

