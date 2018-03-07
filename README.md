SORVA (Significance Of Rare VAriants)
=====================================
SORVA is a stand-alone program for querying the SORVA database and calculating the significance of your NGS sequencing findings. The web-based version can be found at https://sorva.genome.ucla.edu and is highly recommended for most applications. This stand-alone version is recommended for bulk queries, e.g. calculating statistics for every gene.

The SORVA dataset contains calculations on the number of individuals who have a rare variant in a given gene for numerous filtering threshold scenarios, which may be used for calculating the significance of an observed rare variant being causal for disease. Run SORVA to answer the question: How often do individuals carry a mutation in a given gene? You can specify the type of variant, e.g. only count an individual if they are homozygous for a loss-of-function variant, after filtering out common variants with a minor allele frequency >= 5%, by specifying options while running the script.

## Requirements

* Python. Python version 2.7.12 and 3.5.1 have been tested.
* SciPy. Version 0.13.3 has been tested.

## Reference

Aliz R Rao et al. Calculating the statistical significance of rare variants causal for Mendelian disorders
doi: https://doi.org/10.1101/103218

## Installation

After extracting the downloaded zip file, navigate into the directory and for options, run:

    python sorva.py --help
    python sorva-stats.py --help


## Usage

### Tutorial

Let's go through the following example:

We had sequenced 10 unrelated individuals with a rare, Mendelian disorder and 3 of them were heterozygous for rare, loss-of-function variants in DMD. How likely is this to occur by random chance?

First, we need to find out: what fraction of the general ('ALL') population has a heterozygous LOF variant anywhere in the gene DMD, after filtering out all variants with MAF > 0.05 ?

    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD

Just out of curiosity, how about anywhere within described protein domains in DMD?

    python sorva.py --genelist ensembl75 --genomebuild hg19 -c lof -p ALL -m 0.05 -b binary -z het --gene DMD --protdomains

To calculate our P-value, the output of the first command becomes the -f parameter of the following command:

    python sorva-stats.py -f 0.0023961661341853034 -n 10 -s 3 -c 1
    
The -c parameter = 1 indicates that all of the individuals were unrelated to each other. The resulting number is the nominal P-value, which we have to correct for multiple testing. A simple but conservative method is the Bonferroni correction, where we multiply the P-value by the number of genes we had sequenced:

    P-value = 1.63029317901e-06 * 24000 = 0.0391

Our results are statistically significant. This gene is highly suspicious of being causal for disease and we follow up with functional studies to confirm our results.

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
  -h, --help            show this help message and exit
  -f F                  background frequency / fraction of individuals in
                        population who have a mutation in the gene of interest
                        (output of sorva.py)
  --num_families NUM_FAMILIES, -n NUM_FAMILIES
                        The number of families sequenced, grouped by which
                        families have identical family structures. E.g. if we
                        sequenced 5 singletons and 3 sib pairs, this should be
                        [5, 3]
  --num_successes NUM_SUCCESSES, -s NUM_SUCCESSES
                        The number of families sequenced who have a variant in
                        the gene of interest (assumed to be shared IBD within
                        the family), in the order that corresponds to the sets
                        listed in the --num_families parameter. E.g. if we
                        found that 3 singletons and 2 sib pairs had variants
                        in a given gene, this would be [3, 2]
  --f_coefficients F_COEFFICIENTS, -c F_COEFFICIENTS
                        The coefficients for f. In case of an autosomal
                        dominant (AD) disorder, this will usually be the
                        coefficient of relationship. For autosomal recessive
                        disorders, this will be 1 for unrelated individuals
                        and 0.25 for sib pairs. For consanguineous pedigrees,
                        see manuscript for method to calculate this
                        coefficient. E.g. if we sequenced singletons and sib
                        pairs and we're looking for het variants (disorder is
                        AD), this would be [1, 0.25]
  --denovo DENOVO       Out of the singletons who have a variant in the gene
                        of interest, how many are de novo variants?
  --gene GENE, -g GENE  gene name; req'd if denovo > 0. Use standard gene name
                        from HGNC or specify Ensembl gene ID starting with
                        ENSG.
  --consequence {nonsyn,lof}
                        variant consequence filtering threshold; req'd if
                        denovo > 0. nonsyn = missense or more severe, also
                        includes LOF. lof = potential loss-of-function (LOF),
                        which includes stop loss gain, splice site mutation,
                        and frameshift indel.
  --verbose, -v         output additional information on how P-value was
                        calculated
```

## Files

```[distribution]
	data		[Directory contains 840 mutational burden data files]
	COPYING.txt	[Copy of GNU General Public License v3]
	LICENSE.md	[Copyright and license information]
	README.md	[This file]
	sorva.py	[Python script to query data]
	sorvastats.py	[Python script to calculate significance P-value]
```
## Feedback

If you have any questions or comments, please contact alizrrao (at) gmail . com.
