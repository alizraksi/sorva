#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import csv
from os.path import isfile, join
import argparse

class SorvaError(Exception): pass

#--genelist ensembl75 --genomebuild hg19 --consequence nonsyn --population ALL --maf 0.05 --binarity binary --zygosity hom_or_compoundhet --gene DMD
#--genelist ensembl75 --genomebuild hg19 --consequence nonsyn --population ALL --maf 0.05 --binarity binary --zygosity hom_or_compoundhet --gene DMD --verbose
#--genelist ensembl75 --genomebuild hg19 --consequence nonsyn --population ALL --maf 0.05 --binarity binary --zygosity hom_or_compoundhet --gene DMD --protdomains

def file_exists(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        try:
            open(file_path)
            return True
        except IOError as e: return False
    else: return False

def run(consequence, population, maf, binarity, zygosity, gene, protdomains, genelist, genomebuild, verbose):
    
    analysisType = 'genes' if protdomains is None else 'proteindomains'
    mafstr = 'MAF' + str(maf)
    binarity = 'Binary' if binarity=='binary' else 'CountVariants'
    if zygosity == 'het':
        zygosity = 'HeteroOnly'
    elif zygosity == 'hom':
        zygosity = 'HomoOnly'
    elif zygosity == 'hom_or_compoundhet':
        if binarity == 'Binary':
            zygosity = 'HomOrCompoundHet'
        else:
            raise SorvaError('Cannot find hom_or_compoundhet variants in countvariants mode. Set binarity to binary.')
    elif zygosity == 'both':
        zygosity ='Both'
    else:
        raise SorvaError('Value of zygosity parameter is invalid: ' + zygosity)
    filename = '_'.join([genomebuild, genelist, consequence, population, analysisType, mafstr, binarity, zygosity, 'sum']) + '.txt'
    dir_path = os.path.dirname(os.path.realpath(__file__))
    filepath = dir_path + '/data/' + filename

    if not file_exists(filepath):
        raise SorvaError('Cannot find data file: ' + filepath + '\nPlease check that the data directory exists and contains 841 files.')
    
    with open(filepath, 'r') as in_file:
        reader = csv.reader(in_file, delimiter='\t')
    
        count = None
        if analysisType == 'genes':
            
            for row in reader:
                if ((row[0] == gene) or (row[3] == gene)):
                    count = row[2]
                    
            if count is None:
                sys.exit('Cannot find gene. Please check spelling, capitalization, and whether you are using the official HGNC gene symbol or Ensembl gene ID.')
                
            # We've found the count. Now convert the value to a fraction.    
            populations = {'ALL' : all,
                    'EAS' : eas,
                    'EUR' : eur,
                    'AFR' : afr,
                    'AMR' : amr,
                    'SAS' : sas,
            }
            total = populations[population]()
            proportion = float(count) / total
            
            # Generate output    
            if verbose is None:
                print(proportion)
            else:   
                percentage = 100.0 * int(count) / int(total)
                consqstring = "nonsynonymous or worse" if consequence == 'nonsyn' else "loss-of-function (LOF)"
                if zygosity == 'HeteroOnly':
                    zygositystr = "are heterozygous for a"
                elif zygosity == 'HomoOnly':
                    zygositystr = "are homozygous for a"
                else:
                    zygositystr = "have a"
                print(str(count) + ' out of ' + str(total) + ' individuals (' + str(percentage) + ' %) in the ' + population + ' superpopulation ' + zygositystr + ' rare, ' + consqstring + ' variant in the gene ' + gene + '.')
                
        else: # analysisType == 'proteindomains':

            # Check that gene name is valid
            filename = '_'.join([genomebuild, genelist, consequence, population, 'genes', mafstr, binarity, zygosity, 'sum']) + '.txt'
            dir_path = os.path.dirname(os.path.realpath(__file__))
            filepath = dir_path + '/data/' + filename
            if not file_exists(filepath):
                raise SorvaError('Cannot find data file: ' + filepath + '\nPlease check that the data directory exists and contains 841 files.')
            count = None
            with open(filepath, 'r') as in_file:
                gene_reader = csv.reader(in_file, delimiter='\t')
                for row in gene_reader:
                    if ((row[0] == gene) or (row[3] == gene)):
                        count = row[2]
                if count is None:
                    sys.exit('Cannot find gene. Please check spelling, capitalization, and whether you are using the official HGNC gene symbol or Ensembl gene ID.')
                
            # Find protein domains
            count = 0
            for row in reader:
                if ((row[0] == gene) or (row[11] == gene)):
                    count = count+1
                    if count == 1:
                        print('gene_name\tcount\tchrom\tpos\tinterpro_id\tinterpro_descr\tdb\tdb_acc\taa_start\taa_stop\tnum_exons_spanned\tgene_id')
                    print('\t'.join(row))
                    
            if count == 0:
                print('No proteins domains in the gene have any mutations.')
             

def all():
    return 2504
def eas():
    return 504
def eur():
    return 503
def afr():
    return 661
def amr():
    return 347
def sas():
    return 489


def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = '''SORVA (Significance Of Rare VAriants)
        
Run SORVA to answer the question: How often do individuals carry a mutation in a given gene?
You can specify the type of variant, e.g. only count an individual if they are homozygous for a loss-of-function variant, after filtering out common variants with a minor allele frequency >= 5%, by specifying the options below.

IMPORTANT NOTE: If intend to use the output to calculate statistics with, using statistics.py, use the following rules: If sequenced individuals have homozygous mutations in the gene, use the option --zygosity hom. If sequenced individuals have heterozygous mutations in the gene, use --zygosity both.''',
        epilog = '''SORVA version 1.0 (c)2015-2016 Aliz R. Rao. All rights reserved.''',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--consequence', '-c', required=True, choices=['nonsyn', 'lof'],
        help = 'Variant consequence filtering threshold. nonsyn = missense or more severe, also includes LOF. lof = potential loss-of-function (LOF), which includes stop loss gain, splice site mutation, and frameshift indel.')
    parser.add_argument('--population', '-p', required=True, choices=['ALL', 'EAS', 'EUR', 'AFR', 'AMR', 'SAS'],
        help='Superpopulation as specified by 1000 Genomes Project. EAS=East Asian. EUR=European. AFR=African. AMR=ad-mixed American. SAS=South Asian. ALL=all populations.')
    parser.add_argument('--maf', '-m', required=True, type=float, choices=[0.05, 0.01, 0.005, 0.001, 0.0005],
        help='Minor allele frequency (MAF) threshold.')
    parser.add_argument('--binarity', '-b', required=True, choices=['binary', 'countvariants'],
        help='Count the number of individuals (binary), or the total number of variants seen in individuals (countvariants). The latter option will count an individual multiple times if they have multiple variants in the gene.')
    parser.add_argument('--zygosity', '-z', required=True, choices=['het', 'hom', 'both', 'hom_or_compoundhet'],
        help='Zygosity of variants to include. het=heterozygous variants only. hom=homozygous variants only. both=both types of variants. hom_or_compoundhet=homozygous variants or potential compound heterozygous variants (two heterozygous variants at different loci within the gene)')
    parser.add_argument('--gene', '-g', required=True,
        help='Gene name. Use standard gene name from HGNC or specify Ensembl gene ID starting with ENSG.')
    parser.add_argument('--protdomains', required=False, action='count',
        help='Output protein domain info.')
    parser.add_argument('--genelist', required=False, default='ensembl75', choices=['ensembl75'],
        help='Gene list and version.')
    parser.add_argument('--genomebuild', required=False, default='hg19', choices=['hg19'],
        help='Human genome build.')
    parser.add_argument('--verbose', '-v', required=False, action='count',
        help='Gene mutational burden output should be in the form of a full sentence.')

    args = parser.parse_args()
    
    run(consequence=args.consequence, population=args.population, maf=args.maf, binarity=args.binarity, zygosity=args.zygosity, gene=args.gene, protdomains=args.protdomains, genelist=args.genelist, genomebuild=args.genomebuild, verbose=args.verbose)

    return 0

if __name__ == "__main__": sys.exit(main())
