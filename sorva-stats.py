#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import csv
import argparse
from scipy.stats import binom
import math

class StatisticsError(Exception): pass

#-f 0.01 --n1 2 --s1 2

def check_params(n, s):
    if s > n:
        sys.exit('The number of families with a variant in the gene cannot be greater than the number of families sequenced. Please check the input parameters and try again. (' + str(s) + '>' + str(n) + ')')

def check_params_denovo(s, d):
    if d > s:
        sys.exit('The number of individuals with a de novo variant in the gene cannot be greater than the number of individuals with a variant in the gene. Please check the input parameters and try again. (' + str(d) + '>' + str(s) + ')')

def file_exists(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        try:
            open(file_path)
            return True
        except IOError as e: return False
    else: return False

def run(f, n1, s1, denovo, gene, consequence, n1_2, s1_2, n1_4, s1_4, n1_8, s1_8, n1_16, s1_16, n1_32, s1_32, n_custom, s_custom, custom, verbose):
    
    # check that correct number of parameters were passed
    check_params(n1, s1)
    check_params(n1_2, s1_2)
    check_params(n1_4, s1_4)
    check_params(n1_8, s1_8)
    check_params(n1_16, s1_16)
    check_params(n1_32, s1_32)
    check_params(n_custom, s_custom)
    check_params_denovo(s1, denovo)

    #f= 0.00838658146965
    #print(str(8232 * 1.2 * math.pow(10,(-8))))
    #print(str(1 - binom.cdf(1 - 1, 1, 8232 * 1.2 * math.pow(10,(-8)))))
    #print(str(1 - binom.cdf(2 - 1, 5, 8232 * 1.2 * math.pow(10,(-8)))))
    #print(str(1 - binom.cdf(100 - 1, 100, f)))

    # do we have de novo variants?

    pval_dn = 1
    if denovo > 0:
        if gene == '':
            raise StatisticsError('Gene symbol is a required input when calculating with de novo variants. Please specify the --gene parameter.')
        if consequence == '':
            raise StatisticsError('Variant consequence is a required input when calculating with de novo variants. Please specify the --consequence parameter.')
        # get transcript length
        dir_path = os.path.dirname(os.path.realpath(__file__))
        filepath = dir_path + '/data/refseq_tx_lengths.csv'
        if not file_exists(filepath):
            raise StatisticsError(
                'Cannot find data file: ' + filepath + '\nPlease check that the data directory exists and contains 841 files.')
        count = None
        with open(filepath, 'r') as in_file:
            reader = csv.reader(in_file, delimiter=',')
            l = None
            for row in reader:
                if ((row[1] == gene) or (row[2] == gene)):
                    l = int(row[4])
            if l is None:
                sys.exit(
                    'Cannot find gene ' + gene + '. Please check spelling, capitalization, and whether you are using the official HGNC gene symbol or Ensembl gene ID.')

            d = 1.2 * math.pow(10, (-8))
            if consequence == 'nonsyn':
                c = 0.7064
            else:
                c = 0.0285  # LOF variants
            pval_dn = l * d * c
            if n1 > 1:
                pval_dn = 1 - binom.cdf(denovo - 1, n1, (l * d * c))

            # Do we have additional indivs w var who aren't de novo? Later count p-val for seeing these and later multiply the p-vals together
            s1 = s1 - denovo
            n1 = n1 - denovo
    
    # calculate p-values considering family structures of individuals sequenced
    
    if f==0:
        # Equation would always return significant with f==0. Therefore, set f to a very small number. Assume that if we
        # had sequenced twice the number of samples, we would have seen 1 individual with a variant in the gene.
        f = 1.0 / (2 * 2504)

    # The probability of having at least k out of n exomes having at least 1 var
    pb1 = 1 - binom.cdf(s1 - 1, n1, f)          # Interpretation: the prob of seeing 0 or more singletons with at least 1 variant is 1
    pb1_2 = 1 - binom.cdf(s1_2 - 1, n1_2, f / 2)  #          the prob of seeing s or more families sharing 1/2 of genome with at least 1 var shared is ...
    pb1_4 = 1 - binom.cdf(s1_4 - 1, n1_4, f / 4)
    pb1_8 = 1 - binom.cdf(s1_8 - 1, n1_8, f / 8)
    pb1_16 = 1 - binom.cdf(s1_16 - 1, n1_16, f / 16)
    pb1_32 = 1 - binom.cdf(s1_32 - 1, n1_32, f / 32)
    pb_custom = 1 - binom.cdf(s_custom - 1, n_custom, f * custom)

    # Probability of having all of the above events occurring together:
    pval = pval_dn * pb1 * pb1_2 * pb1_4 * pb1_8 * pb1_16 * pb1_32 * pb_custom

    if verbose == None:
        print(pval)
    elif verbose >= 1:
        # Set n1 and s1 back to their original values, in case they were changed
        s1 = s1 + denovo
        n1 = n1 + denovo

        signifstr = ''
        if pval < 0.05:
            signifstr = 'reaches'
        else:
            signifstr = 'does not reach'
        pvalstr = ', which ' + signifstr + ' genome-wide significance.'

        #TODO: incorporate information about follow-up sequencing
        #if ($f_n1 > 1):
        #$pval_followup_str = "Afterwards, the results of the follow-up sequencing were also taken into account, and this reduced the P-value from $pval_before_followup to $pval.";

        percentage = 100.0 * f

        singletonstr = ''
        halfsharedstr = ''
        fourthsharedstr = ''
        eightsharedstr = ''
        sixteenthsharedstr = ''
        thirtytwothsharedstr = ''
        customsharedstr = ''
        andstr = ''

        if (n1 > 0):
            singletonstr = str(s1) + ' out of ' + str(n1) + ' unrelated individuals'
            if denovo > 0:
                singletonstr = singletonstr + ' (' + str(denovo) + ' of which have de novo variants)'

        if (n1_2 > 0):
            if (n1 > 0):
                andstr = ' and '
            else:
                andstr = ''
            halfsharedstr = andstr + str(s1_2) + ' out of ' + str(n1_2) + ' families where individuals share half of their genome'

        if (n1_4 > 0):
            if ((n1 > 0) or (n1_2 > 0)):
                andstr = ' and '
            else:
                andstr = ''
            fourthsharedstr = andstr + str(s1_4) + ' out of ' + str(n1_4) + ' families where individuals share 1/4 of their genome'

        if (n1_8 > 0):
            if ((n1 > 0) or (n1_2 > 0) or (n1_4 > 0)):
                andstr = ' and '
            else:
                andstr = ''
            eightsharedstr = andstr + str(s1_8) + ' out of ' + str(n1_8) + ' families where individuals share 1/8 of their genome'

        if (n1_16 > 0):
            if ((n1 > 0) or (n1_2 > 0) or (n1_4 > 0) or (n1_8 > 0)):
                andstr = ' and '
            else:
                andstr = ''
            sixteenthsharedstr = andstr + str(s1_16) + ' out of ' + str(n1_16) + ' families where individuals share 1/16 of their genome'

        if (n1_32 > 0):
            if ((n1 > 0) or (n1_2 > 0) or (n1_4 > 0) or (n1_8 > 0) or (n1_16 > 0)):
                andstr = ' and '
            else:
                andstr = ''
                thirtytwothsharedstr = andstr + str(s1_32) + ' out of ' + str(n1_32) + ' families where individuals share 1/32 of their genome'

        if (n_custom > 0):
            if ((n1 > 0) or (n1_2 > 0) or (n1_4 > 0) or (n1_8 > 0) or (n1_16 > 0) or (n1_32 > 0)):
                andstr = ' and '
            else:
                andstr = ''
            customsharedstr = andstr + str(s_custom) + ' out of ' + str(n_custom) + ' families where individuals share ' + str(custom) + ' of their genome'

        print('Given that ' + "{0:.2f}".format(percentage) + ' % (f = ' + str(f) + ') of individuals in the population have a rare variant in the given gene, the statistical significance of seeing such variants in ' + singletonstr + halfsharedstr + fourthsharedstr + eightsharedstr + sixteenthsharedstr + thirtytwothsharedstr + customsharedstr + ' is P = ' + str(pval))

        if (verbose >= 2):
            print('\nNote: this P-value must be corrected for the number of genes targeted, e.g. using the Bonferroni adjustment.')

def main():
   
    #command line arguments
    parser = argparse.ArgumentParser(
        description = '''SORVA (Significance Of Rare VAriants) - Statistics package
        
Run SORVA to answer the question: How often do individuals carry a mutation in a given gene?
Next, calculate whether seeing such variants in a certain number of sequenced individuals is significant, given the background rate, and considering the family structures of the individuals sequenced.

The input parameter -p is the output of sorva.py. When running sorva.py, if sequenced individuals have homozygous mutations in the gene, use the option --zygosity hom. If sequenced individuals have heterozygous mutations in the gene, use --zygosity both.
 
IMPORTANT NOTE: Results must be corrected for multiple testing, by the number of genes sequenced (~20,000).''',
        epilog = '''SORVA version 1.0 (c)2015-2016 Aliz R. Rao. All rights reserved.''',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', required=True, type=float,
        help = 'background frequency / fraction of individuals in population who have a mutation in the gene of interest (output of sorva.py)')
    parser.add_argument('--n1', required=False, type=int, default=0,
        help='number of singletons (unrelated individuals) sequenced')
    parser.add_argument('--s1', required=False, type=int, default=0,
        help='number of singletons who have a variant in the gene of interest')
    parser.add_argument('--denovo', required=False, type=int, default=0,
        help='out of the singletons who have a variant in the gene of interest, how many are de novo variants?')
    parser.add_argument('--gene', '-g', required=False, default='',
        help='gene name; req\'d if denovo > 0. Use standard gene name from HGNC or specify Ensembl gene ID starting with ENSG.')
    parser.add_argument('--consequence', '-c', required=False, choices=['nonsyn', 'lof'], default='',
        help='variant consequence filtering threshold; req\'d if denovo > 0. nonsyn = missense or more severe, also includes LOF. lof = potential loss-of-function (LOF), which includes stop loss gain, splice site mutation, and frameshift indel.')
    parser.add_argument('--n1_2', required=False, type=int, default=0,
        help='number of families where individuals share 1/2 of their genome (e.g. parent-child pair or full siblings)')
    parser.add_argument('--s1_2', required=False, type=int, default=0,
        help='number of families sharing 1/2 of their genome who share variants in the gene of interest')
    parser.add_argument('--n1_4', required=False, type=int, default=0,
        help='number of families where individuals share 1/4 of their genome (e.g. grandparent-grandchild pair, aunt/uncle-niece/nephew pair, half siblings, two first cousins or three siblings)')
    parser.add_argument('--s1_4', required=False, type=int, default=0,
        help='number of families sharing 1/4 of their genome who share variants in the gene of interest')
    parser.add_argument('--n1_8', required=False, type=int, default=0,
        help='number of families where individuals share 1/8 of their genome (e.g. second cousins)')
    parser.add_argument('--s1_8', required=False, type=int, default=0,
        help='number of families sharing 1/8 of their genome who share variants in the gene of interest')
    parser.add_argument('--n1_16', required=False, type=int, default=0,
        help='number of families where individuals share 1/16 of their genome (e.g. third cousins)')
    parser.add_argument('--s1_16', required=False, type=int, default=0,
        help='number of families sharing 1/16 of their genome who share variants in the gene of interest')
    parser.add_argument('--n1_32', required=False, type=int, default=0,
        help='number of families where individuals share 1/32 of their genome')
    parser.add_argument('--s1_32', required=False, type=int, default=0,
        help='number of families sharing 1/32 of their genome who share variants in the gene of interest')
    parser.add_argument('--n_custom', required=False, type=int, default=0,
        help='number of families where individuals share a fraction of their genome that is a custom value, specified by --custom')
    parser.add_argument('--s_custom', required=False, type=int, default=0,
        help='number of families sharing user-specified fraction of their genome who share variants in the gene of interest')
    parser.add_argument('--custom', required=False, type=float, default=0,
        help='fraction of genome shared by family members sequenced')
    parser.add_argument('--verbose', '-v', required=False, action='count',
        help='output additional information on how P-value was calculated')

    args = parser.parse_args()
    
    run(f=args.f, n1=args.n1, s1=args.s1, denovo=args.denovo, gene=args.gene, consequence=args.consequence, n1_2=args.n1_2, s1_2=args.s1_2, n1_4=args.n1_4, s1_4=args.s1_4, n1_8=args.n1_8, s1_8=args.s1_8, n1_16=args.n1_16, s1_16=args.s1_16, n1_32=args.n1_32, s1_32=args.s1_32, n_custom=args.n_custom, s_custom=args.s_custom, custom=args.custom, verbose=args.verbose)

    return 0

if __name__ == "__main__": sys.exit(main())
