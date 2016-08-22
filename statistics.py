#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from scipy.stats import binom

class StatisticsError(Exception): pass

#-f 0.01 --n1 2 --s1 2

def check_params(n, s):
    if s > n:
        sys.exit('The number of families with a variant in the gene cannot be greater than the number of families sequenced. Please check the input parameters and try again. (' + str(s) + '>' + str(n) + ')')

def run(f, n1, s1, n1_2, s1_2, n1_4, s1_4, n1_8, s1_8, n1_16, s1_16, n1_32, s1_32, n_custom, s_custom, custom, verbose):
    
    # check that correct number of parameters were passed
    check_params(n1, s1)
    check_params(n1_2, s1_2)
    check_params(n1_4, s1_4)
    check_params(n1_8, s1_8)
    check_params(n1_16, s1_16)
    check_params(n1_32, s1_32)
    check_params(n_custom, s_custom)
    
    #f= 0.00838658146965
    #print(str(8232 * 1.2 * math.pow(10,(-8))))
    #print(str(1 - binom.cdf(1 - 1, 1, 8232 * 1.2 * math.pow(10,(-8)))))
    #print(str(1 - binom.cdf(2 - 1, 5, 8232 * 1.2 * math.pow(10,(-8)))))
    #print(str(1 - binom.cdf(100 - 1, 100, f)))

    
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

    # Probability of having all 9 of the above events occurring together:
    pval = pb1 * pb1_2 * pb1_4 * pb1_8 * pb1_16 * pb1_32 * pb_custom

    print(pval)

    if verbose >= 1:
        print("Hello " + str(pval))

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
        help='number of singletons who have a mutation in the gene of interest')
    parser.add_argument('--n1_2', required=False, type=int, default=0,
        help='number of families where individuals share 1/2 of their genome (e.g. parent-child pair or full siblings)')
    parser.add_argument('--s1_2', required=False, type=int, default=0,
        help='number of families sharing 1/2 of their genome who share mutations in the gene of interest')
    parser.add_argument('--n1_4', required=False, type=int, default=0,
        help='number of families where individuals share 1/4 of their genome (e.g. grandparent-grandchild pair, aunt/uncle-niece/nephew pair, half siblings, two first cousins or three siblings)')
    parser.add_argument('--s1_4', required=False, type=int, default=0,
        help='number of families sharing 1/4 of their genome who share mutations in the gene of interest')
    parser.add_argument('--n1_8', required=False, type=int, default=0,
        help='number of families where individuals share 1/8 of their genome (e.g. second cousins)')
    parser.add_argument('--s1_8', required=False, type=int, default=0,
        help='number of families sharing 1/8 of their genome who share mutations in the gene of interest')
    parser.add_argument('--n1_16', required=False, type=int, default=0,
        help='number of families where individuals share 1/16 of their genome (e.g. third cousins)')
    parser.add_argument('--s1_16', required=False, type=int, default=0,
        help='number of families sharing 1/16 of their genome who share mutations in the gene of interest')
    parser.add_argument('--n1_32', required=False, type=int, default=0,
        help='number of families where individuals share 1/32 of their genome')
    parser.add_argument('--s1_32', required=False, type=int, default=0,
        help='number of families sharing 1/32 of their genome who share mutations in the gene of interest')
    parser.add_argument('--n_custom', required=False, type=int, default=0,
        help='number of families where individuals share a fraction of their genome that is a custom value, specified by --custom')
    parser.add_argument('--s_custom', required=False, type=int, default=0,
        help='number of families sharing user-specified fraction of their genome who share mutations in the gene of interest')
    parser.add_argument('--custom', required=False, type=float, default=0,
        help='fraction of genome shared by family members sequenced')
    parser.add_argument('--verbose', '-v', required=False, action='count',
        help='output additional information on how P-value was calculated')

    args = parser.parse_args()
    
    run(f=args.f, n1=args.n1, s1=args.s1, n1_2=args.n1_2, s1_2=args.s1_2, n1_4=args.n1_4, s1_4=args.s1_4, n1_8=args.n1_8, s1_8=args.s1_8, n1_16=args.n1_16, s1_16=args.s1_16, n1_32=args.n1_32, s1_32=args.s1_32, n_custom=args.n_custom, s_custom=args.s_custom, custom=args.custom, verbose=args.verbose)

    return 0

if __name__ == "__main__": sys.exit(main())
