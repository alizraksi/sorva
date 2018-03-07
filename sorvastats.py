#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import csv
import argparse
from scipy.stats import binom
import math
import scipy


class StatisticsError(Exception): pass

## The following are example inputs used to calculate P-values in Table 2 of Rao & Nelson et al 2018:
# -f 0 -n [5,1] -s [1,1] -c [1,0.25]
# -f 0.00439297124601 -n [4,1] -s [4,1] -c [1,0.0078125]
# -f 0.00199680511182 -n 10 -s 7 -c 1
# -f 0.00119808306709 -n [5,1,1] -s [2,0,1] -c [1,0.25,0.0625]


def parse_list_args(s):
    s = str(s).strip('[]').replace(' ', '')
    if '.' in s:
        l = [float(item) for item in s.split(',')]
    else:
        l = [int(item) for item in s.split(',')]
    return l


def check_params(num_families, num_successes, f_coefficients, denovo):
    for n, s, c in zip(num_families, num_successes, f_coefficients):
        if s > n:
            sys.exit(
                'The number of families with a variant in the gene cannot be greater than the number of families sequenced. Please check the input parameters and try again. (' + str(
                    s) + '>' + str(n) + ')')
        if c < 0 or c > 1:
            sys.exit(
                'The coefficient of f is out of bounds, must be between 0 and 1. Please check the input parameters and try again. (f_coefficient = '
                    + str(c) + ')')
        # if we have singletons, check if any of them are denovo and if the number of denovo indivs is within bounds
        if c == 1 and denovo is not None and denovo > s:
            sys.exit(
                'The number of individuals with a de novo variant in the gene cannot be greater than the number of individuals with a variant in the gene. Please check the input parameters and try again. (' + str(
                    denovo) + '>' + str(s) + ')')


def check_params_zyg_inh(z, i):
    if not z and not i:
        sys.exit(
            'Either --zygosity or --inheritance parameter is required.')
    if (z == 'het' and i == 'AR') or (z == 'hom' and i == 'AD'):
        sys.exit(
            'Conflict between --zygosity or --inheritance parameters. For autosomal recessive disorders, you should be looking for homozygous or compound het variants. For autosomal dominant disorders, you should be looking for heterozygous variants.')


def check_params_custom(n, s, r, e, z, i):
    if s > n:
        sys.exit(
            'The number of families with a variant in the gene cannot be greater than the number of families sequenced. Please check the input parameters and try again. (' + str(
                s) + '>' + str(n) + ')')
    if n > 0:
        if (z == 'het' or i == 'AD'):
            if not r and not e:
                sys.exit(
                    'Either --r_custom or --e_custom parameter is required.')
        if (z == 'hom' or i == 'AR'):
            if not e:
                sys.exit(
                    '--e_custom parameter is required.')


def file_exists(file_path):
    if isinstance(file_path, str) and os.path.isfile(file_path):
        try:
            open(file_path)
            return True
        except IOError as e:
            return False
    else:
        return False


def calc_P_alt(f, n0, n1, s0, s1, r0, r1):
# Implements calculating exact probability for when we sequence families with 2 types of pedigree structures.
# This will return a number equal to that returned by calc_P and can be useful for debugging and sanity checks.

    P0_X = n0 + n1 - s0 - s1
    P0_probsuccess = 1 - (1.0 * n0 / (n0 + n1)) * r0 * f - (1.0 * n1 / (n0 + n1)) * r1 * f
    P0_ntrials = n0 + n1

    P1_X = s0
    P1_probsuccess = 1.0 * (r0 * n0) / ((r0 * n0) + (r1 * n1))
    P1_ntrials = s0 + s1

    ## Gives same result:
    # P1_X = s1
    # P1_probsuccess = 1 - 1.0 * (r0 * n0) / ((r0 * n0) + (r1 * n1))
    # P1_ntrials = s0 + s1

    P = scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess) * scipy.stats.binom.pmf(P1_X, P1_ntrials, P1_probsuccess)
    #print 'P1 = ' + str(scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess))
    #print 'P2 = ' + str(scipy.stats.binom.pmf(P1_X, P1_ntrials, P1_probsuccess))

    return P


def calc_P_3(f, num_families, num_successes, f_coefficients):
# Implements calculating exact probability for when we sequence families with 3 types of pedigree structures.

    n0 = num_families[0]
    n1 = num_families[1]
    n2 = num_families[2]
    s0 = num_successes[0]
    s1 = num_successes[1]
    s2 = num_successes[2]
    r0 = f_coefficients[0]
    r1 = f_coefficients[1]
    r2 = f_coefficients[2]

    P0_X = n0 + n1 + n2 - s0 - s1 - s2
    P0_probsuccess = 1 - (1.0 * n0 / (n0 + n1 + n2)) * r0 * f - (1.0 * n1 / (n0 + n1 + n2)) * r1 * f - (1.0 * n2 / (n0 + n1 + n2)) * r2 * f
    P0_ntrials = n0 + n1 + n2

    P1_X = s0
    P1_probsuccess = 1.0 * (r0 * n0) / ((r0 * n0) + (r1 * n1) + (r2 * n2))
    P1_ntrials = s0 + s1 + s2

    P2_X = s1
    P2_probsuccess = 1.0 * (r1 * n1) / ((r1 * n1) + (r2 * n2))
    P2_ntrials = s1 + s2

    P = scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess) * scipy.stats.binom.pmf(P1_X, P1_ntrials, P1_probsuccess) * scipy.stats.binom.pmf(P2_X, P2_ntrials, P2_probsuccess)

    return P


def calc_P(f, num_families, num_successes, f_coefficients):
# Calculating exact probability for when we sequence families with any number of pedigree structures.
# If we only have a single type of family structure, it will return P = scipy.stats.binom.pmf(n0 - s0, n0, 1 - coeff0 * f)
# If we have two types of family structures, it will return the same value as calc_P_alt(f, n0, n1, s0, s1, coeff0, coeff1)

    # We will be destroying these lists, keep a copy so the original remains unchanged
    n_families = list(num_families)
    n_successes = list(num_successes)
    f_coeffs = list(f_coefficients)

    P = 1
    n_total_orig = sum(n_families)
    n_total_prev = n_total_orig
    n_failures_prev = 0
    remainingP = 1

    while len(n_families) >= 1:

        n_total = sum(n_families)
        s_total = sum(n_successes)

        P_ntrials = n_total_prev - n_failures_prev
        P_X = P_ntrials - s_total
        P_probsuccess = 1
        for n, s, c in zip(n_families, n_successes, f_coeffs):
            P_probsuccess = P_probsuccess - ((1.0 * n / (n_total_orig) * c * f) / remainingP)
        #print 'Px_X, Px_ntrials, Px_probsuccess = ' + str([P_X, P_ntrials, P_probsuccess])
        #print 'Px = ' + str(scipy.stats.binom.pmf(P_X, P_ntrials, P_probsuccess))
        P = P * scipy.stats.binom.pmf(P_X, P_ntrials, P_probsuccess)

        n_failures_prev = P_X
        n_total_prev = n_total
        remainingP = 0
        for n, s, c in zip(n_families, n_successes, f_coeffs):
            remainingP = remainingP + (1.0 * n / (n_total_orig)) * c * f
        del n_families[0]
        del n_successes[0]
        del f_coeffs[0]

    return P


def calc_P_denovo(f, l, n_sequenced, n_lof_dn, n_missense_dn, n_lof_inh):
    # Implements calculating exact probability for when we sequence families with 3 types of pedigree structures.

    n = n_sequenced
    s0 = n_lof_dn
    s1 = n_missense_dn
    s2 = n_lof_inh

    d = 1.2 * math.pow(10, (-8))
    c_missense = 0.7064
    c_lof = 0.0285

    P0_X = n - s0 - s1 - s2
    P0_probsuccess = 1 - (c_missense * l * d) - (c_lof * l * d) - (f - (c_lof * l * d))
    P0_ntrials = n

    P1_X = s0
    P1_probsuccess = (c_lof * l * d) / ((c_lof * l * d) + (c_missense * l * d) + (f - (c_lof * l * d)))
    P1_ntrials = s0 + s1 + s2

    P2_X = s1
    P2_probsuccess = 1.0 * (c_missense * l * d) / ((c_missense * l * d) + (f - (c_lof * l * d)))
    P2_ntrials = s1 + s2

    P = scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess) * scipy.stats.binom.pmf(P1_X, P1_ntrials,
                                                                                        P1_probsuccess) * scipy.stats.binom.pmf(
        P2_X, P2_ntrials, P2_probsuccess)

    return P


def calc_P_denovo_lofdn_missdn(f, l, n_sequenced, n_lof_dn, n_missense_dn):
    # Implements calculating exact probability for when we sequence families with 3 types of pedigree structures.

    n = n_sequenced
    s0 = n_lof_dn
    s1 = n_missense_dn

    d = 1.2 * math.pow(10, (-8))
    c_missense = 0.7064
    c_lof = 0.0285

    P0_X = n - s0 - s1
    P0_probsuccess = 1 - (c_missense * l * d) - (c_lof * l * d)
    P0_ntrials = n

    P1_X = s0
    P1_probsuccess = (c_lof * l * d) / ((c_lof * l * d) + (c_missense * l * d))
    P1_ntrials = s0 + s1

    P = scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess) * scipy.stats.binom.pmf(P1_X, P1_ntrials,
                                                                                        P1_probsuccess)
    #print 'P0 = ' + str(scipy.stats.binom.pmf(P0_X, P0_ntrials, P0_probsuccess))
    #print 'P1 = ' + str(scipy.stats.binom.pmf(P1_X, P1_ntrials, P1_probsuccess))

    return P


def run(f, num_families, num_successes, f_coefficients, denovo=None, gene=None, consequence=None, verbose=None):

    num_families = parse_list_args(num_families)
    num_successes = parse_list_args(num_successes)
    f_coefficients = parse_list_args(f_coefficients)

    # check that correct number of parameters were passed
    check_params(num_families, num_successes, f_coefficients, denovo)

    # if we only have a single set of families, convert parameters into a list with a single element
    if type(num_families) is not list: num_families = [num_families]
    if type(num_successes) is not list: num_successes = [num_successes]
    if type(f_coefficients) is not list: f_coefficients = [f_coefficients]

    # do we have de novo variants?

    pval_dn = 1
    if denovo > 0:
        if gene == '':
            raise StatisticsError(
                'Gene symbol is a required input when calculating with de novo variants. Please specify the --gene parameter.')
        if consequence == '':
            raise StatisticsError(
                'Variant consequence is a required input when calculating with de novo variants. Please specify the --consequence parameter.')
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

            # Find the number of singletons sequenced, split this into denovo fams and non-denovo fams
            dn_ix = n1 = s1 = None
            ix = 0
            for n, s, f_coeff in zip(num_families, num_successes, f_coefficients):
                if f_coeff == 1:
                    n1 = n
                    s1 = s
                    dn_ix = ix
                ix = ix + 1
            # Do we have additional indivs w var who aren't de novo?
            num_families[dn_ix] = num_families[dn_ix] - denovo
            num_successes[dn_ix] = num_successes[dn_ix] - denovo  #TODO: Calculate how to deal with some denovo fams, some not

            d = 1.2 * math.pow(10, (-8))
            if consequence == 'nonsyn':
                c = 0.7064
            else:
                c = 0.0285  # LOF variants
            pval_dn = l * d * c
            if n1 > 1:
                pval_dn = 1 - binom.cdf(denovo - 1, n1, (l * d * c))



    # calculate p-values considering family structures of individuals sequenced

    if f == 0:
        # Equation would always return significant with f==0. Therefore, set f to a very small number. Assume that if we
        # had sequenced twice the number of samples, we would have seen 1 individual with a variant in the gene. In this case,
        # N = the number of individuals in the 1000 Genomes Project dataset that the SORVA dataset was derived from.
        f = 1.0 / (2 * 2504)

    # If we have a single type of pedigree structure:
    if len(num_successes) == 1:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            P_other = calc_P(f, num_families, [s0_other], f_coefficients)
            if P_other <= P_exact:
                if verbose >= 2:
                    print 'Observing ' + str(s0_other) + ' successes is at least as significant (P = ' + str(P_other) + ')'
                P_combined = P_combined + P_other

    # If we have two types of pedigree structures:
    elif len(num_successes) == 2:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            for s1_other in range(0, num_families[1] + 1):
                P_other = calc_P(f, num_families, [s0_other, s1_other], f_coefficients)
                if P_other <= P_exact:
                    if verbose >= 2:
                        print 'Observing ' + str(s0_other) + ' and ' + str(
                        s1_other) + ' successes is at least as significant (P = ' + str(P_other) + ')'
                    P_combined = P_combined + P_other

    # If we have three types of pedigree structures:
    elif len(num_successes) == 3:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P_3(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            for s1_other in range(0, num_families[1] + 1):
                for s2_other in range(0, num_families[2] + 1):
                    P_other = calc_P_3(f, num_families, [s0_other, s1_other, s2_other], f_coefficients)
                    if P_other <= P_exact:
                        if verbose >= 2:
                            print 'Observing ' + str(s0_other) + ' and ' + str(
                                s1_other) + ' and ' + str(s2_other) + ' successes is at least as significant (P = ' + str(P_other) + ')'
                        P_combined = P_combined + P_other

    else:
        sys.exit('Calculating P-value when more than 3 types of pedigree structures sequenced has not yet been implemented.')
    #TODO: use recursion to make algorithm work for any number of types of family structures

    if verbose == None:
        print(P_combined)
    elif verbose >= 1:

        signifstr = ''
        if P_combined < 0.05:
            signifstr = 'reaches'
        else:
            signifstr = 'does not reach'
        pvalstr = ', which ' + signifstr + ' genome-wide significance.'

        # TODO: incorporate information about follow-up sequencing
        # if ($f_n1 > 1):
        # $pval_followup_str = "Afterwards, the results of the follow-up sequencing were also taken into account, and this reduced the P-value from $pval_before_followup to $pval.";

        percentage = 100.0 * f

        print('Given that ' + "{0:.2f}".format(percentage) + ' % (f = ' + str(f)
            + ') of individuals in the population have a rare variant in the given gene, the statistical significance of seeing such variants in the constellation of families that were sequenced is P = ' + str(
            P_combined) + pvalstr)

        print('\nNote: this P-value must be corrected for the number of genes targeted, e.g. using the Bonferroni adjustment.')


def run_denovo(f, num_families, num_successes, f_coefficients, denovo=None, gene=None, consequence=None, verbose=None, n_lof_dn=None, n_missense_dn=None, n_lof_inh=None):
    num_families = parse_list_args(num_families)
    num_successes = parse_list_args(num_successes)
    f_coefficients = parse_list_args(f_coefficients)

    # check that correct number of parameters were passed
    check_params(num_families, num_successes, f_coefficients, denovo)

    # if we only have a single set of families, convert parameters into a list with a single element
    if type(num_families) is not list: num_families = [num_families]
    if type(num_successes) is not list: num_successes = [num_successes]
    if type(f_coefficients) is not list: f_coefficients = [f_coefficients]

    if f == 0:
        # Equation would always return significant with f==0. Therefore, set f to a very small number. Assume that if we
        # had sequenced twice the number of samples, we would have seen 1 individual with a variant in the gene. In this case,
        # N = the number of individuals in the 1000 Genomes Project dataset that the SORVA dataset was derived from.
        f = 1.0 / (2 * 2504)

    # Calculate p-values considering family structures of individuals sequenced

    # Do we have de novo variants?
    if n_lof_dn > 0 or n_missense_dn > 0:
        if gene == '':
            raise StatisticsError(
                'Gene symbol is a required input when calculating with de novo variants. Please specify the --gene parameter.')
        if consequence == '':
            raise StatisticsError(
                'Variant consequence is a required input when calculating with de novo variants. Please specify the --consequence parameter.')
        # get transcript length
        dir_path = os.path.dirname(os.path.realpath(__file__))
        filepath = dir_path + '\\data\\refseq_tx_lengths.csv'
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

            # Find the number of singletons sequenced, split this into denovo fams and non-denovo fams
            dn_ix = n1 = s1 = None
            ix = 0
            for n, s, f_coeff in zip(num_families, num_successes, f_coefficients):
                if f_coeff == 1:
                    n1 = n
                    s1 = s
                    dn_ix = ix
                ix = ix + 1

            if n_lof_inh > 0:
                P_exact = calc_P_denovo(f, l, n1, n_lof_dn, n_missense_dn, n_lof_inh)
            else:
                P_exact = calc_P_denovo_lofdn_missdn(f, l, n1, n_lof_dn, n_missense_dn)
            if verbose >= 2:
                print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                    P_exact)
            P_combined = 0
            iter_max = min(n1 + 1, 30)  #in case n1 is large, assume that probability of finding >= 30 de novo or inherited variants is negligible, don't iterate those possible cases
            for s0_other in range(0, iter_max):
                for s1_other in range(0, iter_max):
                    for s2_other in range(0, iter_max):
                        P_other = calc_P_denovo(f, l, n1, s0_other, s1_other, s2_other)
                        if (n_lof_inh == 0 and s2_other == 0):
                            P_other = calc_P_denovo_lofdn_missdn(f, l, n1, s0_other, s1_other)
                        if P_other <= P_exact:
                            if n_lof_inh > 0 or (n_lof_inh == 0 and s2_other == 0):
                                if verbose >= 2:
                                    print 'Observing ' + str(s0_other) + ', ' + str(s1_other) + ', ' + str(s2_other) + ' successes is at least as significant (P = ' + str(
                                        P_other) + ')'
                                P_combined = P_combined + P_other
        print(P_combined)
        return 0

    # If we have a single type of pedigree structure:
    if len(num_successes) == 1:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            P_other = calc_P(f, num_families, [s0_other], f_coefficients)
            if P_other <= P_exact:
                if verbose >= 2:
                    print 'Observing ' + str(s0_other) + ' successes is at least as significant (P = ' + str(
                        P_other) + ')'
                P_combined = P_combined + P_other

    # If we have two types of pedigree structures:
    elif len(num_successes) == 2:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            for s1_other in range(0, num_families[1] + 1):
                P_other = calc_P(f, num_families, [s0_other, s1_other], f_coefficients)
                if P_other <= P_exact:
                    if verbose >= 2:
                        print 'Observing ' + str(s0_other) + ' and ' + str(
                            s1_other) + ' successes is at least as significant (P = ' + str(P_other) + ')'
                    P_combined = P_combined + P_other

    # If we have three types of pedigree structures:
    elif len(num_successes) == 3:
        # The probability of observing the exact same number of successes in the constellation of families that were sequenced:
        P_exact = calc_P_3(f, num_families, num_successes, f_coefficients)
        if verbose >= 2:
            print 'Probability of observing the exact same number of successes in the constellation of families that were sequenced: P = ' + str(
                P_exact)
        P_combined = 0
        for s0_other in range(0, num_families[0] + 1):
            for s1_other in range(0, num_families[1] + 1):
                for s2_other in range(0, num_families[2] + 1):
                    P_other = calc_P_3(f, num_families, [s0_other, s1_other, s2_other], f_coefficients)
                    if P_other <= P_exact:
                        if verbose >= 2:
                            print 'Observing ' + str(s0_other) + ' and ' + str(
                                s1_other) + ' and ' + str(
                                s2_other) + ' successes is at least as significant (P = ' + str(P_other) + ')'
                        P_combined = P_combined + P_other

    else:
        sys.exit(
            'Calculating P-value when more than 3 types of pedigree structures sequenced has not yet been implemented.')
    # TODO: use recursion to make algorithm work for any number of types of family structures

    if verbose == None:
        print(P_combined)
    elif verbose >= 1:

        signifstr = ''
        if P_combined < 0.05:
            signifstr = 'reaches'
        else:
            signifstr = 'does not reach'
        pvalstr = ', which ' + signifstr + ' genome-wide significance.'

        # TODO: incorporate information about follow-up sequencing
        # if ($f_n1 > 1):
        # $pval_followup_str = "Afterwards, the results of the follow-up sequencing were also taken into account, and this reduced the P-value from $pval_before_followup to $pval.";

        percentage = 100.0 * f

        print('Given that ' + "{0:.2f}".format(percentage) + ' % (f = ' + str(f)
              + ') of individuals in the population have a rare variant in the given gene, the statistical significance of seeing such variants in the constellation of families that were sequenced is P = ' + str(
                    P_combined) + pvalstr)

        print(
            '\nNote: this P-value must be corrected for the number of genes targeted, e.g. using the Bonferroni adjustment.')



def main():
    # command line arguments
    parser = argparse.ArgumentParser(
        description='''SORVA (Significance Of Rare VAriants) - Statistics package

Run SORVA to answer the question: How often do individuals carry a mutation in a given gene?
Next, calculate whether seeing such variants in a certain number of sequenced individuals is significant, given the background rate, and considering the family structures of the individuals sequenced.

The input parameter -p is the output of sorva.py. When running sorva.py, if sequenced individuals have homozygous mutations in the gene, use the option --zygosity hom. If sequenced individuals have heterozygous mutations in the gene, use --zygosity both.

IMPORTANT NOTE: Results must be corrected for multiple testing, by the number of genes sequenced (~20,000).''',
        epilog='''SORVA version 2.0 (c)2017-2018 Aliz R. Rao. All rights reserved.''',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', required=True, type=float,
                        help='background frequency / fraction of individuals in population who have a mutation in the gene of interest (output of sorva.py)')
    parser.add_argument('--num_families', '-n', required=True, type=str, default='',
                        help='The number of families sequenced, grouped by which families have identical family structures. E.g. if we sequenced 5 singletons and 3 sib pairs, this should be [5, 3]')
    parser.add_argument('--num_successes', '-s', required=True, type=str, default='',
                        help='The number of families sequenced who have a variant in the gene of interest (assumed to be shared IBD within the family), in the order that corresponds to the sets listed in the --num_families parameter. E.g. if we found that 3 singletons and 2 sib pairs had variants in a given gene, this would be [3, 2]')
    parser.add_argument('--f_coefficients', '-c', required=True, type=str, default='',
                        help='The coefficients for f. In case of an autosomal dominant (AD) disorder, this will usually be the coefficient of relationship. For autosomal recessive disorders, this will be 1 for unrelated individuals and 0.25 for sib pairs. For consanguineous pedigrees, see manuscript for method to calculate this coefficient. E.g. if we sequenced singletons and sib pairs and we\'re looking for het variants (disorder is AD), this would be [1, 0.25]')
    parser.add_argument('--denovo', required=False, type=int, default=0,
                        help='Out of the singletons who have a variant in the gene of interest, how many are de novo variants?')
    parser.add_argument('--gene', '-g', required=False, default='',
                        help='gene name; req\'d if denovo > 0. Use standard gene name from HGNC or specify Ensembl gene ID starting with ENSG.')
    parser.add_argument('--consequence', required=False, choices=['nonsyn', 'lof'], default='',
                        help='variant consequence filtering threshold; req\'d if denovo > 0. nonsyn = missense or more severe, also includes LOF. lof = potential loss-of-function (LOF), which includes stop loss gain, splice site mutation, and frameshift indel.')
    #parser.add_argument('--e_custom', required=False, type=float, default=0,
    #                    help='For consanguinious pedigrees, the coefficient is pow(0.5, E), where E is the number of independent edges in the paths connecting the sequenced individuals through a single common ancestor, in a path diagram. In other words, this many edges needed to have passed down the variant in order for all sequenced individuals to share the variants observed, assuming the variants are shared IBD.')
    parser.add_argument('--verbose', '-v', required=False, action='count',
                        help='output additional information on how P-value was calculated')

    args = parser.parse_args()

    run(f=args.f, num_families=args.num_families, num_successes=args.num_successes, f_coefficients=args.f_coefficients,
        denovo=args.denovo, gene=args.gene, consequence=args.consequence, verbose=args.verbose)

    return 0


if __name__ == "__main__": sys.exit(main())
