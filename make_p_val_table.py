#!/usr/bin/env python3

import sys
import csv
import tqdm

import scipy.special as sp
from scipy.optimize import minimize
from scipy.stats import chi2
import numpy as np


def L_fun(k, n, a, b):
    if a < 0 or b < 0:
        return -1000000
    result = 0
    for pos_n, pos_k in zip(n, k):
        result += -np.log(pos_n + 1) - sp.betaln(pos_n - pos_k + 1, pos_k + 1) - sp.betaln(a, b) \
            + sp.betaln(a + pos_k, b + pos_n - pos_k)
    return result


def L(k, n):
    return lambda ab: -L_fun(k, n, ab[0], ab[1])


def write_empty(chrom, start, end, writer):
    out_row = dict(chrom=chrom, start=start, end=end, p_val=1.0)
    writer.writerow(out_row)


def write_p_value(chrom, start, end, k_1, n_1, k_2, n_2, writer, method='nelder-mead'):
    if not n_1:
        assert not n_2
        write_empty(chrom, start, end, writer)
        return

    reads1 = ','.join('{}/{}'.format(x, y) for x, y in zip(k_1, n_1))
    reads2 = ','.join('{}/{}'.format(x, y) for x, y in zip(k_2, n_2))

    x0 = np.array([2, 2])
    min_l_1 = minimize(L(k_1, n_1), x0, method=method)
    a1, b1 = min_l_1.x
    l1 = -min_l_1.fun

    min_l_2 = minimize(L(k_2, n_2), x0, method=method)
    a2, b2 = min_l_2.x
    l2 = -min_l_2.fun

    min_l_0 = minimize(L(k_1 + k_2, n_1 + n_2), x0, method=method)
    a0, b0 = min_l_0.x
    l0 = -min_l_0.fun

    D = 2 * (l1 + l2 - l0)
    p_val = 1 - chi2.cdf(D, 2)

    out_row = dict(chrom=chrom, start=start, end=end, reads1=reads1, reads2=reads2,
        a1=a1, b1=b1, l1=l1, a2=a2, b2=b2, l2=l2, a0=a0, b0=b0, l0=l0, D=D, p_val=p_val, n_pos=len(k_1))
    writer.writerow(out_row)


def interval_table(reader, writer, step=300):
    chromosome = None
    start = None
    end = None
    k_1 = []
    n_1 = []
    k_2 = []
    n_2 = []
    for row in tqdm.tqdm(reader):
        row_start = int(row['start'])

        if start is None or row['chromosome'] != chromosome or row_start >= end:
            if k_1:
                write_p_value(chromosome, start, end, k_1, n_1, k_2, n_2, writer)
                start += step
                end = start + step
            else:
                chromosome = row['chromosome']
                start = row_start
                end = start + step

            if row['chromosome'] != chromosome:
                chromosome = row['chromosome']
                start = row_start
            else:
                # Now:     end <= row_start
                # We need: start <= row_start < end
                while end <= row_start:
                    write_empty(chromosome, start, end, writer)
                    start += step
                    end = start + step

            k_1.clear()
            n_1.clear()
            k_2.clear()
            n_2.clear()

        curr_k1 = int(row['met_hp_1'])
        curr_k2 = int(row['met_hp_2'])
        k_1.append(curr_k1)
        k_2.append(curr_k2)
        n_1.append(int(row['unm_hp_1']) + curr_k1)
        n_2.append(int(row['unm_hp_2']) + curr_k2)

    if k_1:
        write_p_value(chromosome, start, end, k_1, n_1, k_2, n_2, writer)


def main():
    with open(sys.argv[1]) as inp, open(sys.argv[2], 'w') as outp:
        reader = csv.DictReader(inp, delimiter='\t')
        outp.write('# %s\n' % ' '.join(sys.argv))
        writer = csv.DictWriter(outp, delimiter='\t', quoting=csv.QUOTE_NONE, dialect='unix', restval='NA',
            fieldnames=["chrom", "start", "end", "reads1", "reads2", "a1", "b1", "l1", "a2", "b2", "l2",
            "a0", "b0", "l0", "D", "n_pos", "p_val"])
        writer.writeheader()

        if len(sys.argv) > 3:
            step = int(sys.argv[3])
        else:
            step = 300
        interval_table(reader, writer, step=step)


if __name__ == '__main__':
    main()