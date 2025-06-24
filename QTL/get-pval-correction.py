#!/usr/bin/env python3

import subprocess

# Get total number of tests 
m = int(subprocess.check_output(['wc', '-l', 'chrALL.sorted.txt']).split()[0])

Q = 0
prev_val = 0

# Step through sorted nominal p-values
# Incrementing rank counter by 1
# Calculate adjusted p-val as 
# p_adj = p_nominal * {m, # tests} / {q, rank)
# Break when over alpha of 0.05
with open('chrALL.sorted.txt', 'r') as infile:
    with open('pval-correction-mapping.txt', 'w') as outfile:
        outfile.write('p_nominal\tp_bh\tp_Bonferroni\n')
        for line in infile:
            Q += 1
            p = float(line.strip())
            p_adj = (p*m)/Q
            p_Bonferroni = p * m
            if p_adj < prev_val:
                p_adj = prev_val
            outfile.write(f'{p}\t{p_adj}\t{p_Bonferroni}\n')
            if p_adj > 0.05:
                break
            prev_val = p_adj


