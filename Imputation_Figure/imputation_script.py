### you will need to install msprime, typically easiest with pip IMHO

#!/usr/bin/env python3
import msprime as msp
import tskit
import numpy as np
import random

rng = np.random.default_rng(20251010)

# Human-like Parameters for a 200kb region
# I bet I could go small, like 100kb, maybe even 50kb
# this would give stronger haplotypes
n_hap = 20 # number of haplotypes (rows)
L = 200000 # 200 kb
Ne = 10000
mu = 1e-8
recomb = 1e-8

# Simulate tree sequence with mutations
# IMPORTANT: simulate haploid samples by setting ploidy=1
# This yields exactly n_hap genotypes per variant

ts = msp.sim_ancestry(
    samples=n_hap, ploidy=1, recombination_rate=recomb, sequence_length=L,
    population_size=Ne, random_seed=20251010
)
ts = msp.sim_mutations(ts, rate=mu, random_seed=20251011)

# Collect biallelic 0/1 variants only
variants = []
for var in ts.variants():
    if len(var.alleles) >= 2:
        g = var.genotypes.astype(np.int8)
        unique = np.unique(g)
        if set(unique.tolist()).issubset({0, 1}):
            variants.append(g)

if len(variants) < 20:
    raise SystemExit(f"Only {len(variants)} segregating biallelic SNPs; increase parameters and rerun.")

G = np.vstack(variants) # shape: num_sites x n_hap
# Sanity: enforce expected width
if G.shape[1] != n_hap:
    raise SystemExit(f"Unexpected number of genotypes per site: {G.shape[1]} (expected {n_hap}).")

# Compute MAF per site
allele_counts = G.sum(axis=1)
maf = np.minimum(allele_counts, n_hap - allele_counts) / n_hap

# Select top-20 MAF SNPs
idx = np.argsort(maf)[::-1][:20]
G20 = G[idx, :].T # shape: n_hap x 20 (rows=hyp, cols=sites)
maf20 = maf[idx]

# Map 0/1 to allele letters per SNP using allowed pairs
snp_pairs_list = [
    ("A", "T"), ("A", "C"), ("A", "G"), ("C", "G"), ("C", "T"), ("G", "T")
]
assignments = rng.integers(0, len(snp_pairs_list), size=G20.shape[1])
allele_pairs = [snp_pairs_list[k] for k in assignments]

# Build Case 0 matrix of letters
case0 = np.empty_like(G20, dtype='<U1')
for j in range(G20.shape[1]):
    a0, a1 = allele_pairs[j]
    case0[:, j] = np.where(G20[:, j] == 0, a0, a1)

# Case 1: lower panel rows 10..19 keep 6 columns, others '?'
case1 = case0.copy()
keep_cols = rng.choice(G20.shape[1], size=6, replace=False)
for i in range(10, 20):
    for j in range(G20.shape[1]):
        if j not in keep_cols:
            case1[i, j] = '?'

# Case 2: lower panel 80% missing at random (per cell)
case2 = case0.copy()
for i in range(10, 20):
    missing = rng.choice(G20.shape[1], size=int(np.ceil(0.8 * G20.shape[1])), replace=False)
    case2[i, missing] = '?'

# Case 3: all rows 80% missing at random (per cell)
case3 = case0.copy()
for i in range(0, 20):
    missing = rng.choice(G20.shape[1], size=int(np.ceil(0.8 * G20.shape[1])), replace=False)
    case3[i, missing] = '?'


def format_panel(mat: np.ndarray) -> list[str]:
    lines = []
    for i in range(mat.shape[0]):
        row = ''.join(mat[i, :])
        if i == 9:
            lines.append(row + "\n")
        else:
            lines.append(row)
    return lines

out = []
out.append("IMPUTATION METHODS FIGURE (msprime-based)")
out.append("==========================================")
out.append("")

out.append("CASE 0: Complete Genotype Data (Top-20 MAF SNPs)")
out.append("===============================================")
out.extend(format_panel(case0))
out.append("")

out.append("CASE 1: Reference Panel + SNP Chip (6 columns known in lower panel)")
out.append("==================================================================")
out.extend(format_panel(case1))
out.append("")

out.append("CASE 2: Reference Panel + Low Coverage Sequencing (0.2X)")
out.append("========================================================")
out.extend(format_panel(case2))
out.append("")

out.append("CASE 3: All Individuals with 80% Missing Data (no reference)")
out.append("==========================================================")
out.extend(format_panel(case3))

with open('formatted.txt', 'w') as f:
    f.write('\n'.join(out) + '\n')

# Print a brief summary to stdout
for j, (a0, a1) in enumerate(allele_pairs, start=1):
    print(f"Col {j:02d}: {a0} vs {a1} (MAF={maf20[j-1]:.2f})")
