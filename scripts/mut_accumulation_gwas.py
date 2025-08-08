import numpy as np, pandas as pd, click, ete3, gzip, collections
import scipy.stats as stats
from scipy.stats import f_oneway, levene, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# identify branches leading to each genome
# load in mutations.gz
# remove those in recombination
# count mut frequency for each genome
@click.command()
@click.option('-t', '--tree')
@click.option('-m', '--mutation')
@click.option('-r', '--recombination')
@click.option('-d', '--diverged')
@click.option('-g', '--group')
def main(tree, mutation, recombination, diverged, group) :
    tre = ete3.Tree(tree, format=1)
    path_to_leaves = []
    ignored_nodes = [tre] + tre.get_children()
    for leaf in tre.get_leaves() :
        path_to_leaves.append([leaf.name, [leaf.name]])
        for n in leaf.get_ancestors() :
            if n not in ignored_nodes :
                path_to_leaves[-1][1].append(n.name)

    rec_regions = {}
    with open(recombination, 'rt') as fin :
        for line in fin :
            if line.startswith('\tImportation') :
                p = line.strip().split('\t')
                if p[1] not in rec_regions :
                    rec_regions[p[1]] = {}
                if p[2] not in rec_regions[p[1]] :
                    rec_regions[p[1]][p[2]] = []
                rec_regions[p[1]][p[2]].append([int(p[3]), int(p[4])])

    div_regions = []
    with open(diverged, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            div_regions.append([p[0], p[1], int(p[2]), int(p[3])])

    regions = {}
    for d in div_regions :
        if d[1] not in regions or regions[d[1]] < int(d[3]) :
            regions[d[1]] = int(d[3])

    for n, r in regions.items() :
        regions[n] = np.empty(r, dtype=int)
        regions[n][:] = -1
    
    for idx, d in enumerate(div_regions) :
        regions[d[1]][d[2]-1:d[3]] = idx
    
    branches = {}
    with gzip.open(mutation, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            if p[0] not in branches :
                branches[p[0]] = []
            s = int(p[2])
            if p[1] not in regions or s >= regions[p[1]].size or regions[p[1]][s-1] < 0 :
                continue
            inRec = False
            if p[0] in rec_regions and p[1] in rec_regions[p[0]] :
                rec = rec_regions[p[0]][p[1]]
                for r_s, r_e in rec :
                    if r_s <= s and r_e >= s :
                        inRec =True
                        break
                    if r_s > s :
                        break
            if not inRec :
                branches[p[0]].append(regions[p[1]][s-1])
    
    leaf_snps = []
    for idx, (leaf, path) in enumerate(path_to_leaves) :
        snps = np.bincount([reg for n in path for reg in branches.get(n, [])], minlength=len(div_regions))
        leaf_snps.append(snps)
    leaf_snps = np.array(leaf_snps)
    leaf_snps = pd.DataFrame(leaf_snps, index=[p[0] for p in path_to_leaves], columns=[d[0] for d in div_regions])
    
    # with open('tmp.out', 'wt') as fout :
    #     fout.write(leaf_snps.to_csv())
    
    groups = {}
    with open(group, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            if p[1] not in ('', '0') and p[0] in leaf_snps.index :
                if p[1] not in groups :
                    groups[p[1]] = []
                groups[p[1]].append(p[0])
    
    groups['others'] = []
    for grp, members in list(groups.items()) :
        if len(members) <= 5 and grp != 'others' :
            groups['others'].extend(members)
            groups.pop(grp)

    country = {}
    for grp, members in groups.items() :
        for m in members :
            country[m] = grp

    groups = list(groups.items())
    for region in leaf_snps.columns :
        stat = leaf_snps[region]
        stat = stat.loc[[s in country for s in stat.index]]
        sgrp = [country[s] for s in stat.index]
        stat2 = stat.to_dict()
        grp_stat = [[stat2.get(mm) for mm in m if mm in stat2] for g, m in groups]
    
        max_idx = np.argmax([np.mean(g) for g in grp_stat])
        statistic, pvalue = f_oneway(*grp_stat)
        print(f"{region}\tANOVA F-statistic = {statistic:.2f}, p-value = {pvalue:.4f}; Bonferroni = {pvalue*leaf_snps.columns.size:.4f}; Mean = {np.mean(stat):.4f}; Max = {groups[max_idx][0]}; Max_val = {[np.mean(g) for g in grp_stat][max_idx]}")


if __name__ == '__main__' :
        main()