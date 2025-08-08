import numpy as np, pandas as pd, click, ete3, gzip, collections, re
import scipy.stats as stats
from scipy.stats import f_oneway, levene, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd

@click.command()
@click.option('-t', '--tree')
@click.option('-m', '--mutation')
@click.option('-r', '--recombination')
@click.option('-G', '--gff')
@click.option('-g', '--group')
def main(tree, mutation, recombination, gff, group) :
    tre = ete3.Tree(tree, format=1)
    path_to_leaves = []
    ignored_nodes = [tre] + tre.get_children()
    for leaf in tre.get_leaves() :
        path_to_leaves.append([leaf.name, [leaf.name]])
        for n in leaf.get_ancestors() :
            if n not in ignored_nodes :
                path_to_leaves[-1][1].append(n.name)

    regions = []
    prev = ['', 0, '']
    with gzip.open(gff, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            p[3], p[4] = int(p[3]), int(p[4])
            gene_id = p[-1].split(';', 2)[1].split('=')[-1]
            if p[2] == 'CDS' :
                if p[0] != prev[0] :
                    prev = [p[0], 0, '']
                if prev[1] + 1 < p[3] :
                    regions.append([p[0], prev[1]+1, p[3]-1, f'{prev[2]}__{gene_id}'])
                regions.append([p[0], p[3], p[4], gene_id])
                if p[4] >= prev[1] :
                    prev = [prev[0], p[4], gene_id]
    regions.sort()

    rec_regions = {}
    # with open(recombination, 'rt') as fin :
    #     for line in fin :
    #         if line.startswith('\tImportation') :
    #             p = line.strip().split('\t')
    #             if p[1] not in rec_regions :
    #                 rec_regions[p[1]] = {}
    #             if p[2] not in rec_regions[p[1]] :
    #                 rec_regions[p[1]][p[2]] = []
    #             rec_regions[p[1]][p[2]].append([int(p[3]), int(p[4])])

    branches = {}
    with gzip.open(mutation, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            if p[0] not in branches :
                branches[p[0]] = []
            s = int(p[2])

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
                branches[p[0]].append([(p[1], s)] + re.split(r'->', p[4]))

    leaf_snps = {}
    for idx, (leaf, path) in enumerate(path_to_leaves) :
        if idx % 100 == 0 :
            print(idx)
        leaf_snps[leaf] = {}
        sites = {}
        for br in path[::-1] :
            for site, anc, desc in branches.get(br, []) :
                if site not in sites :
                    sites[site] = [anc, desc]
                else :
                    sites[site][1] = desc
        x = 0
        for site, snvs in sorted(sites.items()) :
            if snvs[0] == snvs[1] :
                # print()
                continue
            while x < len(regions) :
                reg = regions[x]
                if reg[0] < site[0] or ((reg[0] == site[0]) and reg[2] < site[1]) :
                    x += 1
                else :
                    break
            for i in range(x, len(regions)) :
                reg = regions[i]
                if reg[0] > site[0] or ((reg[0] == site[0]) and reg[1] > site[1]) :
                    break
                if (reg[0] == site[0]) and (reg[2] >= site[1]) and  (reg[1] <= site[1]) :
                    leaf_snps[leaf][reg[3]] = leaf_snps[leaf].get(reg[3], 0) + 1

    groups = {}
    with open(group, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            if p[1] not in ('', '0') and p[0] in leaf_snps.keys() :
                if p[1] not in groups :
                    groups[p[1]] = []
                groups[p[1]].append(p[0])

    groups['others'] = []
    for grp, members in list(groups.items()) :
        if len(members) <= 5 and grp != 'others' :
            groups['others'].extend(members)
            groups.pop(grp)
    if len(groups['others']) < 1 :
        groups.pop('others')

    # states = []
    # for grp, members in groups.items() :
    #     for m in members :
    #         states.append([m, grp])
    # states = np.array(states)

    sorted_groups = sorted(groups.items())
    with open('tmp', 'wt') as fout :
        for reg in regions :
            dat = []
            for grp, members in sorted(groups.items()) :
                dat.append([ leaf_snps.get(m, {}).get(reg[3], 0) for m in members ])
            means = [np.mean(g) for g in dat]
            if max(means) == min(means) :
                continue
            max_idx = int(np.argmax(means))
            statistic, pvalue = f_oneway(*dat)
            means = [np.mean(g) for g in dat]
            means = '\t'.join([f'{g[0]}: {m:.4f}' for g, m in zip(sorted_groups, means)])
            print(f"{reg[0]}\t{reg[1]}\t{reg[2]}\t{reg[3]}\tANOVA F-statistic = {statistic:.2f}, p-value = {pvalue:.4g}; Max = {sorted_groups[max_idx][0]};\t{means}")
            fout.write(f"{reg[0]}\t{reg[1]}\t{reg[2]}\t{reg[3]}\tANOVA F-statistic = {statistic:.2f}, p-value = {pvalue:.4g}; Max = {sorted_groups[max_idx][0]};\t{means}\n")



if __name__ == '__main__' :
        main()
