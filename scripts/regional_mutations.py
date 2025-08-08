import click, gzip, re, numpy as np


@click.command()
@click.option('-g', '--gff', help='gff file')
@click.option('-m', '--mutations', help='mutations.gz file')
@click.option('-r', '--rechmm', help='RecHMM file')
def main(gff, mutations, rechmm) :
    rec = {}
    with open(rechmm, 'rt') as fin :
        for line in fin :
            if line.startswith('\tImportation') :
                p = line.strip().split('\t')
                if p[1] not in rec :
                    rec[p[1]] = []
                rec[p[1]].append([p[2], int(p[3]), int(p[4])])
    rec = { n:sorted(r) for n, r in rec.items() }
    genes = []
    seqLen = {}
    with gzip.open(gff, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                continue
            p = line.strip().split('\t')
            if p[2] == 'region' :
                seqLen[p[0]] = int(p[4])
            if p[2] == 'CDS' :
                n = re.findall(r'Parent=([^;]+)', p[8])[0]
                genes.append([p[0], int(p[3]), int(p[4]), n])
    intergenic = []
    p = ['', 0, 0, '']
    for g in genes :
        if g[0] != p[0] :
            if p[0] in seqLen and p[2] < seqLen[p[0]] :
                intergenic.append([p[0], p[2]+1, seqLen[p[0]], f'{p[3]}__'])
            if g[1] > 1 :
                intergenic.append([g[0], 1, g[1]-1, f'__{g[3]}'])
            p = g[:]
        elif g[1] - p[2] > 1 :
            intergenic.append([g[0], p[2]+1, g[1]-1, f'{p[3]}__{g[3]}'])
        if g[2] > p[2] :
            p = g[:]

    regions = sorted(genes + intergenic)
    seqs = {n:np.repeat([-1], s).astype(np.int16) for n, s in seqLen.items() }
    for i, (n, s, e, g) in enumerate(regions) :
        if n in seqs :
            seqs[n][s-1:e] = i
    
    region_mutations = np.zeros(len(regions)*6, dtype=int).reshape([len(regions), 6])
    with gzip.open(mutations, 'rt') as fin :
        for line in fin :
            if line.startswith('#') :
                if line.startswith('## Missing_region:') :
                    _, _, n, s, e = line.strip().split()
                    seqs[n][int(s)-1:int(e)] = -1
                continue
            p = line.strip().split('\t')
            i = seqs[p[1]][int(p[2])-1]
            if i < 0 :
                continue
            inRec = 0
            for r in rec.get(p[0], []) :
                if (p[1] < r[0]) or ((p[1] == r[0]) and (int(p[2]) > r[2])) :
                    break
                if (p[1] == r[0]) and (int(p[2]) >= r[1]) :
                    inRec = 1
            if inRec :
                continue
            t = 4
            if len(p) > 5 :
                if (p[5].find('Nonsyn') > 0) :
                    t = 1
                elif (p[5].find('Syn') > 0) :
                    t = 0
                elif (p[5].find('Nonsense') > 0) or (p[5].find('Frameshift') > 0) :
                    t = 2
                elif (p[5].find('Indel') > 0) :
                    t = 3
            region_mutations[i][t] += 1
    reg, size = np.unique([ss for s in seqs.values() for ss in s], return_counts=True)
    for r, s in zip(reg, size) :
        if r >= 0 :
            region_mutations[r][5] = s
    print('#Gene\tS\tNS\tDisruption\tINDEL\tIntergenic\tSize')
    for r, m in zip(regions, region_mutations) :
        if m[5] > 0 :
            print(f'{r[3]}\t{m[0]}\t{m[1]}\t{m[2]}\t{m[3]}\t{m[4]}\t{m[5]}')



if __name__ == '__main__' :
    main()