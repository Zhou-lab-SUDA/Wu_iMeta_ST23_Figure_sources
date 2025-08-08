# import seq
# import group
# seq pairwise distrance
# fisher

import click, numpy as np
from scipy.stats import fisher_exact


def pairwise(sequence) :
    seq = list(sequence.items())
    pairs = []
    for i, (n0, ss0) in enumerate(seq) :
        for j, (n1, ss1) in enumerate(seq[:i]) :
            tot_size = (ss0 != '-') & (ss1 != '-')
            if np.sum(tot_size) < 10 :
                continue
            snp = (ss0 != ss1) & tot_size
            dist = np.sum(snp)/np.sum(tot_size)
            pairs.append([n0, n1, dist])
    return pairs



@click.command()
@click.option('-s', '--seq')
@click.option('-g', '--group')
@click.option('-d', '--dists', default='0.001,0.005,0.01,0.015,0.02,0.04')
def main(seq, group, dists) :
    groups = {}
    individuals = {}
    with open(group, 'rt') as fin :
        for line in fin :
            p = line.strip().split()
            groups[p[0]] = p[2]
            individuals[p[0]] = p[1].split('_')[0]

    sequence = {}
    with open(seq, 'rt') as fin :
        for line in fin :
            if line.startswith('>') :
                n = line[1:].strip().split()[0]
                sequence[n] = []
            else :
                sequence[n].extend(line.strip().split())
    for n, s in sequence.items() :
        sequence[n] = ''.join(s).upper()
    sequence = {n:np.array(list(s)) for n, s in sequence.items() if n.find('|') and n.split('|')[0] in groups}
    print(len(sequence))
    sample_pair = pairwise(sequence)

    for dist in dists.split(',') :
        dist = float(dist)
        cross_table = [[0, 0], [0, 0]]
        for n0, n1, d in sample_pair :
            try :
                if individuals[n0.split('|')[0]] == individuals[n1.split('|')[0]] :
                    continue
                if groups[n0.split('|')[0]] == groups[n1.split('|')[0]] :
                    j = 0
                else :
                    j = 1
                if d <= dist :
                    i = 0
                else :
                    i = 1
                cross_table[i][j] += 1
            except :
                continue
        
        p = fisher_exact(cross_table)[1]
        print(f'{seq}\t{dist}\t{p}\t{cross_table[0][0]}\t{cross_table[0][1]}\t{cross_table[1][0]}\t{cross_table[1][1]}')



if __name__ == '__main__' :
    main()