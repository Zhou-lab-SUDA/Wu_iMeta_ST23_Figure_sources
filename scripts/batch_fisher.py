import os, sys, numpy as np, pandas as pd, click, scipy.stats


@click.command()
@click.argument('d_file')
def main(d_file) :
    data = pd.read_csv(d_file, header=None, sep='\t').values
    n_test = data.shape[1]
    
    n_size, n_mut = data.T[1].sum(), data.T[2].sum()
    probs = []
    for gene, s, m in data :
        s0, m0 = n_size - s, n_mut - m
        p = scipy.stats.fisher_exact([[m, max(0, s-m)], [m0, s0-m0]])
        probs.append(p.pvalue)
    
    probs = np.array(probs)
    ps2 = probs*probs.size #scipy.stats.false_discovery_control(probs)
    
    with open(d_file +'.res', 'wt') as fout :
        for (gene, s, m), p1, p2 in zip(data, probs, ps2) :
            fout.write(f'{gene}\t{s}\t{m}\t{p1}\t{p2}\n')





if __name__ == '__main__' :
    main()
