import click, numpy as np, pandas as pd, gzip, ete3, sys

@click.command()
@click.argument('region')
def main(region) :
    regions = pd.read_csv(region, sep='\t').values
    blocks = []
    seqs = []
    p = 0
    for r in regions :
        seqs.extend([0]*r[1] + [1]*r[2])
        blocks.append([p, len(seqs)])
        p = len(seqs)
    seqs = np.array(seqs)
    blocks = np.array(blocks)
    
    n_iter = 10000
    for x in np.arange(100) :
        simulations = np.zeros([3, regions.shape[0]])
        for ite in np.arange(n_iter) :
            sys.stderr.write(f'{x} {ite}\n')
            seqs = np.random.permutation(seqs)
            sim = np.array([np.sum(seqs[blk[0]:blk[1]]) for blk in blocks])
            simulations[0] += (regions.T[2] >= sim)
            simulations[1] += (regions.T[2] <= sim)
            simulations[2] += sim
                
        for i, r in enumerate(regions) :
            p0, p1, e = simulations.T[i]/n_iter
            print(f'{x}\t{r[0]}\t{r[1]}\t{r[2]}\t|\t{p0}\t{p1}\t|\t{e}')


if __name__ == '__main__' :
    main()
