import ete3_extensions, click, pandas as pd, numpy as np, os, collections


@click.command()
@click.option('-d', '--outdir', help='output folder')
@click.option('-n', '--nexus', help='nexus file')
@click.option('-t', '--trait', help='trait file')
@click.option('-m', '--at_most', help='default: 10', type=int, default=10)
def main(outdir, nexus, trait, at_most) :
    tre = ete3_extensions.read_nexus(nexus)[0]
    
    traits = collections.defaultdict(list)
    for n, t in pd.read_csv(trait, sep=',').values :
        traits[t].append(n)
    
    subsamples = []
    for t, names in traits.items() :
        n = names if len(names) <= at_most else np.random.choice(names, at_most, replace=False).tolist()
        subsamples.extend(n)
    
    subtre = ete3_extensions.prune(tre, subsamples)
    
    if not os.path.isdir(outdir) :
        os.makedirs(outdir)
    with open(os.path.join(outdir, 'dating.out.nex'), 'wt') as fout :
        fout.write(ete3_extensions.write_nexus([subtre]))



if __name__ == '__main__' :
    main()
