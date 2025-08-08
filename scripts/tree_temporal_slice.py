import ete3_extensions, click

@click.command()
@click.option('-n','--nexus')
@click.option('-y', '--year', type=float)
def main(nexus, year) :
    tre = ete3_extensions.read_nexus(nexus)[0]
    
    distributions = {}
    for node in tre.iter_descendants('postorder') :
        parent = node.up
        y1 = node.annotations['date']
        y2 = parent.annotations['date']
        if y1 >= year and y2 < year :
            state = node.annotations['state']
            prop = node.annotations['state.prop']
            distributions[state] = distributions.get(state, 0) + prop
    
    for n, p in sorted(distributions.items()) :
        print(f'{n}\t{p}')




if __name__ == '__main__' :
    main()