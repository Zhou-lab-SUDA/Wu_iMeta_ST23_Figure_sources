import pandas as pd, numpy as np, click
from sklearn.metrics import mutual_info_score
from scipy.stats import entropy


def Theil_U(x, y):
    # Mutual information
    mi = mutual_info_score(x, y)
    # Entropy of x (A)
    p_x = pd.Series(x).value_counts(normalize=True)
    h_a = entropy(p_x, base=2)  # Base 2 for bits
    u = mi / h_a
    
    simulation = np.zeros(3000)
    for i in np.arange(3000) :
        idx = np.random.choice(x.size, x.size)
        mi_ = mutual_info_score(x[idx], y[idx])
        p_x_ = pd.Series(x[idx]).value_counts(normalize=True)
        h_a_ = entropy(p_x_, base=2)
        simulation[i] = mi_ / h_a_
    return u, np.std(simulation)


def Goodman_Kruskal_lambda(x, y):
    # Create contingency table
    cont_table = pd.crosstab(x, y)
    # Get row and column maxima
    row_maxes = cont_table.max(axis=1)
    col_maxes = cont_table.max(axis=0)
    # Calculate error without predictor
    E1 = len(x) - col_maxes.sum()
    # Calculate error with predictor
    E2 = len(x) - row_maxes.sum()
    # Calculate lambda
    if E1 == 0:
        return 0
    lambda_val = (E1 - E2) / E1
    return lambda_val


def pseudo_R(x, y) :
    import statsmodels.api as sm
    _, x = np.unique(x, return_inverse=True)
    _, y = np.unique(y, return_inverse=True)
    model = sm.MNLogit(x, y).fit()
    null_model = sm.MNLogit(x, pd.DataFrame([1]*x.size)).fit()
    return 1 - model.llf / null_model.llf

@click.command()
@click.argument('dat_file')
def main(dat_file) :
    data = pd.read_csv(dat_file, sep='\t', na_filter=False)
    for j, cj in enumerate(data.columns) :
        for i, ci in enumerate(data.columns[:j]) :
            di, dj = data.values.T[i], data.values.T[j]
            idx = (di != '-') & (dj != '-') & (di != '') & (dj != '')
            u, s = Theil_U(di[idx], dj[idx])
            print(f"{u*100.:.1f}% +- {s*100.:.1f}% of the uncertainty in {ci} is explained by {cj}")
            u, s = Theil_U(dj[idx], di[idx])
            print(f"{u*100.:.1f}% +- {s*100.:.1f}% of the uncertainty in {cj} is explained by {ci}")



if __name__ == '__main__' :
    main()
