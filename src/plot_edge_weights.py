import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import powerlaw
from scipy.stats import binned_statistic

def drop_zeros(a_list):
    return [i for i in a_list if i>0]
    
# https://stackoverflow.com/questions/16489655/plotting-log-binned-network-degree-distributions
def log_binning(counter_dict,bin_count=35):
    max_x = np.log10(max(counter_dict.keys()))
    max_y = np.log10(max(counter_dict.values()))

    max_base = max([max_x,max_y])

    min_x = np.log10(min(drop_zeros(counter_dict.keys())))

    bins = np.logspace(min_x,max_base,num=bin_count)
    
    # Based off of: http://stackoverflow.com/questions/6163334/binning-data-in-python-with-scipy-numpy
    # bin_means_y = (np.histogram(list(counter_dict.keys()),bins,weights=list(counter_dict.values()))[0] / np.histogram(list(counter_dict.keys()),bins)[0])
    # bin_means_x = (np.histogram(list(counter_dict.keys()),bins,weights=list(counter_dict.keys()))[0] / np.histogram(list(counter_dict.keys()),bins)[0])

    bin_means = binned_statistic(list(counter_dict.keys()),
                                 [list(counter_dict.keys()), list(counter_dict.values())],
                                 bins=bins)[0]
    
    return bin_means

def fit_powerlaw(x, y, xmin):
    l = [([xx]*yy) for (xx, yy) in zip(x, y)]
    l = [item for sublist in l for item in sublist]
    dist = powerlaw.Fit(l, xmin=xmin, discrete=True)
    return dist

datasets = ['tags-stack-overflow', 'wikipedia']
for dataset in datasets:
    path = '../output/edge_weights_{}'.format(dataset)
    with open(path, 'r') as f:
        lines = f.readlines()
        x = [int(line.split()[0]) for line in lines]
        y = [int(line.split()[1]) for line in lines]
        m = np.sum(y)
        xmin = None
        if dataset == 'wikipedia':
            xmin = 11
        dist = fit_powerlaw(x, y, xmin)
        print('Datset: {}'.format(dataset))
        print('alpha: {}, xmin: {}'.format(dist.power_law.alpha, dist.power_law.xmin))
        ba_c2 = {xx : yy / m for (xx,yy) in zip(x,y)}

    ba_x,ba_y = log_binning(ba_c2,50)

    fig, ax = plt.subplots()
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y/m), max(y/m))
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(ba_x, ba_y, color='red', marker='s', s=50)
    alpha, xmin = dist.power_law.alpha, dist.power_law.xmin
    p_x = [t for t in ba_x if t >= xmin]
    p_y = [(alpha-1)/xmin * (t/xmin)**(-alpha) for t in p_x]
    ax.plot(p_x, p_y)
    ax.get_lines()[0].set_color('blue')
    plt.xlabel('edge weight', fontsize=20)
    plt.ylabel('fraction', fontsize=20)
    plt.title('Dataset: {}'.format(dataset), fontsize=20)

    red_patch = mpatches.Patch(color='red', label='empirical fraction')
    blue_patch = mpatches.Patch(color='blue', label='best fit powerlaw')
    plt.legend(handles=[blue_patch, red_patch], fontsize=20)

    fig.savefig('../figs/edge_weights_{}.pdf'.format(dataset), bbox_inches='tight')
    plt.close()