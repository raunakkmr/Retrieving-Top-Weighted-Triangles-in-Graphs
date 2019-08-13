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

def fit_powerlaw(x, y):
    l = [([xx]*yy) for (xx, yy) in zip(x, y)]
    l = [item for sublist in l for item in sublist]
    dist = powerlaw.Fit(l, discrete=True)
    return dist

datasets = ['tags-stack-overflow', 'threads-stack-overflow'] #, 'MAG']
for dataset in datasets:
    path = '../output/edge_weights_{}'.format(dataset)
    with open(path, 'r') as f:
        lines = f.readlines()
        x = [int(line.split()[0]) for line in lines]
        y = [int(line.split()[1]) for line in lines]
        m = np.sum(y)
        dist = fit_powerlaw(x, y)
        print('Datset: {}'.format(dataset))
        print('alpha: {}, xmin: {}'.format(dist.power_law.alpha, dist.power_law.xmin))
        ba_c2 = {xx : yy / m for (xx,yy) in zip(x,y)}

    ba_x,ba_y = log_binning(ba_c2,50)

    _, ax = plt.subplots()
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y/m), max(y/m))
    plt.xscale('log')
    plt.yscale('log')
    ax.scatter(ba_x, ba_y, color='red', marker='s', s=50)
    dist.plot_pdf(ax)
    ax.get_lines()[0].set_color('blue')
    plt.xlabel('edge weight', fontsize=12)
    plt.ylabel('fraction', fontsize=12)
    plt.title('Dataset: {}'.format(dataset))

    red_patch = mpatches.Patch(color='red', label='empirical fraction')
    blue_patch = mpatches.Patch(color='blue', label='best fit powerlaw')
    plt.legend(handles=[blue_patch, red_patch])

    plt.savefig('../figs/edge_weights_{}'.format(dataset))
    plt.close()

# hist, bins, _ = plt.hist(x_a, weights=y_a, bins=10)
# plt.close()
# logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
# plt.hist(x_a, weights=y_a, bins=logbins, log=True, color=['red','green','blue'], label=datasets)
# plt.xscale('log')
# plt.xlabel('edge weight', fontsize=12)
# plt.ylabel('count', fontsize=12)
# plt.legend()
# plt.savefig('../figs/edge_weights')
# plt.close()