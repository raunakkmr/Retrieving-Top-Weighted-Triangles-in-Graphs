import matplotlib.pyplot as plt
import numpy as np

datasets = ['tags-stack-overflow', 'threads-stack-overflow', 'MAG']
for dataset in datasets:
    path = '../output/edge_weights_{}'.format(dataset)
    with open(path, 'r') as f:
        lines = f.readlines()
        x = [np.log(int(line.split()[0])) for line in lines]
        y = [np.log(int(line.split()[1])) for line in lines]
        xmax = int(np.max(x))
        for y_ in np.arange(0, int(np.max(y))+1, 2):
            plt.plot(range(0, xmax+1), [y_] * len(range(0, xmax+1)), lw=-.5, color="black", alpha=0.3)
        plt.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="on", left="off", right="off", labelleft="on") 
        plt.ylim(0, np.max(y))
        plt.xlim(0, np.max(x))
        plt.plot(x, y, color='red')
        plt.xlabel('edge weight')
        plt.ylabel('count')
        plt.title('Dataset: {}'.format(dataset))
        plt.savefig('../figs/edge_weights_{}'.format(dataset))
        plt.close()