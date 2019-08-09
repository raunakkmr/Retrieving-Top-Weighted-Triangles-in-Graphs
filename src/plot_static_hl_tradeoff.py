import matplotlib.pyplot as plt
import numpy as np

def make_plot(thresholds, accs, times, k, dataset):

  thresholds = [0] + thresholds
  times = [0] + times
  accs = [0] + accs

  plt.subplot(2, 1, 1)
  for y in np.arange(0, 1, 0.2):
    plt.plot(range(0, 100), [y] * len(range(0, 100)), "--", lw=0.5, color="black", alpha=0.3)
  plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on") 
  plt.ylim(0, accs[-1])
  plt.xlim(0, thresholds[-1])
  plt.plot(thresholds, accs, label='accuracy', color='green')
  plt.ylabel('accuracy')
  # plt.title('accuracy/time vs threshold trade-off\n for {} and k = {}'.format(dataset, k))

  plt.subplot(2, 1, 2)
  if dataset == 'eth':
    step = 20
  else:
    step = 2.5
  for y in np.arange(0, int(times[-1])+1, step):
    plt.plot(range(0, 100), [y] * len(range(0, 100)), "--", lw=0.5, color="black", alpha=0.3)
  plt.tick_params(axis="both", which="both", bottom="off", top="off",
                labelbottom="on", left="off", right="off", labelleft="on") 
  plt.ylim(0, times[-1])
  plt.xlim(0, thresholds[-1])
  plt.plot(thresholds, times, label='time (s)', color='red')
  plt.ylabel('time (s)')
  plt.xlabel('thresholds: percentage of edges labelled heavy')

  plt.savefig('../figs/static_hl_tradeoff_{}_{}'.format(dataset, k))
  plt.close()

def get_info_and_plot(dataset, lines):
  sort_time = float(lines[7].split()[-1])
  idx = -1
  kvals = [25, 1000, 40000]
  times, accs = [], []
  for line in lines:
    if 'brute' in line:
      idx += 1
      k = kvals[idx]
      times, accs = [], []
    if 'Total Time' in line:
      times.append(float(line.split()[-1])+sort_time)
    if 'Accuracy' in line:
      accs.append(float(line.split()[-1]))
    if len(accs) == len(thresholds):
      make_plot(thresholds, accs, times[1:], k, dataset)
      times, accs = [], []

thresholds = list(range(5, 105, 5))
with open('../output/static_hl') as f:
  lines = f.readlines()
  get_info_and_plot('tags-stack-overflow', lines[:1055])
  get_info_and_plot('eth', lines[1055:])