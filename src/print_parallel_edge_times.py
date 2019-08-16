import argparse
import os

files = [
  'tags-stack-overflow',
  'threads-stack-overflow',
  'wikipedia',
  'eth',
  'aminer',
  'temporal-reddit-reply',
  'MAG',
  'spotify'
]

parser = argparse.ArgumentParser()
parser.add_argument('-k')
args = parser.parse_args()

for fname in files:
  path = os.path.join('../output/compare_parallel_edge_{}'.format(args.k), fname)
  with open(path, 'r') as f:
    lines = f.readlines()
    sort_time = 0
    pre_times = []
    times = []
    total_times = []
    accuracies = []
    for line in lines:
      if "Sort time" in line:
        sort_time = float(line.split()[-1])
      if "Pre-processing" in line:
        pre_times.append(float(line.split()[-1]))
      if "Total Time" in line:
        times.append(float(line.split()[-1]))
      if "Accuracy" in line:
        accuracies.append(float(line.split()[-1]))
    for (x, y) in zip(pre_times, times):
      total_times.append(sort_time+x+y)
    print("Dataset: {}".format(fname))
    for (t, a) in zip(total_times, accuracies):
      print("Time: {:.3f}, accuracy: {:.3f}".format(t, a))
    print("---------------------------")
