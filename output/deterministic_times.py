import argparse
import os

files = [
  'congress-bills',
  'threads-math-sx',
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
  path = os.path.join(os.getcwd(), 'compare_deterministic_{}'.format(args.k), fname)
  with open(path, 'r') as f:
    lines = f.readlines()
    sort_time = float(lines[7].split()[-1])
    static_alg_time = float(lines[27].split()[-1])
    static_time = sort_time + static_alg_time
    dynamic_pre_time = float(lines[32].split()[-1])
    dynamic_alg_time = float(lines[36].split()[-1])
    dynamic_time = sort_time + dynamic_pre_time + dynamic_alg_time
    auto_pre_time = float(lines[41].split()[-1])
    auto_alg_time = float(lines[45].split()[-1])
    auto_time = sort_time + auto_pre_time + auto_alg_time
    if fname != 'spotify':
      brute_time = float(lines[52].split()[-1])
      static_acc = float(lines[59].split()[-1])
      dynamic_acc = float(lines[67].split()[-1])
      auto_acc = float(lines[75].split()[-1])
    else:
      brute_time = '>86400'
      static_acc = '\\xmark'

    row = ''
    if fname != 'congress-bills':
      row += '& '
    else:
      row += '  '
    if fname != 'spotify':
      row += '{:>25} & {:10.3f} & {:7.3f} & {:7.3f} & {:7.3f} & {:7.3f} \\\\'.format(fname, brute_time, static_time, dynamic_time, auto_time, static_acc)
    else:
      row += '{:>25} & {:>10} & {:7.3f} & {:7.3f} & {:7.3f} & {:7} \\\\'.format(fname, brute_time, static_time, dynamic_time, auto_time, static_acc)

    print(row)
