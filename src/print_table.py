import os

import numpy as np

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

datasets = [
  'tags-stack-overflow',
  'threads-stack-overflow',
  'Wikipedia-clickstream',
  'Ethereum',
  'Aminer',
  'reddit-reply',
  'MAG',
  'Spotify'
]

l_sort = 7
l_shl_alg = 27
l_dhl_pre = 32
l_dhl_alg = 36
l_auto_pre = 41
l_auto_alg = 45
l_brute = 52
l_shl_acc = 59
spot_shl_acc = 52

edge_sort = 7

for k in [1000, 100000]:
    for (fname, dataset) in zip(files, datasets):
        brute_avg, shl_avg, dhl_avg, auto_avg = 0, 0, 0, 0
        edge_time_avg, edge_acc_avg = 0, 0
        wedge_time_avg = 0

        for run in range(1,11):

            # Deterministic stuff.
            path = os.path.join('../output/compare_deterministic_{}_{}'.format(k, run), fname)
            with open(path, 'r') as f:
                lines = f.readlines()
            sort_time = float(lines[l_sort].split()[-1])
            shl_alg_time = float(lines[l_shl_alg].split()[-1])
            dhl_pre_time = float(lines[l_dhl_pre].split()[-1])
            dhl_alg_time = float(lines[l_dhl_alg].split()[-1])
            auto_pre_time = float(lines[l_auto_pre].split()[-1])
            auto_alg_time = float(lines[l_auto_alg].split()[-1])
            if fname != 'spotify':
                shl_acc = float(lines[l_shl_acc].split()[-1])
                brute_time = float(lines[l_brute].split()[-1])
            else:
                shl_acc = float(lines[spot_shl_acc].split()[-1])
                brute_time = '>86400'
            shl_time = sort_time + shl_alg_time
            dhl_time = sort_time + dhl_pre_time + dhl_alg_time
            auto_time = sort_time + auto_pre_time + auto_alg_time
            if fname != 'spotify':
                brute_avg += brute_time
            else:
                brute_avg = brute_time
            shl_avg += shl_time
            dhl_avg += dhl_time
            auto_avg += auto_time

            # Edge sampling stuff.
            path = os.path.join('../output/compare_parallel_edge_{}_{}'.format(k, run), fname)
            with open(path, 'r') as f:
                lines = f.readlines()
            pre_times, times, accs = [], [], []
            edge_sort_time = float(lines[edge_sort].split()[-1])
            for line in lines:
                if 'Pre-processing' in line:
                    pre_times.append(float(line.split()[-1]))
                if 'Total Time' in line:
                    times.append(float(line.split()[-1]))
                if 'Accuracy' in line:
                    accs.append(float(line.split()[-1]))
            times = times[:-1]  # don't need time taken by adaptive heavy light
            tot_times = [edge_sort_time + p + t for (p,t) in zip(pre_times, times)]

            if k == 1000:
                thresh = 0.98
            else:
                thresh = 0.49
            flag = False
            for (i, (acc, t)) in enumerate(zip(accs, tot_times)):
                if fname == 'wikipedia' and i == 0:
                    continue
                if acc > thresh:
                    edge_acc_avg += acc
                    edge_time_avg += t
                    flag = True
                    break
            if not flag:
                edge_time_avg += times[np.argmax(accs)]

            # Wedge sampling stuff.
            if k == 100000:
                pass
            else:
                if fname not in ['MAG', 'tags-stack-overflow', 'temporal-reddit-reply', 'spotify']:
                    pass
                else:
                    if fname != 'spotify':
                        path = os.path.join('../output/compare_parallel_wedge_{}_{}'.format(k, run), fname)
                        with open(path, 'r') as f:
                            lines = f.readlines()
                        pre_times, times, accs = [], [], []
                        wedge_sort_time = float(lines[edge_sort].split()[-1])
                        for line in lines:
                            if 'Pre-processing' in line:
                                pre_times.append(float(line.split()[-1]))
                            if 'Total Time' in line:
                                times.append(float(line.split()[-1]))
                            if 'Accuracy' in line:
                                accs.append(float(line.split()[-1]))
                        times = times[:-1]  # don't need time taken by adaptive heavy light
                        tot_times = [wedge_sort_time + p + t for (p,t) in zip(pre_times, times)]

                        thresh = 0.49
                        flag = False
                        for (i, (acc, t)) in enumerate(zip(accs, tot_times)):
                            if fname == 'MAG' and i == 0:
                                continue
                            if acc > thresh:
                                wedge_time_avg += t
                                flag = True
                                break
                        if not flag:
                            wedge_time_avg += times[np.argmax(accs)]
                    else:
                        wedge_time_avg = 0

        if fname != 'spotify':
            brute_avg /= 10
        shl_avg /= 10
        dhl_avg /= 10
        auto_avg /= 10
        edge_acc_avg /= 10
        edge_time_avg /= 10
        wedge_time_avg /= 10

        row = '& {:>25}'.format(dataset)
        if fname != 'spotify':
            row += ' & {:10.2f}'.format(brute_avg)
        else:
            row += ' & {:>10}'.format(brute_time)
        row += ' & {:8.2f}'.format(edge_time_avg)
        if k == 100000 or fname not in ['MAG', 'tags-stack-overflow', 'temporal-reddit-reply', 'spotify']:
            if fname != 'spotify':
                tmp = '{:.2f}'.format(brute_avg)
                row += ' & {:>8}'.format('>{}'.format(tmp))
            else:
                row += ' & {:8.2f}'.format(0)
        else:
            row += ' & {:8.2f}'.format(wedge_time_avg)

        row += ' & {:7.2f} & {:7.2f} & {:7.2f} & {:7.2f} \\\\'.format(
            dhl_avg, auto_avg, shl_avg, shl_acc
        )
        print(row)
    print('\n')