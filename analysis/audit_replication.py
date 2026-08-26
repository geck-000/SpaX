# -*- coding: utf-8 -*-
"""Replication audit: for every campaign the paper draws a claim from, how many
independent packings actually stand behind it?

A null result at n=4 and a null result at n=40 are very different claims. This
lists n per condition so the text can be checked against what the data support.
"""
import os, glob, re
import pandas as pd

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')

# campaign -> (file, how to group into conditions)
CAMPAIGNS = [
    ('column ensemble',   'results_colseeds_all.csv',   r'(z\d\d)'),
    ('base size sweep',   'results_basesweep_new.csv',  r'(L\d+)'),
    ('full tensor',       'results_coltensor.csv',      r'(z\d\d)'),
    ('tilt 0',            'results_tilt00.csv',         None),
    ('tilt 15',           'results_tilt15.csv',         None),
    ('tilt 30',           'results_tilt30.csv',         None),
    ('percolation',       'results_perc.csv',           None),
    ('morphology',        'results_morph.csv',          None),
    ('orientation',       'results_orient.csv',         None),
    ('channel geometry',  'results_channel.csv',        None),
    ('size/anisotropy',   'results_sizechan.csv',       r'^(\w+?_L\d+)'),
    ('bending (eringen)', 'results_eringen.csv',        r'_(L\d+)'),
    ('bending control',   'results_eringen_homog.csv',  None),
    ('skeletal sweep',    'results_skeletal.csv',       r'SKEL_([cp]\d+)'),
    ('basal laminae',     'results_skeletal_laminae.csv', r'(z\d+)'),
    ('gas sweep',         'results_gas.csv',            None),
    ('brine K/G',         'results_brine.csv',          None),
    ('seasonal',          'results_seas.csv',           r'SEAS_(w\d+)'),
    ('FY/MY',             'results_fymy.csv',           r'FYMY_(fy|my)'),
    ('salinity family',   'results_salfamily.csv',      r'SAL_(\w+?)_'),
    ('failure sweep',     'results_failure.csv',        r'(z\d\d)'),
    ('weibull SCF',       'results_weibull.csv',        None),
    ('steep column',      'results_steep_column.csv',   r'(z\d+)'),
    ('low base',          'results_lowbase.csv',        r'(z\d\d)'),
    ('cantilever col',    'results_gogo_column.csv',    r'(z\d+)'),
]

print('%-20s %6s %8s %s' % ('campaign', 'rows', 'n/cond', 'note'))
print('-' * 72)
thin = []
for name, f, pat in CAMPAIGNS:
    if not os.path.exists(f):
        print('%-20s %6s %8s %s' % (name, '--', '--', 'MISSING'))
        continue
    d = pd.read_csv(f)
    n = len(d)
    if 'n_packings' in d.columns:
        per = sorted(set(d['n_packings']))
        note = 'replicates recorded in file'
        s = ','.join(str(int(x)) for x in per)
    elif pat:
        g = d['run_id'].astype(str).str.extract(pat)[0]
        counts = g.value_counts()
        s = '%d-%d' % (counts.min(), counts.max())
        note = '%d conditions' % counts.nunique()
        if counts.min() < 5:
            thin.append((name, counts.min()))
    else:
        s = '1'
        note = 'one cell per condition'
        thin.append((name, 1))
    print('%-20s %6d %8s %s' % (name, n, s, note))

print()
print('CAMPAIGNS WITH FEWER THAN 5 REPLICATES PER CONDITION')
print('(a null or a small effect from these is weakly supported)')
for name, k in thin:
    print('   %-20s min n = %d' % (name, k))
