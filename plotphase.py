#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from getMetadata import getMetadata
from tqdm import tqdm

for arg in sys.argv[1:]:
    fname = arg
    maxsize = 100e6
    if not arg.endswith('.dat'):
        if arg.startswith('--size='):
            maxsize = int(float(arg.split('=')[1])*1e6)
        else:
            print(f'Argument `{arg}` not understood, ignoring...')
        continue

    fsize = os.path.getsize(fname)
    md = getMetadata(fname)
    if fsize > maxsize:
        print(f'Large file detected {fsize/1e9:.2f} GB')
        if not any(['--size=' in x for x in sys.argv]):
            tmp = input('Enter max filesize in MB: ')
            maxsize = min(int(float(tmp)*1e6), maxsize)
        d = round(fsize / maxsize)  #Plot 100 MB files at most
        data = []
        t = []
        with open(fname, 'r') as f:
            next(f)
            next(f)
            lines = int(md['iters']) // int(md['log-interval'])
            count = 0
            for i in tqdm(range(lines)):
                if count == d:
                    line = next(f)
                    tmp = [float(x) for x in line.split(',')]
                    data.append(tmp[:-1])
                    t.append(tmp[-1])
                    count = 0
                else:
                    next(f)
                count += 1
    else:
        data = pd.read_csv(fname, header=None, skiprows=2)
        t = data[data.columns[-1]].values
        data = data[data.columns[:-1]].values

    print('File loaded successfully!')
    figsize = (16, 9)
    plt.figure(figsize=figsize)
    plt.imshow(
            data,
            aspect='auto',
            origin='lower',
            cmap='viridis',
            extent=(0, 1, t[0], t[-1])
            )

    md = getMetadata(fname)
    for k in md:
        print(k, md[k])
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title(fr'N={md["N"]}  K={md["K"]}  a={md["a"]}  ic={md["ic"]}  $g \sim N({md["gmean"]},{md["gstddev"]})$')

plt.show()
