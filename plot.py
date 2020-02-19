#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from getMetadata import getMetadata

def plttrial(f, figsize=(16,16)):
    md = getMetadata(f)
    data = pd.read_csv(f, header=1)
    data = data.groupby(np.arange(len(data))//50).mean()
    N = data.N0.iloc[0] + data.N1.iloc[0] + data.N2.iloc[0]
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=figsize)
    t = ''
    for k in md:
        t += k+'='+md[k]+'  '
    fig.suptitle(t)

    axs[0].plot(data.t, data.N0/N, 'r-', lw=1, label='N0')
    axs[0].plot(data.t, data.N1/N, 'b-', lw=1, label='N1')
    axs[0].plot(data.t, data.N2/N, 'g-', lw=1, label='N2')
    axs[0].set_ylabel(r'$N_0$ $N_1$ $N_2$', rotation=90, fontsize=20)
    axs[0].legend(loc=0)
    axs[0].set_ylim(-0.1, 1.1)

    axs[1].plot(data.t, data.r, label='r')
    axs[1].plot(data.t, data.psi, label=r'$\psi$')
    axs[1].set_ylabel(r'$r$ $\psi$', rotation=90, fontsize=20)
    axs[1].legend(loc=0)
    axs[1].set_ylim(-0.1, 1.1)

    axs[2].plot(data.t, data.cycles/3, label='num_waves')
    axs[2].legend(loc=0)

if __name__ == "__main__":
    for f in sys.argv[1:]:
        plttrial(f)

        md = getMetadata(f)
        pad =  max([len(s) for s in md])
        pad2 = max([len(s) for s in md.values()])
        for k in md:
            print(f'{k:{pad}} = {md[k]:{pad2}}')
    plt.tight_layout(rect=(0, 0, 1, .95))
    plt.show()
