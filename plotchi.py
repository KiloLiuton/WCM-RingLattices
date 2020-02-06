#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for f in sys.argv[1:]:
    data = pd.read_csv(f, header=1)
    #data = data.groupby(np.arange(len(data))//50).mean()
    print(data.head())
    print(data.shape)
    figsize = (16, 9)
    plt.figure(figsize=figsize)

    plt.subplot(311)
    plt.plot(data.a, data.chir, 'r-', lw=.4, label=r'$\chi_r$')
    plt.yscale('log')
    plt.legend(loc='best')

    plt.subplot(312)
    plt.plot(data.a, data.chipsi, 'r-', lw=.4, label=r'$\chi_\psi$')
    plt.yscale('log')
    plt.legend(loc='best')

    plt.subplot(313)
    plt.plot(data.a, data.r, 'r-', lw=.4, label=r'$r$')
    plt.legend(loc='best')

plt.show()
