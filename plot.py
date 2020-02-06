#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for f in sys.argv[1:]:
    data = pd.read_csv(f, header=1)
    data = data.groupby(np.arange(len(data))//50).mean()
    print(data.head())
    print(data.shape)
    N = data.N0.iloc[0] + data.N1.iloc[0] + data.N2.iloc[0]
    figsize = (16, 9)
    plt.figure(figsize=figsize)
    plt.plot(data.t, data.N0/N, 'r-', lw=.4, label='N0')
    plt.plot(data.t, data.N1/N, 'b-', lw=.4, label='N1')
    plt.plot(data.t, data.N2/N, 'g-', lw=.4, label='N2')
    plt.legend(loc=0)
    plt.ylim(0, 1)

plt.show()
