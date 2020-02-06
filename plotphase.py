#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for f in sys.argv[1:]:
    data = pd.read_csv(f, header=None, skiprows=2)
    t = data[data.columns[-1]]
    data = data[data.columns[:-1]]
    data = data.groupby(np.arange(len(data))//50).mean()
    print(data.shape)
    figsize = (16, 9)
    plt.figure(figsize=figsize)
    plt.imshow(
            data,
            aspect='auto',
            origin='lower',
            cmap='viridis',
            extent=(0, data.shape[1], t.iloc[0], t.iloc[-1])
            )

plt.show()
