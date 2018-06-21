import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def read_csvs():
    chars = ['u', 'v', 'p']
    dfs = []
    for char in chars:
        df = pd.read_csv(f"../log/log{char}.csv")
        dfs.append(df)
    return dfs


if __name__ == '__main__':
    dfu, dfv, dfp = read_csvs()
    print(dfu.shape)
