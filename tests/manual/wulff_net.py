import numpy as np
import matplotlib.pyplot as plt

from crystalipf import cubic


if __name__ == "__main__":
    hs = ks = ls = np.arange(-10, 11)

    pxs, pys = [], []

    for h in hs:
        for k in ks:
            for l in ls:
                px, py = cubic.IPF.hkl_to_ipf_xy(h, k, l)

                pxs.append(px)
                pys.append(py)

    plt.scatter(pxs, pys, s=1)
    plt.show()
