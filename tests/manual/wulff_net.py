import numpy as np
import matplotlib.pyplot as plt

from crystalipf import cubic


if __name__ == "__main__":
    hs = ks = ls = np.arange(-10, 11)

    pxs, pys = [], []

    for h in hs:
        for k in ks:
            for l in ls:
                if h == k == l == 0:
                    continue

                px, py = cubic.IPF.hkl_to_ipf_xy(h, k, l)

                pxs.append(px)
                pys.append(py)

    plt.scatter(pxs, pys, s=1)

    sst_xs, sst_ys = cubic.edges()
    ipf = cubic.IPF()
    plt.plot(sst_xs, sst_ys, c='black', lw=2)
    plt.scatter(ipf.bary_x, ipf.bary_y, c='red', s=50)

    plt.show()
