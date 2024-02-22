import numpy as np
import matplotlib.pyplot as plt

from ipf.cubic import image, image_plt, rgb, edges, position


def example_orientation_matrices():
    from scipy.spatial.transform import Rotation

    alphas = np.linspace(0, np.pi / 4.0, 20)
    betas = np.linspace(0, np.pi / 4.0, 20)
    gammas = np.linspace(0, np.pi / 4.0, 20)

    Us = []
    for alpha in alphas:
        for beta in betas:
            for gamma in gammas:
                U = Rotation.from_euler('zxz', [alpha, beta, gamma]).as_matrix()
                Us.append(U)

    return Us


if __name__ == "__main__":
    # image_plt(showimage=True)

    fig, ax = plt.subplots(1, 2)

    ax[0].imshow(image())

    Us = example_orientation_matrices()
    rgbs = []
    pxs, pys = [], []
    h = 0, 0, 1
    for U in Us:
        px, py = position(U, h)
        pxs.append(px)
        pys.append(py)
        rgbs.append(rgb(U))

    ax[1].scatter(pxs, pys, c=rgbs, s=1)
    sst_xs, sst_ys = edges()
    ax[1].plot(sst_xs, sst_ys, color='black')

    plt.show()
