import numpy as np

from crystalipf import cubic
from scipy.spatial.transform import Rotation

import matplotlib.pyplot as plt


def orientation_matrix_from_euler(alpha, beta, gamma):
    """
    This function returns the orientation matrix from the euler angles.
    :param alpha: float
    :param beta: float
    :param gamma: float
    :return: numpy.ndarray
    """
    return Rotation.from_euler('zxz', [alpha, beta, gamma]).as_matrix()


if __name__ == "__main__":
    sample_direction = 0, 0, 1

    # Let's define the euler angles.
    # alpha and gamma can be in the range [0, 2 * pi].
    # beta can be in the range [0, pi].
    n_alpha = 50
    n_beta = 50
    n_gamma = 50
    alphas = np.linspace(0, 2 * np.pi, n_alpha)
    betas = np.linspace(0, np.pi, n_beta)
    gammas = np.linspace(0, 2 * np.pi, n_gamma)

    rotated_vectors = []
    for alpha in alphas:
        for beta in betas:
            for gamma in gammas:
                # We calculate the orientation matrix from the euler angles.
                orientation_matrix = orientation_matrix_from_euler(alpha, beta, gamma)

                rotated_vector = np.dot(orientation_matrix, sample_direction)
                rotated_vectors.append(rotated_vector)

    rotated_vectors = np.array(rotated_vectors)

    # Allows to check if the rotated vectors are reasonable.
    # If the rotated vectors are reasonable, the plot should show a spherish object.
    check_rotated_vector = False
    if check_rotated_vector:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(rotated_vectors[:, 0], rotated_vectors[:, 1], rotated_vectors[:, 2])
        plt.show()
        plt.close()

    # Now we calculate the position of the rotated vectors in the standard stereographic triangle.
    cubic_ipf = cubic.IPF()

    sst_positions = []
    for idx, rotated_vector in enumerate(rotated_vectors):
        print("Doing vector {} of {}: {}".format(idx + 1, len(rotated_vectors), rotated_vector))

        if idx == 408:
            a = 5
            print("BP")

        sst_position = cubic_ipf.unit_vector_to_ipf_xy(rotated_vector)
        sst_positions.append(sst_position)

    sst_positions = np.array(sst_positions)
    plt.scatter(sst_positions[:, 0], sst_positions[:, 1])
    sst_xs, sst_ys = cubic.edges()
    plt.plot(sst_xs, sst_ys, c='black', lw=2)

    plt.show()
