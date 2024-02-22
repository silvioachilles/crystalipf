import numpy as np

from ipf import cubic


def test_Math_linear_slope():
    assert cubic.Math.linear_slope(0, 0, 1, 1) == 1
    assert cubic.Math.linear_slope(0, 0, 1, 2) == 2
    assert cubic.Math.linear_slope(0, 0, 2, 1) == 0.5
    assert cubic.Math.linear_slope(0, 0, 1, -1) == -1


def test_Math_linear_y_intercept():
    assert cubic.Math.linear_y_intercept(1, 0, 0) == 0
    assert cubic.Math.linear_y_intercept(1, 1, 2) == 1
    assert cubic.Math.linear_y_intercept(0.5, 1, 1) == 0.5
    assert cubic.Math.linear_y_intercept(0, 0, 0) == 0


def test_Math_linear_intersection():
    assert cubic.Math.linear_intersection(1, -1, 0, 0) == (0, 0)
    assert cubic.Math.linear_intersection(1, -1, 1, 1) == (0, 1)


def test_Math_angle_between_vectors():
    precision = 6

    phi = cubic.Math.angle_between_vectors((0, 1), (1, 0))
    assert np.round(phi, precision) == np.round(np.pi / 2.0, precision)
    phi = cubic.Math.angle_between_vectors((1, 0), (1, 0))
    assert np.round(phi, precision) == 0
    phi = cubic.Math.angle_between_vectors((1, 0), (-1, 0))
    assert np.round(phi, precision) == np.round(np.pi, precision)
    phi = cubic.Math.angle_between_vectors((1, 0), (1, 1))
    assert np.round(phi, precision) == np.round(np.pi / 4.0, precision)


def test_Math_cartesian_to_spherical():
    pass


def test_Math_project_spherical_on_plane():
    pass


def test_Math_project_cartesian_on_plane():
    pass


def test_Crystallography_reciprocal_vector():
    pass


def test_IPF_filter_by_xy():
    ipf = cubic.IPF()

    # These are all inside the ipf.
    ps_inside = np.array([
        [0, 0],
        [0.29, 0.16],
        [0.15, 0.05],
        [0.39, 0.135],
        [0.4038, 0.0667],
        [0.3767, 0.2989],
    ])
    assert ipf.filter_by_xy(ps_inside[:, 0], ps_inside[:, 1]) == (ps_inside[:, 0], ps_inside[:, 1])

    # These are all outside the ipf.
    ps_outside = np.array([
        [0.5, 0.5],
        [0.5, 0],
        [0.6, 0.2],
        [0.1, 0.3],
        [-0.1, 1.0],
        [-0.2, 0.0],
        [0.0, -0.1]
    ])
    assert ipf.filter_by_xy(ps_outside[:, 0], ps_outside[:, 1]) == ([], [])


def test_IPF_hkl_to_ipf_xy():
    pass


def test_IPF_g_h_to_ipf_xy():
    pass


def test_IPF_unit_vector_to_ipf_xy():
    pass


def test_IPF_hkl_intersection():
    pass


def test_IPF_bary_azimuthal():
    pass


def test_IPF_hsl_hue():
    pass


def test_IPF_hsl_lightness():
    pass


def test_IPF_h1l_from_pxy():
    pass


def test_IPF_rgb_from_pxy():
    pass


def test_IPF_rgb_from_U():
    pass


def test_IPF_is_in_sst():
    pass


def test_Plots_sst_axes():
    pass


def test_Plots_get_colored_ipf_img():
    pass


def test_plot_colored_ipf():
    pass


def test_symmetry_operations():
    pass


def test_image():
    pass


def test_image_plt():
    pass


def test_rgb():
    pass


def test_position():
    pass


def test_edges():
    pass
