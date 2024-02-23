import numpy as np

from crystalipf import cubic


class TestCubic:
    ipf = cubic.IPF()
    plots = cubic.Plots()

    precision = 6

    def test_Math_linear_slope(self):
        assert cubic.Math.linear_slope(0, 0, 1, 1) == 1
        assert cubic.Math.linear_slope(0, 0, 1, 2) == 2
        assert cubic.Math.linear_slope(0, 0, 2, 1) == 0.5
        assert cubic.Math.linear_slope(0, 0, 1, -1) == -1

    def test_Math_linear_y_intercept(self):
        assert cubic.Math.linear_y_intercept(1, 0, 0) == 0
        assert cubic.Math.linear_y_intercept(1, 1, 2) == 1
        assert cubic.Math.linear_y_intercept(0.5, 1, 1) == 0.5
        assert cubic.Math.linear_y_intercept(0, 0, 0) == 0

    def test_Math_linear_intersection(self):
        assert cubic.Math.linear_intersection(1, -1, 0, 0) == (0, 0)
        assert cubic.Math.linear_intersection(1, -1, 1, 1) == (0, 1)

    def test_Math_angle_between_vectors(self):
        phi = cubic.Math.angle_between_vectors((0, 1), (1, 0))
        assert np.round(phi, self.precision) == np.round(np.pi / 2.0, self.precision)

        phi = cubic.Math.angle_between_vectors((1, 0), (1, 0))
        assert np.round(phi, self.precision) == 0

        phi = cubic.Math.angle_between_vectors((1, 0), (-1, 0))
        assert np.round(phi, self.precision) == np.round(np.pi, self.precision)

        phi = cubic.Math.angle_between_vectors((1, 0), (1, 1))
        assert np.round(phi, self.precision) == np.round(np.pi / 4.0, self.precision)

    def test_Math_cartesian_to_spherical(self):
        x, y, z = 0, 0, 1
        theta, phi = cubic.Math.cartesian_to_spherical(x, y, z)

        assert theta == 0.0
        assert phi == 0.0

        x, y, z = 1, 0, 0
        theta, phi = cubic.Math.cartesian_to_spherical(x, y, z)

        assert theta == np.pi / 2.0
        assert phi == 0.0

    def test_Math_project_spherical_on_plane(self):
        theta = 0.0
        phi = 0.0
        x, y = cubic.Math.project_spherical_on_plane(theta, phi)

        assert np.round(x, self.precision) == 0.0
        assert np.round(y, self.precision) == 0.0

        theta = np.pi / 2.0
        phi = 0.0
        x, y = cubic.Math.project_spherical_on_plane(theta, phi)

        assert np.round(x, self.precision) == 1.0
        assert np.round(y, self.precision) == 0.0

        theta = np.pi / 2.0
        phi = np.pi / 2.0
        x, y = cubic.Math.project_spherical_on_plane(theta, phi)

        assert np.round(x, self.precision) == 0.0
        assert np.round(y, self.precision) == 1.0

    def test_Math_project_cartesian_on_plane(self):
        cx, cy, cz = 0, 0, 1
        px, py = cubic.Math.project_cartesian_on_plane(cx, cy, cz)

        assert np.round(px, self.precision) == 0.0
        assert np.round(py, self.precision) == 0.0

        cx, cy, cz = 1, 0, 0
        px, py = cubic.Math.project_cartesian_on_plane(cx, cy, cz)

        assert np.round(px, self.precision) == 1.0
        assert np.round(py, self.precision) == 0.0

        cx, cy, cz = 0, 1, 0
        px, py = cubic.Math.project_cartesian_on_plane(cx, cy, cz)

        assert np.round(px, self.precision) == 0.0
        assert np.round(py, self.precision) == 1.0

        cx, cy, cz = 0, -1, 0
        px, py = cubic.Math.project_cartesian_on_plane(cx, cy, cz)

        assert np.round(px, self.precision) == 0.0
        assert np.round(py, self.precision) == -1.0

    def test_Crystallography_reciprocal_vector(self):
        h, k, l = 1, 0, 0
        g = cubic.Crystallography.reciprocal_vector(h, k, l)

        assert np.round(g[0], self.precision) == 1.0
        assert np.round(g[1], self.precision) == 0.0
        assert np.round(g[2], self.precision) == 0.0

        h, k, l = 0, 1, 0
        g = cubic.Crystallography.reciprocal_vector(h, k, l)

        assert np.round(g[0], self.precision) == 0.0
        assert np.round(g[1], self.precision) == 1.0
        assert np.round(g[2], self.precision) == 0.0

        h, k, l = 0, 0, 1
        g = cubic.Crystallography.reciprocal_vector(h, k, l)

        assert np.round(g[0], self.precision) == 0.0
        assert np.round(g[1], self.precision) == 0.0
        assert np.round(g[2], self.precision) == 1.0

        h, k, l = 1, 1, 0
        g = cubic.Crystallography.reciprocal_vector(h, k, l)

        assert np.round(g[0] - 1 / np.sqrt(2), self.precision) == 0.0
        assert np.round(g[1] - 1 / np.sqrt(2), self.precision) == 0.0
        assert np.round(g[2], self.precision) == 0.0

    def test_IPF_filter_by_xy(self):
        # These are all inside the ipf.
        ps_inside = np.array([
            [0, 0],
            [0.29, 0.16],
            [0.15, 0.05],
            [0.39, 0.135],
            [0.4038, 0.0667],
            [0.3767, 0.2989],
        ])

        p_xs_inside_filtered, p_ys_inside_filtered = self.ipf.filter_by_xy(ps_inside[:, 0], ps_inside[:, 1])
        assert len(p_xs_inside_filtered) == len(ps_inside)
        assert len(p_ys_inside_filtered) == len(ps_inside)

        for ip, p in enumerate(ps_inside):
            assert p[0] == p_xs_inside_filtered[ip]
            assert p[1] == p_ys_inside_filtered[ip]

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

        p_xs_outside_filtered, p_ys_outside_filtered = self.ipf.filter_by_xy(ps_outside[:, 0], ps_outside[:, 1])
        assert len(p_xs_outside_filtered) == 0
        assert len(p_ys_outside_filtered) == 0

    def test_IPF_hkl_to_ipf_xy(self):
        h, k, l = 0, 0, 1
        x, y = cubic.IPF.hkl_to_ipf_xy(h, k, l)

        assert np.round(x, self.precision) == 0.0
        assert np.round(y, self.precision) == 0.0

    def test_IPF_g_h_to_ipf_xy(self):
        g = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        h = [0, 0, 1]

        x, y = self.ipf.g_h_to_ipf_xy(g, h)

        assert np.round(x, self.precision) == 0.0
        assert np.round(y, self.precision) == 0.0

    def test_IPF_unit_vector_to_ipf_xy(self):
        g = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        h = [0, 0, 1]

        v = np.dot(g, h)

        x, y = self.ipf.unit_vector_to_ipf_xy(v)

        assert np.round(x, self.precision) == 0.0
        assert np.round(y, self.precision) == 0.0

    def test_IPF_hkl_intersection(self):
        px, py = self.ipf.bary_x + 0.4, self.ipf.bary_y + 0.05

        m = cubic.Math.linear_slope(self.ipf.bary_x, self.ipf.bary_y, px, py)
        b = cubic.Math.linear_y_intercept(m, px, py)

        intersect_x, intersect_y = self.ipf.hkl_intersection(px, py, m, b)

    def test_IPF_bary_azimuthal(self):
        px, py = self.ipf.bary_x, self.ipf.bary_y + 0.5
        azimuthal = self.ipf.bary_azimuthal(px, py)
        assert np.round(azimuthal, self.precision) == 0.0

    def test_IPF_hsl_hue(self):
        px, py = self.ipf.bary_x, self.ipf.bary_y + 0.5
        hue = self.ipf.hsl_hue(px, py)

        assert np.round(hue, self.precision) == 0.0

    def test_IPF_hsl_lightness(self):
        px, py = self.ipf.bary_x, self.ipf.bary_y
        lightness = self.ipf.hsl_lightness(px, py)

        assert np.round(lightness, self.precision) == 0.0

        px, py = self.ipf.p001_x, self.ipf.p001_y
        lightness = self.ipf.hsl_lightness(px, py)

        assert np.round(lightness, self.precision) == 1.0

        px, py = self.ipf.p101_x, self.ipf.p101_y
        lightness = self.ipf.hsl_lightness(px, py)

        assert np.round(lightness, self.precision) == 1.0

        px, py = 0.0, 0.0
        lightness = self.ipf.hsl_lightness(px, py)

        assert np.round(lightness, self.precision) == 1.0

        px, py = self.ipf.p111_x, self.ipf.p111_y
        lightness = self.ipf.hsl_lightness(px, py)

        assert np.round(lightness, self.precision) == 1.0

    def test_IPF_h1l_from_pxy(self):
        px, py = self.ipf.bary_x, self.ipf.bary_y
        h, s, l = self.ipf.h1l_from_pxy(px, py)
        assert np.round(l, self.precision) == 0.0

    def test_IPF_rgb_from_pxy(self):
        px, py = self.ipf.bary_x, self.ipf.bary_y
        r, g, b = self.ipf.rgb_from_pxy(px, py)

        assert np.round(r, self.precision) == 0.0
        assert np.round(g, self.precision) == 0.0
        assert np.round(b, self.precision) == 0.0

    def test_IPF_rgb_from_U(self):
        U = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        g = np.array([0, 0, 1])

        r, g, b = self.ipf.rgb_from_U(U, g)

    def test_IPF_is_in_sst(self):
        ps_inside = np.array([
            [0, 0],
            [0.29, 0.16],
            [0.15, 0.05],
            [0.39, 0.135],
            [0.4038, 0.0667],
            [0.3767, 0.2989],
        ])

        for p in ps_inside:
            assert self.ipf.is_in_sst(p[0], p[1])

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

        for p in ps_outside:
            assert not self.ipf.is_in_sst(p[0], p[1])

    def test_Plots_sst_axes(self):
        cubic.Plots.sst_axes()

    def test_Plots_get_colored_ipf_img(self):
        cubic.Plots.get_colored_ipf_img()

    def test_plot_colored_ipf(self):
        self.plots.plot_colored_ipf()

    def test_symmetry_operations(self):
        symmetry_ops = cubic.symmetry_operations()

    def test_image(self):
        cubic.image()

    def test_image_plt(self):
        cubic.image_plt()

    def test_rgb(self):
        g = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        h = [0, 0, 1]

        rgb = cubic.rgb(g, h)

    def test_position(self):
        g = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])

        h = [0, 0, 1]

        px, py = cubic.position(g, h)

    def test_edges(self):
        sst_xs, sst_ys = cubic.edges()
