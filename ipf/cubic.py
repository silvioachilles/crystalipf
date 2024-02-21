import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors


def symmetry_operations():
    # Applying 24 Crystal Symmetry elements O[432] to the orientation matrix
    Sym_crys = np.zeros((24, 3, 3))
    # 1
    Sym_crys[0, 0, 0] = 1
    Sym_crys[0, 1, 1] = 1
    Sym_crys[0, 2, 2] = 1
    # 2
    Sym_crys[2, 0, 0] = 1
    Sym_crys[2, 1, 1] = -1
    Sym_crys[2, 2, 2] = -1
    # 3
    Sym_crys[6, 0, 0] = -1
    Sym_crys[6, 1, 1] = 1
    Sym_crys[6, 2, 2] = -1
    # 4
    Sym_crys[7, 0, 0] = -1
    Sym_crys[7, 1, 1] = -1
    Sym_crys[7, 2, 2] = 1
    # 5
    Sym_crys[1, 0, 0] = 1
    Sym_crys[1, 1, 2] = -1
    Sym_crys[1, 2, 1] = 1
    # 6
    Sym_crys[9, 0, 1] = -1
    Sym_crys[9, 1, 0] = 1
    Sym_crys[9, 2, 2] = 1
    # 7
    Sym_crys[4, 0, 2] = -1
    Sym_crys[4, 1, 1] = 1
    Sym_crys[4, 2, 0] = 1
    # 8
    Sym_crys[17, 0, 2] = 1
    Sym_crys[17, 1, 0] = 1
    Sym_crys[17, 2, 1] = 1
    # 9
    Sym_crys[16, 0, 1] = 1
    Sym_crys[16, 1, 2] = 1
    Sym_crys[16, 2, 0] = 1
    # 10
    Sym_crys[20, 0, 2] = 1
    Sym_crys[20, 1, 1] = -1
    Sym_crys[20, 2, 0] = 1
    # 11
    Sym_crys[3, 0, 0] = 1
    Sym_crys[3, 1, 2] = 1
    Sym_crys[3, 2, 1] = -1
    # 12
    Sym_crys[5, 0, 2] = 1
    Sym_crys[5, 1, 1] = 1
    Sym_crys[5, 2, 0] = -1
    # 13
    Sym_crys[8, 0, 1] = 1
    Sym_crys[8, 1, 0] = -1
    Sym_crys[8, 2, 2] = 1
    # 14
    Sym_crys[23, 0, 1] = -1
    Sym_crys[23, 1, 0] = -1
    Sym_crys[23, 2, 2] = -1
    # 15
    Sym_crys[21, 0, 0] = -1
    Sym_crys[21, 1, 2] = -1
    Sym_crys[21, 2, 1] = -1
    # 16
    Sym_crys[22, 0, 2] = -1
    Sym_crys[22, 1, 1] = -1
    Sym_crys[22, 2, 0] = -1
    # 17
    Sym_crys[18, 0, 1] = 1
    Sym_crys[18, 1, 0] = 1
    Sym_crys[18, 2, 2] = -1
    # 18
    Sym_crys[14, 0, 1] = 1
    Sym_crys[14, 1, 2] = -1
    Sym_crys[14, 2, 0] = -1
    # 19
    Sym_crys[12, 0, 1] = -1
    Sym_crys[12, 1, 2] = -1
    Sym_crys[12, 2, 0] = 1
    # 20
    Sym_crys[10, 0, 1] = -1
    Sym_crys[10, 1, 2] = 1
    Sym_crys[10, 2, 0] = -1
    # 21
    Sym_crys[13, 0, 2] = -1
    Sym_crys[13, 1, 0] = 1
    Sym_crys[13, 2, 1] = -1
    # 22
    Sym_crys[15, 0, 2] = -1
    Sym_crys[15, 1, 0] = -1
    Sym_crys[15, 2, 1] = 1
    # 23
    Sym_crys[11, 0, 2] = 1
    Sym_crys[11, 1, 0] = -1
    Sym_crys[11, 2, 1] = -1
    # 24
    Sym_crys[19, 0, 0] = -1
    Sym_crys[19, 1, 2] = 1
    Sym_crys[19, 2, 1] = 1

    return Sym_crys


class Math:
    @staticmethod
    def project_spherical_on_plane(theta, phi):
        """
        Stereographic projection of spherical coordinates onto a plane at z = 0.

        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        """
        p_x = np.tan(theta / 2) * np.cos(phi)
        p_y = np.tan(theta / 2) * np.sin(phi)
        return p_x, p_y

    @staticmethod
    def project_cartesian_on_plane(hx, hy, hz):
        """
        Stereographic projection of a cartesian vector onto a plane at z = 0.

        :param hx: x component of the vector.
        :param hy: y component of the vector.
        :param hz: z component of the vector.
        """
        theta, phi = Math.cartesian_to_spherical(hx, hy, hz)
        px, py = Math.project_spherical_on_plane(theta, phi)
        return px, py

    @staticmethod
    def cartesian_to_spherical(h_x, h_y, h_z):
        """
        Converts cartesian coordinates (x, y, z) to spherical coordinates (theta, phi).

        :param h_x: x component of the vector
        :param h_y: y component of the vector
        :param h_z: z component of the vector
        """
        if h_z < 0:
            h_z = -h_z

        if np.round(h_z, 4) == 1.0:
            h_z = 1.0

        theta = np.arccos(h_z)
        phi = np.arctan2(h_y, h_x)
        return theta, phi

    @staticmethod
    def linear_slope(x1, y1, x2, y2):
        """
        Calculates the slope of a linear function from two points.

        :param x1: x-coordinate of the first point.
        :param y1: y-coordinate of the first point.
        :param x2: x-coordinate of the second point.
        :param y2: y-coordinate of the second point.

        :return: Slope of the function.
        """
        return (y2 - y1) / (x2 - x1)

    @staticmethod
    def linear_y_intercept(m, x, y):
        """
        Calculates the y-intercept of a linear function.

        :param m: Slope of the function.
        :param x: x-coordinate of a point on the function.
        :param y: y-coordinate of a point on the function.

        :return: y-intercept of the function.
        """
        return y - m * x

    @staticmethod
    def linear_intersection(m1, m2, b1, b2):
        """
        Calculates the intersection of linear functions.

        :param m1: Slope of the first function.
        :param m2: Slope of the second function.
        :param b1: y-intercept of the first function.
        :param b2: y-intercept of the second function.

        :return: x, y coordinates of the intersection.
        """
        x = (b2 - b1) / (m1 - m2)
        y = m1 * x + b1
        return x, y

    @staticmethod
    def angle_between_vectors(v1, v2):
        """
        Calculates the angle between two vectors.

        :param v1: First vector.
        :param v2: Second vector.

        :return: Angle between the two vectors.
        """
        a_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)
        v1v2 = np.dot(v1, v2)
        phi = np.arccos(v1v2 / (a_norm * v2_norm))
        return phi


class Crystallography:
    @staticmethod
    def reciprocal_vector(
            h, k, l,
            a1=np.array([1, 0, 0], dtype=np.float32),
            a2=np.array([0, 1, 0], dtype=np.float32),
            a3=np.array([0, 0, 1], dtype=np.float32),
    ):
        """
        Calculates the reciprocal vector for a given set of Miller indices and a primitive lattice.

        :param h: Miller index
        :param k: Miller index
        :param l: Miller index
        :param a1: Primitive lattice vector
        :param a2: Primitive lattice vector
        :param a3: Primitive lattice vector

        :return: reciprocal vector
        """
        norm = np.dot(a1, np.cross(a2, a3))

        b1 = np.cross(a2, a3) / norm
        b2 = np.cross(a3, a1) / norm
        b3 = np.cross(a1, a2) / norm

        Ghkl = h * b1 + k * b2 + l * b3
        Ghkl /= np.linalg.norm(Ghkl)
        return Ghkl


class IPF:
    def __init__(self):
        self.hkls = (
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
        )
        self.p001_x, self.p001_y = IPF.hkl_to_ipf_xy(*self.hkls[0])
        self.p101_x, self.p101_y = IPF.hkl_to_ipf_xy(*self.hkls[1])
        self.p111_x, self.p111_y = IPF.hkl_to_ipf_xy(*self.hkls[2])

        # The barycenter of the standard stereographic triangle.
        self.bary_x = np.mean([self.p001_x, self.p101_x, self.p111_x])
        self.bary_y = np.mean([self.p001_y, self.p101_y, self.p111_y])

        self.m_001_111 = Math.linear_slope(self.p001_x, self.p001_y, self.p111_x, self.p111_y)
        self.b_001_111 = Math.linear_y_intercept(self.m_001_111, self.p111_x, self.p111_y)

        curve_xs = []
        curve_ys = []
        hs = np.linspace(1.0, 0.0, 3000)
        for h in hs:
            sample_normal = Crystallography.reciprocal_vector(1, h, 1)
            px, py = Math.project_cartesian_on_plane(*sample_normal)

            curve_xs.append(px)
            curve_ys.append(py)

        self.curve_xs = np.array(curve_xs)
        self.curve_ys = np.array(curve_ys)

        self.m_001_101 = Math.linear_slope(self.p001_x, self.p001_y, self.p101_x, self.p101_y)
        self.b_001_101 = Math.linear_y_intercept(self.m_001_101, self.p101_x, self.p101_y)

        self.eta_p001 = self.calc_phi(self.p001_x, self.p001_y, np.array(([0, 1])))
        self.eta_p101 = self.calc_phi(self.p101_x, self.p101_y, np.array(([0, 1])))
        self.eta_p111 = self.calc_phi(self.p111_x, self.p111_y, np.array(([0, 1])))

        phi = np.array([0, 1])
        self.phi_reference = phi / np.linalg.norm(phi)

        self.symmetry_operations = symmetry_operations()

    def filter_by_xy(self, p_xs, p_ys):
        px_in_sst, py_in_sst = [], []
        for p_x, p_y in zip(p_xs, p_ys):
            if self.is_in_sst(p_x, p_y):
                px_in_sst.append(p_x)
                py_in_sst.append(p_y)
            else:
                self.is_in_sst(p_x, p_y)

        if len(px_in_sst) == 0:
            raise Exception("Grain had multiple contributions before filtering and 0 afterwards. Cant be!")

        if len(px_in_sst) == 1:
            return [px_in_sst[0]], [py_in_sst[0]]

        if len(px_in_sst) > 1:
            px_in_sst_round = np.round(px_in_sst, 3)
            py_in_sst_round = np.round(py_in_sst, 3)
            if len(np.unique(px_in_sst_round)) == 1 and len(np.unique(py_in_sst_round)) == 1:
                return [px_in_sst[0]], [py_in_sst[0]]
            else:
                raise Exception("Grain has multiple contributions to sst")

    @staticmethod
    def hkl_to_ipf_xy(h, k, l):
        """
        Converts a set of Miller indices to the corresponding stereographic projection coordinates.

        :param h: Miller index
        :param k: Miller index
        :param l: Miller index
        """
        plane_normal = Crystallography.reciprocal_vector(h, k, l)
        theta, phi = Math.cartesian_to_spherical(*plane_normal)
        px, py = Math.project_spherical_on_plane(theta, phi)
        return px, py

    def g_h_to_ipf_xy(self, g, h):
        transformed_direction = np.dot(g, h)

        return self.unit_vector_to_ipf_xy(transformed_direction)

    def unit_vector_to_ipf_xy(self, v):
        """
        Determines the inverse pole-figure position for a given orientation matrix g and sample direction h.

        :param g: Orientation matrix
        :param h: Sample direction
        """

        if np.round(np.linalg.norm(v), 3) != 1.0:
            raise Exception("Vector is not a unit vector")

        # Apply symmetry
        n_sym = self.symmetry_operations.shape[0]
        vs_applied_symmetries = np.zeros([n_sym, 3])
        thetas = np.zeros(n_sym)
        phis = np.zeros(n_sym)

        for sym_idx, sym_op in enumerate(self.symmetry_operations):
            v_applied_symmetry = np.dot(sym_op, v)

            # if z is negative, project to upper hemisphere
            if v_applied_symmetry[2] < 0.:
                v_applied_symmetry[2] *= -1

            vs_applied_symmetries[sym_idx] = v_applied_symmetry

            theta, phi = Math.cartesian_to_spherical(*v_applied_symmetry)
            thetas[sym_idx] = theta
            phis[sym_idx] = phi

        pxs, pys = Math.project_spherical_on_plane(thetas, phis)

        sst_xs, sst_ys = self.filter_by_xy(pxs, pys)

        if len(sst_xs) == 0:
            raise Exception("No valid IPF coordinates found.")

        elif len(sst_xs) == 1:
            return sst_xs[0], sst_ys[0]

        else:
            raise Exception("More than one contribution")

    def hkl_intersection(self, px, py, m, b):
        dx = 0.001
        x_temp = px
        y_temp = py
        while True:
            x_temp += dx
            y_temp += m * dx
            dx_array = self.curve_xs - x_temp
            x_temp_idx = np.argmin(np.abs(self.curve_xs - x_temp))
            y_hkl = self.curve_ys[x_temp_idx]

            if y_temp >= y_hkl:
                return x_temp, y_temp

            if x_temp > self.p101_x*1.1:
                # TODO: This condition is entered in rare cases, probably due to bad rounding at some places
                print("Cannot find intersection of sst and linear curve. x: {}\ty: {}, returning x, y from p101".format(px, py))
                return self.p101_x, self.p101_y

    def calc_phi(self, px, py, reference):
        v1_temp = reference
        x_temp = px - self.bary_x
        y_temp = py - self.bary_y
        v2_temp = np.array([x_temp, y_temp])

        phi = Math.angle_between_vectors(v1_temp, v2_temp)
        if px > self.bary_x:
            phi = 2 * np.pi - phi

        return phi

    def calc_color_phi(self, px, py):
        phi = self.calc_phi(px, py, self.phi_reference)
        # phi += np.pi / 4
        if phi > (2 * np.pi):
            phi -= 2 * np.pi
        phi /= (2 * np.pi)
        return phi

    def calc_color_theta(self, px, py, intersect_x, intersect_y):
        dx_temp = intersect_x - self.bary_x
        dy_temp = intersect_y - self.bary_y
        d_full = np.hypot(dx_temp, dy_temp)

        dx_temp = px - self.bary_x
        dy_temp = py - self.bary_y
        d_partial = np.hypot(dx_temp, dy_temp)

        # Normalization still unclear, theta should be in range [0, pi/2]
        # But matplotlib wants [0, 1]
        # theta = d_partial / d_full * (np.pi / 2)
        theta = d_partial / d_full
        return theta

    def h1l_from_pxy(self, px, py):
        """"
        Calculates the hue and lightning for a given position in the inverse pole-figure.

        :param px: x-coordinate of the point in the inverse pole-figure.
        :param py: y-coordinate of the point in the inverse pole-figure.
        """
        m = Math.linear_slope(self.bary_x, self.bary_y, px, py)
        b = Math.linear_y_intercept(m, px, py)

        angle = self.calc_phi(px, py, self.phi_reference)
        angle_round = round(angle, 2)

        eta_p001 = round(self.eta_p001, 2)
        eta_p101 = round(self.eta_p101, 2)
        eta_p111 = round(self.eta_p111, 2)

        if (angle_round >= eta_p001) and (angle_round <= eta_p101):
            intersect_x, intersect_y = Math.linear_intersection(m, self.m_001_101, b, self.b_001_101)
        elif (angle_round > eta_p101) and (angle_round < eta_p111):
            intersect_x, intersect_y = self.hkl_intersection(px, py, m, b)
        else:
            intersect_x, intersect_y = Math.linear_intersection(m, self.m_001_111, b, self.b_001_111)

        theta = self.calc_color_theta(px, py, intersect_x, intersect_y)
        phi = self.calc_color_phi(px, py)
        return phi, 1., theta

    def rgb_from_pxy(self, px, py, check_value_space=True):
        """
        Calculates the RGB color for a given position in the inverse pole-figure.

        :param px: x-coordinate of the point in the inverse pole-figure.
        :param py: y-coordinate of the point in the inverse pole-figure.
        :param check_value_space: If True, the RGB values are checked to be in the range [0, 1].
        """
        h1l = self.h1l_from_pxy(px, py)

        rgb = matplotlib.colors.hsv_to_rgb(h1l)

        if check_value_space:
            if rgb[0] > 1.:
                rgb[0] = 1.
            if rgb[1] > 1.:
                rgb[1] = 1.
            if rgb[2] > 1.:
                rgb[2] = 1.

        return rgb

    def rgb_from_U(self, U, sample_direction):
        """
        Calculates the inverse pole-figure RGB color for a given orientation and sample direction.

        :param U: Orientation matrix.
        :param sample_direction: Sample direction.
        """
        px, py = self.g_h_to_ipf_xy(U, sample_direction)
        rgb = self.rgb_from_pxy(px, py)
        return rgb

    def is_in_sst(self, x, y):
        """
        Checks if a given point is inside the standard stereographic triangle.

        :param x: x-coordinate of the point.
        :param y: y-coordinate of the point.
        """
        y_min = self.p001_y
        if (x >= self.p001_x) and (x <= self.p111_x):
            y_max = self.m_001_111 * x + self.b_001_111
        elif (x >= self.p111_x) and (x <= self.p101_x):
            x_idx = np.argmin(np.abs(self.curve_xs - x))
            y_max = self.curve_ys[x_idx]
        else:
            return False

        if (y >= y_min) and (y <= y_max):
            return True
        else:
            return False

    @staticmethod
    def hsl_to_rgb(hsl):
        """
        Converts a color from HSL to RGB.

        :param hsl: (hue, saturation, lightning).
        """
        h, s, l = hsl[0], hsl[1], hsl[2]
        c = (1 - np.abs(2 * l - 1)) * s
        x = c * (1 - np.abs(h / 60))

        m = l - (c / 2)

        if (h >= 0) and (h < 60):
            rgb_tilde = (c, x, 0)
        elif (h >= 60) and (h < 120):
            rgb_tilde = (x, c, 0)
        elif (h >= 120) and (h < 180):
            rgb_tilde = (0, c, x)
        elif (h >= 180) and (h < 240):
            rgb_tilde = (0, x, c)
        elif (h >= 240) and (h < 300):
            rgb_tilde = (x, 0, c)
        elif (h >= 300) and (h < 360):
            rgb_tilde = (c, 0, x)
        else:
            raise Exception("Unhandled value for h: {}".format(h))

        rgb = (np.array(rgb_tilde) + m)

        return rgb


class Plots:
    def __init__(self):
        self.image, self.x_max, self.y_max = self.get_colored_ipf_img()

    @staticmethod
    def sst_xyz_equal_area():
        points = []
        for h in np.linspace(0.0, 1.0, 50):
            reciprocal = Crystallography.reciprocal_vector(1, h, 1)
            points.append(reciprocal)
        for h in np.linspace(1.0, 0.0, 50):
            reciprocal = Crystallography.reciprocal_vector(h, h, 1)
            points.append(reciprocal)
        for h in np.linspace(0.0, 1.0, 50):
            reciprocal = Crystallography.reciprocal_vector(h, 0, 1)
            points.append(reciprocal)

        return np.array(points)

    @staticmethod
    def sst_xy_equal_angle():
        points = Plots.sst_xyz_equal_area()
        pxs, pys = [], []
        for point in points:
            px, py = Math.project_cartesian_on_plane(*point)
            pxs.append(px)
            pys.append(py)

        return pxs, pys

    @staticmethod
    def plot_points_in_sst(pxs, pys):
        plt.figure()
        sst_xs, sst_ys = Plots.sst_xy_equal_angle()
        plt.plot(sst_xs, sst_ys)
        plt.scatter(pxs, pys)
        plt.show()

    @staticmethod
    def get_colored_ipf_img():
        # print("Creating IPF with colorcode now. This may cause some warnings.")
        x_max = 0.45
        y_max = 0.4
        rows = np.linspace(0.0, y_max, 450)
        cols = np.linspace(0.0, x_max, 450)

        image = np.zeros((len(rows), len(cols), 3))
        image += 1.0
        image_transposed = np.zeros((len(cols), len(rows), 3))
        image_transposed += 1.0

        ipf = IPF()
        for row_idx, row in enumerate(rows):
            for col_idx, col in enumerate(cols):
                if ipf.is_in_sst(col, row):
                    try:
                        h1l = ipf.h1l_from_pxy(col, row)
                        rgb = matplotlib.colors.hsv_to_rgb(h1l)
                        rgb = np.round(rgb, 4)
                        image[-row_idx-1, col_idx, :] = rgb
                    except Exception as ex:
                        print(ex)
                        # image[row_idx, col_idx, :] = [1.0, 0.0, 0.0]

                    image_transposed[-col_idx, row_idx, :] = image[row_idx, col_idx, :]

        # print("IPF color done. Warnings should end here.")
        # return image_transposed, x_max, y_max
        return image, x_max, y_max

    def plot_colored_ipf(self, axis='off', path=None, showimage=False):
        plt.imshow(self.image, extent=(0.0, self.x_max, 0.0, self.y_max))

        sst_xs, sst_ys = Plots.sst_xy_equal_angle()
        plt.plot(sst_xs, sst_ys, c='black', lw=2)

        plt.axis(axis)
        plt.xlim(-0.2, self.x_max+0.2)
        plt.ylim(-0.2, self.y_max+0.2)
        plt.annotate("(0 0 1)", (-0.2, -0.075), fontsize=20)
        plt.annotate("(1 0 1)", (0.42, -0.075), fontsize=20)
        plt.annotate("(1 1 1)", (0.38, 0.4), fontsize=20)

        if path is not None:
            plt.savefig(path)
        if showimage:
            plt.show()

    def plot_colored_ipf_raw(self, axis='on', path=None, showimage=False):
        plt.imshow(self.image, extent=(0.0, self.x_max, 0.0, self.y_max))

        sst_xs, sst_ys = Plots.sst_xy_equal_angle()
        plt.plot(sst_xs, sst_ys, c='black', lw=2)

        plt.axis(axis)

        if path is not None:
            plt.savefig(path)
        if showimage:
            plt.show()
