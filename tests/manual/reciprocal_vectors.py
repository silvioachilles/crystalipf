from ipf.cubic import Crystallography


if __name__ == "__main__":
    # The reciprocal vector G of the 001 Miller index is [0, 0, 1].
    hkl = 0, 0, 1
    G_must = 0, 0, 1
    G_calculated = Crystallography.reciprocal_vector(*hkl)

    assert G_must[0] == G_calculated[0]
    assert G_must[1] == G_calculated[1]
    assert G_must[2] == G_calculated[2]
