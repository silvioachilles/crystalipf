### Let's summarize here what is happening
The purpose of the inverse pole figure (ipf) is to give a neat visualization of one or multiple orientation. 
An orientation must be given as an orientation matrix.
The ipf will always be specific for a certain direction, which is often [0, 0, 1].

The ipf can be caluated by the following:
1. The orientation matrix is multiplied by the sample direction.
2. The resulting vector is normalized and will then point onto the surface of a unit sphere.
3. Due to crystal symmetry, each vector pointing onto the lower surface can be projected onto the upper surface by giving the z component a positive sign.
4. A stereographic projection is applied to the upper half of the unit sphere.
5. Again, due to symmetry, the stereographic projection can be divied into an amount of triangles that contain equal information. The inverse pole figure is exactly one of these triangles. The symmetry operations of the space group are applied to each projected vector, such that all projected vectors lie inside the chosen triangle.
6. A HSL coloring is applied to the triangle.

## Purpose
An ipf can be used to visualize correlations of the orientations of multiple grains. This is often useful when doing grain resolved diffraction, such as 3DXRD.