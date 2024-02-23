import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from ipf.cubic import position, rgb, image, image_plt, edges

# orientation matrix
grain = Rotation.from_euler(
    'zxz',
    [np.pi / 8.0, np.pi / 8.0, np.pi / 8.0]
).as_matrix()

print(grain)

# sample direction
g = np.array([0, 0, 1])

# ipf position
grain_ipf_x, grain_ipf_y = position(grain, g)

# ipf color
grain_ipf_rgb = rgb(grain, g)

# Visualization of the results
fig, ax = plt.subplots()

# The edges of the ipf
ipf_edges_xs, ipf_edges_ys = edges()
ax.plot(ipf_edges_xs, ipf_edges_ys, c='black')
# The colored position of the grain's orientation
ax.scatter([grain_ipf_x], [grain_ipf_y], c=[grain_ipf_rgb])
ax.axis('off')

plt.savefig('ipf_grain.png')
plt.close()


# A fully colored ipf with annotations
fig, ax = plt.subplots()

image_plt(ax)

plt.savefig('ipf_colored_annotations.png')
plt.close()


# A fully colored ipf
fig, ax = plt.subplots()

ipf_image, xmax, ymax = image()
ax.imshow(ipf_image, extent=(0, xmax, 0, ymax))
ax.plot(ipf_edges_xs, ipf_edges_ys, c='black')
ax.axis('off')

plt.savefig('ipf_colored.png')
plt.close()
