# crystalipf
This is an easy-to-use library to create inverse
pole-figures (ipf) for crystal lattices with
cubic symmetry.

## Installation
The library can be installed in three different ways:
1. Install it from the Python Package Index.
    ```bash
    pip install crystalipf
    ```
2. Clone the repository and install it.
    ```bash
    git clone https://github.com/silvioachilles/crystalipf.git
    cd crystalipf
    pip install .
    ``` 
3. Copy `crystalipf/cubic.py` to your project and include it in your code.


## Usage
An orientation must be expressed as an orientation matrix.
In the code snippet below it is shown how to retrieve
the rgb color of the orientation and its position in the ipf.

```python
import numpy as np
import matplotlib.pyplot as plt

from crystalipf.cubic import position, rgb, image, image_plt, edges

# example orientation matrix of a grain
grain = np.array([
    [ 0.71825437, -0.68019413,  0.14644661],
    [ 0.68019413,  0.6421339,  -0.35355339],
    [ 0.14644661,  0.35355339,  0.92387953]
])

# sample direction
g = np.array([0, 0, 1])

# ipf position
grain_ipf_x, grain_ipf_y = position(grain, g)

# ipf color
grain_ipf_rgb = rgb(grain, g)

# to illustrate, the grain's orientation is plotted
# in the ipf with respective coloring
fig, ax = plt.subplots()
ax.scatter([grain_ipf_x], [grain_ipf_y], c=[grain_ipf_rgb])

# the edges of the ipf
ipf_edges_xs, ipf_edges_ys = edges()
ax.plot(ipf_edges_xs, ipf_edges_ys)

ax.axis('off')
plt.show()
```
<img src="https://github.com/silvioachilles/crystalipf/blob/main/doc/ipf_grain.png?raw=true" width="300" height="200" />

Using the `image_plt` function, you can create a matplotlib plot of a fully
colored ipf with annotations for the corners.

```python
fig, ax = plt.subplots()
image_plt(ax)
plt.show()
```

<img src="https://github.com/silvioachilles/crystalipf/blob/main/doc/ipf_colored_annotations.png?raw=true" width="300" height="200" />

An raw image of a fully colored ipf can be retrieved with the `image` function,
allowing to create a customized ipf plot.

```python
fig, ax = plt.subplots()

ipf_image, xmax, ymax = image()
ax.imshow(ipf_image, extent=0, xmax, 0, ymax)
ax.plot(ipf_edges_xs, ipf_edges_ys, c='black')

plt.show()
```

<img src="https://github.com/silvioachilles/crystalipf/blob/main/doc/ipf_colored.png?raw=true" width="300" height="200" />

## Troubleshooting
Please open a github issue if anything appears not correct. 
