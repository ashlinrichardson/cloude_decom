# cloude_decom
High-performance interactive implementation of **Cloude Decomposition** [1] by [Shane R Cloude](https://scholar.google.ca/citations?hl=en&user=h-ZWMcUAAAAJ&view_op=list_works&sortby=pubdate). 

### Abstract
"In this paper [1] we show how to develop the idea of orthogonality in radar polarimetry for enhanced target detection and land-use classification of POLSAR data. It is well known that every elliptical polarization has its orthogonal partner, uniquely defined as the antipodal point on the Poincaré sphere. Here we show that for scatterers, we can extend this idea so that every scattering matrix has a corresponding multi-dimensional orthogonal space [2]. We give a geometrical interpretation of this approach using a generalized Poincaré sphere representation. We then derive an algorithm for finding the peak signal in this ortho-space. We derive this optimum for both monostatic and bistatic radar systems, to illustrate how bistatic polarimetry offers great potential benefits in future POLSAR studies."

[1] **"Generalized Poincaré Orthogonality: A New Approach to POLSAR Data Analysis"** SR Cloude, A Richardson, proceedings of The 7th Asia-Pacific Conference on Synthetic Aperture Radar (2021)

[2] **"Target Detection Using Rank-1 Polarimetric Processing"** SR Cloude, IEEE Geoscience and Remote Sensing Letters (2020)
### Demo over forested area
[![IMAGE ALT TEXT](http://img.youtube.com/vi/03ddjowiCyI/0.jpg)](http://www.youtube.com/watch?v=03ddjowiCyI "Video Title")

## How to run
Tested on **Ubuntu 20 LTS** OS. Windows, MacOS, other Linux to be supported imminently
### Building and running the interactive version on the provided test data:
```
git clone git@github.com:ashlinrichardson/cloude_decom.git
cd cloude_decom
python3 cpp/compile.py
cd T3
cloude_view pauli.bin
```
Note: **compile.py** installs:
* dependencies g++ and freeglut3-dev
* binaries to **/usr/bin/cloude_decom** and **/usr/bin/cloude_view**

### Building and running the non-interactive version:

### Changing the target window size

### Acknowledgements
Thanks to [Eric Pottier](https://scholar.google.it/citations?hl=en&user=wObZqM0AAAAJ&view_op=list_works&sortby=pubdate) and JAXA for providing ALOS-1 quad-pol data over San Fransisco, California.
