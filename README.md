# cloude_decom
Interactive implementation of **Cloude Decomposition** [1] by [Shane R Cloude](https://scholar.google.ca/citations?hl=en&user=h-ZWMcUAAAAJ&view_op=list_works&sortby=pubdate)

[![IMAGE ALT TEXT](http://img.youtube.com/vi/03ddjowiCyI/0.jpg)](http://www.youtube.com/watch?v=03ddjowiCyI "Video Title")

## Abstract
"In this paper [1] we show how to develop the idea of orthogonality in radar polarimetry for enhanced target detection and land-use classification of POLSAR data. It is well known that every elliptical polarization has its orthogonal partner, uniquely defined as the antipodal point on the Poincaré sphere. Here we show that for scatterers, we can extend this idea so that every scattering matrix has a corresponding multi-dimensional orthogonal space [2]. We give a geometrical interpretation of this approach using a generalized Poincaré sphere representation. We then derive an algorithm for finding the peak signal in this ortho-space. We derive this optimum for both monostatic and bistatic radar systems, to illustrate how bistatic polarimetry offers great potential benefits in future POLSAR studies."

## Biblio
[1] **"Generalized Poincaré Orthogonality: A New Approach to POLSAR Data Analysis"** SR Cloude, A Richardson, proceedings of The 7th Asia-Pacific Conference on Synthetic Aperture Radar (2021)

[2] **"Target Detection Using Rank-1 Polarimetric Processing"** SR Cloude, IEEE Geoscience and Remote Sensing Letters (2020)

## Usage
### ubuntu 20 deps:
```sudo apt install g++ freeglut3-dev```

### Compile and run:
**compile.py** installs binaries to **/usr/bin/cloude_decom** and **/usr/bin/cloude_view**
```cd cpp
python3 compile.py
cd ..
cd T3
cloude_view
```
