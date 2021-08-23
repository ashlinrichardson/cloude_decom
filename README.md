# cloude_decom
Near-real-time interactive implementation of **Cloude Decomposition** [1] by [Shane R Cloude](https://scholar.google.ca/citations?hl=en&user=h-ZWMcUAAAAJ&view_op=list_works&sortby=pubdate) supporting rapid exploration of fully-polarimetric [SAR](https://earthdata.nasa.gov/learn/backgrounders/what-is-sar) data such as from [nisar](https://nisar.jpl.nasa.gov/), [uavsar](https://uavsar.jpl.nasa.gov/), [alos palsar](https://asf.alaska.edu/data-sets/sar-data-sets/alos-palsar/), or [radarsat2](https://www.asc-csa.gc.ca/eng/satellites/radarsat2/Default.asp)

## Abstract
"In this paper [1] we show how to develop the idea of orthogonality in radar polarimetry for enhanced target detection and land-use classification of POLSAR data. It is well known that every elliptical polarization has its orthogonal partner, uniquely defined as the antipodal point on the Poincaré sphere. Here we show that for scatterers, we can extend this idea so that every scattering matrix has a corresponding multi-dimensional orthogonal space [2]. We give a geometrical interpretation of this approach using a generalized Poincaré sphere representation. We then derive an algorithm for finding the peak signal in this ortho-space. We derive this optimum for both monostatic and bistatic radar systems, to illustrate how bistatic polarimetry offers great potential benefits in future POLSAR studies."
### Biblio
[1] **"Generalized Poincaré Orthogonality: A New Approach to POLSAR Data Analysis"** SR Cloude, A Richardson, proceedings of The 7th Asia-Pacific Conference on Synthetic Aperture Radar (2021)

[2] **"Target Detection Using Rank-1 Polarimetric Processing"** SR Cloude, IEEE Geoscience and Remote Sensing Letters (2020)
## Demo over forested area
[![IMAGE ALT TEXT](http://img.youtube.com/vi/03ddjowiCyI/0.jpg)](http://www.youtube.com/watch?v=03ddjowiCyI "Video Title")

## Data
Inputs: A) fully-polarimetric SAR data in standard [PolSARPro](https://ietr-lab.univ-rennes1.fr/polsarpro-bio/) "T3 coherency matrix" format with **config.txt** file B) an RGB encoding such as the well-known Pauli encoding, used for selection of a target location C) image coordinates of user-selected target to be "cancelled". Output: *optimised radar cross section* (a greyscale image)

The interactive version supports exploring fully-polarimetric SAR data quickly, by selecting a target location using the mouse, then viewing the resulting **optimised radar cross section** promptly

## Running
Tested on **Ubuntu 20 LTS** OS. Windows, MacOS, other Linux to be supported imminently via revised compilation and/or binaries
### Building and running the interactive method, using the test data provided:
```
git clone git@github.com:ashlinrichardson/cloude_decom.git        # 1) download the project..
cd cloude_decom                                                   # 2) enter the project folder 
python3 cpp/compile.py                                            # 3) build the project codes..
cd T3                                                             # 4) enter the test data folder.. 
cloude_view pauli.bin                                             # 5) run the interactive program! 
```
Notes:
1) **compile.py** installs:
* required dependencies g++ and freeglut3-dev
* binaries to **/usr/bin/cloude_decom** and **/usr/bin/cloude_view**

2) Alternate downloading method, please follow this instead of step 1) above if you're not already [connected to github by ssh](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh):
```
wget -O cloude_decom.zip https://github.com/ashlinrichardson/cloude_decom/archive/refs/heads/master.zip
unzip cloude_decom.zip
```

### Building and running the non-interactive version:

### Changing the target window size

### Running the method on other data

## Thanks
Thanks to [Eric Pottier](https://scholar.google.it/citations?hl=en&user=wObZqM0AAAAJ&view_op=list_works&sortby=pubdate) and JAXA for providing ALOS-1 quad-pol data over San Fransisco, California. Please [click here] to see other sample data generously provided by Dr. Pottier.

# Contributing
At your convenience, please be welcome to:
* [open an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue) on this repository,
* [submit a pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request),

or provide feedback by email
