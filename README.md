# cloude_decom
Near-real-time interactive implementation of **Cloude Decomposition** [1] by [Shane R Cloude](https://scholar.google.ca/citations?hl=en&user=h-ZWMcUAAAAJ&view_op=list_works&sortby=pubdate) supporting rapid exploration of [quad-polarized](https://www.nrcan.gc.ca/maps-tools-and-publications/satellite-imagery-and-air-photos/satellite-imagery-products/educational-resources/tutorial-radar-polarimetry/polarization-radar-systems/9567) [SAR](https://earthdata.nasa.gov/learn/backgrounders/what-is-sar) data such as from [nisar](https://nisar.jpl.nasa.gov/), [uavsar](https://uavsar.jpl.nasa.gov/), [palsar](https://asf.alaska.edu/data-sets/sar-data-sets/alos-palsar/), or [radarsat2](https://www.asc-csa.gc.ca/eng/satellites/radarsat2/Default.asp)

## Demo over forested area
[![IMAGE ALT TEXT](http://img.youtube.com/vi/03ddjowiCyI/0.jpg)](http://www.youtube.com/watch?v=03ddjowiCyI "Video Title")

## Abstract
"In this paper [1] we show how to develop the idea of orthogonality in radar polarimetry for enhanced target detection and land-use classification of POLSAR data. It is well known that every elliptical polarization has its orthogonal partner, uniquely defined as the antipodal point on the Poincaré sphere. Here we show that for scatterers, we can extend this idea so that every scattering matrix has a corresponding multi-dimensional orthogonal space [2]. We give a geometrical interpretation of this approach using a generalized Poincaré sphere representation. We then derive an algorithm for finding the peak signal in this ortho-space. We derive this optimum for both monostatic and bistatic radar systems, to illustrate how bistatic polarimetry offers great potential benefits in future POLSAR studies."
### Biblio
[1] "Generalized Poincaré Orthogonality: A New Approach to POLSAR Data Analysis" SR Cloude, A Richardson, proceedings of The 7th Asia-Pacific Conference on Synthetic Aperture Radar (2021)

[2] "Target Detection Using Rank-1 Polarimetric Processing" SR Cloude, IEEE Geoscience and Remote Sensing Letters (2020)

## Data
Inputs: A) fully-polarimetric SAR data in standard [PolSARPro](https://ietr-lab.univ-rennes1.fr/polsarpro-bio/) "T3 coherency matrix" format with **config.txt** file B) an RGB encoding such as the well-known Pauli encoding, used for selection of a target location C) image coordinates of user-selected target to be "cancelled". Output: *optimised radar cross section* (a greyscale image)
## Running
The interactive version supports exploring fully-polarimetric SAR data quickly, by selecting a target location using the mouse, then viewing the resulting optimised radar cross section promptly

* Tested on **Ubuntu 20 LTS** and MacOS. Windows to be supported soon
### Instructions for building and running the interactive method, using the test data provided:
Please open your terminal and run the following commands:
```
# 1) download the code:
curl -o cloude_decom.zip https://codeload.github.com/ashlinrichardson/cloude_decom/zip/refs/heads/master
unzip cloude_decom.zip
mv cloude_decom-master cloude_decom
cd cloude_decom                                                   # 2) enter the project folder 
python3 cpp/compile.py                                            # 3) build the project codes..
cd T3                                                             # 4) enter the test data folder.. 
cloude_view pauli.bin                                             # 5) run the interactive program! 
```
If all goes well, you should see an interactive visualization of the test data (as below). 

Notes:
1) **[compile.py](https://github.com/ashlinrichardson/cloude_decom/blob/master/cpp/compile.py)** installs, provided you enter your super-user password:
* required dependencies g++ and freeglut3-dev (on ubuntu). Mac users need to install **xcode** "command line" development tools
* binaries to **/usr/bin/cloude_decom** and **/usr/bin/cloude_view**

2) Alternate downloading method, you could use this instead of step "1)" above if you're already [connected to github by ssh](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh):

```
git clone git@github.com:ashlinrichardson/cloude_decom.git        # 1) download the project..
```

and then proceed to "2) enter the project folder". This approach would make sense if you're experienced with GitHub and considering making code revisions

### Using the mouse
Assuming the mouse pointer is positioned somewhere over the image display:
* engaging it, restores the default visualization used (e.g. the pauli encoding)
* **releasing it, runs the decomposition** and displays the **optimized radar cross section** associated with the target area under the cursor, when the button was released

So it's necessary to engage and then release the left mouse button, to generate an output. Again the location where the mouse button is released, becomes the target area for processing

### Sample output
ALOS PALSAR data over SanFransisco, California in pauli encoding **(r, g, b) = (T22, T33, T11)**:
<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/T3/pauli.png" width="800">

For example, triggering the processing by pointing at water, demonstrates (by "cancelling" water) the ability to pinpoint ships: 

<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/T3/opt_cancel_water.png" width="800">

whereas cancelling an urban area highlights oceanographic and topographic features:

<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/T3/opt_cancel_urban.png" width="800">


### Uninstalling
At your terminal and from within the cloude_decom folder:
```
python3 cpp/uninstall.py
```
### Building and running the non-interactive version:
### Changing the target window size
### Running the method on other data
## Thanks
Thanks to [Eric Pottier](https://scholar.google.it/citations?hl=en&user=wObZqM0AAAAJ&view_op=list_works&sortby=pubdate) and [JAXA](https://global.jaxa.jp/) for providing ALOS-1 quad-pol data over San Fransisco, California. Please [click here](https://ietr-lab.univ-rennes1.fr/polsarpro-bio/sample_datasets/) to see other sample data generously provided by Dr. Pottier.
# Contributing
At your convenience, please be welcome to:
* [open an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue) on this repository,
* [submit a pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request), or
* provide feedback by email
