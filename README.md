# cloude_decom
***NEW: for instructions on new Python version*** please [click here](https://github.com/ashlinrichardson/cloude_decom?tab=readme-ov-file#new-running-the-python-version)! 

Near-real-time, interactive, self-contained, non-proprietary implementation of the **Cloude Decomposition** [1] by [Shane R Cloude](https://scholar.google.ca/citations?hl=en&user=h-ZWMcUAAAAJ&view_op=list_works&sortby=pubdate) for:
* rapid exploration of [quad-polarized](https://www.nrcan.gc.ca/maps-tools-and-publications/satellite-imagery-and-air-photos/satellite-imagery-products/educational-resources/tutorial-radar-polarimetry/polarization-radar-systems/9567) [SAR](https://earthdata.nasa.gov/learn/backgrounders/what-is-sar) data such as from [nisar](https://nisar.jpl.nasa.gov/), [uavsar](https://uavsar.jpl.nasa.gov/), [palsar](https://asf.alaska.edu/data-sets/sar-data-sets/alos-palsar/), or [radarsat2](https://www.asc-csa.gc.ca/eng/satellites/radarsat2/Default.asp) and
* inspiring new PolSAR applications.

Sample data included.

### Poincaré's radar by Dr. Cloude
This July 26 2002 presentation is a companion resource accompanying the article below [1]. 

[![IMAGE ALT TEXT](https://i.ytimg.com/vi/fvjGcp0XKNA/hqdefault.jpg)](https://www.youtube.com/watch?v=fvjGcp0XKNA)

### Cloude_decom demonstration video (over forested area)
[![IMAGE ALT TEXT](http://img.youtube.com/vi/03ddjowiCyI/0.jpg)](http://www.youtube.com/watch?v=03ddjowiCyI)

## Abstract
"In this paper [1] we show how to develop the idea of orthogonality in radar polarimetry for enhanced target detection and land-use classification of POLSAR data. It is well known that every elliptical polarization has its orthogonal partner, uniquely defined as the antipodal point on the Poincaré sphere. Here we show that for scatterers, we can extend this idea so that every scattering matrix has a corresponding multi-dimensional orthogonal space [2]. We give a geometrical interpretation of this approach using a generalized Poincaré sphere representation. We then derive an algorithm for finding the peak signal in this ortho-space. We derive this optimum for both monostatic and bistatic radar systems, to illustrate how bistatic polarimetry offers great potential benefits in future POLSAR studies."
### Biblio
[1] "Generalized Poincaré Orthogonality: A New Approach to POLSAR Data Analysis" SR Cloude, A Richardson, https://arxiv.org/abs/2109.09093

[2] "Target Detection Using Rank-1 Polarimetric Processing" SR Cloude, IEEE Geoscience and Remote Sensing Letters (2020)

## Data
Inputs: A) fully-polarimetric SAR data in standard [PolSARPro](https://ietr-lab.univ-rennes1.fr/polsarpro-bio/) "T3 coherency matrix" format with **config.txt** file B) an RGB encoding such as the well-known Pauli encoding, used for selection of a target location C) image coordinates of user-selected target to be "cancelled". Output: *optimised radar cross section* (a greyscale image)
## Running
The interactive version supports exploring fully-polarimetric SAR data quickly, by selecting a target location using the mouse, then viewing the resulting optimised radar cross section promptly

* Tested on **Ubuntu 20 LTS** and MacOS. Windows is supported via WSL ( Windows Subsystem: Linux )
## Instructions for starting WSL ( on Windows ):
1. Open your windows terminal (AKA "cmd.exe")
2. If you didn't install WSL yet:
```
wsl --install
```
3. Enter the wsl environment:
```
wsl
```

4. Install "unzip" command (needed to open code) 
At the WSL prompt:
```
sudo apt install unzip
```
## Instructions for building and running the interactive method, using the test data provided:
Please open your terminal and run the following commands:

### 1) download the code
(by the way, make sure to use the "copy" icon to the right of the code fragment, followed by pasting in your terminal, in order to avoid re-typing manually). All code fragments require pressing "RETURN" key afterwards to run the command. 
```
curl -o cloude_decom.zip https://codeload.github.com/ashlinrichardson/cloude_decom/zip/refs/heads/master; unzip cloude_decom.zip; mv cloude_decom-master cloude_decom
```

### 2) enter the project folder 
```
cd cloude_decom
```

### 3) compile the code
```
python3 cpp/compile.py
```

### 4) enter the sample data folder
```
cd T3
```

### 5) run the interactive program
```
cloude_view
```
If all goes well, you should see an interactive visualization of the test data (as below). 

Notes:
1) **[compile.py](https://github.com/ashlinrichardson/cloude_decom/blob/master/cpp/compile.py)** installs, provided you enter your super-user password:
* required dependencies g++ and freeglut3-dev (on ubuntu). Mac users need to install **xcode** "command line" development tools
* binaries to **/usr/bin/cloude_decom** and **/usr/bin/cloude_view**

2) Alternate downloading method, you could use this instead of step "1)" above if you're already [connected to github by ssh](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh):

```
git clone git@github.com:ashlinrichardson/cloude_decom.git
```

and then proceed to "2) enter the project folder". This approach would make sense if you're experienced with GitHub and considering making code revisions

3) Another alternative downloading method (assuming you have **wget** installed)
```
wget https://github.com/ashlinrichardson/cloude_decom/archive/refs/heads/master.zip; unzip master.zip; mv cloude_decom-master cloude_decom
```
### Using the mouse
Assuming the mouse pointer is positioned somewhere over the image display:
* **depressing the left button**, restores the default visualization used (e.g. the pauli encoding)
* **releasing the left button, runs the decomposition** and displays the **optimized radar cross section** associated with the target area under the cursor, when the button was released

So it's necessary to engage and then release the left mouse button, to generate an output. Again the location where the mouse button is released, becomes the target area for processing

### Closing the program
Depressing the right mouse botton, or the "esc" key on the keyboard, both close the program 

### Sample output
ALOS PALSAR data over SanFransisco, California in pauli encoding **(r, g, b) = (T22, T33, T11)**:
<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/pauli.png" width="800">

For example, triggering the processing by pointing at water, demonstrates (by "cancelling" water) the ability to pinpoint ships: 

e.g. optimum radar cross section produced by selecting target at row/col index: (x,y)=863,739


<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/opt_cancel_water.png" width="800">

whereas cancelling an urban area highlights oceanographic and topographic features: e.g. optimum radar cross section produced by selecting target at point (x,y) = (688, 663) 

<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/opt_cancel_urban.png" width="800">


### Uninstalling
At your terminal and from within the cloude_decom folder:
```
python3 cpp/uninstall.py
```
### Changing the target window size
Default target window size is 3. To run with a different window size, e.g. 5:
```
cd T3
cloude_view stack.bin 5
```
### NEW: running the Python version
If you already downloaded the code, proceed to step 3)
### 1) download the code
```
curl -o cloude_decom.zip https://codeload.github.com/ashlinrichardson/cloude_decom/zip/refs/heads/master; unzip cloude_decom.zip; mv cloude_decom-master cloude_decom
```

### 2) enter the project folder 
```
cd cloude_decom
```

### 3) make sure python dependencies are installed ( quite standard ):
```
python3 -m pip install numpy matplotlib
```

### 4) run the python code on the sample data provided:
```
python3 py/cloude_decom.py T3
```

This gui works "the same" as the C++ gui:
* Depressing the mouse button will restore the "default" visualization e.g. (r,g,b) = (T22, T33, T11).
* Releasing the mouse button will run the decom ( and generate and display the product opt.bin).

### 5) Polygon / area target
Using the right mouse button draws an area target. Click the left mouse button once to close the polygon

<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/poly_target.gif" width="800">

### 6) Polygon / area target: specified from shapefile + georeferenced image
Here is the sample command for running the decom using a "shapefile" target:

```
python3 py/cloude_decom.py T3 --shapefile=T3/shapefiles/water.shp
```

Note that T11.bin must include georeferencing information in the T11.hdr file. For this option, the GUI is suppressed and outputs ( opt.bin, etc ) are made available in the T3 data folder.

### 7) Extra options
To run the python version at a specific target location (column, row index) examine the top of the python file for examples. For example, the flag ```--special_rgb``` can be used to display a special color-encoding:

<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/rgb_special.png" width="800">

### 8) Extra parameters
Quite a few parameters are output by the program, such as the diagonal elements of the rank-1 T3 matrix:
<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/rank_1_t3.png" width="800">

Comparing the pauli representation with the eigenvalues ( sorted ) and also with the rank-1 matrix elements:
<img src="https://raw.githubusercontent.com/ashlinrichardson/cloude_decom/master/plots/compare.gif" width="800">

### Running the method on other data
* Make sure your dataset is in T3 matrix format ( PolSARPro standard ). ENVI header files ( one for each .bin file) are expected, as is the PolSARPro "config.txt" file.

### Running on T3 data produced by SNAP
## Thanks
Thanks to [Eric Pottier](https://scholar.google.it/citations?hl=en&user=wObZqM0AAAAJ&view_op=list_works&sortby=pubdate) and [JAXA](https://global.jaxa.jp/) for providing ALOS-1 quad-pol data over San Fransisco, California. Please [click here](https://ietr-lab.univ-rennes1.fr/polsarpro-bio/sample_datasets/) to see other sample data generously provided by Dr. Pottier.

Thanks to:
* Hao Chen, Canadian Forest Service
* Subhadip Dey, Indian Institute of Technology Bombay
* Shane Cloude, AELc
* Adithi Balaji (U.Vic)
for initial beta testing

# Contributing
At your convenience, please be welcome to:
* [open an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue) on this repository,
* [submit a pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request), or
* provide feedback by email

# Notes: steps performed on SAN_FRANCISCO_ALOS1.zip:
1. Radar --> Radiometric --> Calibrate (write complex output)
2. Radar --> Geometric --> ALOS Deskewing
3. Radar --> Polarimetric --> Faraday Rotation Correction
4. Radar -> Polarimetrc --> Polarimetric Matrix Generation ( T3 matrix ) 
5. Radar --> Radiometric --> Radiometric Terrain Flattening ( Small et al RTF )
6. Radar --> Utilities --> Multilooking ( to approx square pixel ) 
7. Radar -> Polarimetric --> Speckle filter (Boxcar filter, 7x7) 
8. Radar -> Geometric --> Terrain Correction --> Range Doppler Terrain Correction ( to put in geographic coordinates )
