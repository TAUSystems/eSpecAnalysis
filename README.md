# eSpecAnalysis 
v2023.10.08.01
modes (NOTE: Currently only Modes 0,1,4,5,6 are implemented.):
	0: Calibrate Perspective Correction and semi-automatic calibration for pixel scaling.
	1: Retrieve pointing from pointing screen and generate electron spectrum for electron
		spectrometer screens.
	2: Retrieve pointing from electron spectrometer screens and generate electron spectrum
		for spectrometer screens.
	3: User Input Mode. Allow to selection of day and shot(s) to analyze.
	4: View pointing screen with physical scaling.
	5: View eScreen A with physical scaling.
	6: View eScreen B with physical scaling.

	Mode 0:
	-Use points defined in perspective.cal to generate a perspective transformation matrix.
	-Use information in perspective.cal to generate physical scale for pixels using a
		reference image.
	-Transformation matrix and physical scales cached into perspective.cache.
	-Reference Images must be located in each screen (Pointing, eScreenA, eScreenB)
		directory for current date.
	-Reference Images must be tiff format (support both .tif, .tiff extension) and be named
		RefImage or RefIMG.
	Mode 1:
	-Loads perspective.cache and use spectrometer (x,y,z,phi,theta) screen coordinates and
		orientation to generate mrad scaling.
	-Loads eScreenX.map and use pointing information from pointing screen to generate
		energy spectrum on eScreenX.
	-Axis used is (z,x,y) where z is propagation and y is vertical axis.
		-Detector screen is assume to have normal pointing in +z direction.
		-Angle theta is rotation about x-axis.
		-Angle phi is rotation about y-axis.
	-Scans eScreenA directory for new file at specified refresh rate. Retrieve file name and
		locate similar file name in eScreenB and Pointing directories.
	-File structure is assumed to be in the form of:
		[dirRoot]\[Date:YYYYMMDD]\[dirDetector]
	-Generate an Analysis Folder in day folder and saves generated analysis image with a 
		similar file name.
	-Blue crosshair on pointing screen specifies peak in region defined by maximum acceptance
		angles.
	-Black grosshair on pointing screen specifies peak in window.
	-Red region on pointing screen specifies maximum acceptance angle.
	-Automatically close after inactivity for specified timeout.
	Mode 2:
	N/A
	Mode 3:
	N/A
	Mode 4:
	-Scans pointing directory for new file at specified refresh rate. Draw image with 5 mm
		axis scale.
	Mode 5:
	-Scans eScreenA directory for new file at specified refresh rate. Draw image with 5 mm
		axis scale.
	Mode 6:
	-Scans eScreenB directory for new file at specified refresh rate. Draw image with 5 mm
		axis scale.
refreshRate (Hz):
	Specify how often to scan for new files in specified directories. Note that this is used
	to calculate the pause time. The true refresh rate will at be less due to compute time.
timeout (s):
	Automatically terminate program after inactivity in all modes except 0.
dirRoot:
	Root of file structure:
		[dirRoot]\[Date:YYYYMMDD]\[dirDetector]
dirPointing:
	folder name for pointing screen data
dirEScreenA:
	folder name for upstream electron spectrometer spectrum screen data
dirEScreenB:
	folder name for downstream electron spectrometer spectrum screen data
x,y,z,phi,theta (mm, deg):
	Associated coordinate and orientation of screens. Source is assumed to be at (0,0,0).
	Information is only used to generate mrad scaling for axis. Energy scaling requires .map
xMaxAngle,yMaxAngle (mrad):
	Maximum acceptance angles for spectrometer.

perspective.cal Notes:
x specifies horizontal pixel location and y specifies vertical pixel location.
xTEPointing, yTEPointing (px):
	pixel value for top point for perspective correction of pointing screen.
xLEPointing, yLEPointing (px):
	pixel value for left point for perspective correction of pointing screen.
xBEPointing, yBEPointing (px):
	pixel value for bottom point for perspective correction of pointing screen.
xREPointing, yREPointing (px):
	pixel value for right point for perspective correction of pointing screen.
xPadEPointing, yPadEPointing (px):
	extend window in x,y by pad amount. Minimum window size is defined by min(x,y), max(x,y)
	of points specified above.
xTLEScreenX, yTLEScreenX (px):
	pixel value for top left point for perspective correction of eScreenX.
xBLEScreenX, yBLEScreenX (px):
	pixel value for bottom left point for perspective correction of eScreenX.
xBREScreenX, yBREScreenX (px):
	pixel value for bottom right point for perspective correction of eScreenX.
xTREScreenX, yTREScreenX (px):
	pixel value for top right point for perspective correction of eScreenX.
xPadEScreenX yPadEScreenX (px):
	extend window in x,y by pad amount. Minimum window size is defined by min(x,y), max(x,y)
	of points specified above. Note that xPad only extends the right edge of the window.
	yPad is used to extend the left edge, along with top and bottom edges.
xContrast, yContrast:
	value used to get more defined edges of the physical scale in reference images.
xThreshold, yThreshold:
	value to specify physical scale tick marks. Above threshold is set to True, and below
	threshold is set to False.
perspective.cache Notes:
	Cache from calibration mode.
	Each screen occupy 4 lines.
		-Window Size in pixels (x,y), and pixel value for zeros (x,y).
		-Homographay matrix row(0), row(1), row(2).
		-Incremental pixel value for each mm tick found by calibration mode for x.
		-Incremental pixel value for each mm tick found by calibration mode for y.
eScreenX.map Notes:
	Precomputed deflection of various electron energies and pointings at spectrometer screen.
	First line specifies energy and pointing axis for the following array.
		-Lowest Energy (MeV), dE (MeV), Start Pointing (mrad), dpointing (mrad)
		-Array of screen positions for electron trajectory of a given energy (row) and
			pointing (column)

## History:
v2023.10.08.01
	-Reworked source code structure to be cmake compatible.
	-Removed windows only dependencies. Code should be unix compatible.
	-Recast some functions to allow for better integration with openCV data structure.
	-Restructure certain functions to improve stability associated with missing files. 

## Authors
This repository contains software developed through a collaboration involving personnel 
affiliated with TAU Systems and The University of Texas at Austin. Individual contributions 
were made by contributors acting in their respective employment or collaborative 
capacities.

This project was developed by Thanh Ha <thanh.ha@tausystems.com> at [TAU Systems Inc.](https://www.tausystems.com) 
in collaboration with the UT3 group at The University of Texas at Austin.
