This folder contains my functions to perform tracking and segmentation in PhyloCell.

Useful files:
	trackYeastCells.m:
		the Matlab implementation of my tracking algorithm
	
	extractMapping.m + applyMapping.m + countMappingDifferences.m:
		functions that can be used to compare and evaluate mapping algorithms
	
	exportMontage.m + javitools.jar:
		function to export a sequence of images directly into an MJPEG AVI
	
	phy_segmentCellsWatershedAG.m:
		several segmentation algorithms using watershed
	
	mapCellsUsingYTracker.m + YTracker.jar:
		the Java implementation of my tracking algorithm, but most of it is unusable development code
		the static method ytracker.ContourPredictionTools.borderGeodesicDistances() can be used to compute a geodesic distance transform from an image file

The other files are support files that are necessary to some functions above or that I used during development and testing and may still come in handy.
