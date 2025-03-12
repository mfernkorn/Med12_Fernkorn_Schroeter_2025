// ImageJ Macro for nuclei segmentation and intensity measurement of three channeled images. 
// These have to be in a single folder and every channel is saved as a sepereated tif
// Naming convention currecntly expects _ch0x.tif with x beeing 0 for nuclei, 1 or 2.
// Measurements are saved for every channel seperatly as a csv with the same name as the corresponding image

// Define Directory of images
dir = getDirectory("Choose Source Directory ");
//dir = "/Users/fernkorn/Documents/Immunostaining_Analysis/20230621_Epi:PrE_wt_vs_Med12KOA11/Images/";
list = getFileList(dir);
run("Set Measurements...", "area mean standard integrated redirect=None decimal=3");

setBatchMode(true)
// Loop through all files (all threee channels in every loop)
for (i=0; i<list.length; i+=4) {
	open(dir+list[i]);
	rename('nuclei');
	
	// Segmentation
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'nuclei', 'modelChoice':'Model (.zip) from File', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.4', 'nmsThresh':'0.4', 'outputType':'Both', 'modelFile':'/Users/fernkorn/Documents/Immunostaining_Analysis/2101_StardistModel_mes.zip', 'nTiles':'1', 'excludeBoundary':'15', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectImage('nuclei');
	roiManager("Deselect");
	
	// Measure and save measurements as csv files with same name as image
	roiManager("multi-measure measure_all");
	saveAs("Results", dir+replace(list[i],".tif",".csv"));
	
	open(dir+list[i+3]);
	roiManager("multi-measure measure_all");
	saveAs("Results", dir+replace(list[i+3],".tif",".csv"));
	
	open(dir+list[i+1]);
	roiManager("multi-measure measure_all");
	saveAs("Results", dir+replace(list[i+1],".tif",".csv"));
	
	open(dir+list[i+2]);
	roiManager("multi-measure measure_all");
	saveAs("Results", dir+replace(list[i+2],".tif",".csv"));

	close();
	close();
	close();
	close();
	close();
}
setBatchMode(false)