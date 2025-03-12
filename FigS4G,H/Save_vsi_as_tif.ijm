//Define Directory of images to be saved as tif
dir = getDirectory("Choose Source Directory ");
//dir = "/Volumes/abt2/group/agschroeter/Max/Microscopy Data/Clonogenicity_Assays_CellR2/20221013_50hN2B27_ESL/";
list = getFileList(dir);

dir2 = getDirectory("Select Target Directory ");

setBatchMode(true);
// Make an array of files ending ".vsi"
vsilist = newArray(0);
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".vsi")) {
    	vsilist = Array.concat(vsilist,list[i]);
        //vsilist = append(vsilist, list[i]);
    }
}

// Save image as Tiff in target directory as originalname.tif
for (i=0; i<vsilist.length; i++) {
        showProgress(i+1, vsilist.length);
        // open file using Bio-Formats, you may need to edit these two lines
        s = "open=["+dir+vsilist[i]+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";
        run("Bio-Formats Importer", s);

        saveAs("Tiff", dir2+replace(vsilist[i],".vsi",".tif"));
        close();
}

setBatchMode(false);

