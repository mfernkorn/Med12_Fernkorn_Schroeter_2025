// Macro to quickly measure defined ROIs in all channels
// Images and ROIs have to be opened
// Specify Output path

run("Set Measurements...", "area mean standard integrated median display redirect=None decimal=3");

run("Split Channels");

list = getList("image.titles");

selectImage(list[0]);
roiManager("Deselect");
roiManager("Measure");

selectImage(list[1]);
roiManager("Measure");

selectImage(list[2]);
roiManager("Measure");

saveAs("Results", "");