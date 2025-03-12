// Input directory
dir = getDirectory("Choose a Directory");

// Define channel of interest for segmentation
searchString_BF = "BF Penta";
searchString_channel1 = "RFP Cube";
searchString_channel2 = "A647 Penta";
Media = "2iL"

// Get the list of files in the directory
fileList = getFileList(dir);

// function to manually segment cells from BF image
function manual_segment(dir, searchstring) {
	// Loop through each file in the directory
	for (i = 0; i < fileList.length; i++) {
    
    	// Check if the file name contains the search string and open iterativly through all Brightfiled images
    	if (indexOf(fileList[i], searchstring) >= 0) {
    			open(dir + fileList[i]);
    			setTool("freehand");
    			
    			// manual segmentation and save ROIs
    			waitForUser("Manually select cells (with freehand selection) and add ROIs to ROI manager (press t)");
    			
    			// Define the name of the Segmentations for saving containing only the position information
    			// Define the fixed substrings
				beforeText = Media + "_";
				afterText = "_" + searchString_BF;

				// Find the start and end positions of the desired part
				startPos = indexOf(fileList[i], beforeText) + lengthOf(beforeText);
				endPos = indexOf(fileList[i], afterText);

				// Extract the part between the fixed substrings
				Position = substring(fileList[i], startPos, endPos);
    			
    			// Save Segmentations with Position name
    			roiManager("Save", dir + Position + "_ROIs.zip");
    			close();
    			close("ROI Manager");
    	}
    }
}

// Detect spots in one channel
function detect_spots(dir, searchstring, threshold, maximum) {
	for (i = 0; i < fileList.length; i++) {
		// Check if the file name contains the search string and open iterativly through all Brightfiled images
		if (indexOf(fileList[i], searchstring) >= 0) {
			// Open image
			open(dir + fileList[i]);
			
			// Open Segmentations, named by corresponding Position
    		// Define the fixed substrings
			beforeText = Media + "_";
			afterText = "_" + searchstring;

			// Find the start and end positions of the desired part
			startPos = indexOf(fileList[i], beforeText) + lengthOf(beforeText);
			endPos = indexOf(fileList[i], afterText);

			// Extract the part between the fixed substrings
			Position = substring(fileList[i], startPos, endPos);
    			
			open(dir + Position + "_ROIs.zip"); // Check that this is really the BF segmentation name
			
			// Perform spot detection for each cell
    		cellnumber = roiManager("count");
    		for(j=0; j<cellnumber; j++) {
				roiManager("select", j);
				run("Duplicate...", "duplicate");
				rename("current");
				run("Clear Outside", "stack");
				run("RS-FISH", "image=[current] mode=Advanced anisotropy=1.2000 robust_fitting=[No RANSAC] use_anisotropy spot_intensity=[Linear Interpolation] image_min=0 image_max=" + maximum +" add sigma=0.99400 threshold=0.00910 support=3 min_inlier_ratio=0.10 max_error=1.50 spot_intensity_threshold=" + threshold + " background=[No background subtraction] background_subtraction_max_error=0.05 background_subtraction_min_inlier_ratio=0.10 results_file=[] num_threads=8 block_size_x=128 block_size_y=128 block_size_z=16");
				selectWindow("smFISH localizations");
				saveAs("text", dir + fileList[i] + "_" + j + ".csv");
				close(fileList[i] + "_" + j + ".csv");
				close("current");
    		}
    		close();
		}
		close("ROI Manager");
	}
}

// run Segmentation function
//manual_segment(dir, searchString_BF)

// Rund Spot detection function
// Regarding paramter setting: Use RS-FISh manually on a few images and manually set suitable threshold. For maximum use an approximate for the average maximum intensity per image for each channel.
detect_spots(dir, searchString_channel1, 1000, 28000)
detect_spots(dir, searchString_channel2, 500, 10000)

