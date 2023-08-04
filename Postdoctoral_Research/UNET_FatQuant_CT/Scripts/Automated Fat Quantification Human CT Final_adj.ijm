//Open up DICOM Image and Set Image Properties for the scan
run("Bio-Formats", "display_metadata");
imageTitle=getTitle();
selectWindow(imageTitle);
run("Properties...");
close();
//set the directorites for the original iamges and segmentation images
OriginalImageDir=getDirectory("Select the directory");
SegmentDir=getDirectory("Select the directory"); 
OriginalImageList=getFileList(OriginalImageDir); 
SegmentList=getFileList(SegmentDir); 
// set scan number, first image number, last image number, SAT and VAT prediction image identifiers
//Number = getString("Scan Number?", "");
//first = getString("First image number?","");
//last = getString("Last image number?","");
Origid= getString("Original images suffix?", ".tiff");
SATid = getString("SAT images suffix?", "_Probability1.nii.gz");
VATid = getString("VAT images suffix?", "_Probability2.nii.gz");

//Open + Stack Original Images
// a loop through your original images... (first ... last)
for(i=0; i < OriginalImageList.length; i++){
	if (endsWith(OriginalImageList[i], Origid)){
		Orig = OriginalImageDir + OriginalImageList[i];
	// opens the images and stacks them
	open(Orig);
	}
}
run("Images to Stack", "name=Stack title=[] use");

//Open + Stack SAT Prediction Images
// a loop through your SAT images... (string searching for SATid)
for(i=0; i<SegmentList.length; i++){
	if (endsWith(SegmentList[i], SATid)){
		SAT = SegmentDir + SegmentList[i]; //SegmentDir + Number + "-" + i + SATid;
			// opens the images and stacks them
			run("Bio-Formats (Windowless)", "open=SAT");
	}
}
run("Images to Stack", "name=Stack2 title=[] use");


//Open + Stack VAT Prediction Images
// a loop through your SAT images... (first ... last)
for(i=0; i<SegmentList.length; i++){
	if (endsWith(SegmentList[i], VATid)){
		VAT = SegmentDir + SegmentList[i];
		// opens the images and stacks them
		run("Bio-Formats (Windowless)", "open=VAT");
	}
}
run("Images to Stack", "name=Stack3 title=[] use");

//Quantify SAT By Using Segmentation image to generate a selection on the 
//original image, threshold on given intensity values, and quantify area
selectWindow("Stack2");
setOption("BlackBackground", true);
run("Make Binary", "method=Default background=Dark calculate black");
for (i=1; i<=nSlices; i++) { 
selectWindow("Stack2");
setSlice(i); 
run("Create Selection");
selectWindow("Stack");
setSlice(i);
run("Restore Selection");
run("8-bit");
run("Threshold...");
setThreshold(82, 97);
run("Analyze Particles...", "summarize slice");
run("Restore Selection");
}

//Quantify VAT By Using Segmentation image to generate a selection on the 
//original image, threshold on given intensity values, and quantify area
selectWindow("Stack3");
setOption("BlackBackground", true);
run("Make Binary", "method=Default background=Dark calculate black");
for (i=1; i<=nSlices; i++) { 
selectWindow("Stack3");
setSlice(i); 
run("Create Selection");
run("Make Inverse");
selectWindow("Stack");
setSlice(i);
run("Restore Selection");
run("8-bit");
run("Threshold...");
setThreshold(82, 97);
run("Analyze Particles...", "summarize slice");
run("Restore Selection");
}

//Quantify Nonfat by thresholding the original images such that all tissue
//except fat is selected, then quantifying the area
for (i=1; i<=nSlices; i++) {
selectWindow("Stack");
setSlice(i);
run("8-bit");
run("Threshold...");
setThreshold(97, 255);
run("Analyze Particles...", "summarize slice");
}

//Disable Global Calibration
run("Properties...");

