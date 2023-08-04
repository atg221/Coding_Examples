OriginalImageDir=getDirectory("Select the directory");
OriginalImageList=getFileList(OriginalImageDir); 
Origid= getString("Original images prefix?", "I");
for(i=0; i < OriginalImageList.length; i++){
	if (startsWith(OriginalImageList[i], Origid)){
		Orig = OriginalImageDir + OriginalImageList[i];
	// opens the images and stacks them
	open(Orig);
	title = getTitle;
	selectWindow(title);
	setOption("ScaleConversions", true);
	run("8-bit");
		title = getTitle;
        print(title);
        saveAs("tiff", OriginalImageDir+title);
    run("Close All");
	}
}