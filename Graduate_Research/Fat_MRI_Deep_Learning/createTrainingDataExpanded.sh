##########
#
# Need to change the following two directories:
#   * baseDir
#   * utilitiesDir
#
##########

baseDir=/Users/ntustison/Data/unet_fat_mri/
utilitiesDir=/Users/ntustison/Pkg/Utilities/bin/

inputImageDir=${baseDir}/Data/PreprocessedImages/
inputSegmentationDir=${baseDir}/Data/Segmentations/
templateDir=${baseDir}/Data/Template/
outputImageDir=${baseDir}/Data/TrainingDataExpanded/Images/
outputSegmentationDir=${baseDir}/Data/TrainingDataExpanded/Segmentations/

mkdir -p $outputImageDir
mkdir -p $outputSegmentationDir

images=( `ls ${inputImageDir}/*main_N4Denoised.nii.gz` )
segmentations=( `ls ${inputSegmentationDir}/*segmentation.nii.gz` ) 

for (( i=0; i < ${#images[@]}; i++ ));
  do
    imageSource=${images[$i]}
    imageSourceRoot=`basename $imageSource`
    imageSourceRoot=${imageSourceRoot/_main_N4Denoised\.nii\.gz/}

    maskSource=${segmentations[$i]}

    for (( j=0; j < ${#images[@]}; j++ ));
      do
        imageTarget=${images[$j]}
        imageTargetRoot=`basename $imageTarget`
        imageTargetRoot=${imageTargetRoot/_main_N4Denoised\.nii\.gz/}

        echo "Processing ${imageTargetRoot}x${imageSourceRoot}"

        imageSourceWarped=${outputImageDir}/${imageTargetRoot}x${imageSourceRoot}_main_N4Denoised_Rotation0.nii.gz;
        maskSourceWarped=${outputSegmentationDir}/${imageTargetRoot}x${imageSourceRoot}_segmentation_Rotation0.nii.gz;

        if [ ! -f $imageSourceWarped ]
          then 
            antsApplyTransforms -d 2 -i $imageSource \
                                    -r $imageTarget \
                                    -o $imageSourceWarped \
                                    -n Linear \
                                    -t [${templateDir}/T_${imageTargetRoot}_segmentation${j}Affine.txt,1] \
                                    -t ${templateDir}/T_${imageTargetRoot}_segmentation${j}InverseWarp.nii.gz \
                                    -t ${templateDir}/T_${imageSourceRoot}_segmentation${i}Warp.nii.gz \
                                    -t ${templateDir}/T_${imageSourceRoot}_segmentation${i}Affine.txt 
          fi                          

        if [ ! -f $maskSourceWarped ]
          then 
            antsApplyTransforms -d 2 -i $maskSource \
                                    -r $imageTarget \
                                    -o $maskSourceWarped \
                                    -n GenericLabel[Linear] \
                                    -t [${templateDir}/T_${imageTargetRoot}_segmentation${j}Affine.txt,1] \
                                    -t ${templateDir}/T_${imageTargetRoot}_segmentation${j}InverseWarp.nii.gz \
                                    -t ${templateDir}/T_${imageSourceRoot}_segmentation${i}Warp.nii.gz \
                                    -t ${templateDir}/T_${imageSourceRoot}_segmentation${i}Affine.txt 
          fi                          

        imageSourceWarpedFlipped=${imageSourceWarped/N4Denoised/N4Denoised_Flipped}
        if [ ! -f $imageSourceWarpedFlipped ]
          then 
            PermuteFlipImageOrientationAxes 2 ${imageSourceWarped} ${imageSourceWarpedFlipped} 0 1 1 0
            CopyImageHeaderInformation ${imageSourceWarped} ${imageSourceWarpedFlipped} ${imageSourceWarpedFlipped} 1 1 1
          fi  
        
        maskSourceWarpedFlipped=${maskSourceWarped/segmentation/segmentation_Flipped}
        if [ ! -f $maskSourceWarpedFlipped ]
          then 
            PermuteFlipImageOrientationAxes 2 ${maskSourceWarped} ${maskSourceWarpedFlipped} 0 1 1 0
            CopyImageHeaderInformation ${maskSourceWarped} ${maskSourceWarpedFlipped} ${maskSourceWarpedFlipped} 1 1 1
          fi  

        deltaAngle=20
        for(( angle=$deltaAngle; angle < 360; angle+=deltaAngle ));
          do
            echo "  Calculating images for angle ${angle}"

            inputImage=${imageSourceWarped}
            outputImage=${inputImage/Rotation0/Rotation${angle}}
            if [ ! -f $outputImage ]
              then 
                ${utilitiesDir}/RigidTransformImage 2 $inputImage $outputImage ${angle}
              fi  

            inputMask=${maskSourceWarped}
            outputMask=${inputMask/Rotation0/Rotation${angle}}
            if [ ! -f $outputMask ]
              then 
                ${utilitiesDir}/RigidTransformImage 2 $inputMask $outputMask ${angle}
                ThresholdImage 2 $outputMask $outputMask 0.5 10 1 0
              fi  

            inputImage=${imageSourceWarpedFlipped}
            outputImage=${inputImage/Rotation0/Rotation${angle}}
            if [ ! -f $outputImage ]
              then 
                ${utilitiesDir}/RigidTransformImage 2 $inputImage $outputImage ${angle}
              fi  

            inputMask=${maskSourceWarpedFlipped}
            outputMask=${inputMask/Rotation0/Rotation${angle}}
            if [ ! -f $outputMask ]
              then 
                ${utilitiesDir}/RigidTransformImage 2 $inputMask $outputMask ${angle}
                ThresholdImage 2 $outputMask $outputMask 0.5 10 1 0
              fi  
          done
      done  
  done
