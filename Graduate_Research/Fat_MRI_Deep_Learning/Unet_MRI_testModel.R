#load required packages
library( ANTsR )
library( keras )
library( abind )
library( ggplot2 )

keras::backend()$clear_session()

#set directories where data, ANTsRNET models/scripts, and testing images are
setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Rcn2/Rcn2 Fat MRI/Female/Rcn2 KO/')
baseDirectory <- './'
dataDirectory <- paste0( baseDirectory, '8013RP/' )
modelDirectory <- '/Users/owc/Downloads/Programs/ANTsRNet-master/Models/'

testingDirectory <- paste0( dataDirectory, '/Axial - Water-Filtered/' )
testingImagesDirectory <- '/Users/owc/Documents/UVA Grad Labs/Weibin/Rcn2/Rcn2 Fat MRI/Female/Rcn2 KO/8013RP/Axial - Water-Filtered/'
#set where predicted images are going
predictedDirectory <- '/Users/owc/Documents/UVA Grad Labs/Weibin/Rcn2/Rcn2 Fat MRI/Female/Rcn2 KO/8013RP/Automatic_Segmentations/Axial_Water-Filtered/'
if( !dir.exists( predictedDirectory ) )
  {
  dir.create( predictedDirectory, recursive = TRUE )
  }

#load ANTsR scripts
source( paste0( modelDirectory, 'createUnetModel.R' ) )
source( paste0( modelDirectory, 'unetUtilities.R' ) )

#for unprocessed Images (Use whatever identifier you used that is shared among all images)
testingImageFiles <- list.files( path = testingImagesDirectory, 
   pattern = "main", recursive = TRUE, full.names = TRUE )

#for processed images (Use whatever identifier you used that is shared among all images; used N4 for bias corrected/de-noised images)
testingImageFiles <- list.files( path = testingImagesDirectory, 
  pattern = "_N4", recursive = TRUE, full.names = TRUE )

#run for all images (predictions) once you have set up 'testingImageFiles'
testingImage <- antsImageRead( testingImageFiles[1], dimension = 2 )    
imageDimension <- dim( testingImage )

testingImageArrays <- list()

batchSize <- length( testingImageFiles )
testingDataArray <- array( NA, dim = c( batchSize, imageDimension, 1 ) )

for ( i in 1:batchSize )
  {
  testingImage <- antsImageRead( testingImageFiles[i], dimension = 2 ) 
  testingArray <- as.array( testingImage )
  testingArray <- ( testingArray - mean( as.vector( testingArray ) ) ) / 
    sd( as.vector( testingArray ) )
  testingDataArray[i,,,1] <- testingArray
  }
X_test <- testingDataArray  

segmentationLabels <- c( 0, 1 )
numberOfLabels <- length( segmentationLabels )

unetModelTest <- createUnetModel2D( c( imageDimension, 1 ), 
  numberOfClassificationLabels = numberOfLabels, layers = 1:4 )
load_model_weights_hdf5( unetModelTest, 
  filepath = paste0( baseDirectory, 'unetWeights.h5' ) )
unetModelTest %>% compile( loss = loss_multilabel_dice_coefficient_error,
  optimizer = optimizer_adam( lr = 0.0001 ),  
  metrics = c( multilabel_dice_coefficient ) )

predictedData <- unetModelTest %>% predict( X_test, verbose = 1 )

probabilityImages <- decodeY( predictedData, testingImage )

for( i in 1:length( probabilityImages ) )
  {
  for( j in 1:length( probabilityImages[[i]] ) )
    {
    imageFileName <- gsub( ".nii.gz", 
      paste0( "_Probability", segmentationLabels[j], ".nii.gz" ), 
      testingImageFiles[[i]] )
    imageFileName <- 
      gsub( testingImagesDirectory, predictedDirectory, imageFileName )
    outputDirectory <- dirname( imageFileName )

    if( !dir.exists( outputDirectory ) )
      {
      dir.create( outputDirectory, recursive = TRUE )
      }
    antsImageWrite( probabilityImages[[i]][[j]], imageFileName ) 
    }  
  }
