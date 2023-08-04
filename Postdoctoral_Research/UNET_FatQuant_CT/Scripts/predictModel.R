library( ANTsR )
library( ANTsRNet )
library( keras )
library( abind )
library( ggplot2 )

keras::backend()$clear_session()

baseDirectory <- '../'
dataDirectory <- paste0( baseDirectory, 'Data/' )
testingDirectory <- paste0( dataDirectory, 'TestingData/' )
predictedDirectory <- paste0( testingDirectory, 'PredictedSegmentations/' )
dir.create( predictedDirectory, showWarnings = FALSE, recursive = TRUE )

# classes <- c( "Background", "SubqFat", "AbdominalMuscle", "AbdominalCavity" )
classes <- c( "Background", "SubqFat", "Abdominal" )
segmentationLabels <- 0:( length( classes ) - 1 )

testingImageFiles <- list.files( 
  path = paste0( testingDirectory, '/OriginalImages/Tiff/' ), 
  pattern = "*.tiff", recursive = TRUE, full.names = TRUE )

testingImages <- list()
testingMasks <- list()
testingImageArrays <- list()
testingMaskArrays <- list()

for ( i in 1:length( testingImageFiles ) )
  {
  testingImages[[i]] <- antsImageRead( testingImageFiles[i], dimension = 2 )    

  id <- basename( testingImageFiles[i] ) 
  id <- gsub( ".tiff", '', id )

  testingImageArrays[[i]] <- as.array( testingImages[[i]] )
  }

testingData <- abind( testingImageArrays, along = 3 )  
testingData <- aperm( testingData, c( 3, 1, 2 ) )
testingData <- ( testingData - mean( testingData ) ) / sd( testingData )

X_test <- array( testingData, dim = c( dim( testingData ), 1 ) )

unetModelTest <- createUnetModel2D( c( dim( testingImages[[1]] ), 1 ), 
  numberOfClassificationLabels = length( classes ),
  convolutionKernelSize = c( 5, 5 ), deconvolutionKernelSize = c( 5, 5 ),
  numberOfLayers = 4, numberOfFiltersAtBaseLayer = 32 )
if( file.exists( paste0( baseDirectory, "/ModelWeights/unetWeights.h5" ) ) )
  {
  load_model_weights_hdf5( unetModelTest, 
    filepath = paste0( baseDirectory, "/ModelWeights/unetWeights.h5" ) )
  }
unetModelTest %>% compile( loss = "categorical_crossentropy",
 optimizer = optimizer_adam( lr = 0.0001 ),  
 metrics = c( "acc", multilabel_dice_coefficient ) )

predictedData <- unetModelTest %>% predict( X_test, verbose = 1 )

probabilityImages <- decodeUnet( predictedData, testingImages[[1]] )

for( i in 1:length( probabilityImages ) )
  {
  cat( "Writing probability segmentation images for", testingImageFiles[[i]], "\n" )  
  for( j in 1:length( probabilityImages[[i]] ) )
    {
    imageFileName <- gsub( ".tiff", 
      paste0( "_Probability", segmentationLabels[j], ".nii.gz" ), 
      testingImageFiles[[i]] )
    imageFileName <- 
      gsub( testingDirectory, predictedDirectory, imageFileName )

    dir.create( dirname( imageFileName ), showWarnings = FALSE, recursive = TRUE )

    probabilityArray <- as.array( probabilityImages[[i]][[j]] )

    antsImageWrite( 
      as.antsImage( probabilityArray, reference = testingImages[[i]] ), 
      imageFileName ) 
    }  
  }
