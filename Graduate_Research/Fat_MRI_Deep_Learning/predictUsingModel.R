library( ANTsR )
library( keras )
library( abind )
library( ggplot2 )

setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr9/Chr9 Congenics/B6 vs Chr9 Congenic Fat MRI/Unet Fat MRI/')
baseDirectory <- './'
dataDirectory <- paste0( baseDirectory, 'Data/' )
modelDirectory <- '/Users/owc/Downloads/Programs/ANTsRNet-master/Models/'
trainingDirectory <- paste0( dataDirectory, 'TestData/Water-Filtered/7853-LC_WHITE/' )

numberOfLabels <- 2

predictedDirectory <- paste0( dataDirectory, 'TestPredictedData_New/Water-Filtered/7853-LC_WHITE_All_Tail_Added_2/' )
dir.create( predictedDirectory, showWarnings = FALSE )

source( paste0( modelDirectory, 'createUnetModel.R' ) )

# I know we're predicting from the same cohort we're training
# but we can easily change this when we get the actual prediction
# image data

predictionImageFiles <- list.files( path = trainingDirectory, 
  pattern = "main", full.names = TRUE )

predictionImage <- antsImageRead( predictionImageFiles[1], dimension = 2 )   
predictionData <- array( dim = c( length( predictionImageFiles ), dim( predictionImage ), 1 ) )

cat( "\n\nReading prediction data (n = ", length( predictionImageFiles ), ")\n", sep = '' )
pb <- txtProgressBar( min = 0, max = length( predictionImageFiles ), style = 3 )
for ( i in 1:length( predictionImageFiles ) )
  {
  predictionImage <- antsImageRead( predictionImageFiles[i], dimension = 2 )   
  predictionArray <- as.array( predictionImage )
  predictionArray <- ( predictionArray - mean( predictionArray ) ) / sd( predictionArray )
  predictionData[i,,,1] <- as.array( predictionImage )  

  setTxtProgressBar( pb, i )
  }
cat( "\n") 

X_prediction <- predictionData

cat( "\nCreating Unet model from existing weights and doing prediction\n", sep = '' )
unetModelPrediction <- createUnetModel2D( c( dim( predictionImage ), 1 ), 
  numberOfClassificationLabels = numberOfLabels, layers = 1:3 )
load_model_weights_hdf5( unetModelPrediction, 
  filepath = paste0( baseDirectory, 'Models/Tail_Added_2/unetWeights.h5' ) )
unetModelPrediction %>% compile( loss = loss_multilabel_dice_coefficient_error,
  optimizer = optimizer_adam( lr = 0.0001 ),  
  metrics = c( multilabel_dice_coefficient ) )

predictedData <- unetModelPrediction %>% predict( X_prediction, verbose = 1 )

cat( "\nConverting prediction data (n = ", length( predictionImageFiles ), ")\n", sep = '' )
pb <- txtProgressBar( min = 0, max = length( predictionImageFiles ), style = 3 )
for( i in 1:length( predictionImageFiles ) )
  {
  for( j in 1:numberOfLabels )
    {
    imageArray <- predictedData[i,,,j]  
    image <- as.antsImage( imageArray, reference = predictionImage )

    imageFileName <- gsub( ".nii.gz", paste0( "_Probability", j, ".nii.gz" ), 
      predictionImageFiles[[i]] )
    imageFileName <- gsub( trainingDirectory, predictedDirectory, imageFileName )
    imageFileName <- gsub( 'Images/', '', imageFileName )

    antsImageWrite( image, imageFileName ) 
    }  
  setTxtProgressBar( pb, i )
  }
cat( "\n\n\n") 

