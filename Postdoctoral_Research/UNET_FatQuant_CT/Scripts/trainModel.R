library( ANTsRNet )
library( keras )
library( ggplot2 )
library( tensorflow )


keras::backend()$clear_session()

baseDirectory <- '../'
scriptsDirectory <- './'
dataDirectory <- paste0( baseDirectory, 'Data/' )

source( paste0( scriptsDirectory, 'unetBatchGenerator.R' ) )

# classes <- c( "Background", "SubqFat", "AbdominalMuscle", "AbdominalCavity" )
classes <- c( "Background", "SubqFat", "Abdominal" )

trainingImageDirectory <- paste0( dataDirectory, 'TrainingData/' )
trainingSegmentationFiles <- list.files( 
  path = paste0( trainingImageDirectory, 'Segmentations/' ), 
  pattern = "segmentation.nii.gz", full.names = TRUE )

trainingTransformDirectory <- paste0( dataDirectory, 'Template/' )

trainingTransforms <- list()
trainingImages <- list()
trainingSegmentations <- list()

pb <- txtProgressBar( min = 0, max = length( trainingSegmentationFiles ), style = 3 )
for( i in 1:length( trainingSegmentationFiles ) )
  {
  segmentation <- antsImageRead( trainingSegmentationFiles[i], dimension = 2 )
  trainingSegmentations[[i]] <- thresholdImage( segmentation, 1, 1, 1, 0 ) + 
    thresholdImage( segmentation, 2, 3, 2, 0 )

  id <- basename( trainingSegmentationFiles[i] ) 
  id <- gsub( "_segmentation.nii.gz", '', id )

  trainingImageFile <- paste0( trainingImageDirectory, "OriginalImages/Nifti/", id, ".tiff_main.nii.gz" )
  trainingImages[[i]] <- antsImageRead( trainingImageFile, dimension = 2 )

  xfrmPrefix <- paste0( trainingTransformDirectory, 'T_', id, '_segmentation', i-1 )

  fwdtransforms <- c()
  fwdtransforms[1] <- paste0( xfrmPrefix, 'Warp.nii.gz' )
  fwdtransforms[2] <- paste0( xfrmPrefix, 'Affine.txt' )
  invtransforms <- c()
  invtransforms[1] <- paste0( xfrmPrefix, 'Affine.txt' )
  invtransforms[2] <- paste0( xfrmPrefix, 'InverseWarp.nii.gz' )

  trainingTransforms[[i]] <- list( 
    fwdtransforms = fwdtransforms, invtransforms = invtransforms )

  setTxtProgressBar( pb, i )  
  }

# unetModel <- createUnetModel2D( c( dim( trainingImages[[1]] ), 1 ), 
#  numberOfClassificationLabels = length( classes ),
#  convolutionKernelSize = c( 3, 3 ), deconvolutionKernelSize = c( 2, 2 ),
#  numberOfLayers = 4, numberOfFiltersAtBaseLayer = 16 )
# if( file.exists( paste0( baseDirectory, "/unetWeights.h5" ) ) )
#  {
#  load_model_weights_hdf5( unetModel, 
#    filepath = paste0( baseDirectory, "/unetWeights.h5" ) )
#  }

# multilabel Dice loss function
# unetModel %>% compile( loss = loss_multilabel_dice_coefficient_error,
#   optimizer = optimizer_adam( lr = 0.00001 ),  
#   metrics = c( multilabel_dice_coefficient ) )

# with( tf$device( "/cpu:0" ), {
  unetModel <- createUnetModel2D( c( dim( trainingImages[[1]] ), 1 ),
    numberOfClassificationLabels = length( classes ),
    convolutionKernelSize = c( 5, 5 ), deconvolutionKernelSize = c( 5, 5 ),
    numberOfLayers = 4, numberOfFiltersAtBaseLayer = 32 )
  # } )

# parallel_unetModel <- multi_gpu_model( unetModel, gpus = 1 )
parallel_unetModel <- unetModel

parallel_unetModel %>% compile( loss = "categorical_crossentropy",
 optimizer = optimizer_adam( lr = 0.0001 ),
 metrics = c( "acc", multilabel_dice_coefficient ) )

###
#
# Set up the training generator
#
batchSize <- 16L

# Split trainingData into "training" and "validation" componets for
# training the model.

numberOfData <- length( trainingImages )

sampleIndices <- sample( numberOfData )

trainingIndices <- sampleIndices[1:round( 0.8 * numberOfData )]
numberOfTrainingData <- length( trainingIndices )
validationIndices <- sampleIndices[length( trainingIndices ):numberOfData]
numberOfValidationData <- length( validationIndices )

trainingData <- unetImageBatchGenerator$new( 
  imageList = trainingImages[trainingIndices], 
  segmentationList = trainingSegmentations[trainingIndices], 
  transformList = trainingTransforms[trainingIndices], 
  referenceImageList = trainingImages, 
  referenceTransformList = trainingTransforms
  )

trainingDataGenerator <- trainingData$generate( batchSize = batchSize )

validationData <- unetImageBatchGenerator$new( 
  imageList = trainingImages[validationIndices], 
  segmentationList = trainingSegmentations[validationIndices], 
  transformList = trainingTransforms[validationIndices],
  referenceImageList = trainingImages, 
  referenceTransformList = trainingTransforms
  )

validationDataGenerator <- validationData$generate( batchSize = batchSize )

###
#
# Run training
#

track <- parallel_unetModel$fit_generator( 
  generator = reticulate::py_iterator( trainingDataGenerator ), 
  steps_per_epoch = ceiling( numberOfTrainingData / batchSize ),
  epochs = 50,
  validation_data = reticulate::py_iterator( validationDataGenerator ),
  validation_steps = ceiling( numberOfValidationData / batchSize ),
  callbacks = list( 
    callback_model_checkpoint( paste0( baseDirectory, "ModelWeights/unetWeights.h5" ), 
      monitor = 'val_loss', save_best_only = TRUE, save_weights_only = TRUE,
      verbose = 1, mode = 'auto', period = 1 ),
     callback_reduce_lr_on_plateau( monitor = 'val_loss', factor = 0.1,
       verbose = 1, patience = 10, mode = 'auto' ),
     callback_early_stopping( monitor = 'val_loss', min_delta = 0.0001, 
       patience = 20 )
    )
  )

save_model_weights_hdf5( unetModel, paste0( baseDirectory, "ModelWeights/unetWeights.h5" ), overwrite = TRUE )

## Plot the model fitting

epochs <- track$epoch

unetModelDataFrame <- data.frame( 
  Epoch = rep( epochs, 2 ), 
  Type = c( rep( 'Training', length( epochs ) ), rep( 'Validation', length( epochs ) ) ),
  Accuracy = c( unlist( track$history$multilabel_dice_coefficient ), 
                unlist( track$history$val_multilabel_dice_coefficient ) )
  )

unetModelAccuracyPlot <- 
  ggplot( data = unetModelDataFrame, aes( x = Epoch, y = Accuracy, colour = Type ) ) +
  geom_point( shape = 1, size = 0.5 ) +
  geom_line( size = 0.3 ) +
  ggtitle( "Accuracy")

ggsave( paste0( baseDirectory, "unetModelAccuracyPlot.pdf" ), 
  plot = unetModelAccuracyPlot, width = 5, height = 2, units = 'in' )
