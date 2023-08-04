library( ANTsR )
library( keras )
library( abind )
library( magick )

keras::backend()$clear_session()

baseDirectory <- '../'
dataDirectory <- paste0( baseDirectory, 'Data/' )
trainingTransformDirectory <- paste0( dataDirectory, 'Template/' )

modelDirectory <- '/Users/ntustison/Pkg/ANTsRNet/Models/'
source( paste0( modelDirectory, 'createUnetModel.R' ) )
source( paste0( modelDirectory, 'unetUtilities.R' ) )
source( paste0( baseDirectory, '/Scripts/unetBatchGenerator.R' ) )

trainingImageFiles <- list.files( 
  path = paste0( dataDirectory, 'TrainingImages' ), 
  pattern = "main_N4Denoised", full.names = TRUE )
trainingMaskFiles <- list.files( 
  path = paste0( dataDirectory, 'Segmentations' ), 
  pattern = "segmentation", full.names = TRUE )

trainingTransforms <- list()
trainingImages <- list()
trainingSegmentations <- list()

pb <- txtProgressBar( min = 0, max = length( trainingImageFiles ), style = 3 )
for( i in 1:length( trainingImageFiles ) )
  {
  trainingImages[[i]] <- antsImageRead( trainingImageFiles[i], dimension = 2 )
  trainingSegmentations[[i]] <- antsImageRead( trainingMaskFiles[i], dimension = 2 )

  id <- basename( trainingMaskFiles[i] ) 
  id <- gsub( "_segmentation.nii.gz", '', id )

  xfrmPrefix <- paste0( trainingTransformDirectory, 'T_', id, '_segmentation', i - 1 )

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

unetModel <- createUnetModel2D( c( dim( trainingImages[[1]] ), 1 ), 
  numberOfClassificationLabels = 2, layers = 1:4 )

unetModel %>% compile( loss = loss_multilabel_dice_coefficient_error,
  optimizer = optimizer_adam( lr = 0.0001 ),  
  metrics = c( multilabel_dice_coefficient ) )

###
#
# Set up the training generator
#
batchSize <- 100

# Split trainingData into "training" and "validation" componets for
# training the model.

numberOfTrainingData <- length( trainingImageFiles )

sampleIndices <- sample( numberOfTrainingData )

validationSplit <- 30
trainingIndices <- sampleIndices[1:validationSplit]
validationIndices <- sampleIndices[( validationSplit + 1 ):numberOfTrainingData]

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

track <- unetModel$fit_generator( 
  generator = reticulate::py_iterator( trainingDataGenerator ), 
#  steps_per_epoch = ceiling( 400 / batchSize ),
  steps_per_epoch = 30,
  epochs = 100,
  validation_data = reticulate::py_iterator( validationDataGenerator ),
  validation_steps = ceiling( 100 / batchSize ),
  callbacks = list( 
    callback_model_checkpoint( paste0( baseDirectory, "unetWeights.h5" ), 
      monitor = 'val_loss', save_best_only = TRUE, save_weights_only = TRUE,
      verbose = 1, mode = 'auto', period = 1 )
    # callback_early_stopping( monitor = 'val_loss', min_delta = 0.001, 
    #   patience = 10 ),
    # callback_reduce_lr_on_plateau( monitor = 'val_loss', factor = 0.5,
    #   patience = 0, epsilon = 0.001, cooldown = 0 )
                  # callback_early_stopping( patience = 2, monitor = 'loss' ),
    )
  )







