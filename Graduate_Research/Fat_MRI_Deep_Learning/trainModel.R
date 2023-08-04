library( ANTsR )
library( keras )
library( ggplot2 )

baseDirectory <- '../'
dataDirectory <- paste0( baseDirectory, 'Data/' )
modelDirectory <- '/Users/ntustison/Pkg/ANTsRNet/Models/'
trainingDirectory <- paste0( dataDirectory, 'TrainingDataExpanded/' )

source( paste0( modelDirectory, 'createUnetModel.R' ) )
source( paste0( modelDirectory, 'unetUtilities.R' ) )

trainingImageFiles <- list.files( 
  path = paste0( trainingDirectory, 'Images' ), 
  pattern = "main_N4Denoised", full.names = TRUE )
trainingMaskFiles <- list.files( 
  path = paste0( trainingDirectory, 'Segmentations' ), 
  pattern = "segmentation", full.names = TRUE )

trainingProportion <- 0.3
set.seed( 1234 )
trainingIndices <- sample.int( 
  length( trainingMaskFiles ), size = length( trainingMaskFiles ) * trainingProportion )

trainingImage <- antsImageRead( trainingImageFiles[trainingIndices[1]], dimension = 2 )   
trainingData <- array( dim = c( length( trainingIndices ), dim( trainingImage ), 1 ) )

trainingMask <- antsImageRead( trainingMaskFiles[trainingIndices[1]], dimension = 2 )    
segmentationLabels <- sort( unique( as.vector( as.array( trainingMask ) ) ) )
numberOfLabels <- length( segmentationLabels )
trainingLabelData <- array( dim = c( length( trainingIndices ), dim( trainingImage ) ) )

cat( "Segmentation with ", numberOfLabels, 
  " labels: ", segmentationLabels, ".\n", sep = "" )

cat( "Reading training data (n = ", length( trainingIndices ), ")\n", sep = '' )
pb <- txtProgressBar( min = 0, max = length( trainingIndices ), style = 3 )
for ( i in 1:length( trainingIndices ) )
  {
  trainingImage <- antsImageRead( trainingImageFiles[trainingIndices[i]], dimension = 2 )   
  trainingArray <- as.array( trainingImage )
  trainingArray <- ( trainingArray - mean( trainingArray ) ) / sd( trainingArray )
  trainingData[i,,,1] <- trainingArray
  trainingMask <- as.array( antsImageRead( trainingMaskFiles[trainingIndices[i]], dimension = 2 ) )
  trainingLabelData[i,,] <- trainingMask

  if( i %% 100 == 0 )
    {
    gc( verbose = FALSE )
    }

  setTxtProgressBar( pb, i )
  }
cat( "\n")

X_train <- trainingData
Y_train <- encodeY( trainingLabelData, segmentationLabels )

unetModel <- createUnetModel2D( c( dim( trainingImage ), 1 ), 
  numberOfClassificationLabels = numberOfLabels, layers = 1:4 )

unetModel %>% compile( loss = loss_multilabel_dice_coefficient_error,
  optimizer = optimizer_adam( lr = 0.0001 ),  
  metrics = c( multilabel_dice_coefficient ) )

track <- unetModel %>% fit( X_train, Y_train, 
                 epochs = 100, batch_size = 32, verbose = 1, shuffle = TRUE,
                 callbacks = list( 
                   callback_model_checkpoint( paste0( baseDirectory, "unetWeights.h5" ), 
                     monitor = 'val_loss', save_best_only = TRUE, save_weights_only = TRUE,
                     verbose = 1, mode = 'auto', period = 1 ),
                   callback_early_stopping( monitor = 'val_loss', min_delta = 0.001, 
                     patience = 10 ),
                   callback_reduce_lr_on_plateau( monitor = 'val_loss', factor = 0.5,
                     patience = 0, epsilon = 0.001, cooldown = 0 )
                 ), 
                 validation_split = 0.2 )
## Save the model

save_model_weights_hdf5( 
  unetModel, filepath = paste0( baseDirectory, 'unetWeights.h5' ) )
save_model_hdf5( 
  unetModel, filepath = paste0( baseDirectory, 'unetModel.h5' ), overwrite = TRUE )

## Plot the model fitting

epochs <- 1:length( track$metrics$loss )

unetModelDataFrame <- data.frame( 
  Epoch = rep( epochs, 2 ), 
  Type = c( rep( 'Training', length( epochs ) ), rep( 'Validation', length( epochs ) ) ),
  Loss =c( track$metrics$loss, track$metrics$val_loss ), 
  Accuracy = c( track$metrics$multilabel_dice_coefficient, 
                track$metrics$val_multilabel_dice_coefficient )
  )

unetModelLossPlot <- 
  ggplot( data = unetModelDataFrame, aes( x = Epoch, y = Loss, colour = Type ) ) +
  geom_point( shape = 1, size = 0.5 ) +
  geom_line( size = 0.3 ) +
  ggtitle( "Loss" )

unetModelAccuracyPlot <- 
  ggplot( data = unetModelDataFrame, aes( x = Epoch, y = Accuracy, colour = Type ) ) +
  geom_point( shape = 1, size = 0.5 ) +
  geom_line( size = 0.3 ) +
  ggtitle( "Accuracy")

ggsave( paste0( baseDirectory, "unetModelLossPlot.pdf" ), 
  plot = unetModelLossPlot, width = 5, height = 2, units = 'in' )
ggsave( paste0( baseDirectory, "unetModelAccuracyPlot.pdf" ), 
  plot = unetModelAccuracyPlot, width = 5, height = 2, units = 'in' )









