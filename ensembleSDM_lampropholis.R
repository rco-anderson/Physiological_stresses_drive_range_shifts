library(sf)
library(dismo)
library(randomForest)
library(gbm)
library(e1071)
library(maxnet)
library(xgboost)
library(lightgbm)
library(caret)
library(ROCR)
library(dplyr)
library(geosphere)
library(reshape2)
library(ggspatial)


Sys.setenv(JAVA_HOME='/Library/Java/JavaVirtualMachines/zulu-23.jdk/Contents/Home') #rJava

set.seed(123)

# Load data
data <- read.csv("data/similis.csv")

# Load shapefile for species distribution
shp <- read_sf("shapefiles/guichenoti.shp")
                          #%>% st_buffer(dist = 50000)

# Convert to sf
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = st_crs(shp))

# Presence-absence data based on spatial intersection
data_sf$presence <- ifelse(st_intersects(data_sf, shp, sparse = FALSE), 1, 0)

# Filter current scenario
current_data_sf <- data_sf %>% filter(scenario == "current")

# Drop geometry and select relevant columns
current_df <- current_data_sf %>%
  st_drop_geometry()

current_df <- subset(current_df, select = -c(scenario, id))

# Manually specify the column names to select
#selected_cols <- c("TC.mean", "maxTC", "speed", "met", "ewl", "dehy", 
#                   "energ_comsump", "shade", "depth", "TC.act", 
#                   "maxTC.act", "speed.act", "met.act", "ewl.act", 
#                   "dehy.act", "energ_comsump.act","presence")

# Subset the dataframe based on selected column names
#current_df <- current_df[, selected_cols]


# Convert presence to factor
current_df$presence <- as.factor(current_df$presence)

current_df <- na.omit(current_df)


# Split data into training (70%) and testing (30%) sets
trainIndex <- createDataPartition(current_df$presence, p = 0.7, list = FALSE)
trainData <- current_df[trainIndex, ]
testData <- current_df[-trainIndex, ]


### GLM Model ###

glm_model <- glm(presence ~ ., data = trainData, family = binomial)
glm_preds <- predict(glm_model, newdata = testData, type = "response")

# AUC
glm_pred_obj <- prediction(glm_preds, testData$presence)
glm_auc <- performance(glm_pred_obj, measure = "auc")@y.values[[1]]

# TSS
glm_confusion <- table(testData$presence, glm_preds > 0.5)
glm_sensitivity <- glm_confusion[2, 2] / sum(glm_confusion[2, ])
glm_specificity <- glm_confusion[1, 1] / sum(glm_confusion[1, ])
glm_tss <- glm_sensitivity + glm_specificity - 1


### GBM Model ###

trainData$presence <- as.numeric(trainData$presence) - 1

gbm_model <- gbm(presence ~ ., data = trainData, distribution = "bernoulli", n.trees = 100)
gbm_preds <- predict(gbm_model, newdata = testData, type = "response")

# AUC
gbm_pred_obj <- prediction(gbm_preds, testData$presence)
gbm_auc <- performance(gbm_pred_obj, measure = "auc")@y.values[[1]]

# TSS
gbm_confusion <- table(testData$presence, gbm_preds > 0.5)
gbm_sensitivity <- gbm_confusion[2, 2] / sum(gbm_confusion[2, ])
gbm_specificity <- gbm_confusion[1, 1] / sum(gbm_confusion[1, ])
gbm_tss <- gbm_sensitivity + gbm_specificity - 1


### MaxEnt Model ###

# Convert data frame to a matrix of predictors
predictors <- trainData[, -which(names(trainData) == "presence")]

# Train MaxEnt model
maxent_model <- maxent(x = predictors, p = trainData$presence)

# Predict on test data
test_predictors <- testData[, -which(names(testData) == "presence")]
maxent_preds <- predict(maxent_model, test_predictors)

# AUC
maxent_pred_obj <- prediction(maxent_preds, testData$presence)
maxent_auc <- performance(maxent_pred_obj, measure = "auc")@y.values[[1]]

# TSS
maxent_confusion <- table(testData$presence, maxent_preds > 0.5)
maxent_sensitivity <- maxent_confusion[2, 2] / sum(maxent_confusion[2, ])
maxent_specificity <- maxent_confusion[1, 1] / sum(maxent_confusion[1, ])
maxent_tss <- maxent_sensitivity + maxent_specificity - 1



### Random Forest (RF) Model ###

rf_model <- randomForest(presence ~ ., data = trainData, probability = TRUE, ntree = 100)
rf_preds <- predict(rf_model, newdata = testData, type = "response")

# AUC
rf_pred_obj <- prediction(rf_preds, testData$presence)
rf_auc <- performance(rf_pred_obj, measure = "auc")@y.values[[1]]

# TSS
rf_confusion <- table(testData$presence, rf_preds > 0.5)
rf_sensitivity <- rf_confusion[2, 2] / sum(rf_confusion[2, ])
rf_specificity <- rf_confusion[1, 1] / sum(rf_confusion[1, ])
rf_tss <- rf_sensitivity + rf_specificity - 1



### XGBoost Model ###

xgb_train <- xgb.DMatrix(data = as.matrix(trainData[, -which(names(trainData) == "presence")]), 
                         label = trainData$presence)

xgb_model <- xgboost(data = xgb_train, objective = "binary:logistic", nrounds = 100)
xgb_preds <- predict(xgb_model, newdata = as.matrix(testData[, -which(names(testData) == "presence")]))

# AUC
xgb_pred_obj <- prediction(xgb_preds, testData$presence)
xgb_auc <- performance(xgb_pred_obj, measure = "auc")@y.values[[1]]

# TSS
xgb_confusion <- table(testData$presence, xgb_preds > 0.5)
xgb_sensitivity <- xgb_confusion[2, 2] / sum(xgb_confusion[2, ])
xgb_specificity <- xgb_confusion[1, 1] / sum(xgb_confusion[1, ])
xgb_tss <- xgb_sensitivity + xgb_specificity - 1



### LightGBM Model ###

lgb_train <- lgb.Dataset(data = as.matrix(trainData[, -which(names(trainData) == "presence")]), 
                         label = trainData$presence)

lgb_test <- lgb.Dataset(data = as.matrix(testData[, -1]), label = testData$total_tourists, reference = lgb_train)

lgb_model <- lgb.train(obj = "binary", 
                       data = lgb_train, 
                       nrounds = 100,
                       valids = list(test = lgb_test), 
                       verbose = 1)

lgb_preds <- predict(lgb_model, as.matrix(testData[, -which(names(testData) == "presence")]))

# AUC
lgb_pred_obj <- prediction(lgb_preds, testData$presence)
lgb_auc <- performance(lgb_pred_obj, measure = "auc")@y.values[[1]]

# TSS
lgb_confusion <- table(testData$presence, lgb_preds > 0.5)
lgb_sensitivity <- lgb_confusion[2, 2] / sum(lgb_confusion[2, ])
lgb_specificity <- lgb_confusion[1, 1] / sum(lgb_confusion[1, ])
lgb_tss <- lgb_sensitivity + lgb_specificity - 1



### Summary of AUC and TSS for all models ###

model_performance <- data.frame(
  Model = c("GLM", "GBM", "MaxEnt", "RF", "XGBoost", "LightGBM"),
  AUC = c(glm_auc, gbm_auc, maxent_auc, rf_auc, xgb_auc, lgb_auc),
  TSS = c(glm_tss, gbm_tss, maxent_tss, rf_tss, xgb_tss, lgb_tss)
)

print(model_performance)

# Plot the ROC curves for all models
roc_glm <- performance(glm_pred_obj, measure = "tpr", x.measure = "fpr")
roc_gbm <- performance(gbm_pred_obj, measure = "tpr", x.measure = "fpr")
roc_maxent <- performance(maxent_pred_obj, measure = "tpr", x.measure = "fpr")
roc_rf <- performance(rf_pred_obj, measure = "tpr", x.measure = "fpr")
roc_xgb <- performance(xgb_pred_obj, measure = "tpr", x.measure = "fpr")
roc_lgb <- performance(lgb_pred_obj, measure = "tpr", x.measure = "fpr")

plot(roc_glm, col = "#0072B2", main = "ROC Curves", lwd = 2)
plot(roc_gbm, col = "#D55E00", add = TRUE, lwd = 2)
plot(roc_maxent, col = "#009E73", add = TRUE, lwd = 2)
plot(roc_rf, col = "#CC79A7", add = TRUE, lwd = 2)
plot(roc_xgb, col = "#56B4E9", add = TRUE, lwd = 2)
plot(roc_lgb, col = "#E69F00", add = TRUE, lwd = 2)

legend("bottomright", legend = c("GLM", "GBM", "MaxEnt", "RF","XGBoost", "LightGBM"),
       col = c("#0072B2", "#D55E00", "#009E73", "#CC79A7",  "#56B4E9", "#E69F00"), lwd = 2)

# Plot TSS comparison
barplot(model_performance$TSS, names.arg = model_performance$Model, col = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#56B4E9", "#E69F00", "#000000"),
        main = "TSS Comparison Across Models", ylim = c(0, 1), ylab = "TSS")



### Build ensemble model ###


# TSS values for all models
tss_values <- c(glm_tss, gbm_tss, maxent_tss, rf_tss, xgb_tss, lgb_tss)
names(tss_values) <- c("GLM", "GBM", "MaxEnt", "RF", "XGBoost", "LightGBM")

# Filter models with TSS > 0.6
filtered_tss_values <- tss_values[tss_values > 0.6]
print(filtered_tss_values)

# Calculate normalized weights based on filtered TSS
normalized_weights <- filtered_tss_values / sum(filtered_tss_values)
print(normalized_weights)


# Combine predictions from the filtered models
ensemble_preds <- rep(0, length(glm_preds))  # Initialize with zeros


if ("GLM" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["GLM"] * glm_preds
}
if ("GBM" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["GBM"] * gbm_preds
}
if ("MaxEnt" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["MaxEnt"] * maxent_preds
}
if ("RF" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["RF"] * rf_preds
}
if ("SVM" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["SVM"] * svm_preds_prob
}
if ("XGBoost" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["XGBoost"] * xgb_preds
}
if ("LightGBM" %in% names(normalized_weights)) {
  ensemble_preds <- ensemble_preds + normalized_weights["LightGBM"] * lgb_preds
}

# AUC ensemble model
ensemble_pred_obj <- prediction(ensemble_preds, testData$presence)
ensemble_auc <- performance(ensemble_pred_obj, measure = "auc")@y.values[[1]]
print(paste("Ensemble AUC:", ensemble_auc))

# TSS ensemble model
ensemble_confusion <- table(testData$presence, ensemble_preds > 0.6)
ensemble_sensitivity <- ensemble_confusion[2, 2] / sum(ensemble_confusion[2, ])
ensemble_specificity <- ensemble_confusion[1, 1] / sum(ensemble_confusion[1, ])
ensemble_tss <- ensemble_sensitivity + ensemble_specificity - 1
print(paste("Ensemble TSS:", ensemble_tss))



### Mapping species distribution  ###


## current ##

# Filter the dataset for the 'plus2' scenario and drop geometry
current <- data_sf %>% filter(scenario == "current") %>%
                          st_drop_geometry() %>%
                          subset(select = -c(scenario, id))
current$presence <- NA

# Convert new data to a matrix of predictors (excluding the presence column)
new_predictors <- current[, -which(names(current) == "presence")]

# Initialize a vector for combined predictions
combined_preds_new <- rep(0, nrow(current))

# GLM Predictions (if TSS > 0.6)
if ("GLM" %in% names(normalized_weights)) {
  glm_preds_new <- predict(glm_model, newdata = current, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GLM"] * glm_preds_new
}

# GBM Predictions (if TSS > 0.6)
if ("GBM" %in% names(normalized_weights)) {
  gbm_preds_new <- predict(gbm_model, newdata = current, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GBM"] * gbm_preds_new
}

# MaxEnt Predictions (if TSS > 0.6)
if ("MaxEnt" %in% names(normalized_weights)) {
  maxent_preds_new <- predict(maxent_model, current)
  combined_preds_new <- combined_preds_new + normalized_weights["MaxEnt"] * maxent_preds_new
}

# RF Predictions (if TSS > 0.6)
if ("RF" %in% names(normalized_weights)) {
  rf_preds_new <- predict(rf_model, newdata = current, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["RF"] * rf_preds_new
}

# XGBoost Predictions (if TSS > 0.6)
if ("XGBoost" %in% names(normalized_weights)) {
  xgb_preds_new <- predict(xgb_model, as.matrix(new_predictors), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["XGBoost"] * xgb_preds_new
}

# LightGBM Predictions (if TSS > 0.6)
if ("LightGBM" %in% names(normalized_weights)) {
  lgb_preds_new <- predict(lgb_model, as.matrix(new_predictors), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["LightGBM"] * lgb_preds_new
}

# Convert combined predictions to probabilities (ensure values between 0 and 1)
combined_probs_new <- pmin(pmax(combined_preds_new, 0), 1)

# Prepare the map data with the new predictions
current_map <- data_sf %>% filter(scenario == "current")
current_map$prob <- combined_probs_new

# Convert sf object to Spatial object (if not already done)
sp_current_map <- as_Spatial(current_map)

# Create an empty raster with the desired resolution and extent
current.r <- raster(extent(sp_current_map), resolution = 0.25) # Adjust resolution as needed

# Rasterize the spatial data
current.r <- rasterize(sp_current_map, current.r, field = "prob", fun = mean)

# Plot the raster using ggplot
ggplot() +
  geom_raster(data = as.data.frame(current.r, xy = TRUE), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(option = "viridis", name = "Probability") +
  geom_sf(data = shp, fill = NA, color = "red", size = 0.5) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")



## plus2 ##


# Filter the dataset for the 'plus2' scenario and drop geometry
plus2 <- data_sf %>% filter(scenario == "plus2") %>%
  st_drop_geometry() %>%
  subset(select = -c(scenario, id))

plus2$presence <- NULL

# Convert new data to a matrix of predictors (excluding the presence column)
#new_predictors <- current[, -which(names(plus2) == "presence")]

# Initialize a vector for combined predictions
combined_preds_new <- rep(0, nrow(plus2))

# GLM Predictions (if TSS > 0.6)
if ("GLM" %in% names(normalized_weights)) {
  glm_preds_new <- predict(glm_model, newdata = plus2, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GLM"] * glm_preds_new
}

# GBM Predictions (if TSS > 0.6)
if ("GBM" %in% names(normalized_weights)) {
  gbm_preds_new <- predict(gbm_model, newdata = plus2, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GBM"] * gbm_preds_new
}

# MaxEnt Predictions (if TSS > 0.6)
if ("MaxEnt" %in% names(normalized_weights)) {
  maxent_preds_new <- predict(maxent_model, plus2)
  combined_preds_new <- combined_preds_new + normalized_weights["MaxEnt"] * maxent_preds_new
}

# RF Predictions (if TSS > 0.6)
if ("RF" %in% names(normalized_weights)) {
  rf_preds_new <- predict(rf_model, newdata = plus2, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["RF"] * rf_preds_new
}

# XGBoost Predictions (if TSS > 0.6)
if ("XGBoost" %in% names(normalized_weights)) {
  xgb_preds_new <- predict(xgb_model, as.matrix(plus2), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["XGBoost"] * xgb_preds_new
}

# LightGBM Predictions (if TSS > 0.6)
if ("LightGBM" %in% names(normalized_weights)) {
  lgb_preds_new <- predict(lgb_model, as.matrix(plus2), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["LightGBM"] * lgb_preds_new
}

# Convert combined predictions to probabilities (ensure values between 0 and 1)
combined_probs_new <- pmin(pmax(combined_preds_new, 0), 1)

# Prepare the map data with the new predictions
plus2_map <- data_sf %>% filter(scenario == "plus2")
plus2_map$prob <- combined_probs_new

# Convert sf object to Spatial object (if not already done)
sp_plus2_map <- as_Spatial(plus2_map)

# Create an empty raster with the desired resolution and extent
plus2.r <- raster(extent(sp_plus2_map), resolution = 0.25) # Adjust resolution as needed

# Rasterize the spatial data
plus2.r <- rasterize(sp_plus2_map, plus2.r, field = "prob", fun = mean)

# Plot the raster using ggplot
ggplot() +
  geom_raster(data = as.data.frame(plus2.r, xy = TRUE), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(option = "viridis", name = "Probability") +
  geom_sf(data = shp, fill = NA, color = "red", size = 0.5) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")



## plus4 ##


# Filter the dataset for the 'plus2' scenario and drop geometry
plus4 <- data_sf %>% filter(scenario == "plus4") %>%
  st_drop_geometry() %>%
  subset(select = -c(scenario, id))

plus4$presence <- NULL

# Convert new data to a matrix of predictors (excluding the presence column)
#new_predictors <- current[, -which(names(plus2) == "presence")]

# Initialize a vector for combined predictions
combined_preds_new <- rep(0, nrow(plus4))

# GLM Predictions (if TSS > 0.6)
if ("GLM" %in% names(normalized_weights)) {
  glm_preds_new <- predict(glm_model, newdata = plus4, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GLM"] * glm_preds_new
}

# GBM Predictions (if TSS > 0.6)
if ("GBM" %in% names(normalized_weights)) {
  gbm_preds_new <- predict(gbm_model, newdata = plus4, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["GBM"] * gbm_preds_new
}

# MaxEnt Predictions (if TSS > 0.6)
if ("MaxEnt" %in% names(normalized_weights)) {
  maxent_preds_new <- predict(maxent_model, plus4)
  combined_preds_new <- combined_preds_new + normalized_weights["MaxEnt"] * maxent_preds_new
}

# RF Predictions (if TSS > 0.6)
if ("RF" %in% names(normalized_weights)) {
  rf_preds_new <- predict(rf_model, newdata = plus4, type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["RF"] * rf_preds_new
}

# XGBoost Predictions (if TSS > 0.6)
if ("XGBoost" %in% names(normalized_weights)) {
  xgb_preds_new <- predict(xgb_model, as.matrix(plus4), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["XGBoost"] * xgb_preds_new
}

# LightGBM Predictions (if TSS > 0.6)
if ("LightGBM" %in% names(normalized_weights)) {
  lgb_preds_new <- predict(lgb_model, as.matrix(plus4), type = "response")
  combined_preds_new <- combined_preds_new + normalized_weights["LightGBM"] * lgb_preds_new
}

# Convert combined predictions to probabilities (ensure values between 0 and 1)
combined_probs_new <- pmin(pmax(combined_preds_new, 0), 1)

# Prepare the map data with the new predictions
plus4_map <- data_sf %>% filter(scenario == "plus4")
plus4_map$prob <- combined_probs_new

# Convert sf object to Spatial object (if not already done)
sp_plus4_map <- as_Spatial(plus4_map)

# Create an empty raster with the desired resolution and extent
plus4.r <- raster(extent(sp_plus4_map), resolution = 0.25) # Adjust resolution as needed

# Rasterize the spatial data
plus4.r <- rasterize(sp_plus4_map, plus4.r, field = "prob", fun = mean)

# Plot the raster using ggplot
ggplot() +
  geom_raster(data = as.data.frame(plus4.r, xy = TRUE), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(option = "viridis", name = "Probability") +
  geom_sf(data = shp, fill = NA, color = "red", size = 0.5) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")


# rasters #

crs(current.r) <- crs(shp)
crs(plus2.r) <- crs(shp)
crs(plus4.r) <- crs(shp)

## adonis
#adonis_current <- current.r
#adonis_plus2 <- plus2.r
#adonis_plus4 <- plus4.r

## amicula
#amicula_current <- current.r
#amicula_plus2 <- plus2.r
#amicula_plus4 <- plus4.r

## caligula
#caligula_current <- current.r
#caligula_plus2 <- plus2.r
#caligula_plus4 <- plus4.r

## couperi
#couperi_current <- current.r
#couperi_plus2 <- plus2.r
#couperi_plus4 <- plus4.r

## delicata brisb
#delicata.brisb_current <- current.r
#delicata.brisb_plus2 <- plus2.r
#delicata.brisb_plus4 <- plus4.r

## delicata nsw
#delicata.nsw_current <- current.r
#delicata.nsw_plus2 <- plus2.r
#delicata.nsw_plus4 <- plus4.r

## delicata twv
#delicata.twv_current <- current.r
#delicata.twv_plus2 <- plus2.r
#delicata.twv_plus4 <- plus4.r

## guichenoti
#guichenoti_current <- current.r
#guichenoti_plus2 <- plus2.r
#guichenoti_plus4 <- plus4.r

## similis
#similis_current <- current.r
#similis_plus2 <- plus2.r
#similis_plus4 <- plus4.r


### Plot all maps ###

# Load necessary libraries
library(ggplot2)
library(sf)
library(dplyr)

# Convert rasters to data frames
convert_to_df <- function(raster, species, scenario) {
  df <- as.data.frame(raster, xy = TRUE)
  df$species <- species
  df$scenario <- scenario
  return(df)
}

# List of species
species_list <- c("adonis", "amicula", "caligula", "couperi", 
                  "delicata.brisb", "delicata.nsw", "delicata.twv", 
                  "guichenoti", "similis")

# Read shapefiles for each species
shapefiles <- list(
  "adonis" = read_sf("shapefiles/adonis_shp.shp"),
  "amicula" = read_sf("shapefiles/amicula_shp.shp"),
  "caligula" = read_sf("shapefiles/caligula.shp"),
  "couperi" = read_sf("shapefiles/couperi_shp.shp"),
  "delicata.brisb" = read_sf("shapefiles/delicata_shp.shp"),
  "delicata.nsw" = read_sf("shapefiles/delicata_shp.shp"),
  "delicata.twv" = read_sf("shapefiles/delicata_shp.shp"),
  "guichenoti" = read_sf("shapefiles/guichenoti.shp"),
  "similis" = read_sf("shapefiles/similis.shp")
)

# Create a list to store all data frames
all_data <- list()

# Loop through each species and scenario to create data frames
for (species in species_list) {
  current_raster <- get(paste0(species, "_current"))
  plus2_raster <- get(paste0(species, "_plus2"))
  plus4_raster <- get(paste0(species, "_plus4"))
  
  current_df <- convert_to_df(current_raster, species, "current")
  plus2_df <- convert_to_df(plus2_raster, species, "plus2")
  plus4_df <- convert_to_df(plus4_raster, species, "plus4")
  
  all_data[[species]] <- rbind(current_df, plus2_df, plus4_df)
}

# Combine all data frames into one
combined_data <- do.call(rbind, all_data) %>% na.omit()

# Define groups for splitting
delicata_guichenoti <- c("delicata.brisb", "delicata.nsw", "delicata.twv", "guichenoti")
other_species <- setdiff(species_list, delicata_guichenoti)

# Create separate data frames for the two groups
data_dg <- combined_data %>% filter(species %in% delicata_guichenoti)
data_other <- combined_data %>% filter(species %in% other_species)

# Spectral palette
custom_palette <- scale_fill_distiller(palette = "Spectral")

# Function to create plots ensuring each species gets its own shapefile
plot_species_group <- function(data) {
  ggplot() +
    geom_raster(data = data, aes(x = x, y = y, fill = layer)) +
    facet_grid(species ~ scenario) +
    custom_palette +
    theme_void() +
    theme(strip.text = element_blank(), legend.position = "none", panel.spacing = unit(.5, "lines")) +
    # Add shapefile dynamically for each species
    geom_sf(data = bind_rows(lapply(unique(data$species), function(sp) shapefiles[[sp]] %>% mutate(species = sp))),
            aes(geometry = geometry), color = "black", fill = NA, size = 0.25, inherit.aes = FALSE)
}

# Generate plots
plot_dg <- plot_species_group(data_dg)
plot_other <- plot_species_group(data_other)

# Print plots
plot_dg
plot_other


### Estimating range shift  ###


# Define a function to calculate suitable cells and area
calculate_suitable_area <- function(raster_layer, threshold = 0.4) {
  # Convert to binary (1 for suitable, 0 for unsuitable)
  binary_layer <- raster_layer >= threshold
  
  # Count the number of suitable cells
  suitable_cells <- sum(values(binary_layer), na.rm = TRUE)
  
  # Calculate the total area in square kilometers
  cell_area_km2 <- area(raster_layer)  # Area of each cell in km²
  total_area_km2 <- sum(values(cell_area_km2 * binary_layer), na.rm = TRUE)
  
  return(list(suitable_cells = suitable_cells, total_area_km2 = total_area_km2))
}

# List of species and their corresponding raster objects
species_rasters <- list(
  adonis = list(current = adonis_current, plus2 = adonis_plus2, plus4 = adonis_plus4),
  amicula = list(current = amicula_current, plus2 = amicula_plus2, plus4 = amicula_plus4),
  caligula = list(current = caligula_current, plus2 = caligula_plus2, plus4 = caligula_plus4),
  couperi = list(current = couperi_current, plus2 = couperi_plus2, plus4 = couperi_plus4),
  delicata.brisb = list(current = delicata.brisb_current, plus2 = delicata.brisb_plus2, plus4 = delicata.brisb_plus4),
  delicata.nsw = list(current = delicata.nsw_current, plus2 = delicata.nsw_plus2, plus4 = delicata.nsw_plus4),
  delicata.twv = list(current = delicata.twv_current, plus2 = delicata.twv_plus2, plus4 = delicata.twv_plus4),
  guichenoti = list(current = guichenoti_current, plus2 = guichenoti_plus2, plus4 = guichenoti_plus4),
  similis = list(current = similis_current, plus2 = similis_plus2, plus4 = similis_plus4)
)

# Create an empty data frame to store results
results <- data.frame(
  Species = character(),
  Scenario = character(),
  Suitable_Cells = numeric(),
  Total_Area_km2 = numeric(),
  stringsAsFactors = FALSE
)

# Function to calculate suitable area and centroid
calculate_suitable_area_centroid <- function(raster_layer, threshold = 0.4) {
  # Convert to binary (1 for suitable, 0 for unsuitable)
  binary_layer <- raster_layer >= threshold
  
  # Get coordinates of suitable cells
  suitable_coords <- coordinates(raster_layer)[values(binary_layer) == 1, ]
  
  # Calculate centroid
  if (nrow(suitable_coords) > 0) {
    centroid <- colMeans(suitable_coords, na.rm = TRUE)
  } else {
    centroid <- c(NA, NA)  # No suitable habitat
  }
  
  # Count suitable cells
  suitable_cells <- sum(values(binary_layer), na.rm = TRUE)
  
  # Calculate total area in square kilometers
  cell_area_km2 <- raster::area(raster_layer)
  total_area_km2 <- sum(values(cell_area_km2 * binary_layer), na.rm = TRUE)
  
  return(list(suitable_cells = suitable_cells, 
              total_area_km2 = total_area_km2, 
              centroid_lon = centroid[1], 
              centroid_lat = centroid[2]))
}

# Create an empty data frame to store results
results <- data.frame(Species = character(), Scenario = character(), 
                      Suitable_Cells = numeric(), Total_Area_km2 = numeric(),
                      Centroid_Lon = numeric(), Centroid_Lat = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through each species and scenario
for (species in names(species_rasters)) {
  for (scenario in names(species_rasters[[species]])) {
    raster_layer <- species_rasters[[species]][[scenario]]
    
    # Calculate stats and centroid
    stats <- calculate_suitable_area_centroid(raster_layer)
    
    # Append results
    results <- rbind(results, data.frame(
      Species = species,
      Scenario = scenario,
      Suitable_Cells = stats$suitable_cells,
      Total_Area_km2 = stats$total_area_km2,
      Centroid_Lon = stats$centroid_lon,
      Centroid_Lat = stats$centroid_lat
    ))
  }
}

# Compute shifts in centroids
shifts <- data.frame(Species = character(), Shift_2C_km = numeric(), Shift_4C_km = numeric(), stringsAsFactors = FALSE)

for (species in unique(results$Species)) {
  current <- subset(results, Species == species & Scenario == "current", select = c(Centroid_Lon, Centroid_Lat))
  plus2 <- subset(results, Species == species & Scenario == "plus2", select = c(Centroid_Lon, Centroid_Lat))
  plus4 <- subset(results, Species == species & Scenario == "plus4", select = c(Centroid_Lon, Centroid_Lat))
  
  if (!any(is.na(current)) & !any(is.na(plus2))) {
    shift_2C <- distHaversine(as.numeric(current), as.numeric(plus2)) / 1000
  } else {
    shift_2C <- NA
  }
  
  if (!any(is.na(current)) & !any(is.na(plus4))) {
    shift_4C <- distHaversine(as.numeric(current), as.numeric(plus4)) / 1000
  } else {
    shift_4C <- NA
  }
  
  shifts <- rbind(shifts, data.frame(Species = species, Shift_2C_km = shift_2C, Shift_4C_km = shift_4C))
}

# Print results
print(results)
print(shifts)


# Create a lookup table for species names
species_lookup <- c(
  "adonis" = "Lampropholis adonis",
  "amicula" = "Lampropholis amicula",
  "caligula" = "Lampropholis caligula",
  "couperi" = "Lampropholis couperi",
  "delicata.brisb" = "Lampropholis delicata Brisb",
  "delicata.nsw" = "Lampropholis delicata NSW",
  "delicata.twv" = "Lampropholis delicata Towns",
  "guichenoti" = "Lampropholis guichenoti",
  "similis" = "Lampropholis similis"
)

# Replace species names using the lookup table
results <- results %>%
  mutate(Species = species_lookup[Species])

shifts <- shifts %>%
  mutate(Species = species_lookup[Species])

# Reshape data from wide to long format
shifts_long <- shifts %>%
  pivot_longer(cols = c(Shift_2C_km, Shift_4C_km), 
               names_to = "Scenario", 
               values_to = "Shift_km") %>%
  mutate(Distance = recode(Scenario, 
                           "Shift_2C_km" = "plus2", 
                           "Shift_4C_km" = "plus4"))

### plot area reduction ###

# Convert 'Scenario' to factor with ordered levels
results$Scenario <- factor(results$Scenario, levels = c("current", "plus2", "plus4"))

# Create the line plot
ggplot(results, aes(x = Scenario, y = Total_Area_km2, group = Species, colour = Species)) +
  geom_line(linewidth = 1) +  # Line connecting points
  geom_point(size = 2) +      # Points on the line
  labs(x = "Climate Scenario", y = "Total Suitable Area (km²)") +
  scale_color_viridis_d(option = "plasma") +
  theme_bw() +
  theme(aspect.ratio = 1.2)  # Adjust legend position


### plot centroid shifts ###

# Reshape data for arrows
centroids_wide <- results %>%
  dplyr::select(Species, Scenario, Centroid_Lon, Centroid_Lat) %>%
  pivot_wider(names_from = Scenario, values_from = c(Centroid_Lon, Centroid_Lat))

# Load the map of Australia
australia_map <- map_data("world") %>%
  filter(region == "Australia")

# Get species bounding box for zooming (slightly zoomed out)
buffer <- 4.5  # Adjust this value for more/less zoom
lon_min <- min(centroids_wide$Centroid_Lon_current, na.rm = TRUE) - buffer
lon_max <- max(centroids_wide$Centroid_Lon_current, na.rm = TRUE) + buffer
lat_min <- min(centroids_wide$Centroid_Lat_current, na.rm = TRUE) - 3.9
lat_max <- max(centroids_wide$Centroid_Lat_current, na.rm = TRUE) + 6

# Main plot (zoomed-in area)
main_plot <- ggplot() +
  # Australia map (background)
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group), 
               fill = "gray95", colour = "black", linewidth = .25) +
  
  # Current species locations (black dots)
  geom_point(data = centroids_wide, 
             aes(x = Centroid_Lon_current, y = Centroid_Lat_current), 
             colour = "white", fill = "black", size = 2, shape = 21) +  
  
  # Arrows for future shifts (blue for +2°C, red for +4°C)
  geom_segment(data = centroids_wide, 
               aes(x = Centroid_Lon_current, y = Centroid_Lat_current,
                   xend = Centroid_Lon_plus2, yend = Centroid_Lat_plus2), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5, 
               colour = "#F39C12", na.rm = TRUE) +
  
  geom_segment(data = centroids_wide, 
               aes(x = Centroid_Lon_current, y = Centroid_Lat_current,
                   xend = Centroid_Lon_plus4, yend = Centroid_Lat_plus4), 
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5, 
               colour = "#E74C3C", na.rm = TRUE) +
  
  # Labels for species
  geom_text_repel(data = centroids_wide, 
                  aes(x = Centroid_Lon_current, y = Centroid_Lat_current, 
                      label = Species), size = 3)  +
  
  # Labels and theme
  labs(x = "Longitude", y = "Latitude") +
  
  # Zoom into species distribution area
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), crs = 4326) +
  annotation_scale(location = "br") +
  
  theme_bw() +
  theme(legend.position = "none")

main_plot

# Inset map (full Australia with zoomed region)
inset_map <- ggplot() +
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group), 
               fill = "gray90", colour = "black") +
  
  # Highlight the zoomed-in area with a red rectangle
  geom_rect(aes(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max), 
            fill = NA, colour = "red", linewidth = 1) +
  
  theme_void()  # Remove all axes, titles, etc.

# Combine main plot and inset
grid.newpage()
print(main_plot, vp = viewport(width = 0.9, height = 0.9))  # Main plot
print(inset_map, vp = viewport(x = 0.2, y = 0.8, width = 0.3, height = 0.3))  # Inset map position
