# --- Libraries ---
# Load necessary libraries
library(sf)
library(dplyr)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)
library(MuMIn)
library(phyr)
library(rr2)
library(Matrix)
library(ggplot2)
library(tidyr)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)

setwd("~/Library/CloudStorage/Dropbox/Lampropholis range shift")

# --- Data Input ---
# Read data files for different species
dat.ad <- read.csv("data/adonis.csv")            # Lampropholis adonis
dat.am <- read.csv("data/amicula.csv")           # Lampropholis amicula
dat.ca <- read.csv("data/caligula.csv")          # Lampropholis caligula
dat.co <- read.csv("data/couperi.csv")           # Lampropholis couperi
dat.del.b <- read.csv("data/delicata.brisb.csv") # Lampropholis delicata brisbane
dat.del.n <- read.csv("data/delicata.nsw.csv")   # Lampropholis delicata NSW
dat.del.t <- read.csv("data/delicata.twv.csv")   # Lampropholis delicata Townsville
dat.gu <- read.csv("data/guichenoti.csv")        # Lampropholis guichenoti
dat.si <- read.csv("data/similis.csv")        # Lampropholis guichenoti

# Combine species data into a list for processing
data_files <- list(
  "Lampropholis adonis" = dat.ad,
  "Lampropholis amicula" = dat.am,
  "Lampropholis caligula" = dat.ca,
  "Lampropholis couperi" = dat.co,
  "Lampropholis delicata Brisb" = dat.del.b,
  "Lampropholis delicata NSW" = dat.del.n,
  "Lampropholis delicata Towns" = dat.del.t,
  "Lampropholis guichenoti" = dat.gu,
  "Lampropholis similis" = dat.si
)

# Read shapefiles for each species
shapefiles <- list(
  "Lampropholis adonis" = read_sf("shapefiles/adonis_shp.shp"),
  "Lampropholis amicula" = read_sf("shapefiles/amicula_shp.shp"),
  "Lampropholis caligula" = read_sf("shapefiles/caligula.shp"),
  "Lampropholis couperi" = read_sf("shapefiles/couperi_shp.shp"),
  "Lampropholis delicata Brisb" = read_sf("shapefiles/delicata_shp.shp"),
  "Lampropholis delicata NSW" = read_sf("shapefiles/delicata_shp.shp"),
  "Lampropholis delicata Towns" = read_sf("shapefiles/delicata_shp.shp"),
  "Lampropholis guichenoti" = read_sf("shapefiles/guichenoti.shp"),
  "Lampropholis similis" = read_sf("shapefiles/similis.shp")
)

# --- Data Merging and Preparation ---

# Combine data from all species and check presence-absence based on spatial intersection
merged_data <- do.call(rbind, lapply(names(data_files), function(species) {
  df <- data_files[[species]]
  shp <- shapefiles[[species]]
  
  # Convert dataframe to spatial format (sf)
  df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = st_crs(shp))
  
  # Check for presence-absence (1 for presence, 0 for absence)
  df_sf$presence <- ifelse(st_intersects(df_sf, shp, sparse = FALSE), 1, 0)
  
  # Add species column
  df_sf$species <- species
  return(df_sf)
}))

# Clean merged data and drop geometry
merged_data <- merged_data %>% na.omit() %>% st_drop_geometry()

# Filter presence data for analysis
presence_data <- merged_data %>% filter(presence == 1)
#presence_data <- as.data.table(presence_data)


# --- Phylogenetic Tree ---

# Load the phylogenetic tree
filename <- "Lampropholis_tree.nex"
my_tree <- read.nexus(filename)

# Rename tree tip labels 
new_tiplabels <- c("Lampropholis guichenoti","Lampropholis miriabilis", 
                   "Lampropholis caligula","Lampropholis amicula",
                   "Lampropholis delicata Towns","Lampropholis delicata NSW",
                   "Lampropholis delicata Brisb","Lampropholis similis",
                   "Lampropholis adonis", "Lampropholis couperi",
                   "Carinascinus pretiosus")

my_tree$tip.label <- new_tiplabels

## for phylogenetic signal

# Define species to remove
species_to_remove <- c("Carinascincus pretiosus", "Lampropholis mirabilis")

# Remove them from the tree
pruned_tree <- drop.tip(my_tree, species_to_remove)

# Filter the dataset to include only species in the pruned tree
filtered_data <- presence_data[presence_data$species %in% pruned_tree$tip.label, ]

# Create vectors
act_h <- setNames(filtered_data$total_activity_hours, filtered_data$species)
dehydration_risk <- setNames(filtered_data$dehydrated.10, filtered_data$species)
ctmax_h <- setNames(filtered_data$ctmax, filtered_data$species)
energ.c <- setNames(filtered_data$energ_comsump, filtered_data$species)


# --- PGLMM Models ---

# - Total Activity Hours - #

# Full model
act.model.full  <- pglmm(total_activity_hours ~ species + scenario + (1|id), 
                   data = presence_data, 
                   cov_ranef = list(species = my_tree),
                   family = "gaussian", REML = FALSE)

act.model.full

# Reduced model
act.model.red  <- pglmm(total_activity_hours ~ scenario + (1|id), 
                         data = presence_data, 
                         cov_ranef = list(species = my_tree),
                         family = "gaussian", REML = FALSE)

act.model.red

# Likelihood ratio test
ll_full <- act.model.full$logLik
ll_reduced <- act.model.red$logLik

# Degrees of freedom
df_full <- length(act.model.full$B) + length(act.model.full$theta)
df_reduced <- length(act.model.red$B) + length(act.model.red$theta)

# Compute LRT
LRT <- 2 * (ll_full - ll_reduced)
df_diff <- df_full - df_reduced
pval <- pchisq(LRT, df = df_diff, lower.tail = FALSE)

format.pval(pval, digits = 3, eps = .Machine$double.eps)


# - Dehydration (hours above 10%) - #

# Full model
dehy.model.full <- pglmm(dehydrated.10 ~ species + scenario + (1|id), 
                    data = presence_data, 
                    cov_ranef = list(species = my_tree),
                    family = "gaussian", REML = FALSE)
dehy.model.full

#Reduced model
dehy.model.red <- pglmm(dehydrated.10 ~ scenario + (1|id), 
                       data = presence_data, 
                       cov_ranef = list(species = my_tree),
                       family = "gaussian", REML = FALSE)
dehy.model.red

# Likelihood ratio test
ll_full <- dehy.model.full$logLik
ll_reduced <- dehy.model.red$logLik

# Degrees of freedom (number of parameters)
df_full <- length(dehy.model.full$B) + length(dehy.model.full$theta)
df_reduced <- length(dehy.model.red$B) + length(dehy.model.red$theta)

# Compute LRT
LRT <- 2 * (ll_full - ll_reduced)
df_diff <- df_full - df_reduced
pval <- pchisq(LRT, df = df_diff, lower.tail = FALSE)

format.pval(pval, digits = 3, eps = .Machine$double.eps)


# - Hours above CTmax - #

# Full model
ctmax.model.full <- pglmm(ctmax ~ species + scenario + (1|id), 
                     data = presence_data, 
                     cov_ranef = list(species = my_tree),
                     family = "gaussian", REML = FALSE)

ctmax.model.full

# Reduced model
ctmax.model.red <- pglmm(ctmax ~ scenario + (1|id), 
                          data = presence_data, 
                          cov_ranef = list(species = my_tree),
                          family = "gaussian", REML = FALSE)

ctmax.model.red

# Likelihood ratio test
ll_full <- ctmax.model.full$logLik
ll_reduced <- ctmax.model.red$logLik

# Degrees of freedom (number of parameters)
df_full <- length(ctmax.model.full$B) + length(ctmax.model.full$theta)
df_reduced <- length(ctmax.model.red$B) + length(ctmax.model.red$theta)

# Compute LRT
LRT <- 2 * (ll_full - ll_reduced)
df_diff <- df_full - df_reduced
pval <- pchisq(LRT, df = df_diff, lower.tail = FALSE)

format.pval(pval, digits = 3, eps = .Machine$double.eps)


# - Energy consumption - #

#Full model
energy.model.full <- pglmm(energ_comsump ~ species + scenario + (1|id), 
                      data = presence_data, 
                      cov_ranef = list(species = my_tree),
                      family = "gaussian", REML = FALSE)
energy.model.full

#Reduced model
energy.model.red <- pglmm(energ_comsump ~ scenario + (1|id), 
                           data = presence_data, 
                           cov_ranef = list(species = my_tree),
                           family = "gaussian", REML = FALSE)
energy.model.red


# Likelihood ratio test
ll_full <- energy.model.full$logLik
ll_reduced <- energy.model.red$logLik

# Degrees of freedom (number of parameters)
df_full <- length(energy.model.full$B) + length(energy.model.full$theta)
df_reduced <- length(energy.model.red$B) + length(energy.model.red$theta)

# Compute LRT
LRT <- 2 * (ll_full - ll_reduced)
df_diff <- df_full - df_reduced
pval <- pchisq(LRT, df = df_diff, lower.tail = FALSE)

format.pval(pval, digits = 3, eps = .Machine$double.eps)


# Model for activity-based energy consumption (energ_comsump.act)
#energy.act.model <- pglmm(energ_comsump.act ~ species + scenario + (1|id), 
#                          data = presence_data, 
#                          cov_ranef = list(species = my_tree),
#                          family = "gaussian", REML = FALSE)
#energy.act.model


# Model for shade activity (shade.act)
#shade.model <- pglmm(shade.act ~ species*scenario + (1|id), 
#                     data = presence_data, 
#                     cov_ranef = list(species = my_tree),
#                     family = "gaussian", REML = FALSE)
#shade.model


## range shitfs ###

### Total area change ###

area.model <- pglmm(Total_Area_km2 ~ Scenario + (1 | Species), 
                    data = results, 
                    cov_ranef = list(Species = my_tree),
                    family = "gaussian", REML = FALSE)

area.model

# Phylogenetic signal
filtered_results <- results[results$Species %in% pruned_tree$tip.label, ]
total_area <- setNames(filtered_results$Total_Area_km2, filtered_results$Species)

lambda_total_area <- phylosig(pruned_tree, total_area, method = "lambda")

lambda_total_area


### range shift distance ###

shift.model <- pglmm(Shift_km ~ Scenario + (1 | Species), 
                     data = shifts_long, 
                     cov_ranef = list(Species = my_tree),
                     family = "gaussian", REML = FALSE)
shift.model


# Phylogenetic signal
filtered_shifts <- shifts_long[shifts_long$Species %in% pruned_tree$tip.label, ]
range_shift <- setNames(filtered_shifts$Shift_km, filtered_shifts$Species)

lambda_range_shift <- phylosig(pruned_tree, range_shift, method = "lambda")

lambda_range_shift



### Plot it! ###



# Function to create box plots for different variables
plot_boxplot <- function(data, variable, y_label) {
  ggplot(data, aes(x = species, y = .data[[variable]], fill = scenario)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = c("current" = "#F9E79F", 
                                 "plus2" = "#F39C12", 
                                 "plus4" = "#E74C3C")) +
    theme_bw() +
    theme(panel.grid = element_blank(),  
          strip.text = element_text(size = 10),  
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = y_label, x = "Species") +
    theme(legend.position = "none")
}

# Generate plots for each variable
p1 <- plot_boxplot(presence_data, "dehydrated.10", "Dehydration Level")
p2 <- plot_boxplot(presence_data, "total_activity_hours", "Total Activity Hours")
p3 <- plot_boxplot(presence_data, "ctmax", "CTmax")
p4 <- plot_boxplot(presence_data, "energ_comsump", "Energy Consumption")

# Display plots
print(p1)
print(p2)
print(p3)
print(p4)



#ex <- merged_data %>%
#  group_by(species, scenario) %>%
#  summarize(mean = mean(dehydrated.10), .groups = "drop") %>%
#  pivot_wider(names_from = scenario, values_from = mean)

#ex <- as.data.frame(ex)
#rownames(ex) <- ex$species

#tree_data <- treedata(my_tree, ex)$phy

#tree_data <- chronos(tree_data)


# Remove species column and convert data to matrix
#ex_matrix <- as.matrix(ex[, -1])  # Remove the first column (species)

# Ensure the rownames of the matrix are the species names
#rownames(ex_matrix) <- ex$species

# Function to create tree plots for each variable
plot_tree_with_heatmap <- function(data, tree, variable, color_palette, color_name, limits = NULL, reverse = FALSE) {
  ex <- data %>%
    group_by(species, scenario) %>%
    summarize(mean = mean(.data[[variable]]), .groups = "drop") %>%
    pivot_wider(names_from = scenario, values_from = mean)
  
  ex <- as.data.frame(ex)
  rownames(ex) <- ex$species
  
  tree_data <- treedata(tree, ex)$phy
  tree_data <- chronos(tree_data)  # Converts to ultrametric tree
  
  ex_matrix <- as.matrix(ex[, -1])  # Remove species column
  rownames(ex_matrix) <- ex$species
  
  # Create tree
  p <- ggtree(tree_data) + geom_tiplab(size=3)
  
  # Define color scale
  color_values <- brewer.pal(5, color_palette)
  if (reverse) {
    color_values <- rev(color_values)  # Reverse the color scale
  }
  
  # Create heatmap
  gheatmap(p, ex_matrix, offset = 0.5, width = 0.3, font.size = 3,
           colnames_angle = -45, hjust = 0) +
    scale_fill_gradientn(colors = color_values, name = color_name, limits = limits,
                         guide = guide_colorbar(direction = "horizontal")) +
    theme(legend.position = "bottom")
}


# Plot each tree with the respective color scale
h1 <- plot_tree_with_heatmap(merged_data, my_tree, "dehydrated.10", "Blues", "Dehydration", limits = c(0, 9000), reverse = TRUE)
h2 <- plot_tree_with_heatmap(merged_data, my_tree, "total_activity_hours", "Oranges", "Activity Hours")
h3 <- plot_tree_with_heatmap(merged_data, my_tree, "ctmax", "Reds", "CTmax")
h4 <- plot_tree_with_heatmap(merged_data, my_tree, "energ_comsump", "Greens", "Energy Consumption")

# Display plots
print(h1)
print(h2)
print(h3)
print(h4)

