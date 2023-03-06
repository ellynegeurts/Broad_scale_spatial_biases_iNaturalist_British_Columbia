#### Set up ####

# Load libraries
library(bcmaps)
library(dplyr)
library(dismo)
library(ENMeval)
library(ENMTools)
library(ggplot2)
library(ggnewscale)
library(ggsn) #allows customize scale bar and arrow for ggplot2
library(ggspatial)
library(here)
library(legendMap)
library(parallel) # Optional, ENMeval does allow for parallel processing, but the model sometimes fails
library(raster)
library(rJava) # Needed to run "maxent.jar" algorithm in ENMeval::evaluate
library(sf)
library(terra)

# Set the JVM heap size in rJava's options:
options(java.parameters = "-Xmx60g") # Allows up to 60 gb of ram usage
# This was needed to run ENMeval later on on my personal computer


#### Read in raster layers that match resolution of MODIS (277 m)####

# All layers are masked in the same projection - BC Albers Projection

bec_277 <- raster(here("Data", "rasters", "bec_277")) # Biogeoclimatic ecosystem classification variable - categorical - 16 levels
park_277 <- raster(here("Data", "rasters", "park_277")) # National and provincial parks - categorical - binary (park vs. not park land)
modis_277 <- raster(here("Data", "rasters", "modis_277")) # MODIS Land Cover - categorical
rd_277 <- raster(here("Data", "rasters", "rd_277")) # Distance to roads - continuous
dem_277 <- raster(here("Data", "rasters", "dem_277")) # Digital elevation model - continuous 
pop_277 <- raster(here("Data", "rasters", "pop_277")) # Human population density - continuous

#### Raster prep ####

# Logging human population density raster
pop_log <- pop_277
values(pop_log) <- log(getValues(pop_277) + 1)

# Create rasterstack
stack <- raster::stack(bec_277, park_277, modis_277, rd_277, dem_277, pop_log)

# Produces a heat map to check for collinearity 
raster.cor.plot(stack)$cor.heatmap

# Gives collinearity values
ras_cor_matrix <- raster.cor.matrix(stack)


# Make nominal (categorical) layers as factor for ENMeval-Maxent models
stack$layer.1 <- raster::as.factor(stack$layer.1) # BEC variable
stack$layer.2 <- raster::as.factor(stack$layer.2) # Park variable
stack$MCD12Q1_LC1_2020_001 <- raster::as.factor(stack$MCD12Q1_LC1_2020_001) # MODIS variable


#### Reading csv files of iNaturalist observations ####

# Read in British Columbia (BC) terrestrial iNaturalist observation. Data downloaded from iNat in Dec. 2021
bc_obs <- read.csv(here("Data", "Prepared_iNat_obs", "terrestrial_obs_only_duplicates_removed.csv"))


## Need to convert df for spatial analysis (an sf object)

# Create a sf object and attach a geo ref system of WGS 1964
bc_obs_sf <- st_as_sf(bc_obs, coords = c("longitude", "latitude"), crs = 4326) 

# Check the geo ref system
st_crs(bc_obs_sf) 

# Project points to NAD83/BC Albers, matching the raster stack
bc_obs_sf_proj <- st_transform(bc_obs_sf, crs = st_crs(stack)) 

# Checking structure (looking to see that the projection worked)
st_crs(bc_obs_sf_proj)

# Now confirming thaT crs matches raster stack
st_crs(bc_obs_sf_proj) == st_crs(stack) # Should say TRUE

plot(bc_obs_sf_proj$geometry)

#### Selecting one presence point per grid cell ####

# For this part, we followed this vignette:
# https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html

# Removing duplicates (i.e. multiple points per grid cell)
occs_cell <- raster::extract(stack[[1]], bc_obs_sf_proj, cellnumbers = TRUE)
occs_cellDups <- duplicated(occs_cell[, 1])
occs <- bc_obs_sf_proj[!occs_cellDups, ]

class(occs) # sf object
st_crs(occs) == st_crs(stack) # TRUE


# Lastly, make occs a df for ENMeval
occs_df <- occs %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])
occs_df <- as.data.frame(cbind(occs_df$lon, occs_df$lat))


#### Running ENMeval - ENMevaluate ####

set.seed(789)

# Run model - note it will likely take many hours
eval_LQHP <- ENMevaluate(
  occ = occs_df,
  envs = stack,
  bg = NULL,
  n.bg = 300000, 
  tune.args = list(fc = c("L", "Q", "H", "P", 
                          "LQ", 
                          "LQH", "LQP", 
                          "LQHP")), # ENMeval doesn't do threshold (T didn't work)
  partitions = "randomkfold",
  algorithm = "maxent.jar",
  categoricals = c("layer.1",
                   "layer.2",
                   "MCD12Q1_LC1_2020_001"),
  RMvalues = seq(0.5, 2, 0.5)
  # numCores = 4, # Optional to speed process, though sometimes the model fails to run if too many cores used
  # parallel = TRUE,
)


##  The messages we received after running the model for your comparison
# * Found 68243 raster cells that were NA for one or more, but not all, predictor variables. Converting these cells to NA for all predictor variables.
# * Randomly sampling 3e+05 background points ...
# * Removed 4664 occurrence points with NA predictor variable values.

# Examine results for eval_LQHP
results <- eval_LQHP@results

# Get best model based on AICc value
bestmod = which(eval_LQHP@results$AICc == min(eval_LQHP@results$AICc))
best <- eval_LQHP@results[bestmod, ]
best

# Get variable importance values for the top model
var.impt <- eval_LQHP@variable.importance 
var.impt.df <- as.data.frame(var.impt$fc.LQHP_rm.0.5)
# We followed instructions from here:
# https://www.rdocumentation.org/packages/ENMeval/versions/0.2.2/topics/Extract%20percent%20contribution%20and%20permutation%20importance%20from%20a%20Maxent%20model


# Create response curve for top model - note our manuscript used the Maxent GUI to help create the response curves
response(eval_LQHP@models[[3]])


# Calculate niche overlap between models
# Schoener's D
overlap_stats_D <- ENMeval::calc.niche.overlap(eval_LQHP@predictions, overlapStat = "D") # Schoeners D. Took about 20 minutes
overlap_stats_D

# Gives matrix locations of the top 10 and bottom 10 prediction overlap ("Schoeners D)
min10 <- which(overlap_stats_D<=sort(overlap_stats_D)[10], arr.ind = TRUE)
max10 <- which(-overlap_stats_D<=sort(-overlap_stats_D)[10], arr.ind = TRUE)

# Create a table of Schoener's D values
tbl_overlap_D <- as.data.frame(as.table(overlap_stats_D))
tbl_overlap_D
# Lowest prediction overlap was fc.Q_rm.2 / fc.LQHP_rm.0.5 = 0.8183719
# Highest prediction overlap was fc.LQH_rm.1.5 / fc.H_rm.1.5 = 0.9982565

#### Make prediction figure ####

# Make predictions
pr <- predict(stack, eval_LQHP@models[[bestmod]], type = "cloglog") # Took several minutes
pr_df <- as.data.frame(pr, xy = T)

# Get neighbouring provinces and country for prediction map
neighbour <- bcmaps::bc_neighbours()

# Create figure
par(mar = c(5, 5, 5, 5))
ggplot(data = neighbour) +
  geom_sf(mapping = aes(fill = iso_a2), show.legend = F) +
  coord_sf() + 
  scale_fill_manual(values = c("azure3", "azure1", "azure2")) + # Canada, ocean, US
  ggspatial::annotation_scale(location = "bl") +
  theme_bw() + 
  north(neighbour, symbol = 10) +
  new_scale_fill() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = layer)) +
  theme_classic() + 
  scale_fill_gradientn(colours = viridis::viridis(99),
                       na.value = "transparent") +
  labs(x = "Longitude", y = "Latitude", fill = "Probability")  +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0)) +
  theme(axis.text = element_text(colour = "black"))


#### Run null model ####

# Run model
null_LQHP_0.5 <- ENMnulls(eval_LQHP, 
                          mod.settings = list(fc = "LQHP", rm = 0.5), 
                          no.iter = 100,
                          parallel = TRUE,
                          numCores = 4
)

# Check results
null_LQHP_0.5@null.results
null.emp.results(null_LQHP_0.5) 

# Visualize between null and empirical models
par(mfrow = c(2, 1))

# Histogram
evalplot.nulls(null_LQHP_0.5, stats = "or.10p", plot.type = "histogram")
evalplot.nulls(null_LQHP_0.5, stats = "auc.val", plot.type = "histogram")

# Violin plot
evalplot.nulls(null_LQHP_0.5, stats = c("or.10p", "auc.val"), plot.type = "violin")
evalplot.nulls(null_LQHP_0.5, stats = "auc.val", plot.type = "violin")
evalplot.nulls(null_LQHP_0.5, stats = "or.10p", plot.type = "violin")
