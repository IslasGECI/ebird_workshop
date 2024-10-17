library(dplyr)
library(ebirdst)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(terra)



# ├ Application 1: multi-species migration chronology ----

# Goal: plot migration chronologies with uncertainty estimates for the set of
# six grassland species below, comparing changes in the proportion of the
# population in Montana throughout the year.

# species list
grassland_species <- c("Baird's Sparrow",
                       "Bobolink",
                       "Chestnut-collared Longspur",
                       "Sprague's Pipit",
                       "Upland Sandpiper",
                       "Western Meadowlark")

# Montana boundary polygon from Natural Earth
# note: you could use any region here, e.g. a shapefile, read in using read_sf
mt_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Montana") |>
  # transform coordinate system to match the raster data
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# loop over each species, generate a migration chronology for each
chronology <- NULL
for (species in grassland_species) {
  # download the data products for this species
  # only weekly 27km relative abundance, median and confidence limits
  ebirdst_download_status(species,
                          pattern = "abundance_(median|lower|upper)_27km")

  # load the median weekly relative abundance and lower/upper confidence limits
  abd_median <- load_raster(species, resolution = "27km")
  abd_lower <- load_raster(species, metric = "lower", resolution = "27km")
  abd_upper <- load_raster(species, metric = "upper", resolution = "27km")

  # convert from relative abundance to proportion of population
  abd_total <- global(abd_median, fun = sum, na.rm = TRUE)$sum
  pop_median <- abd_median / abd_total
  pop_lower <- abd_lower / abd_total
  pop_upper <- abd_upper / abd_total

  # estimate the proportion of population in Montana with confidence limits
  # note: you could also mask and sum here
  pop_median_region <- extract(pop_median, mt_boundary,
                               fun = sum, na.rm = TRUE, ID = FALSE)
  pop_lower_region <- extract(pop_lower, mt_boundary,
                              fun = sum, na.rm = TRUE, ID = FALSE)
  pop_upper_region <- extract(pop_upper, mt_boundary,
                              fun = sum, na.rm = TRUE, ID = FALSE)

  # convert to a data frame in long format (one row per week)
  chronology <- data.frame(species = species,
                           week = as.Date(names(pop_median_region)),
                           median = as.numeric(pop_median_region),
                           lower = as.numeric(pop_lower_region),
                           upper = as.numeric(pop_upper_region)) |>
    bind_rows(chronology)
}

# plot the migration chronologies for each species
ggplot(chronology) +
  aes(x = week, y = median, color = species, fill = species) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Week",
       y = "% of population in Montana",
       title = "Migration chronologies for grassland birds in Montana",
       color = NULL, fill = NULL) +
  theme(legend.position = "bottom")
ggsave("chronology.png")

# ├ Application 2: areas of importance ----

# Goal: map areas of importance during the breeding season for the set of six
# grassland species in Montana.

# start by producing a map of richness for these species
# produce binary range rasters for each species in Montana
range_mt <- list()
for (species in grassland_species) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species, pattern = "abundance_seasonal_mean_3km")

  # load breeding season relative abundance
  abd <- load_raster(species, period = "seasonal") |>
    subset("breeding")
  # crop and mask to Montana
  abd_masked <- mask(crop(abd, mt_boundary), mt_boundary)
  # convert to binary, presence-absence
  range_mt[[species]] <- abd_masked > 0
}
# sum across species to calculate richness
richness <- sum(rast(range_mt), na.rm = TRUE)
# make a simple map
png(filename="plot_2.png")
plot(richness, axes = FALSE)
dev.off()

# for more granularity, use mean proportion of population
# produce proportion of population layers for each species masked to Montana
prop_pop_mt <- list()
for (species in grassland_species) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species,
                          pattern = "proportion-population_seasonal_mean_3km")

  # load breeding season proportion of population
  prop_pop <- load_raster(species,
                          product = "proportion-population",
                          period = "seasonal") |>
    subset("breeding")
  # crop and mask to Montana
  prop_pop_mt[[species]] <- mask(crop(prop_pop, mt_boundary), mt_boundary)
}
# take mean across species
importance <- mean(rast(prop_pop_mt), na.rm = TRUE)
# drop zeros
importance <- ifel(importance == 0, NA, importance)
# drop anything below the median
cutoff <- global(importance, quantile, probs = 0.5, na.rm = TRUE) |>
  as.numeric()
importance <- ifel(importance > cutoff, importance, NA)
# make a simple map
plot(importance, axes = FALSE)
plot(mt_boundary, col = "grey", axes = FALSE, add = TRUE)
plot(importance, axes = FALSE, legend = FALSE, add = TRUE)

# make a slightly nicer map
# reproject
importance_proj <- trim(project(importance, "ESRI:102003"))
mt_boundary_proj <- project(mt_boundary, "ESRI:102003")
# basemap
par(mar = c(0, 0, 0, 0))
plot(mt_boundary_proj, col = "grey", axes = FALSE,
     main = "Areas of importance for grassland birds in Montana")
# add importance raster
plot(importance_proj, legend = FALSE, add = TRUE)
# add legend
fields::image.plot(zlim = c(0, 1), legend.only = TRUE,
                   col = viridis::viridis(100),
                   breaks = seq(0, 1, length.out = 101),
                   smallplot = c(0.15, 0.85, 0.12, 0.15),
                   horizontal = TRUE,
                   axis.args = list(at = c(0, 0.5, 1),
                                    labels = c("Low", "Medium", "High"),
                                    fg = "black", col.axis = "black",
                                    cex.axis = 0.75, lwd.ticks = 0.5,
                                    padj = -1.5),
                   legend.args = list(text = "Relative Importance",
                                      side = 3, col = "black",
                                      cex = 1, line = 0))


# ├ Exercises ----

# Exercises 1: repeat the migration chronology application demonstrated in the
# webinar, but try plotting the proportion of the population relative to the
# contiguous United States (i.e. all US states except Alaska and Hawaii) rather
# than the proportion of the global population. Hint: this will require cropping
# and masking the relative abundance rasters to a boundary of the unites states,
# which is provided below.

# boundary of the contiguous united states
us_boundary <- ne_states(iso_a2 = "US") |>
  filter(!name %in% c("Alaska", "Hawaii")) |>
  st_union() |>
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# Exercise 2: select a group of 5-10 species and a region of interest to you.
# Generate a map identifying areas of importance for these species by finding
# the mean proportion of population across the species in your region of
# interest for either the breeding or non-breeding season. Experiment with
# different quantile cutoffs (e.g. median, 70th quantile, 80th quantile) to see
# how that impacts the areas of importance identified.
