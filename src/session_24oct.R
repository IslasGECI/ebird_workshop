library(dplyr)
library(ebirdst)
library(geodata)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(terra)


# ├ Introduction to the eBird Trends Data Products ----

# ebirdst runs that have trends results
trends_runs <- filter(ebirdst_runs, has_trends)
# download trends for Sage Thrasher
ebirdst_download_trends("Sage Thrasher")
# load trends into R session
trends <- load_trends("Sage Thrasher")

# convert percent per year trends from data frame format to raster
trends_ppy_raster <- rasterize_trends(trends)
# convert cumulative trends from data frame format to raster
trends_cumulative_raster <- rasterize_trends(
        trends,
        c(
                "abd_trend",
                "abd_trend_lower",
                "abd_trend_upper"
        )
)
# save to geotiff for use in QGIS or ArcGIS
writeRaster(trends_ppy_raster, "sagthr_abd-ppy_2012-2022.tif")

# convert to spatial points for use with sf package
trends_sf <- st_as_sf(trends,
        coords = c("longitude", "latitude"),
        crs = 4326
)
# save to GeoPackage for use in QGIS or ArcGIS, could also save to shapefile
write_sf(trends_sf, "sagthr_trends_2012-2022.gpkg", layer = "sagethr_trends")


# conversion from annual to cumulative trend
ppy_trend <- 2
start_year <- 2012
end_year <- 2022
cumulative_trend <- 100 * ((1 + ppy_trend / 100)^(end_year - start_year) - 1)

# ├ Application 1: regional trends ----

# Goal: estimate the % per year trend for Sage Thrasher in the Great Basin Bird
# Conservation Region (BCR 9). BCRs are are ecologically distinct regions in
# North America with similar bird communities, habitats, and resource management
# issues. Polygons defining the BCR boundaries can be downloaded from
# https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions

# load BCR polygons that you've downloaded and unzipped, then subset to BCR 9
# https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions
bcr <- read_sf("src/maps/BCR_Terrestrial_master_International.shp") |>
        filter(BCR == 9) |>
        st_transform(crs = 4326)

# load trend estimates for Sage Thrasher
trends <- load_trends("Sage Thrasher")
# convert to spatial sf format
trends_sf <- st_as_sf(trends, coords = c("longitude", "latitude"), crs = 4326)
# subset to just cells within BCR 9
trends_bcr <- trends_sf[bcr, ]
# calculate abundance-weighted regional trend
bcr_trend_ppt <- sum(trends_bcr$abd * trends_bcr$abd_ppy) / sum(trends_bcr$abd)

# ├ Application 2: multi-region trends with uncertainty ----

# Goal: estimate the % per year trend with 80% confidence limits for Sage
# Thrasher for each state in the contiguous United States.

# polygon boundaries of each state in the contiguous US
states <- gadm(country = "USA", level = 1, path = tempdir()) |>
  st_as_sf() |>
  select(state = ISO_1) |>
  filter(!state %in% c("US-AK", "US-HI"))

# load fold-level trend estimates for Sage Thrasher
trends_folds <- load_trends("Sage Thrasher", fold_estimates = TRUE)
# convert to spatial sf format
trends_folds_sf <- st_as_sf(trends_folds,
                            coords = c("longitude", "latitude"),
                            crs = 4326)
# attach state to the fold-level trends data
trends_folds <- st_join(trends_folds_sf, states, left = FALSE) |>
  st_drop_geometry()

# abundance-weighted average trend by region and fold
trends_states_folds <- trends_folds |>
  group_by(state, fold) |>
  summarize(abd_ppy = sum(abd * abd_ppy) / sum(abd),
            .groups = "drop")

# summarize across folds for each state
trends_states <- trends_states_folds |>
  group_by(state) |>
  summarise(abd_ppy_median = median(abd_ppy, na.rm = TRUE),
            abd_ppy_lower = quantile(abd_ppy, 0.10, na.rm = TRUE),
            abd_ppy_upper = quantile(abd_ppy, 0.90, na.rm = TRUE),
            .groups = "drop") |>
  arrange(abd_ppy_median)

# join trends to state polygons and make a map
trends_states_sf <- left_join(states, trends_states, by = "state")
ggplot(trends_states_sf) +
  geom_sf(aes(fill = abd_ppy_median)) +
  scale_fill_distiller(palette = "Reds",
                       limits = c(NA, 0),
                       na.value = "grey80") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
  labs(title = "Sage Thrasher state-level breeding trends 2012-2022",
       fill = "Relative abundance trend [% change / year]") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("trends_by_state.png")

# ├ Application 3: multi-species trends ----

# Goal: estimate the mean trend for a set three species representing the
# sagebrush bird community.

# species list
sagebrush_species <- c("Brewer's Sparrow", "Sagebrush Sparrow", "Sage Thrasher")
# ensure that all species are for the same region and season
filter(trends_runs, common_name %in% sagebrush_species)

# download trends for these species and load them into R
ebirdst_download_trends(sagebrush_species)
trends_sagebrush_species <- load_trends(sagebrush_species)

# calculate cell-wise mean trend
trends_sagebrush <- trends_sagebrush_species |>
  group_by(srd_id, latitude, longitude) |>
  summarize(n_species = n(),
            abd_ppy = mean(abd_ppy, na.rm = TRUE),
            .groups = "drop")

# convert the points to sf format
# only consider cells where all three species occur
all_species <- trends_sagebrush |>
  filter(n_species == length(sagebrush_species)) |>
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326)

# make a map
ggplot(all_species) +
  geom_sf(aes(color = abd_ppy), size = 2) +
  scale_color_gradient2(low = "#CB181D", high = "#2171B5",
                        limits = c(-4, 4),
                        oob = scales::oob_squish) +
  guides(color = guide_colorbar(title.position = "left", barheight = 15)) +
  labs(title = "Sagebrush species breeding trends (2012-2022)",
       color = "Relative abundance trend [% change / year]") +
  theme_bw() +
  theme(legend.title = element_text(angle = 90))
ggsave("3-sage-trends.png")
