
#' Find which spatially-smeared cluster different points are in
#'
#' This takes the raster brick that comes out of `tess3Q_map_rasters()`
#' and defines the zones where each particular layer is maximal.  Then,
#' it takes the lat-longs of a bunch of points and identifies which
#' layer at that location has the maximal location, thus effectively
#' telling us which cluster/color that location belongs to.
#' @param TRB the raster brick that comes out of tess3Q_map_rasters.  Typically you
#' will want to name the layers in it before you toss it into this function.
#' @param Locations a tibble.  It needs to have, at least the columns
#' `id`, `Lat`, and `Long`.  It will be returned with a new column
#' called spatial_group which tells which spatial area it is in.
assign_points_to_raster_regions <- function(TRB, Locations) {

  # Now, make a stack that has values for each group, only where
  # that group has maximal values.
  max_groups <- TRB
  maxes <- max(max_groups)
  keep_cells <- max_groups >= maxes
  max_groups[keep_cells <= 0] <- NA
  max_groups[keep_cells > 0] <- 1

  # now, max_groups is a raster brick with 1's corresponding
  # the cells of each layer, and 0's everywhere else.

  # now, get the long-lat locations
  longlats <- Locations %>%
    select(Long, Lat) %>%
    as.data.frame()

  assig_mat <- raster::extract(max_groups, longlats, method = "simple")

  # make sure these is only one non-NA value per row
  stopifnot(all(rowSums(!is.na(assig_mat)) == 1))

  assig_tib <- as_tibble(assig_mat)
  assig_names <- names(assig_tib)

  # now, keep just the cluster each one is in
  bind_cols(
    Locations,
    assig_tib
  ) %>%
    pivot_longer(
      cols = all_of(assig_names),
      names_to = "spatial_group",
      values_to = "sgggg_values"
    ) %>%
    filter(!is.na(sgggg_values)) %>%
    select(-sgggg_values)

}
