# DATA.R
#
# This script populates data/.
#
# Ben Davies
# July 2020



# INITIALISATION
# --------------

# Load packages
library(dplyr)
library(haven)
library(igraph)
library(readr)
library(sf)
library(tidyr)
library(yaml)

# Import raw commuting flows
raw_flows_list <- list()
raw_flows_list[[1]] <- read_dta('data/_raw/TTWA_flows_01.dta')
raw_flows_list[[2]] <- read_dta('data/_raw/TTWA_flows_06.dta')
raw_flows_list[[3]] <- read_dta('data/_raw/TTWA_flows_13.dta')

# Identify census years
years <- c(2001, 2006, 2013)
for (i in 1 : length(years)) raw_flows_list[[i]]$year <- years[i]

# Import lookup data
au_boundaries <- read_sf('data/_raw', 'AU2013_GV_Full')
au_boundaries_clipped <- read_sf('data/_raw', 'AU2013_GV_Clipped')
credentials <- read_yaml('data/_credentials.yaml')
raw_geographies <- read_csv('data/_raw/2013_Areas_Table.txt')

# Create concordance map between area units and regional councils
geographies <- raw_geographies %>%
  select(starts_with('AU'), starts_with('REGC'), starts_with('TA2013')) %>%
  `colnames<-`(c('au13', 'au13_desc', 'rc13', 'rc13_desc', 'ta13', 'ta13_desc')) %>%
  mutate_at(c('au13', 'rc13', 'ta13'), as.numeric) %>%
  distinct() %>%
  arrange(au13)

# Compute area unit centroids
au_centroids <- au_boundaries_clipped %>%
  st_centroid() %>%
  mutate(au13 = as.numeric(AU2013),
         long = sapply(geometry, function(x) x[1]),
         lat = sapply(geometry, function(x) x[2])) %>%
  as.data.frame() %>%
  select(au13, long, lat) %>%
  as_tibble()

# Identify area units with zero land areas
landless_aus <- au_boundaries %>%
  as.data.frame() %>%
  as_tibble() %>%
  select(AU2013, LAND_AREA_) %>%
  filter(LAND_AREA_ == 0) %>%
  {as.numeric(.$AU2013)}

# Define sample filters
data <- raw_flows_list %>%
  bind_rows() %>%
  left_join(au_centroids, by = c('ur_au13' = 'au13')) %>%
  left_join(au_centroids, by = c('wp_au13' = 'au13')) %>%
  mutate(distance_km = sqrt((long.y - long.x) ^ 2 + (lat.y - lat.x) ^ 2) / 1e3) %>%
  select(year, ur_au13, wp_au13, flow, distance_km) %>%
  left_join(geographies, by = c('ur_au13' = 'au13')) %>%
  left_join(geographies, by = c('wp_au13' = 'au13')) %>%
  mutate(ur_is_unknown = !ur_au13 %in% geographies$au13,
         wp_is_unknown = !wp_au13 %in% geographies$au13,
         ur_is_outside_ta = ta13.x == 999 | is.na(ta13.x),
         wp_is_outside_ta = ta13.y == 999 | is.na(ta13.y),
         ur_is_landless = ur_au13 %in% landless_aus,
         wp_is_landless = wp_au13 %in% landless_aus,
         has_bad_au13 = ur_is_unknown | wp_is_unknown | ur_is_outside_ta | wp_is_outside_ta | ur_is_landless | wp_is_landless,
         is_long_commute = distance_km >= 150)

# Check sample identification
stopifnot(sum(!data$has_bad_au13 & is.na(data$is_long_commute)) == 0)



# COMMUNITY DETECTION
# -------------------

# Define function for detecting ensemble of communities
get_ensemble <- function(flow_mat, n_runs = 250) {
  
  # Define network
  adj_mat <- flow_mat + t(flow_mat)
  diag(adj_mat) <- diag(adj_mat) / 2  # Undo self-flow double-counting
  net <- graph.adjacency(adj_mat, mode = 'undirected', weighted = TRUE)
  
  # Get ensemble of community allocations
  get_allocations <- function(run) {
    random_net <- permute(net, sample(gorder(net)))  # Randomise vertex IDs
    clust <- cluster_louvain(random_net)
    tibble(au13 = as.numeric(V(random_net)$name),
           community = clust$membership,
           run = rep(run, gorder(net)))
  }
  set.seed(credentials$random_seed)
  bind_rows(lapply(1 : n_runs, get_allocations))
}

# Define function for identifying modal signatures (see Adam et al., 2018)
get_modal_signatures <- function(ensemble) {
  ensemble %>%
    # Identify unique signatures
    group_by(au13) %>%
    mutate(community_seq = paste(community, collapse = '.')) %>%
    ungroup() %>%
    mutate(signature = group_indices(., community_seq)) %>%
    # Compute signature frequencies
    group_by(signature) %>%
    mutate(signature_freq = n() / max(ensemble$run)) %>%
    # Allocate communities to most frequent signatures
    group_by(community, run) %>%
    slice(which.max(signature_freq)) %>%
    ungroup() %>%
    select(community, run, signature) %>%
    right_join(ensemble) %>%
    # Identify modal signatures
    group_by(au13, signature) %>%
    summarise(allocation_rate = n() / max(ensemble$run)) %>%
    group_by(au13) %>%
    slice(which.max(allocation_rate)) %>%
    ungroup()
}

# Initialise storage
flow_matrix_list <- list()
main_ensembles_list <- list()
main_modal_signatures_list <- list()
sub_ensembles_list <- list()
sub_modal_signatures_list <- list()

# Iterate over census years
i <- 1  # Initialise sub-counter
for (t in 1 : length(years)) {
  
  # Identify flows
  flows <- data %>%
    filter(!(has_bad_au13 | is_long_commute),
           year == years[t])
  stopifnot(min(flows$flow) > 0)  # Check data structure
  
  # Construct flow matrix between area units with at least one resident or employee
  nonzero_au13 <- sort(unique(c(flows$ur_au13, flows$wp_au13))) %>%
    intersect(geographies$au13)
  flow_mat_data <- flows %>%
    filter(ur_au13 %in% nonzero_au13 & wp_au13 %in% nonzero_au13) %>%
    mutate(row = sapply(ur_au13, function(x){which(nonzero_au13 == x)}),
           col = sapply(wp_au13, function(x){which(nonzero_au13 == x)})) %>%
    select(ur_au13, wp_au13, row, col, flow) %>%
    as.matrix()
  flow_mat <- matrix(0, length(nonzero_au13), length(nonzero_au13))
  rownames(flow_mat) <- nonzero_au13
  colnames(flow_mat) <- nonzero_au13
  flow_mat[flow_mat_data[, 3:4]] <- flow_mat_data[, 5]
  flow_matrix_list[[t]] <- flow_mat
  
  # Identify main community allocations and modal signatures
  main_ensembles_list[[t]] <- get_ensemble(flow_mat)
  main_modal_signatures_list[[t]] <- get_modal_signatures(main_ensembles_list[[t]])
  
  # Iterate over modal signatures
  modal_signatures <- sort(unique(main_modal_signatures_list[[t]]$signature))
  for (j in 1 : length(modal_signatures)) {
    idx <- which(main_modal_signatures_list[[t]]$signature == modal_signatures[j])
    if (length(idx) > 1) {  # Proceed if top-level community contains at least two area units
      
      # Identify sub-LMA community allocations and modal signatures
      sub_flow_mat <- flow_mat[idx, idx]
      sub_ensembles_list[[i]] <- get_ensemble(sub_flow_mat)
      sub_modal_signatures_list[[i]] <- get_modal_signatures(sub_ensembles_list[[i]])
      
      # Identify census years
      sub_ensembles_list[[i]]$year <- years[t]
      sub_modal_signatures_list[[i]]$year <- years[t]
      
      # Increment counter
      i <- i + 1
    }
  }
  
  # Identify census years
  main_ensembles_list[[t]]$year <- years[t]
  main_modal_signatures_list[[t]]$year <- years[t]
}

# Collate data
main_ensembles <- bind_rows(main_ensembles_list)
main_modal_signatures <- bind_rows(main_modal_signatures_list)
sub_ensembles <- bind_rows(sub_ensembles_list)
sub_modal_signatures <- bind_rows(sub_modal_signatures_list)



# DATA PROCESSING
# ---------------

# Aggregate flows by year and distance threshold
distance_thresholds <- data %>%
  filter(!has_bad_au13) %>%
  mutate(distance_km_threshold = 50 * floor(distance_km / 50) + 50) %>%
  count(year, distance_km_threshold, wt = flow) %>%
  rename(flow_good_au13 = n)

# Aggregate total and sample flows by resident regional council and year
regional_coverage <- data %>%
  mutate(is_in_sample = !(has_bad_au13 | is_long_commute)) %>%
  left_join(geographies, by = c('ur_au13' = 'au13')) %>%
  group_by(year, rc13) %>%
  summarise(flow_total = sum(flow),
            flow_sample = sum(is_in_sample * flow))

# Aggregate flows by year and sample inclusion criterion
national_coverage <- data %>%
  group_by(year) %>%
  summarise(flow_total = sum(flow),
            flow_sample = sum(flow * !(has_bad_au13 | is_long_commute)),
            flow_bad_au13 = sum(flow * has_bad_au13),
            flow_long_commute = sum(flow * is_long_commute, na.rm = TRUE))

# Compute total area unit resident populations
au_residents <- data %>%
  count(ur_au13, year, wt = flow) %>%
  rename(au13 = ur_au13,
         employed_residents_total = n)

# Compute area unit populations
get_au_populations <- function(t) {
  flow_mat <- flow_matrix_list[[t]]
  tibble(au13 = as.numeric(rownames(flow_mat)),
         employed_residents = rowSums(flow_mat),
         employees = colSums(flow_mat),
         resident_employees = diag(flow_mat)) %>%
    mutate(year = years[t])
}
au_populations <- bind_rows(lapply(1 : length(years), get_au_populations)) %>%
  left_join(au_residents) %>%
  select(au13, year, everything()) %>%
  mutate(employed_residents_total = if_else(is.na(employed_residents_total), 0, employed_residents_total))

# Collate LMA and sub-LMA allocations
allocations <- main_modal_signatures %>%
  group_by(year, signature) %>%
  mutate(min_au13 = min(au13)) %>%
  group_by(year) %>%
  mutate(lma = dense_rank(min_au13)) %>%
  ungroup() %>%
  select(au13, year, lma, lma_stability = allocation_rate) %>%
  left_join(sub_modal_signatures) %>%
  group_by(year, lma, signature) %>%
  mutate(min_au13 = min(au13)) %>%
  group_by(year, lma) %>%
  mutate(sub_lma = dense_rank(min_au13)) %>%
  ungroup() %>%
  select(au13, year, lma, lma_stability, sub_lma, sub_lma_stability = allocation_rate) %>%
  mutate(sub_lma_stability = ifelse(is.na(sub_lma_stability), 1, sub_lma_stability))  # Occurs iff LMA is singleton

# Aggregate inter-sub-LMA flows, and compute residual inflows and outflows
tmp <- data %>%
  filter(!(has_bad_au13 | is_long_commute)) %>%
  left_join(allocations, by = c('year', 'ur_au13' = 'au13')) %>%
  left_join(allocations, by = c('year', 'wp_au13' = 'au13')) %>%
  `names<-`(gsub('(.*?)\\.x', 'ur_\\1', names(.))) %>%
  `names<-`(gsub('(.*?)\\.y', 'wp_\\1', names(.))) %>%
  count(year, ur_lma, ur_sub_lma, wp_lma, wp_sub_lma, wt = flow) %>%
  rename(flow = n)
residual_sub_lma_outflows <- tmp %>%
  filter(flow < 6) %>%  # Will be suppressed
  count(year, ur_lma, ur_sub_lma, wt = flow) %>%
  select(year, lma = ur_lma, sub_lma = ur_sub_lma, residual_outflow = n)
residual_sub_lma_inflows <- tmp %>%
  filter(flow < 6) %>%
  count(year, wp_lma, wp_sub_lma, wt = flow) %>%
  select(year, lma = wp_lma, sub_lma = wp_sub_lma, residual_inflow = n)
residual_sub_lma_flows <- residual_sub_lma_inflows %>%
  left_join(residual_sub_lma_outflows)
inter_sub_lma_flows <- tmp %>%
  filter(flow >= 6)

# Widen ensemble tibbles
main_ensembles_wide <- main_ensembles %>%
  mutate(run = paste0('community_', sprintf('%03d', run))) %>%
  spread(run, community)
sub_ensembles_wide <- sub_ensembles %>%
  mutate(run = paste0('community_', sprintf('%03d', run))) %>%
  spread(run, community)



# CONFIDENTIALISATION
# -------------------

# Define function for applying confidentiality rules
confidentialise <- function(x, base = 3) {
  noise <- numeric(length(x))
  for (i in 1 : length(x)) {
    set.seed(credentials$random_seed)
    noise[i] <- runif(1) * base
  }
  above <- x %% base
  res <- x - above + ifelse(noise < above, base, 0)
  res[x < 6] <- 0
  res
}

# Define function for confidentialising selected columns in tibble
confidentialise_columns <- function(x, col_names) {
  mutate_at(x, col_names, confidentialise)
}

# Create confentialised copies of output data
au_populations_conf <- au_populations %>%
  confidentialise_columns(names(.)[grepl('employ', names(.))])
distance_thresholds_conf <- distance_thresholds %>%
  confidentialise_columns(c('flow_good_au13'))
inter_sub_lma_flows_conf <- inter_sub_lma_flows %>%
  confidentialise_columns(c('flow'))
national_coverage_conf <- national_coverage %>%
  confidentialise_columns(names(.)[grepl('flow_', names(.))])
regional_coverage_conf <- regional_coverage %>%
  confidentialise_columns(names(.)[grepl('flow_', names(.))])
residual_sub_lma_flows_conf <- residual_sub_lma_flows %>%
  confidentialise_columns(names(.)[grepl('residual_', names(.))])



# DATA EXPORT
# -----------

# Define function for exporting data
export_data <- function(x, fn) {
  x %>%
    mutate_all(as.character) %>%  # Convert everything to strings
    write_csv(paste0('data/', fn))
}

# Non-sensitive data
export_data(allocations, 'allocations.csv')
export_data(main_ensembles_wide, 'ensembles-main.csv')
export_data(sub_ensembles_wide, 'ensembles-sub.csv')
export_data(geographies, 'geographies.csv')

# Sensitive data
export_data(au_populations, 'NFR_populations.csv')
export_data(distance_thresholds, 'NFR_distance-thresholds.csv')
export_data(inter_sub_lma_flows, 'NFR_inter-sub-lma-flows.csv')
export_data(national_coverage, 'NFR_national-coverage.csv')
export_data(regional_coverage, 'NFR_regional-coverage.csv')
export_data(residual_sub_lma_flows, 'NFR_residual-sub-lma-flows.csv')

# Confidentialised data
export_data(au_populations_conf, 'populations.csv')
export_data(distance_thresholds_conf, 'distance-thresholds.csv')
export_data(inter_sub_lma_flows_conf, 'inter-sub-lma-flows.csv')
export_data(national_coverage_conf, 'national-coverage.csv')
export_data(regional_coverage_conf, 'regional-coverage.csv')
export_data(residual_sub_lma_flows_conf, 'residual-sub-lma-flows.csv')



# SESSION INFO
# ------------

options(width = 80)
write_lines(capture.output(sessioninfo::session_info()), 'logs/data.log')
