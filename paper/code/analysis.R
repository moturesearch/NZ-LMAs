# ANALYSIS.R
#
# This script populates figures/ and tables/.
#
# Ben Davies
# July 2020



# INITIALISATION
# --------------

# Load packages
library(dplyr)
library(ggforce)
library(ggplot2)
library(readr)
library(rmapshaper)
library(scales)
library(sf)
library(spdep)
library(tidyr)

# Import data on area unit populations and LMA allocations
import_data <- function(fn, ...) read_csv(paste0('data/', fn), ...)
allocations <- import_data('allocations.csv')
geographies <- import_data('geographies.csv')
au_populations <- import_data('populations.csv')

# Import data on inter-sub-LMA flows
inter_sub_lma_flows    <- import_data('inter-sub-lma-flows.csv')
residual_sub_lma_flows <- import_data('residual-sub-lma-flows.csv')

# Import descriptive data
distance_thresholds <- import_data('distance-thresholds.csv')
national_coverage   <- import_data('national-coverage.csv')
regional_coverage   <- import_data('regional-coverage.csv')

# Import ensembles
ensembles_main <- import_data('ensembles-main.csv')

# Import boundary files
au_boundaries <- read_sf('data/_raw', 'AU2013_GV_Full')
au_boundaries_clipped <- read_sf('data/_raw', 'AU2013_GV_Clipped')

# Define function for generating fill indices so that adjacent polygons have different colours
# See https://stat.ethz.ch/pipermail/r-sig-geo/2011-April/011400.html
get_fill_index <- function (x) {
  neighbour_list <- st_intersects(x, x)
  n <- length(x)
  res <- numeric(n)
  res[1] <- 1
  colour_seq <- 1 : n
  for (j in 2 : n) {
    res[j] <- which.min(colour_seq %in% res[neighbour_list[[j]]])
  }
  return (res)
}



# POPULATION DENSITY MAP
# ----------------------

au_land_areas <- au_boundaries %>%
  as_tibble() %>%
  select(au13 = AU2013, land_area = LAND_AREA_) %>%
  mutate(au13 = as.integer(au13))
map_data <- au_boundaries_clipped %>%
  mutate(au13 = as.integer(AU2013)) %>%
  filter(au13 != 597000) %>%  # Ignore Chatham Islands
  ms_simplify()
tmp <- map_data %>%
  left_join(au_populations) %>%
  left_join(au_land_areas) %>%
  mutate(resident_density = employed_residents_total / land_area) %>%
  group_by(au13) %>%
  summarise(resident_density = mean(resident_density, na.rm = T)) %>%
  ungroup() %>%
  as.data.frame() %>%
  select(-geometry) %>%
  as_tibble()
regional_boundaries <- map_data %>%
  inner_join(tmp) %>%
  filter(resident_density < Inf) %>%
  left_join(geographies) %>%
  group_by(rc13) %>%
  summarise() %>%
  ungroup()
map_data %>%
  inner_join(tmp) %>%
  ggplot() +
  geom_sf(aes(fill = log(resident_density + 0.01)), col = NA, show.legend = F) +
  geom_sf(data = regional_boundaries, fill = NA, lwd = 0.1) +
  coord_sf(datum = NA) +
  scale_fill_distiller(palette = 'Greys', direction = 1, na.value = NA) +
  theme_minimal() +
  theme(legend.justification = c(1, 0),
        legend.position = c(1, 0))
ggsave('figures/01-density.pdf', width = 6, height = 9)



# EXCLUSION RATES
# ---------------

# Employee exclusion rates by criterion and census year
national_coverage %>%
  mutate(year = paste(year, 'Census'),
         share_excluded = 1 - flow_sample / flow_total,
         share_bad_au13 = flow_bad_au13 / flow_total,
         share_long_commute = flow_long_commute / flow_total) %>%
  select(-starts_with('flow')) %>%
  mutate_at(setdiff(names(.), 'year'), function(x) round(100 * x, 2)) %>%
  `names<-`(c('year',
              'Total',
              'Residence or workplace address belongs to excluded area unit',
              'Estimated commute distance beyond 150km')) %>%
  gather(key, value, -year) %>%
  spread(year, value) %>%
  rename(`Exclusion criterion` = key) %>%
  slice(c(2, 1, 3)) %>%
  write_csv('tables/01-exclusion-rates-national.csv')

# Employee exclusion rates by region and census year
regional_coverage %>%
  filter(rc13 < 99) %>%
  left_join(geographies) %>%
  mutate(`Region of origin` = gsub(' Region$', '', rc13_desc),
         year = paste(year, 'Census'),
         excl_rate = 1 - flow_sample / flow_total) %>%
  select(year, `Region of origin`, excl_rate) %>%
  mutate(excl_rate = round(100 * excl_rate, 2)) %>%
  distinct() %>%
  spread(year, excl_rate) %>%
  write_csv('tables/02-exclusion-rates-regional.csv')



# COMMUTE DISTANCE DISTRIBUTIONS
# ------------------------------

distance_thresholds %>%
  mutate(d_band = paste('d <', distance_km_threshold)) %>%
  count(year, distance_km_threshold, d_band, wt = flow_good_au13) %>%
  group_by(year) %>%
  mutate(perc = round(100 * cumsum(n) / sum(n), 2)) %>%
  ungroup() %>%
  filter(distance_km_threshold <= 400) %>%
  select(year, distance_km_threshold, `Distance d (km) band` = d_band, perc) %>%
  mutate(year = paste(year, 'Census')) %>%
  spread(year, perc) %>%
  select(-distance_km_threshold) %>%
  write_csv('tables/03-distance-cumulative-proportions.csv')



# AREA UNIT, LMA, AND SUB-LMA COUNTS
# ----------------------------------

allocations %>%
  mutate(sub_lma = paste0(lma, '.', sub_lma),
         year = paste(year, 'Census')) %>%
  select(-ends_with('stability')) %>%
  `names<-`(c('Number of area units',
              'year',
              'Number of LMAs',
              'Number of sub-LMAs')) %>%
  gather(`Variable`, value, -year) %>%
  group_by(year, Variable) %>%
  summarise(n = n_distinct(value)) %>%
  ungroup() %>%
  mutate(n = comma(n, accuracy = 1)) %>%
  spread(year, n) %>%
  write_csv('tables/04-counts.csv')



# AREA UNIT, LMA, AND SUB-LMA ATTRIBUTE DISTRIBUTIONS
# ---------------------------------------------------

# Estimate LMA resident employees
lma_resident_employees <- inter_sub_lma_flows %>%
  count(year, ur_lma, wp_lma, wt = flow) %>%
  filter(ur_lma == wp_lma) %>%
  select(year, lma = ur_lma, resident_employees = n)

# Compute LMA populations
lma_populations <- au_populations %>%
  left_join(allocations) %>%
  select(-sub_lma, -resident_employees, -ends_with('stability')) %>%
  gather(key, value, -au13, -year, -lma) %>%
  count(lma, year, key, wt = value) %>%
  spread(key, n) %>%
  left_join(lma_resident_employees)

# Estimate sub-LMA resident employees
sub_lma_resident_employees <- inter_sub_lma_flows %>%
  filter(ur_lma == wp_lma & ur_sub_lma == wp_sub_lma) %>%
  select(year, lma = ur_lma, sub_lma = ur_sub_lma, resident_employees = flow)

# Compute sub-LMA populations
sub_lma_populations <- au_populations %>%
  left_join(allocations) %>%
  select(-resident_employees, -ends_with('stability')) %>%
  gather(key, value, -au13, -year, -lma, -sub_lma) %>%
  count(lma, sub_lma, year, key, wt = value) %>%
  spread(key, n) %>%
  left_join(sub_lma_resident_employees)

# Define function for computing population attribute distributions
get_population_attribute_distributions <- function(tbl) {
  tbl %>%
    mutate(`Supply-side self-containment (%)` = 100 * resident_employees / employed_residents,
           `Demand-side self-containment (%)` = 100 * resident_employees / employees) %>%
    mutate_at(c('employed_residents_total', 'employed_residents', 'employees'), function(x) x / 1e3) %>%
    rename(`Employed residents (000s)` = employed_residents,
           `Employed residents (000s, before exclusions)` = employed_residents_total,
           `Employees (000s)` = employees) %>%
    select(-resident_employees) %>%
    gather(key, value, -year) %>%
    group_by(year, key) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate_at(c('mean', 'sd'), function(x) comma(x, 0.01)) %>%
    mutate(desc = sprintf('%s (%s)', mean, sd),
           year = paste(year, 'Census')) %>%
    select(Attribute = key, year, desc) %>%
    spread(year, desc) %>%
    slice(c(3, 2, 4, 5, 1))
}

# Area unit, LMA and sub-LMA attribute distributions by census year
bind_rows(
  mutate(get_population_attribute_distributions(select(au_populations, -au13)), Geography = 'Area units'),
  mutate(get_population_attribute_distributions(select(lma_populations, -lma)), Geography = 'LMAs'),
  mutate(get_population_attribute_distributions(select(sub_lma_populations, -lma, -sub_lma)), Geography = 'Sub-LMAs')
) %>%
  select(Geography, everything()) %>%
  write_csv('tables/05-attribute-distributions.csv')



# LMA BOUNDARIES
# --------------

# Compute LMA boundaries
map_data <- map_data %>% 
  right_join(allocations) %>%
  left_join(geographies) %>%
  mutate(island = ifelse(rc13 < 10, 'north', 'south'))
lma_boundaries <- map_data %>%
  group_by(year, lma, island) %>%
  summarise() %>%
  ungroup()

# LMA boundaries using 2001, 2006 and 2013 Census data
years <- c(2001, 2006, 2013)
for (y in 1 : 3) {
  tmp <- map_data %>%
    filter(year == years[y]) %>%
    left_join(allocations) %>%
    filter(lma_stability < 1)
  lma_boundaries %>%
    filter(year == years[y]) %>%
    ggplot() +
    geom_sf(data = tmp, fill = 'grey80', lwd = NA) +
    geom_sf(fill = NA, colour = 'grey50', lwd = 0.1) +
    geom_sf_text(aes(label = lma)) +
    coord_sf(datum = NA) +
    scale_fill_grey() +
    theme_void()
  ggsave(sprintf('figures/%02d-boundaries-%s.pdf', 1 + y, years[y]), width = 6, height = 9)
}



# LMA ATTRIBUTES
# --------------

# Compute LMA attributes by census year
lma_attributes <- au_land_areas %>%
  left_join(allocations) %>%
  group_by(lma, year) %>%
  summarise(n_au13 = n_distinct(au13),
            land_area = sum(land_area)) %>%
  ungroup() %>%
  left_join(lma_populations) %>%
  mutate(`Supply-side self-containment (%)` = 100 * resident_employees / employed_residents,
         `Demand-side self-containment (%)` = 100 * resident_employees / employees) %>%
  select(-resident_employees) %>%
  mutate_at(c('employed_residents_total', 'employed_residents', 'employees'), as.integer) %>%
  mutate(land_area = as.integer(round(land_area))) %>%
  select(year,
         LMA = lma,
         `Number of area units` = n_au13,
         `Land area (km2)` = land_area,
         `Employed residents (before exclusions)` = employed_residents_total,
         `Employed residents` = employed_residents,
         `Employees` = employees,
         everything())

# Attributes of LMAs identified using 2001, 2006 and 2013 Census data
years <- c(2001, 2006, 2013)
for (i in 1 : 3) {
  lma_attributes %>%
    filter(year == years[i]) %>%
    select(-year) %>%
    mutate(LMA = as.integer(LMA)) %>%
    mutate_if(is.integer, comma, accuracy = 1) %>%
    mutate_if(is.double, function(x) comma(x, accuracy = 0.01)) %>%
    mutate_all(function(x) ifelse(x %in% c('0', NA), '-', x)) %>%
    write_csv(sprintf('tables/%02d-lma-attributes-%d.csv', 5 + i, years[i]))
}



# LMA RECONFIGURATIONS ACROSS CENSUS YEARS
# ----------------------------------------

alluvial_data <- allocations %>%
  select(au13, year, lma) %>%
  add_count(au13) %>%
  filter(n == 3) %>%
  select(-n) %>%
  left_join(au_populations) %>%
  group_by(au13) %>%
  mutate(mean_employed_residents = mean(employed_residents)) %>%
  ungroup() %>%
  select(au13, mean_employed_residents, year, lma) %>%
  spread(year, lma) %>%
  count(`2001`, `2006`, `2013`, wt = mean_employed_residents)

alluvial_data %>%
  gather_set_data(1:3) %>%
  ggplot(aes(x, id = id, split = factor(y), value = n)) +
  geom_parallel_sets(alpha = 0.4, axis.width = 0.1) +
  geom_parallel_sets_axes(col = NA, alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_labels(angle = 0) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Census year',
       y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(panel.grid = element_blank())
ggsave('figures/05-alluvial.pdf', width = 6, height = 9)



# WELLINGTON SUB-LMA BOUNDARIES
# -----------------------------

# Compute sub-LMA boundaries
sub_lma_boundaries <- map_data %>%
  group_by(year, lma, sub_lma) %>%
  summarise() %>%
  ungroup()

# Identify Wellington sub-LMAs
wellington_sub_lmas <- allocations %>%
  left_join(geographies) %>%
  group_by(year, lma) %>%
  filter(sum(au13_desc == 'Thorndon-Tinakori Road') > 0) %>%
  filter(au13 != 549902) %>%  # Exclude Mara
  ungroup() %>%
  distinct(year, lma, sub_lma)

# Sub-LMA boundaries in Wellington using 2001, 2006 and 2013 Census data
sub_lma_boundaries %>%
  right_join(wellington_sub_lmas) %>%
  group_by(year, lma) %>%
  mutate(sub_lma_rank = dense_rank(sub_lma)) %>%
  ungroup() %>%
  ggplot() +
  geom_sf(aes(fill = factor(sub_lma_rank)), alpha = 0.4, lwd = 0.1, show.legend = FALSE) +
  geom_sf_text(aes(label = sub_lma)) +
  coord_sf(datum = NA) +
  facet_wrap(~year) +
  scale_fill_grey() +
  theme_void()
ggsave('figures/06-wellington.pdf', width = 9, height = 3)

# Alluvial plot for Wellington sub-LMAs
au_populations %>%
  group_by(au13) %>%
  summarise(size = mean(employed_residents_total)) %>%
  left_join(allocations) %>%
  inner_join(wellington_sub_lmas) %>%
  add_count(au13) %>%
  filter(n == 3) %>%
  select(au13, size, year, sub_lma) %>%
  spread(year, sub_lma) %>%
  count(`2001`, `2006`, `2013`, wt = size) %>%
  gather_set_data(1:3) %>%
  ggplot(aes(x, id = id, split = factor(y), value = n)) +
  geom_parallel_sets(alpha = 0.4, axis.width = 0.1) +
  geom_parallel_sets_axes(col = NA, alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_labels(angle = 0) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Census year',
       y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(panel.grid = element_blank())
ggsave('figures/07-wellington-alluvial.pdf', width = 6, height = 4)



# COMMUNITY COUNTS AND LMA STABILITIES
# ------------------------------------

# Distribution of communities identified across runs
ensembles_long <- ensembles_main %>%
  gather(key, value, -au13, -year) %>%
  separate(key, c('key', 'run'), sep = '_') %>%
  mutate(run = as.integer(run)) %>%
  spread(key, value)
ensembles_long %>%
  group_by(year, run) %>%
  summarise(n_communities = max(community)) %>%
  group_by(year) %>%
  summarise(Mean = mean(n_communities),
            `Std. dev.` = sd(n_communities),
            `Min` = min(n_communities),
            `Median` = median(n_communities),
            `Max` = max(n_communities)) %>%
  ungroup() %>%
  left_join(count(lma_attributes, year, name = 'LMAs')) %>%
  mutate_if(is.double, round, 2) %>%
  rename(`Census year` = year) %>%
  write_csv('tables/09-community-counts.csv')

# Weighted mean stability of area units' allocations to LMAs in each census year
stabilities <- allocations %>%
  left_join(au_populations) %>%
  select(-lma, -starts_with('sub_lma'), -resident_employees, -employed_residents_total) %>%
  rename(`Resident-weighted` = employed_residents,
         `Employee-weighted` = employees) %>%
  mutate(Unweighted = 1,
         year = paste(year, 'Census')) %>%
  gather(Weights, wt, -au13, -year, -lma_stability)
stabilities %>%
  group_by(year, Weights) %>%
  summarise(`Stability` = round(100 * sum(wt * lma_stability) / sum(wt), 2),
            `Stability = 100%` = round(100 * sum(wt * (lma_stability == 1)) / sum(wt), 2)) %>%
  ungroup() %>%
  gather(Variable, value, -year, -Weights) %>%
  spread(year, value) %>%
  arrange(Variable, Weights) %>%
  select(Variable, everything()) %>%
  write_csv('tables/10-stabilities.csv')



# SESSION INFO
# ------------

options(width = 80)
write_lines(capture.output(sessioninfo::session_info()), 'logs/analysis.log')
