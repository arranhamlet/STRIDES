options(scipen = 999)

if(!require("pacman")) install.packages("pacman")

#Load packages
pacman::p_load(
  rio,
  usethis,
  purrr,
  janitor,
  here,
  rlang,
  tidyverse
)

#Import data
all_rds <- list.files("data", full.names = T, pattern = ".rds")
rds_names <- gsub("data/|.rds", "", all_rds)

rds_list <- sapply(all_rds, function(x) readRDS(x))

for(i in 1:length(rds_list)){
  assign(rds_names[i], rds_list[[i]])
}

#Rename where appropriate
contact_matricies <- contact_matricies
WHO_disease_reports <- full_disease_df
UN_WPP_fertility <- fertility
UN_WPP_mortality <- mortality
UN_WPP_migration <- migration
UN_WPP_population_total <- population_all
UN_WPP_population_female <- population_female
VIMC_vaccination_sia <- sia_vaccination  %>%
  select(area = country,
         year = year,
         age = age,
         coverage = coverage,
         disease = disease,
         vaccine = vaccine)
WHO_vaccination_routine <- routine_vaccination_data %>%
  select(iso3 = CODE,
         area = NAME,
         vaccine = ANTIGEN,
         vaccine_description = ANTIGEN_DESCRIPTION,
         coverage_category = COVERAGE_CATEGORY,
         coverage_description = COVERAGE_CATEGORY_DESCRIPTION,
         coverage = COVERAGE,
         dose_order,
         disease = Disease)
WHO_vaccination_schedule <- vaccination_schedule %>%
  select(iso3 = ISO_3_CODE,
         area = COUNTRYNAME,
         WHO_region = WHO_REGION,
         year = YEAR,
         vaccine_code = VACCINECODE,
         vaccine_description = VACCINE_DESCRIPTION,
         schedulerounds = SCHEDULEROUNDS,
         target_pop = TARGETPOP,
         target_pop_description = TARGETPOP_DESCRIPTION,
         geoarea = GEOAREA,
         age_administered = AGEADMINISTERED,
         comment = SOURCECOMMENT)
vaccination_pre1980 <- vaccination_pre1980 %>%
  janitor::clean_names()
vaccine_parameters <- vaccine_parameters
disease_parameters <- disease_parameters %>%
  rename(comment = note)

#Export data
use_data(contact_matricies, overwrite = TRUE)
use_data(WHO_disease_reports, overwrite = TRUE)
use_data(UN_WPP_fertility, overwrite = TRUE)
use_data(UN_WPP_mortality, overwrite = TRUE)
use_data(UN_WPP_migration, overwrite = TRUE)
use_data(UN_WPP_population_total, overwrite = TRUE)
use_data(UN_WPP_population_female, overwrite = TRUE)
use_data(VIMC_vaccination_sia, overwrite = TRUE)
use_data(WHO_vaccination_routine, overwrite = TRUE)
use_data(WHO_vaccination_schedule, overwrite = TRUE)
use_data(vaccination_pre1980, overwrite = TRUE)
use_data(vaccine_parameters, overwrite = TRUE)
use_data(disease_parameters, overwrite = TRUE)











