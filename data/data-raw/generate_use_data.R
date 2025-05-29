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
  data.table,
  squire.page,
  tidyverse
)

#Import data
all_rds <- list.files(path = "data/data-raw/", full.names = T, pattern = ".rds", recursive = T)
rds_names <- gsub("data/data-raw/|.rds", "", all_rds)

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
         vaccine = vaccine) %>%
  subset(!disease %in% c("je", "yf", "tet", "tb"))
WHO_vaccination_routine <- routine_vaccination_data %>%
  select(iso3 = CODE,
         area = NAME,
         year = YEAR,
         vaccine = ANTIGEN,
         vaccine_description = ANTIGEN_DESCRIPTION,
         coverage_category = COVERAGE_CATEGORY,
         coverage_description = COVERAGE_CATEGORY_DESCRIPTION,
         coverage = COVERAGE,
         dose_order,
         disease = Disease) %>%
  subset(
    disease %in% c("Diphtheria (childhood combination vaccine)", "Diphtheria, Tetanus, Pertussis (childhood vaccine)",
                   "Measles") | grepl("mening|Poliomyelitis|Rubella", disease, ignore.case = T)
  )
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
         comment = SOURCECOMMENT) %>%
  subset(vaccine_code %in% c("D_S", "DIP", "DT", "DTAP","DTAPHEPBIPV", "DTAPHIB", "DTAPHIBHEPB", "DTAPHIBHEPBIPV", "DTAPHIBIPV", "DTAPIPV", "DTIPV", "DTWP", "DTWPHEPB", "DTWPHIB", "DTWPHIBHEPB", "DTWPHIBHEPBIPV", "MEASLES",
                             "MEN_A_CONJ", "MEN_A_PS", "MEN_AC_CONJ", "MEN_AC_PS", "MEN_ACYW_135CONJ", "MEN_ACYW_135PS", "MEN_B", "MEN_BC", "MEN_C_CONJ", "MM", "MMR", "MMRV", "MR", "OPV", "PCV_15_VALENT", "PCV10", "PCV13",
                             "PCV20", "PPV23", "RUBELLA", "TD_S", "TDAP_S", "TDAP_S_IPV", "TDIPV_S"))

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











