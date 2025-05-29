#' Social Contact Matrices from Prem et al. (2017)
#'
#' These are the contact matrices generated from the paper "Projecting social contact matrices in 152 countries using contact surveys and demographic data" by Prem et al., (2017). Each matrix estimates contact rates across 16 age groups (0–4, 5–9, ..., 75–79) for 152 countries. Data is found in a list where the name of the object in the list is the country/territories ISO3 code.
#'
#' @format A list of 152 numeric 16x16 matrices.
#' @source Prem K, Cook AR, Jit M (2017). PLoS Comput Biol 13(9): e1005697. https://doi.org/10.1371/journal.pcbi.1005697
"contact_matricies"

#' WHO Global Reported Cases of Vaccine-Preventable Diseases (1980–2023)
#'
#' Annual country-level case reports for vaccine-preventable diseases collected by WHO via the Joint Reporting Form.
#'
#' @format A data frame with columns:
#' - iso3: 3-digit ISO country code
#' - name: Country or territory name
#' - year: Year of report
#' - disease_short: Short disease name
#' - disease_description: Full description of disease
#' - cases: Number of reported cases
#' @source WHO Immunization Dashboard: https://immunizationdata.who.int/
"WHO_disease_reports"

#' UN World Population Prospects 2024: Fertility Rates
#'
#' Estimated Total Fertility Rate (live births per 1,000 women) by age (15–49) and country, from the 2024 UN WPP.
#'
#' @format A data frame with columns:
#' - area: Country or territory
#' - iso3: 3-digit ISO code
#' - year: Year
#' - x15 to x49: Live births per 1,000 women, by age
#' @source UN DESA Population Division, World Population Prospects 2024
"UN_WPP_fertility"

#' UN World Population Prospects 2024: Crude Death Rates
#'
#' Crude Death Rates (deaths per 1,000 population) by single year age bands 0–100, from the 2024 UN WPP.
#'
#' @format A data frame with columns:
#' - area: Country or territory
#' - iso3: 3-digit ISO code
#' - year: Year
#' - x0 to x100: Crude death rates
#' @source UN DESA Population Division, World Population Prospects 2024
"UN_WPP_mortality"

#' UN World Population Prospects 2024: Net Migration Rate
#'
#' Net migration rate per 1,000 population, with no age disaggregation, from the 2024 UN WPP.
#'
#' @format A data frame with columns:
#' - area: Country or territory
#' - iso3: 3-digit ISO code
#' - year: Year
#' - migration_rate_1000: Net migration rate
#' @source UN DESA Population Division, World Population Prospects 2024
"UN_WPP_migration"

#' UN World Population Prospects 2024: Total Population
#'
#' Total population (in thousands) by single year of age (0–100), as of 1 July, from the 2024 UN WPP.
#'
#' @format A data frame with columns:
#' - area: Country or territory
#' - iso3: 3-digit ISO code
#' - year: Year
#' - x0 to x100: Total population counts
#' @source UN DESA Population Division, World Population Prospects 2024
"UN_WPP_population_total"

#' UN World Population Prospects 2024: Female Population
#'
#' Female population (in thousands) by single year of age (0–100), as of 1 July, from the 2024 UN WPP.
#'
#' @format A data frame with columns:
#' - area: Country or territory
#' - iso3: 3-digit ISO code
#' - year: Year
#' - x0 to x100: Female population counts
#' @source UN DESA Population Division, World Population Prospects 2024
"UN_WPP_population_female"

#' Supplementary Immunization Activities (SIA) Coverage Estimates
#'
#' Coverage estimates for SIAs from WHO and supplementary sources, as used in Shattock et al. (2024).
#'
#' @format A data frame with columns:
#' - area: Country
#' - year: Year
#' - age: Age
#' - coverage: Coverage
#' - disease: Disease
#' - vaccine: Vaccine
#' @source Shattock et al., 2024, The Lancet. https://doi.org/10.1016/S0140-6736(24)00850-X
"VIMC_vaccination_sia"

#' WHO Routine Infant Immunization Coverage (WUENIC)
#'
#' Routine immunization coverage estimates (WUENIC) by country, vaccine, dose, and year.
#'
#' @format A data frame with columns:
#' - iso3: 3-digit ISO code
#' - area: Country or territory
#' - vaccine: Vaccine name
#' - vaccine_description: Description of vaccine
#' - coverage_category: Type of coverage
#' - coverage_description: Coverage description
#' - coverage: Coverage percent
#' - dose_order: Order of dose
#' - disease: Disease
#' @source WHO Immunization Dashboard: https://immunizationdata.who.int/global
"WHO_vaccination_routine"

#' WHO Reported Vaccination Schedules
#'
#' Current vaccination schedules for children, adolescents, and adults as reported to the WHO.
#'
#' @format A data frame with columns:
#' - iso3: 3-digit ISO code
#' - area: Country or territory
#' - WHO_region: WHO region
#' - year: Year
#' - vaccine_code: Vaccine identifier code
#' - vaccine_description: Vaccine description
#' - schedulerounds: Schedule rounds
#' - target_pop: Target population
#' - target_pop_description: Description of target group
#' - geoarea: Geographic area
#' - age_administered: Age at administration
#' - comment: Comments
#' @source WHO Immunization Dashboard: https://immunizationdata.who.int/global
"WHO_vaccination_schedule"

#' Historical Vaccination Introduction and Coverage (Pre-1980)
#'
#' Estimated vaccine introduction year, starting coverage, and 1970 coverage for diseases, income groups, and WHO regions.
#'
#' @format A data frame with columns:
#' - disease: Disease name
#' - income_group: World Bank income group
#' - who_region: WHO region
#' - introduction_year: Year of introduction
#' - starting_coverage_percent: Starting coverage (%)
#' @source Historical estimates compiled from various national sources
"vaccination_pre1980"

#' Vaccine Efficacy and Duration Parameters
#'
#' Dose-specific efficacy and duration estimates for measles, diphtheria, and pertussis.
#'
#' @format A data frame with columns:
#' - vaccine: Vaccine name
#' - disease: Disease name
#' - parameter: Parameter name
#' - value: Value of parameter
#' - unit: Unit of measurement
#' - dose: Dose number
#' @source Literature estimates from peer-reviewed sources
"vaccine_parameters"

#' Disease Natural History Parameters
#'
#' Epidemiological parameters such as R0, incubation period, infectious period, and duration of immunity for key diseases.
#'
#' @format A data frame with columns:
#' - disease: Disease name
#' - parameter: Parameter name
#' - value: Value of parameter
#' - unit: Unit
#' - comment: Notes or source details
#' @source Literature estimates from peer-reviewed sources
"disease_parameters"
