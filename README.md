# PREVAIL

**PR**ojection of **E**pidemics and **VA**ccination **I**mpact under **L**apses in coverage

# umbrella <img src="img/PREVAIL_logo_upd.png" align="right" width=10% height=10% />

PREVAIL is an R package designed to support dynamic, adaptable, and comprehensive transmission modeling for vaccine-preventable diseases (VPDs), particularly in settings where routine immunization is disrupted or in decline. This work is designed to inform humanitarian response strategies by providing robust projections of disease outbreaks and intervention impact. The development of the package is supported by Community Jameel and builds on prior work funded by the FCDO for mortality rate modeling in Gaza.

# Key Features

- Stochastic: Stochastic models are essential in crisis settings where population sizes, vaccination coverage, and exposure risks are often variable and unpredictable.
	- Reflect realistic outbreak variability in small or displaced populations.
	- Capture the probability of rare but high-impact events, such as a sudden influx of infectious individuals.

- Flexible Modeling Framework: Developed using odin2, the package allows for modular construction of dynamic transmission models, facilitating easy adaptation for different pathogens.

- Comprehensive Structure: The models incorporate multiple dimensions such as age groups, vaccination statuses, and risk factors to reflect realistic population dynamics.

- Robust Vaccination Modeling: Models account for routine immunization programs, vaccination campaigns, and waning immunity to reflect real-world immunization scenarios.

- Crisis-Responsive: The model structure supports dynamic introduction of infectious individuals to simulate outbreaks triggered by conflict, displacement, or cessation of vaccination services.
