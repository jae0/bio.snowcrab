
# NOTE::: this requires 03.snowcrab_carstm.r to be completed
source( file.path( code_root, "bio_startup.R" )  )
loadfunctions("bio.snowcrab")

year.assessment = 2021

# prep data for discrete version
fishery_model_data_inputs( 
    year.assessment=year.assessment, 
    type="biomass_dynamics", 
    for_julia=TRUE   
)

# prep data for continuous version: 
# fishery landings has a weekly time step = 2/52 ~ 0.0385 ~ 0.04  X dt=0.01 seems to work best
fishery_model_data_inputs( 
    year.assessment=year.assessment, 
    type="size_structured_numerical_dynamics", 
    for_julia=TRUE, 
    time_resolution=2/52  
)

# Rdata files are ready load them through julia and model

# using dynamical_models/snowcrab/04.snowcrab_fishery_model.jl


