
get_sppoly_sims = function(){
  
  setwd("C:/home/jae/projects/bstm/docs/movement") 

  require(aegis)
  require(bio.snowcrab)   # loadfunctions("bio.snowcrab") 
  require(terra)

  year_assessment = 2025
  year_start = 1999

  yrs = year_start:year_assessment

  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
   
  carstm_model_label= paste( "default", snowcrab_filter_class, sep="_" )
 
  # params for probability of observation
  pH = snowcrab_parameters( 
    project_class="carstm", 
    yrs=yrs,  
    areal_units_type="tesselation", 
    carstm_model_label= carstm_model_label,  
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=spec_bio ),
      biologicals_using_snowcrab_filter_class=snowcrab_filter_class
    )
  )
 
  # areal units upon which carstm will operate ... this is made in 01.snowcrab.r
  sppoly = areal_units( p=pH )

  sims = carstm_posterior_simulations( pH=pH, pa_threshold=0.05, qmax=0.95 )  # time sliced at fall
 
  sa = units::drop_units(sppoly$au_sa_km2)
  sppoly$au_sa_km2 = sa

  save( sppoly, sims, file="data/sppoly_hsi.RData" )
  return( list(sppoly, sims))
}
