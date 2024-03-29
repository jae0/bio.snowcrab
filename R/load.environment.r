

  # ----------------------------------------------------------------------------------
  # NOTE "year.assessment" must be changed every year before any other run
  #  It cannot be automatically loaded together with the "load.snowcrab.environment". This is because
  #  running in parallel mode requires overrding some parameters in "p" on occasion which cannot be done cleanly
  #  as "load.snowcrab.environment" is sourced with every initialisation of a  new CPU.
  #  Copying the following into each relevent file is not a solution as it is error prone and  repetitive.
  # ----------------------------------------------------------------------------------

load.environment = function( year.assessment=NULL, libs=NULL, p=NULL, ... ) {

  Sys.setlocale("LC_COLLATE", "C")   # turn off locale-specific sorting,

  if (is.null(p)) p = list()
  if (!is.null(libs)) {
    suppressMessages( RLibrary(libs) ) # in case we are trying to bootstrap into something
    p$libs = c(libs, p$libs)
  }
  if (is.null(year.assessment)) {
    if ( exists("year.assessment", p) ) {
      year.assessment=p$year.assessment
    } else {
      warning( "year.assessment was not set .. assuming it is the current year" )
      year.assessment = lubridate::year(Sys.Date())
    }
  }

  require(aegis)
  require(data.table)
  
  p = bio.snowcrab::snowcrab_parameters( p=p, year.assessment=year.assessment, ... ) 

  return(p)
}

