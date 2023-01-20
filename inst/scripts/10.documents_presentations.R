
# Calls to create the documents to be run from within R:
# change as desired / required
 
# --------------
# SAR: Science Advice Report for CSAS 
 
  year.assessment = 2022
  ID = "snowcrab_sar"

  loc = file.path( Sys.getenv("HOME"), "projects", ID )
  media_loc = file.path( loc, "media" ) 

  fn  = file.path( loc, paste( ID, "_", year.assessment, ".rmd", sep="") )
  pm = list( year.assessment=year.assessment, media_loc=media_loc )
  rmarkdown::render( fn, params=pm, output_dir=dirname(fn) )


# --------------
# SR: Science Report for CSAS 
 
  year.assessment = 2022
  ID = "snowcrab_sr"

  loc = file.path( Sys.getenv("HOME"), "projects", ID )
  media_loc = file.path( loc, "media" ) 

  fn  = file.path( loc, paste( ID, "_", year.assessment, ".rmd", sep="") )
  pm = list( year.assessment=year.assessment, media_loc=media_loc )
  rmarkdown::render( fn, params=pm, output_dir=dirname(fn) )


# --------------
# RESDOC: Snow crab Science Research Document for CSAS 
 
  year.assessment = 2022
  ID = "snowcrab_resdoc"

  loc = file.path( Sys.getenv("HOME"), "projects", ID )
  media_loc = file.path( loc, "media" ) 

  fn  = file.path( loc, paste( ID, "_", year.assessment, ".rmd", sep="") )
  pm = list( year.assessment=year.assessment, media_loc=media_loc )
  rmarkdown::render( fn, params=pm, output_dir=dirname(fn) )


# --------------
# Framework: Snow crab Framework Document for CSAS 
 
  year.assessment = 2022
  ID = "snowcrabframework"

  loc = file.path( Sys.getenv("HOME"), "projects", ID )
  media_loc = file.path( loc, "media" ) 

  fn  = file.path( loc, paste( ID, "_", year.assessment, ".rmd", sep="") )
  pm = list( year.assessment=year.assessment, media_loc=media_loc )
  rmarkdown::render( fn, params=pm, output_dir=dirname(fn) )


# --------------
# Snow crab Science Advisory presentations (N and S-ENS):
 


# --------------
# Snow crab Science Advisory presentations:(4X)

    Move items from 11.4X.Advisory.Presentation.R into a markdown and then add the call here




# --------------
# Snow crab Tagging presentation (OTN and spaghetti):
    see 13.otn.tagging.R
    see 13.Emera.RTagging.R
    see Brent tagging projects


