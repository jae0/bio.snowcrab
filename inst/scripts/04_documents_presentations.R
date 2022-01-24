
# Calls to create the documents to be run from within R:

year.assessment = 2021

require(rmarkdown)


# --------------
# Snow crab Science Advice for CSAS 

    # change to where the rmd file lives
    setwd( project.codedirectory("SCReports", "inst", "SAR" ) )   # NOTE: consider moving SCRports to inside bio.snowcrab

    rmarkdown::render( 
      "snowcrab_sar.rmd", 
      params=list(
        year.assessment = year.assessment, 
        bio.data.dir = data_root  # location of your bio.data   
      ),  
      output_dir = work_root  # change as desired
    )

    # MSWord file is produced in output_dir=work_root


# --------------
# Snow crab Science Reaearch Document for CSAS 

    # NOTE this is just a placeholder for the call .  the document is yet to be created 

    # change to where the rmd file lives
    setwd( project.codedirectory("SCReports", "inst", "RESDOC" ) )  

    rmarkdown::render( 
      "snowcrab_resdoc.rmd", 
      params=list(
        year.assessment = year.assessment, 
        bio.data.dir = data_root  # location of your bio.data   
      ),  
      output_dir = work_root  # change as desired
    )

    # MSWord file is produced in output_dir=work_root



# --------------
# Snow crab Science Advisory presentations (N and S-ENS):





# --------------
# Snow crab Science Advisory presentations:(4X)

    Move items from 11.4X.Advisory.Presentation.R into a markdown and then add the call here




# --------------
# Snow crab Tagging presentation (OTN and spahgetti):
    see 13.otn.tagging.R
    see 13.Emera.RTagging.R
    see Brent tagging projects


