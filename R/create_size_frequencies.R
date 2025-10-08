
create_size_frequencies = function(p, region_groups="default", yr_groups=NULL, 
    regions=NULL,  region_titles=NULL, sizedatadir=NULL, outdir=NULL  
  ) {

  if (is.null(yr_groups)) {
    nyrs = length(p$yrs)
    yr_groups = list()
    for (i in 1:ceiling(nyrs/10)) {
      j = i*10 + c(0:-9)
      j = j[ j %in% 1:nyrs ]
      yr_groups[[ paste0("period", i) ]] = as.character( rev(p$yrs)[ j ] )      
    }
  }
 
  if (region_groups=="default") {
    if (is.null(regions)) regions = c("cfanorth", "cfasouth", "cfa4x")
    if (is.null(region_titles)) region_titles=c(cfanorth="NENS", cfasouth="SENS", cfa4x="4X")
    if (is.null(sizedatadir)) sizedatadir = file.path(p$project.outputdir, "size_structure")
    if (is.null(outdir)) outdir =file.path( p$annual.results, "figures", "size.freq", "survey" )
  }

  if (region_groups=="split") {
    # split 23 and 24
    if (is.null(regions)) regions = c("cfanorth", "cfa23", "cfa24", "cfa4x")
    if (is.null(region_titles)) region_titles=c(cfanorth="NENS", cfa23="CFA23", cfa24="CFA24", cfa4x="4X")
    if (is.null(sizedatadir)) sizedatadir = file.path(p$project.outputdir, "size_structure_split")  # alternate save location
    if (is.null(outdir)) outdir =file.path( p$annual.results, "figures", "size.freq", "survey_split" )
  }

  if ( !exists("span", p) ){
    # note ranges in CW will be log transformed later
    # these spans result in dx=2mm
    p$span = function( sexid) {
        switch(sexid,
            male   = c( 5, 155, 50),  #dx=3
            female = c( 5, 95,  30)  #dx=3
        )
    } 
  }

  M = size_distributions(p=p, toget="rawdata", regions=regions, outdir=sizedatadir, redo=TRUE)  # merge det and set with QA/QC
        
  for (yg in 1:length(yr_groups)) {

    years = yr_groups[[yg]]

    M = size_distributions(p=p, toget="crude", Y=years, regions=regions, outdir=sizedatadir, redo=TRUE )  

    for (sx in c("female", "male")) {

      if (sx=="female") sxcol = c("darkorange", "gray95" ) 
      if (sx=="male")   sxcol = c("slategray", "gray95" ) 

      for (yvar in c("den", "denl")) {
        # den = arithmetic mean density,  denl = geometric mean density  

        xd = unique(diff( discretize_data( span=p$span(sx) ) ))[1]
        
        plot_histogram_carapace_width( 
          M=M, 
          years=years, 
          regions=regions, 
          region_titles=region_titles,
          plot_sex=sx,
          Mdelta=xd,
          rdelta = 3,
          yvar=yvar, 
          cols = sxcol,
          outdir=file.path(outdir, names(yr_groups)[yg]) 
        )

      }
    }
  }

  return("done")

}
