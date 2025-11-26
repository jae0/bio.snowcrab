
create_size_frequencies = function(p, region_groups="default", yr_groups=NULL, span=NULL,
    regions=NULL,  region_titles=NULL, sizedatadir=NULL, outdir=NULL  
  ) {
  
  if (0) {
     region_groups="default"
     
     span=NULL
     regions=NULL
     region_titles=NULL
     sizedatadir=NULL
     outdir=NULL  
  }


  if (is.null(yr_groups)) {
    nyrs = length(p$yrs)
    yr_groups = list()
    for (i in 1:ceiling(nyrs/10)) {
      j = i*10 + c(0:-9)
      j = j[ j %in% 1:nyrs ]
      yr_groups[[ paste0("period", i) ]] = as.character( rev(p$yrs)[ j ] )      
    }
  }
 
  if (is.null(regions)) {
    regions = switch( region_groups,
      default = c("cfanorth", "cfasouth", "cfa4x"),
      split   = c("cfanorth", "cfa23", "cfa24", "cfa4x")
    )
  }

  if (is.null(region_titles)) {
    region_titles = switch( region_groups,
      default = c(cfanorth="NENS", cfasouth="SENS", cfa4x="4X"),
      split   = c(cfanorth="NENS", cfa23="CFA23", cfa24="CFA24", cfa4x="4X")
    )
  }

  if (is.null(sizedatadir)) {
    sizedatadir = switch( region_groups,
      default = file.path(p$project.outputdir, "size_structure"),
      split   = file.path(p$project.outputdir, "size_structure_split")
    )
  }
  
  if (is.null(outdir)) {
    outdir = switch( region_groups,
      default = file.path( p$annual.results, "figures", "size.freq", "survey" ),
      split   = file.path( p$annual.results, "figures", "size.freq", "survey_split" )
    )
  }
 
  # note ranges in CW will be log transformed later
  if (is.null(span)) {
    span = function( sexid) {
        switch(sexid,
            male   = c( 5, 155, 50),  
            female = c( 5, 95,  30)   
        )
    } 
  }

  # merge data and add region identifiers
  M = size_distributions(p=p, toget="rawdata", regions=regions, outdir=sizedatadir, redo=TRUE)  # merge det and set with QA/QC
        
  for (yg in 1:length(yr_groups)) {

    years = yr_groups[[yg]]

    # discretize size and compute crude means along factors
    M = size_distributions(p=p, toget="crude", span=span, Y=years, regions=regions, outdir=sizedatadir, redo=TRUE )  

    for (sx in c("female", "male")) {
 
      sxcol = switch( sx,
        female = c("darkorange", "gray95" ) ,
        male =   c("slategray", "gray95" ) 
      )

      for (yvar in c("den", "denl")) {
        # den = arithmetic mean density,  denl = geometric mean density  

        xvals = discretize_data( span=span(sx) )
        xd = 4  # how many ticks to skip
        
        plot_histogram_carapace_width( 
          M=M, 
          years=years, 
          regions=regions, 
          region_titles=region_titles,
          plot_sex=sx,
          Mdelta=xd,  # x-label intervals
          yvar=yvar, 
          cols = sxcol,
          plotoutdir=file.path(outdir, names(yr_groups)[yg]) 
        )

      }
    }
  }

  return("done")

}
