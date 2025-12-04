
create_size_frequencies = function(
    p, 
    mau="region", 
    yr_groups=NULL, 
    span=NULL,
    outdir=NULL  
  ) {
  
  if (0) {
     mau="region"
     span=NULL   
     outdir=NULL  
  }


  if (is.null(yr_groups)) {
    yr_groups = list(
      period1 = as.character( rev(p$yrs)[ 10:1 ] ),
      period2 = as.character( rev(p$yrs)[ 20:11 ] )
    )
  } else if (yr_groups=="all") {
    nyrs = length(p$yrs)
    yr_groups = list()
    for (i in 1:ceiling(nyrs/10)) {
      j = i*10 + c(0:-9)
      j = j[ j %in% 1:nyrs ]
      yr_groups[[ paste0("period", i) ]] = as.character( rev(p$yrs)[ j ] )      
    }
  }
 
 
  maus = management_areal_units( mau=mau )  
  
  # note ranges in CW will be log transformed later
  if (is.null(span)) {
    span = function( sexid) {
        switch(sexid,
            male   = c( 5, 155, 50),  
            female = c( 5, 95,  30)   
        )
    } 
  }


  for (yg in 1:length(yr_groups)) {

    years = yr_groups[[yg]]

    # discretize size and compute crude means along factors
    M = size_distributions(p=p, toget="crude", span=span, Y=years, mau=mau, redo=TRUE )  

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
          mau=mau,
          plot_sex=sx,
          xspan=span(sx),
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
