
create_size_frequencies_single = function(
    p, 
    mau="region", 
    yr=max(p$yrs), 
    span=NULL,
    outdir=NULL  
  ) {
  
  if (0) {
     mau="region"
     span=NULL   
     outdir=NULL  
  }

  plotoutdir = file.path( outdir, mau )
 
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

  Mdelta = 4  # x-label intervals  # how many ticks to skip

  pcex = 1.0

  width = trunc( 2048 / 1.5 )
  height = trunc( 2048 / 2.0 )

  ftype = "png"

  yvar = "denl"  # den = arithmetic mean density,  denl = geometric mean density  

  # discretize size and compute crude means along factors
  M0 = size_distributions(p=p, toget="crude", span=span, Y=yr, mau=mau, redo=FALSE ) 
   
  for (sx in c("female", "male")) {

    sex_code = switch( sx, male = "0",  female = "1" )

    cols = switch( sx,
      female = c("darkorange", "gray95" ) ,
      male =   c("slategray", "gray95" ) 
    )
    
    xvals = discretize_data( span=span(sx) )
    xspan = span(sx)
       
    xlabs = discretize_data(span=xspan) # midpoints

    dir.create( plotoutdir, recursive=TRUE, showWarnings=FALSE )
 
   # pl = layout( matrix( 1, nrow=1, ncol=1, byrow=F ) )
   # par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
   # par(mar=c(0.2, 0.4, 0.4, 1.2))
    M = copy(M0)
    M$Y = M[,..yvar]

    yranges = M[, .(ymax=sum(Y, na.rm=TRUE)), by=.(auid, year, sex, cwd)]
    yran = yranges[, .(yx=max(ymax, na.rm=TRUE)), by=.(auid, sex) ]
    
    # cols = c("blue3", "darkslategray1")
    M$cw = as.numeric( as.character(M$cwd))
    M = M[ cw < xspan[2] & cw > xspan[1], ]

    M$cwd = factor( M$cw )

    # xaxis values and indices:
    rn = as.numeric( as.character(levels(M$cwd) ))

    rl = length(rn)
    rni = 1:rl

    # for printing labels
    rj = 1:length(xlabs)

    rj = rj[seq(Mdelta, length(rj), Mdelta)]
    xlabs = xlabs[seq(Mdelta, length(xlabs), Mdelta)]

    for (a in 1:maus[["n"]]) {

      fn = file.path( plotoutdir, paste( paste(sx, yvar, maus[["internal"]][a], yr, sep="_"), ftype, sep=".") )

      png(filename=fn, width=width, height=height, pointsize=12, res=192, bg='white')  

      ylim=c(0, yran[ auid==maus[["internal"]][a] & sex==sex_code, yx ])

      ms = M[ 
          auid == maus[["internal"]][a] &
          year == yr &
          sex == sex_code & 
          mat %in% c("0", "1"), 
          which = TRUE]  # return only row numbers
      
      skip_plot = FALSE
      if (length(ms)==0) {
          toplot = data.table(cwd=rn, "0"=0, "1"=0)
          skip_plot =TRUE
      } else {
          toplot = dcast( M[ms,],
              formula= cwd ~ mat, 
              value.var="Y",
              fun.aggregate=mean, na.rm=TRUE   # should not be needing this.. in case a larger subset is chosen
          )
          toplot$cwd = as.numeric(as.character( toplot$cwd  ))
          toplot = toplot[ data.table(cwd=rn), on="cwd"]
#            toplot[  is.na(.(0)), "0" ] = 0
      } 

      if ( a==3 & yr %in% as.character(c(1998:2000)) ) {
          skip_plot = TRUE
      } 
      
      axes = TRUE
      if (skip_plot)  axes=FALSE
      
      axisnames = TRUE  # last row

      barplot( t( toplot[, c("1", "0")] ), 
          space=0L, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, 
          xpd=FALSE, lwd=0.001, las=1,  
          xaxt="n", cex.axis=pcex, cex.names=pcex, axis.lty=1
      )

      # last row labels
      axis(1, at =rj, labels=round(xlabs,0) )
      # mtext(maus[["labels"]][a], side=3, line=0, cex=pcex  ) 
      
      # text( length(rn)*0.9, ylim[2]*2/3, yr, cex=pcex )

      if (!skip_plot) {
          if (sx=="female") {
              l = max(which(rn<55)) 
              abline( v=l, lwd=1.5, lty="dotted", col="gray" )
          }  
          if (sx=="male") {
              l = max(which(rn<95)) 
              abline( v=l, lwd=1.5, lty="dotted", col="gray" )
          }  
          axis(2, las=2, cex=pcex*0.75) 
      }  

      mtext("Carapace width (mm) / Largeur de la carapace (mm) ", side=1, line=2.5, cex=pcex)
      mtext(expression(paste("No. / ", km^2)), side=2, line=3, cex=pcex)
      
      # print(fn)
      print(fn)

      dev.off()  

    }  # end region
  
  } # end sex



  return("done")

}
