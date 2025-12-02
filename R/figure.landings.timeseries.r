
figure.landings.timeseries = function( yearmax, outdir=NULL, outfile=NULL, outfile2=NULL, 
    type="line", plotmethod="default",  plottype="png", 
    mau="region"
  ) {

  dir.create( outdir, recursive=T, showWarnings=F  )
  fn = file.path( outdir, paste( outfile, plottype, sep="." ) ) 


  maus = management_areal_units( mau=mau )  
  
  color_map = maus[["color_map"]]
  shapes = maus[["shapes"]]

  FD = NULL
  FD = fishery_data( mau=mau )
  AN = FD[["summary_annual"]]


  if (plotmethod=="withinset") {
    require(ggplot2)
    
    for (i in 1:maus["n"] ) {
      AN[[mau]] = gsub( maus[["internal"]][i], maus[["label"]][i] , AN[[mau]] )
    }

    AN[[mau]]= factor(AN[[mau]], levels=maus[["label"]])
    AN$region = AN[[mau]]
    AN$landings = AN$landings / 1000

    out = ggplot(AN, aes(x=yr, y=landings, fill=region, colour=region)) +
      geom_line( alpha=0.9, linewidth=1.2 ) +
      geom_point(aes(shape=region), size=5, alpha=0.7 )+
      labs(x="Year / Année", y="Landings (t) / Débarquements (t)") +
      theme_light( base_size = 22) + 
      scale_colour_manual(values=color_map) +
      scale_fill_manual(values=color_map) +
      scale_shape_manual(values = c(15, 17, 19)) +
      theme( legend.position="inside", legend.position.inside=c(0.9, 0.9), legend.title=element_blank()) 
 
    out2=out %+% dplyr::filter(out$data, region %in% c("N-ENS", "4X") ) + 
      scale_colour_manual(values=color_map[c(1,3)]) +
      scale_fill_manual(values=color_map[c(1,3)]) +
      scale_shape_manual(values = c(15, 19)) +
      labs(x=NULL, y=NULL) +
      theme_light( base_size = 16) + 
      theme( legend.position="none") 
 
    require(cowplot)
    o = ggdraw( out ) + draw_plot( out2, x=0.116, y=0.58, width=0.4, height=0.37 )

      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=o, device=plottype, width=12, height = 8)
 
    return( fn )
  }


  if (plotmethod=="default") {
    require(ggplot2)

    for (i in 1:maus["n"] ) {
      AN[[mau]] = gsub( maus[["internal"]][i], maus[["label"]][i] , AN[[mau]] )
    }

    AN[[mau]]= factor(AN[[mau]], levels=maus[["label"]])
    AN$region = AN[[mau]]

    AN$landings = AN$landings / 1000

    o = ggplot(AN, aes(x=yr, y=landings, fill=region, colour=region)) +
      geom_line( alpha=0.9, linewidth=1.2 ) +
      geom_point(aes(shape=region), size=5, alpha=0.7 )+
      labs(x="Year / Année", y="Landings (t) / Débarquements (t)") +
      theme_light( base_size = 22) + 
      scale_colour_manual(values=color_map) +
      scale_fill_manual(values=color_map) +
      scale_shape_manual(values = shapes) +
      theme( legend.position="inside", legend.position.inside=c(0.9, 0.9), legend.title=element_blank()) 
   
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=o, device=plottype, width=12, height = 8)
 
    return( fn )
  }



  if (plotmethod=="old") {
        
      #variables = c("effort", "landings", "cpue")
      #for (v in variables) {
        #Extract data for the raster creation
        #M = K[, c("yr", "lon", "lat", v)]
        #M = M[is.finite(M[,v] *M[,"lon"] *M[,"lat"]),]
 
      AN = as.data.frame( AN )
 
      colnames(AN) = maus[["internal"]]
      rownames(AN) = res$yr
      AN$landings = AN$landings / 1000
    
      AN = AN[ which( as.numeric(rownames(AN)) <= yearmax ), ] 
      uyrs = as.numeric(rownames(AN) ) 
  
      fn2 = file.path( outdir, paste(outfile2,"pdf",sep="." ) )

      pdf(file=fn, width=7, height=7, bg='white')
    # png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')

      if (type=="bar") {
        cols = c("grey10", "grey40",  "grey80")
        reverse = c(3,2,1)
        formed.data = t( AN[,c(3,2,1)] ) # re-order for plot
        formed.data[ is.na(formed.data) ] = 0
        barplot( formed.data, space=0, xlab="Year", ylab="Landings (t)", col=cols)
        legend(x=1, y=10000, c("N-ENS", "S-ENS", "4X"), fill=cols[reverse], bty="n")
      }
      if (type=="line") {
        pts = c(19, 22, 24)
        lns = c(1, 1, 1)
        cols = c("grey10", "grey10",  "grey20") 
        yrange = range (AN, na.rm=T)
        yrange[1] = 0
        xrange = range(uyrs)
        xrange[1] = xrange[1]
        xrange[2] = xrange[2]
        xlabels = seq(xrange[1] +1, xrange[2], 2)

        m=1; plot( uyrs, AN[,m],  type="b", ylab="Landings (t)", xlab="Year", col=cols[m], lwd=4, lty=lns[m], pch=pts[m], xaxt="n", xlim=xrange, ylim=yrange)
        m=2; points(uyrs, AN[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
        m=3; points(uyrs, AN[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
        axis(1, at=xlabels, labels=FALSE)   
        text(x=xlabels+1, y=par('usr')[3], labels=xlabels, srt=45, adj=c(1.5,1), xpd=TRUE)
        axis( 2 )
        legend(x=1980, y=8500, c("N-ENS", "S-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2)

        dev.off()
        pdf(file=fn2, width=7, height=7, bg='white')
        
        #png(file=fn2 ,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
        sm = AN[, c(1, 3)]
        pts = c(19, 24)
        lns = c(1, 1)
        cols = c("grey10", "grey20")
        yrange = range (sm, na.rm=T)
        yrange[1] = 0
        xrange = range(uyrs)
        xrange[1] = xrange[1]
        xrange[2] = xrange[2]
        xlabels = seq(xrange[1]+1, xrange[2], 2)

        m=1; plot( uyrs, sm[,m],  type="b", ylab="Landings (t)", xlab="Year", col=cols[m], lwd=4, lty=lns[m], pch=pts[m], xaxt="n", xlim=xrange, ylim=yrange)
        m=2; points(uyrs, sm[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
        axis( 1, at=xlabels, labels=FALSE )
        text(x=xlabels+1, y=par('usr')[3], labels=xlabels, srt=45, adj=c(1.5,1), xpd=TRUE)
        axis( 2 )
        legend(x=1985, y=1200, c("N-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2 )
      }
      
    dev.off()
      #table.view( AN )
    return( fn )
  }

}

