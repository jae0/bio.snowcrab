  figure.cpue.timeseries = function( yearmax, outdir=NULL, outfile=NULL, plotmethod="default", 
    regions = c("cfanorth", "cfasouth", "cfa4x"), region_label = c("N-ENS", "S-ENS", "4X") ) {
     
    dir.create( outdir, recursive=T, showWarnings=F  )
    fn = file.path( outdir, paste( outfile, "pdf", sep="." ) )
    fnpng = file.path( outdir, paste( outfile, "png", sep="." ) )
  
    if (plotmethod=="default") {
      require(ggplot2)
      k = NULL

      for (i in 1:length(regions)) {
        res = get.fishery.stats.by.region(Reg=regions[i])
        res$region = region_label[i]
        k = rbind( k,  res[, c("yr", "cpue", "region" )] )
      }

      k$region = factor(k$region, levels=region_label)
      color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1:length(regions)]

      out = ggplot(k, aes(x=yr, y=cpue, fill=region, colour=region)) +
        geom_line( alpha=0.9, linewidth=1 ) +
        geom_point(aes(shape=region), size=5, alpha=0.7 )+
        labs(x="Year / Année", y="Catch rate (kg/trap) /\n Taux de prise (kg/casier levé)", size = rel(1.5)) +
        # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
        scale_colour_manual(values=color_map) +
        scale_fill_manual(values=color_map) +
        scale_shape_manual(values = c(15, 17, 19)) +
        theme_light( base_size = 22) + 
        theme( legend.position="inside", legend.position.inside=c(0.1, 0.9), legend.title=element_blank()) 
        # scale_y_continuous( limits=c(0, 300) )  
        ggsave(filename=fn, plot=out, device="pdf", width=12, height = 8)
        ggsave(filename=fnpng, plot=out, device="png", width=12, height = 8)
      return( fn )
    }


    if (plotmethod=="old") {
      k = NULL
      for (r in regions) {
        res = get.fishery.stats.by.region(Reg=r)
        k = cbind( k, res$cpue )
      }
      
      k = as.data.frame( k )
      colnames(k) = regions
      rownames(k) = res$yr
      
      k = k[ which( as.numeric(rownames(k)) <= yearmax ), ] 
      uyrs = as.numeric(rownames(k) ) 

      pts = c(19, 22, 24)
      lns = c(1, 1, 1)
      cols = c("grey10", "grey10",  "grey20")

      yrange = range (k, na.rm=T)
      yrange[1] = 0
      xrange = range(uyrs)
      xrange[1] = xrange[1]
      xrange[2] = xrange[2]
  #    xlabels = c(xrange[1], xrange[1]+8, xrange[1]+18, xrange[1]+28, xrange[2])
      xlabels = seq(xrange[1]+1, xrange[2], 2)

      #Cairo( file=fn, type="png", bg="white", , pointsize=30, units="in", width=6, height=4, dpi=300 )
      pdf(file=fn, width=7, height=7, bg='white')
      #png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
        m=1; plot( uyrs, k[,m],  type="b", ylab="Catch rate (kg/trap)", xlab="Year", col=cols[m], lwd=3, lty=lns[m], pch=pts[m], xaxt="n", xlim=xrange, ylim=yrange)
        m=2; points(uyrs, k[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
        m=3; points(uyrs, k[,m], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
        axis( 1, at=xlabels, labels=FALSE )
        text(x=xlabels+1, y=par('usr')[3], labels=xlabels, srt=45, adj=c(1.5,1), xpd=TRUE)

        axis( 2 )
        legend(x=1980, y=100, c("N-ENS", "S-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2 )
      
      dev.off()

      # cmd( "convert -trim -frame 10x10 -mattecolor white ", fn, fn )
      #table.view(k)

      return( fn )
    }
 
  }


