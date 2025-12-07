figure.effort.timeseries = function(yearmax, outdir = NULL, outfile = NULL, outfile2 =NULL, 
  type = "line", plotmethod="default", plottype="png", 
  mau="region" ) {
 
  dir.create( outdir, recursive=T, showWarnings=F  )
  fn = file.path( outdir, paste( outfile, plottype, sep="." ) )
  fn2 = file.path(outdir, paste(outfile2, "pdf",sep = "."))

  maus = management_areal_units( mau=mau )  
  
  color_map = maus[["color_map"]]
  shapes = maus[["shapes"]]

  FD = NULL
  FD = fishery_data( mau=mau )
  AN = FD[["summary_annual"]]


  if (plotmethod=="default") {
    require(ggplot2)
   
    for (i in 1:maus[["n"]] ) {
      AN[[mau]] = gsub( maus[["internal"]][i], maus[["labels"]][i] , AN[[mau]] )
    }

    AN[[mau]]= factor(AN[[mau]], levels=maus[["labels"]])
    AN$area = AN[[mau]]
    AN$effort = AN$effort / 1000
     
    out = ggplot(AN, aes(x=yr, y=effort, fill=area, colour=area)) +
      geom_line( alpha=0.9, linewidth=1 ) +
      geom_point(aes(shape=area), size=5, alpha=0.7 )+
      labs(x="Year / Année", y="Effort (1000 trap haul) /\n Débarquements (1000 casiers levé)" ) +
      scale_colour_manual(values=maus[["color_map"]]) +
      scale_fill_manual(values=maus[["color_map"]]) +
      scale_shape_manual(values = maus[["shapes"]]) +
      theme_light( base_size = 22) + 
      theme( legend.position="inside", legend.position.inside=c(0.1, 0.9), legend.title=element_blank()) 
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=out, device=plottype, width=12, height = 8)

    return( fn )
  }


  if (plotmethod=="old") {

    AN = as.data.frame(AN)
    colnames(AN) = regions
    rownames(AN) = res$yr
    AN$effort = AN$effort / 1000
    
    AN = AN[which(as.numeric(rownames(AN)) <= yearmax),]
    uyrs = as.numeric(rownames(AN))
    
    #png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
    pdf(file = fn, width = 7, height = 7, bg = 'white')
    
    if (type == "bar") {
      AN[is.na(AN)] = 0
      formed.data = t(as.matrix(AN))
      barplot(
        formed.data, space = 0, xlab = "Year", ylab = "Effort (1000 trap hauls)", col =
          cols
      )
      legend(
        x = 1, y = 130, c("N-ENS", "S-ENS", "4X"), fill = cols[reverse], bty = "n"
      )
    }
    if (type == "line") {
      pts = c(19, 22, 24)
      lns = c(1, 1, 1)
      cols = c("grey10", "grey10",  "grey20")
      AN[which(AN == 0)] = NA
      yrange = range (AN, na.rm = T)
      yrange[1] = 0
      xrange = range(uyrs)
      xrange[1] = xrange[1]
      xrange[2] = xrange[2]
      xlabels = seq(xrange[1] + 1, xrange[2], 2)
      
      
      m = 1; plot(
        uyrs, AN[,m],  type = "b", ylab = "Effort (1000 trap hauls/casiers levés)", xlab = "Year", col =
          cols[m], lwd = 3, lty = lns[m], pch = pts[m], xaxt = "n", xlim = xrange, ylim =
          yrange
      )
      m = 2; points(
        uyrs, AN[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
      )
      m = 3; points(
        uyrs, AN[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
      )
      axis(1, at = xlabels, labels = FALSE)
      text(
        x = xlabels + 1, y = par('usr')[3], labels = xlabels, srt = 45, adj = c(1.5,1), xpd =
          TRUE
      )
      axis(2)
      legend(
        x = 1980, y = 100, c("N-ENS", "S-ENS", "4X"), bty = "n", lty = lns, lwd =
          2, pch = pts, col = cols, cex = 1.2
      )
      dev.off()
      
    # png(file = fn2 ,units = 'in', width = 7, height = 7,pointsize = 10, res = 350,type ='cairo')
    pdf(file=fn2, width=7, height=7, bg='white')
      
      sm = AN[, c(1, 3)]
      pts = c(19, 24)
      lns = c(1, 1)
      cols = c("grey10", "grey20")
      yrange = range (sm, na.rm = T)
      yrange[1] = 0
      xrange = range(uyrs)
      xrange[1] = xrange[1]
      xrange[2] = xrange[2]
      xlabels = seq(xrange[1] + 1, xrange[2], 2)
      
      
      m = 1; plot(
        uyrs, sm[,m],  type = "b", ylab = "Effort (1000 trap hauls/casiers levés)", xlab = "Year", col =
          cols[m], lwd = 4, lty = lns[m], pch = pts[m], xaxt = "n", xlim = xrange, ylim =
          yrange
      )
      m = 2; points(
        uyrs, sm[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
      )
      axis(1, at = xlabels, labels = FALSE)
      text(
        x = xlabels + 1, y = par('usr')[3], labels = xlabels, srt = 45, adj = c(1.5,1), xpd =
          TRUE
      )
      axis(2)
      legend(
        x = 1985, y = 45, c("N-ENS", "4X"), bty = "n", lty = lns, lwd = 2, pch =
          pts, col = cols, cex = 1.2
      )
    }
    dev.off()
  # cmd("convert -trim -frame 10x10 -mattecolor white ", fn, fn)
    #table.view( AN)
    return(fn)
  
  }

}
