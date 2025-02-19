figure.effort.timeseries = function(yearmax, outdir = NULL, outfile = NULL, outfile2 =NULL, 
  type = "line", plotmethod="default", plottype="png", 
  regions =  list(region=c("cfanorth", "cfasouth", "cfa4x")), 
  region_label = c("N-ENS", "S-ENS", "4X") ) {
 
  dir.create( outdir, recursive=T, showWarnings=F  )
  fn = file.path( outdir, paste( outfile, plottype, sep="." ) )
  fn2 = file.path(outdir, paste(outfile2, "pdf",sep = "."))


  if (plotmethod=="default") {
    require(ggplot2)
    k = NULL
    k = fishery_data( toget="summary_annual", regions=regions )
    vn = names(regions)
    reg = unlist(regions)
    for (i in 1:length(reg) ) {
      k[[vn]] = gsub( reg[i], region_label[i] , k[[vn]] )
    }

    k[[vn]]= factor(k[[vn]], levels=region_label)
    k$region = k[[vn]]

    k$effort = k$effort / 1000
    
    color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1:length(reg)]
    shapes = c(15, 17, 19, 21, 23)[1:length(reg)]

    out = ggplot(k, aes(x=yr, y=effort, fill=region, colour=region)) +
      geom_line( alpha=0.9, linewidth=1 ) +
      geom_point(aes(shape=region), size=5, alpha=0.7 )+
      labs(x="Year / Année", y="Effort (1000 trap haul) /\n Débarquements (1000 casiers levé)", size = rel(1.5)) +
      scale_colour_manual(values=color_map) +
      scale_fill_manual(values=color_map) +
      scale_shape_manual(values = shapes) +
      theme_light( base_size = 22) + 
      theme( legend.position="inside", legend.position.inside=c(0.1, 0.9), legend.title=element_blank()) 
      # scale_y_continuous( limits=c(0, 300) )  
      ggsave(filename=fn, plot=out, device=plottype, width=12, height = 8)

    return( fn )
  }


  if (plotmethod=="old") {

    e = NULL
    e = fishery_data( toget="summary_annual", regions=regions )
    e = as.data.frame(e)
    colnames(e) = regions
    rownames(e) = res$yr
    e$effort = e$effort / 1000
    
    e = e[which(as.numeric(rownames(e)) <= yearmax),]
    uyrs = as.numeric(rownames(e))
    
    #png( file=fn,units='in', width=7,height=7,pointsize=10, res=350,type='cairo')
    pdf(file = fn, width = 7, height = 7, bg = 'white')
    
    if (type == "bar") {
      e[is.na(e)] = 0
      formed.data = t(as.matrix(e))
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
      e[which(e == 0)] = NA
      yrange = range (e, na.rm = T)
      yrange[1] = 0
      xrange = range(uyrs)
      xrange[1] = xrange[1]
      xrange[2] = xrange[2]
      xlabels = seq(xrange[1] + 1, xrange[2], 2)
      
      
      m = 1; plot(
        uyrs, e[,m],  type = "b", ylab = "Effort (1000 trap hauls/casiers levés)", xlab = "Year", col =
          cols[m], lwd = 3, lty = lns[m], pch = pts[m], xaxt = "n", xlim = xrange, ylim =
          yrange
      )
      m = 2; points(
        uyrs, e[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
      )
      m = 3; points(
        uyrs, e[,m], type = "b", col = cols[m], lwd = 3, lty = lns[m], pch = pts[m]
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
      
      sm = e[, c(1, 3)]
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
    #table.view( e)
    return(fn)
  
  }

}
