
  figure.timeseries.survey = function( p, outdir, variables, plotyears, type="biologicals", all.areas=T, minN=10, u=NULL, 
    graphic='png', bg="white", plotmethod="default",
    regions = c("cfanorth", "cfasouth", "cfa4x"),  region_label = c("N-ENS", "S-ENS", "4X"), backtransform=FALSE  ) {

    
    if (!all.areas) {
      # not for default method .. only for lattice-based stuff
      regions = c("cfasouth", "cfanorth" )
      region_label = c("S-ENS", "N-ENS")
    }

    n.region_label = length(region_label)
    n.regions = length(regions)

    # base data
    tdb = snowcrab.timeseries.db( DS=type, p=p )

    if(missing(variables)){
      variables =  c(
         bio.snowcrab::snowcrab.variablelist("all.to.model"),
         bio.snowcrab::snowcrab.variablelist("snowcrab.cw"),
         bio.snowcrab::snowcrab.variablelist("physical"),
         'sexratio.mat','sexratio.imm','sexratio.all'
      )
      variables = intersect( variables, unique(tdb$variable))
    }
 
    if(missing(plotyears)) plotyears = unique(tdb$year)

    tdb = subset(tdb, variable %in% variables & year %in% plotyears)
    tdb = tdb[ which(tdb$region %in% regions), ]
    tdb$region = factor(tdb$region, levels=regions, labels =region_label)
    tdb = tdb[(which(!is.na(tdb$region))), ]

    #  load transformation tables associated with a given variable

    REPOS = snowcrab.db( DS="data.transforms", p=p)

    tvars = REPOS$varname[which(REPOS$transform=='log10')]

    for ( v in variables ) {
      td = tdb[ which( tdb$variable == v) ,]

      oo = which( td$mean > 0 | is.finite( td$mean ) )
      if (length(oo) < minN ) next()

      if(is.null(u)) u=NULL # for variable specific units, needs a lookup table

      ylim=c(0, max(c(td$ub, td$mean),na.rm=T))
      #browser()
      if(length(grep('ratio',v))==1)ylim=c(0,1)

      xlim=range(td$year)
      if(v %in% tvars){
        ylab = list( v , cex=1)
        ylim[1] = min(c(td$lb, td$mean), na.rm=T)
      } else if (grepl("mass", v)) {
        ylab = list( "Mass per km^2" , cex=1)
      } else if (grepl("no", v)) {
        ylab = list( "Number per km^2" , cex=1)
      } else { 
        ylab = list( v , cex=1)
      }
      xlab = list("Year", cex=1)
      ylim[1] = ylim[1]-diff(ylim)*0.04
      ylim[2] = ylim[2]+diff(ylim)*0.04

      main = capwords(gsub("."," ",v,fixed=T),F,F)

      xlabels = seq(xlim[1], xlim[2], 1)
      ylabels = pretty(ylim,7)

      dir.create( outdir, recursive=T, showWarnings=F )

      fn = file.path( outdir, paste( v, graphic,  sep="." ) )
      if (type=='groundfish.t'){
        fn = file.path( outdir, paste( type, graphic,  sep="." ) )
        main = "Groundfish Survey Temperature"
        xlabels = seq(xlim[1], xlim[2], 2)
      }

      if (v=="R0.mass") {
        main = ""
        ylab = list( "Biomass density log10(kg/km^2) /\n Densité de la biomasse log10(kg/km^2)", cex=1)
      }
    
      dline = ifelse(length(grep('ratio',v))==1,0.5,NA)

      if (plotmethod=="default") {
        require(ggplot2)

        color_map = c("#E69F00", "#56B4E9",  "#CC79A7" , "#D55E00", "#F0E442")[1:length(regions)]
        shapes = c(15, 17, 19, 21, 23)[1:length(regions)]

        if (backtransform) {
          td$mean = 10^(td$mean)
          td$lb = 10^(td$lb)
          td$ub = 10^(td$ub)
        }
        
        out = ggplot(td, aes(x=year, y=mean, fill=region, colour=region)) +
          geom_line( alpha=0.9, linewidth=1.2 ) +
          geom_point(aes(shape=region), size=3, alpha=0.7 ) +
          geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
          labs(x="Year / Année", y=ylab ) +
          scale_colour_manual(values=color_map) +
          scale_fill_manual(values=color_map) +
          scale_shape_manual(values = shapes) +
          theme_light( base_size = 22) + 
          theme( legend.position="inside", legend.position.inside=c(0.1, 0.9), legend.title=element_blank()) 

          # scale_y_continuous( limits=c(0, 300) )  
          ggsave(filename=fn, plot=out, device=graphic, width=12, height = 8)
  
        print( fn )
      }

      if (plotmethod=="lattice") {

        if (graphic=='png') Cairo::Cairo( file=fn, type="png", bg=bg, units="in", dpi=350 )
        if (graphic=='pdf') pdf(file=fn, bg=bg, width=6, height=8 )
        if (graphic=='R') plot.new()
        
        setup.lattice.options()
        pl = xyplot( mean~year|region, data=td, ub=td$ub, lb=td$lb, dline=dline,
              layout=c(1,n.region_label),
              par.strip.text=list(
                plot.symbol=list(col='black', fill='darkgrey', cex=0.75, pch=21),
                axis.text=list(cex=0.7),
                par.main.text=list(cex=1),
                layout.heights=list(strip=1, panel=1, main=0.5),
                strip.background=list(col='lightgrey')),
                #xlim=xlim,
                ylim=ylim,
                scales=list(y=list(at=ylabels, labels=ylabels, cex=0.65, alternating=FALSE), x=list(at=xlabels, labels=xlabels, rot=50, cex=0.65), tck=c(1,0)),
                  main=main, xlab=xlab, ylab=ylab,
                  cex.axis=0.2,
                  cex.main = 1.4,
                  panel = function(x, y, subscripts, ub, lb, ...) {
                larrows(x, lb[subscripts], x, ub[subscripts], angle = 90, code = 3, length=0.05)
                panel.abline(h=median(y,na.rm=T), col="gray", ...)
                panel.abline(h=dline, col="gray", lty=2,...)
                panel.xyplot(x, y, type="b", lty=1, lwd=1.5, pch=21, fill='darkgrey', col="black", ...)
            }
        )

        print(pl)
        print(fn)
        if(graphic!='R') dev.off()
      }

    }  # end for each variable
    # print(fn)
    return("Done")
  }
