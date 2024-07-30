
  figure.timeseries.survey.direct = function( p, outdir, variables, plotyears, type="biologicals", all.areas=T, minN=10, u=NULL, graphic='pdf', bg="white", plotmethod="default",
  regions = c("cfanorth", "cfasouth", "cfa4x"),  region_label = c("N-ENS", "S-ENS", "4X"), backtransform=FALSE  ) {

    if (!all.areas) {
      # not for default method .. only for lattice-based stuff
      regions = c("cfasouth", "cfanorth" )
      region_label = c("S-ENS", "N-ENS")
    }

    n.region_label = length(region_label)
    n.regions = length(regions)

    # base data
    tdb = snowcrab.db( DS="det.georeferenced", p=p)
    tdb = tdb[ which(tdb$sex==0 & tdb$mat==1), ]
  
    cfa4x = polygon_inside(tdb, aegis.polygons::polygon_internal_code("cfa4x"))
    cfanorth = polygon_inside(tdb, aegis.polygons::polygon_internal_code("cfanorth"))
    cfasouth = polygon_inside(tdb, aegis.polygons::polygon_internal_code("cfasouth"))

    tdb$region = NA
    tdb$region[cfa4x] = "cfa4x"
    tdb$region[cfanorth] = "cfanorth"
    tdb$region[cfasouth] = "cfasouth"

    tdb$region = factor(tdb$region, levels=regions, labels =region_label)
    tdb = tdb[(which(!is.na(tdb$region))), ]

    #  load transformation tables associated with a given variable
    tdb = tdb[ which(is.finite(tdb$cw) ), ]
    setDT(tdb)

    td = tdb[, .(mean=mean(cw), lb=quantile(cw, 0.025), ub=quantile(cw, 0.975)), by=.(region, yr) ]
    setDF(td)

    outdir=file.path(p$annual.results, "timeseries", "survey")

      dir.create( outdir, recursive=T, showWarnings=F )

      fn = file.path( outdir, paste( "cw.mat", "pdf",  sep="." ) )
     
      require(ggplot2)
 
    color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
  
    out = ggplot( td, aes(x=yr, y=mean, fill=region, colour=region)) +
        geom_line( alpha=0.9, linewidth=1.2 ) +
        geom_point(aes(shape=region), size=3, alpha=0.7 ) +
        geom_errorbar(aes(ymin=lb,ymax=ub), linewidth=0.8, alpha=0.8, width=0.3)  +
        labs(x=NULL, y=NULL) +
        # labs(x="Year", y="", size = rel(1.5)) +
        scale_colour_manual(values=color_map) +
        scale_fill_manual(values=color_map) +
        scale_shape_manual(values = c(15, 17, 19)) +
        theme_light( base_size = 22) + 
        theme( legend.position="inside", legend.position.inside=c(0.1, 0.9), legend.title=element_blank()) 

        # scale_y_continuous( limits=c(0, 300) )  
        ggsave(filename=fn, plot=out, device="pdf", width=12, height = 8)

    print( fn )

     return("Done")
  }
