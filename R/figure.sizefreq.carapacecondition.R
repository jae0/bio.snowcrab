 
figure.sizefreq.carapacecondition = function( X, cwbr=4, xrange=c(44, 184),
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    outdir=file.path( p$annual.results, "figures", "size.freq", "carapacecondition" )  ) {

    dir.create(outdir, recursive=TRUE)

    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    breaks = seq(xrange[1], xrange[2], by=cwbr)
    mids = breaks[-length(breaks)] + cwbr/2

    # color_map <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # color_map = c("#E69F00", "#56B4E9",  "#CC79A7" )
    # color_map <- c("1"="#999999", "2"="#E69F00", "3"="#56B4E9", "4"="#009E73", "5"="#F0E442" )  #, "#0072B2", "#D55E00", "#CC79A7")[1:5]
    color_map <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442" )  #, "#0072B2", "#D55E00", "#CC79A7")[1:5]

    X$fishyr = X$yr 

    ii = which( X$cw > 50 & X$cw < 170 & X$mat==1 & X$sex==male) 
    X = X[ii,]

    X$region = NA
    for ( region in regions ) {
        r = polygon_inside(x=X, region=aegis.polygons::polygon_internal_code(region), planar=F)
        if (length(r) > 0) X[r, "region"] = region
    }
    X= X[!is.na(X$region), ]

    X$cwd = as.numeric( as.character( cut( X$cw, breaks=breaks, labels =mids ) ) )
    X = X[is.finite(X$shell),]
    X$shell = factor( X$shell )
    shells = levels(X$shell)
    

    yrs = sort( unique( X$fishyr ) )

    for (y in yrs) {
    for (r in regions) {
        Y = X[X$region==r & X$fishyr == y, ]
        setDT(Y)

        N = Y[ CJ( shell, cwd, unique=TRUE ), on=.( shell, cwd ), .(N=.N), by=.EACHI ]
        N$shell = factor(N$shell, levels=shells)

        Ntot = sum(N$N)
        NS = N[CJ(shell, unique=TRUE), on=.(shell), .(Stot=sum(N)), by=.EACHI ]
        NS$p = round( NS$Stot / Ntot*100, 1)
        
        all = data.table(shell=shells)
        NS = NS[all, on="shell"]
        oo = which(is.na(NS$Stot))
        if (length(oo)>0) {
            NS$Stot[oo] = 0
            NS$p[oo] = 0
        }

        llabs = paste( c("1", "2", "3", "4", "5"), ": ", NS$p, "%", sep="" )
        tit = paste(y, ": n=", Ntot, sep="")

        out = ggplot(N, aes(fill=shell, y=N, x=cwd)) + 
            geom_bar(position="stack", stat="identity") +
            scale_fill_manual( name=tit, values=color_map, labels=llabs, drop = FALSE) +
            xlim( xrange[1]*1.1, xrange[2]*0.75 ) +
            theme_light( base_size = 20) + 
            theme( legend.position.inside=c(0.15, 0.8),   axis.title.x=element_blank(), axis.title.y=element_blank()) 
         
        out
        fn = file.path( outdir, paste("sizefreq", r, y, "pdf", sep=".") )
        ggsave(filename=fn, plot=out, device="pdf", width=8, height = 6)

    }}

}

 
