
plot_histogram_carapace_width = function( M, 
    years=as.character(2012:2022), 
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    region_titles=c(cfanorth="NENS", cfasouth="SENS", cfa4x="4X"),
    plot_sex="male",
    yvar="den",
    xlim=NULL,
    pcex = 1.0,
    cols = c("gray40", "gray90" ),
    width=1792, 
    height=2048,
    ftype = "png",
    outdir = NULL ) {
    

    if (plot_sex=="male") {
        if (is.null(xlim)) xlim=c(10, 150)
        sex_code = "0"
    } else if (plot_sex=="female") {
        if (is.null(xlim)) xlim=c(10,  90)
        sex_code = "1"
    }

    if (!is.null(outdir)) {
        dir.create( outdir, recursive=TRUE, showWarnings=FALSE )

        fn = file.path( outdir, paste( plot_sex, yvar, ftype, sep=".") )

        if (grepl("\\.pdf$", fn)) {
            message("warning: PDF options needs to be tweaked ")
            pdf(file=fn, width=6, height=8, bg='white')
        }  else if (grepl("\\.png$", fn)) {
            png(filename=fn, width=width, height=height, pointsize=12, res=192, bg='white')  
        } else{
            message("Format not implemented... ", ftype)
        }

    } else {
        dev.new()
    }    

    nrows = length(years)
    ncols = length(regions)

    if (!exists("region", M)) {
        if (length(regions)==1) {
            M$region = regions
        } else {
            message("regions not found in data")
        }
    }

    pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
    # par(oma=c(5, 4, 4, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
    par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
    par(mar=c(0.2, 0.4, 0.4, 1.2))
 
    M$Y = M[,..yvar]
  
    yranges = M[year %in% years, .(ymax=sum(Y, na.rm=TRUE)), by=.(region, year, sex, cwd)]
    yran = yranges[, .(yx=max(ymax, na.rm=TRUE)), by=.(region,sex) ]
    
    # cols = c("blue3", "darkslategray1")
    
    M$cw = as.numeric( as.character(M$cwd))
    M = M[ cw < xlim[2] & cw > xlim[1], ]

    M$cwd = factor( M$cw )
    rn = as.numeric( as.character(levels(M$cwd) ))

    if (length(region_titles) !=length(regions) ){
        region_titles = region_titles[regions] 
        if (length(region_titles) !=length(regions) ){
            region_titles = toupper(regions) 
        }
    }

    for (a in 1:length(regions)) {
    
    ylim=c(0, yran[ region==regions[a] & sex==sex_code, yx ])
 
    for (b in 1:length(years)) {
        ms = M[ 
            region == regions[a] &
            year == years[b] &
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

        if ( a==3 & years[b] %in% as.character(c(1998:2000)) ) {
            skip_plot = TRUE
        } 
        
        axes = TRUE
        if (skip_plot)  axes=FALSE
        
        axisnames = FALSE
        if ( b == nrows ) axisnames = TRUE  # last row
  
         barplot( t( toplot[, c("1", "0")] ), 
            space=0, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, 
            xpd=FALSE, lwd=0.001, las=1,  
            names.arg=rn, cex.axis=pcex, cex.names=pcex, axis.lty=1
        )
 
        if ( b == 1 ) {
            mtext(region_titles[a], side=3, line=0, cex=pcex  ) 
        }

        if ( a == ncols) {
            text( length(rn)*0.9, ylim[2]*2/3, years[b], cex=pcex )
        }
  
        if (!skip_plot) {
            if (plot_sex=="female") {
                l = max(which(rn<55)) 
                abline( v=l, lwd=1.5, lty="dotted", col="gray" )
            }  
            if (plot_sex=="male") {
                l = max(which(rn<95)) 
                abline( v=l, lwd=1.5, lty="dotted", col="gray" )
            }  
            axis(2, las=2, cex=pcex*0.75) 
        }  

    }}
 
    mtext("Carapace width (mm)", side=1, outer=TRUE, line=2.5, cex=pcex)
    mtext(expression(paste("No. / ", km^2)), side=2, outer=TRUE, line=3, cex=pcex)
    
    # print(fn)

    if (!is.null(outdir)) {
        print(fn)
        dev.off()  
    }

}


