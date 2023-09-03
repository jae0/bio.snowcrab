
plot_histogram_carapace_width = function( M, 
    years=as.character(2012:2022), 
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    plot_sex="male",
    yvar="den",
    xlim=NULL,
    pcex = 1.0,
    width=1792, 
    height=2048,
    fn = NULL ) {
    
    if (plot_sex=="male") {
        if (is.null(xlim)) xlim=c(0, 140)
        sex_code = "0"
    } else if (plot_sex=="female") {
        if (is.null(xlim)) xlim=c(0, 100)
        sex_code = "1"
    }

    if (!is.null(fn)) {

        if (grepl("\\.pdf$", fn)) {
            message("warning: PDF options needs to be tweaked ")
            pdf(file=fn, width=6, height=8, bg='white')
        }  
        if (grepl("\\.png$", fn)) {
            png(filename=fn, width=width, height=height, pointsize=12, res=192, bg='white')  
        }

    } else {
        dev.new()
    }    

    nrows = length(years)
    ncols = length(regions)

    pl = layout( matrix( c(1:(ncols*nrows)), nrow=nrows, ncol=ncols, byrow=F ) )
    # par(oma=c(5, 4, 4, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
    par(oma=c(6, 6, 6, 1)) # outer margins default:  c(0, 1, 0, 1)'c(bottom, left, top, right)'
    par(mar=c(0, 0.4, 0.4, 1.2))
 
    M$Y = M[,..yvar]
    
    yranges = M[, .(ymax=sum(Y, na.rm=TRUE)), by=.(region, year, sex, cwd)]
    yran = yranges[, .(yx=max(ymax)), by=.(region,sex) ]
    
    cols = c("gray10", "gray100" )
    # cols = c("blue3", "darkslategray1")

    
    M$cw = as.numeric( as.character(M$cwd))
    M = M[ cw < xlim[2] & cw > xlim[1], ]

    M$cwd = factor( M$cw )
    rn = as.numeric( as.character(levels(M$cwd) ))
  
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
        }
        

        if (regions[a]==regions[3] & years[b] %in% c(1998:2000) ) {
            skip_plot = TRUE
        } 
        
        axes = TRUE
        #if (regions[a]==regions[1] ) axes=T  # first col

        axisnames = FALSE
        if ( b == nrows ) axisnames = TRUE  # last row

        barplot( t( toplot[, c("1", "0")] ), 
            space=0, axisnames=axisnames, ylim=ylim, axes=axes, col=cols, 
            xpd=FALSE, lwd=0.001, las=1,  
            names.arg=rn, cex.axis=pcex, cex.names=pcex, axis.lty=1
        )

        if ( a == ncols) {
            text( length(rn)*0.9, ylim[2]*2/3, years[b], cex=pcex )
        }
  
        if (!skip_plot) {
            l = max(which(rn<95)) 

            if (plot_sex=="male") {
                abline( v=l, lwd=1.5, lty="dotted", col="gray" )
            }  

            axis(2, las=2, cex=pcex*0.75) 
        }

    }}
     

    mtext("Carapace width (mm)", side=1, outer=TRUE, line=2.5, cex=pcex)
    mtext(expression(paste("No. / ", km^2)), side=2, outer=TRUE, line=3, cex=pcex)
    mtext("NENS", side=3, outer=T, line=1, at=0.15, cex=pcex)
    mtext("SENS", side=3, outer=T, line=1, at=0.5, cex=pcex)
    mtext("4X", side=3, outer=T, line=1, at=0.85, cex=pcex)
    
    print(fn)

    if (!is.null(fn)) dev.off()  
    
}


