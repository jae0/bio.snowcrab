

fisherydata_summary = function( FD=NULL, toget="data", regions = c("cfanorth", "cfasouth", "cfa4x"), y=NULL ) {

    if(is.null(FD)) {
        require(data.table)
        fishery = setDT( snowcrab_landings_db() )
        if (regions[1]=="cfaall") {
            fishery$region = "cfaall"
        } else {
            fishery$region = NA
            for (rg in regions) {
                if ( rg=="cfanorth" ) cfa = c("cfa20", "cfa21", "cfa22", "cfanorth", "north")
                if ( rg=="cfasouth" ) cfa = c("cfa23", "cfa24", "cfasouth", "cfaslope")
                if ( rg=="cfa4x") cfa = "cfa4x"
                ii = which( fishery$cfa %in% cfa )
                if (length(ii) > 0) {
                    fishery$region[ii] = rg 
                } 
            }
        }

        if (is.null(y)) {
            y = sort(unique(fishery$yr) )
        } else {
            fishery = fishery[ which(fishery$yr %in% y) , ]
        }

        fishery$cpue[ which( fishery$cpue > (650*0.454)) ] = NA  # same rule as in snowcrab_landings_db -- 650lbs/trap is a reasonable upper limit 

        # due to landings occuring in GULF region, some tweaks here: 
        ii = which( FD$region=="cfanorth" & FD$year==2013 ) 
        if (length(ii) ==1) FD$cpue[ii ] = 106 
        ii = which( FD$region=="cfanorth" & FD$year==2014 ) 
        if (length(ii) ==1) FD$cpue[ii ] = 104.5 
        
        FD = fishery[, .(landings=sum(landings, na.rm=TRUE), cpue=mean(cpue, na.rm=TRUE)), by=.(region, year) ]
        FD$effort =  FD$landings / FD$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
        FD = FD[ which(!is.na(FD$region)) , ]
        FD = FD[ order(region, year),]
        setDF(FD)

        FD$landings = FD$landings / 1000  # t
        FD$effort = FD$effort / 1000  # 1000 th
    }

    if (toget=="data") {
        return(FD)
    }
 
    xrange = range(FD$year)
    xrange[1] = xrange[1]
    xrange[2] = xrange[2]
    xlabels = pretty(xrange, n=10)

    pts = c(19, 22, 24)
    lns = c(1, 1, 1)
    cols = c("grey10", "grey10",  "grey20") 
   
    if (toget=="timeseries_landings") {
        plot( landings ~ year, data=FD, type="n", ylab="Landings (t)", xlab="Year", xaxt="n", xlim=xrange, ylim=c(0, max(FD$landings, na.rm=T)) ) 
        for (m in 1:length(regions)) {
            ii = which( FD$region == regions[m] )
            if (length(ii) > 0){
                if (regions[m] %in% c("cfanorth", "cfa4x" )) {
                    points( I(landings*5) ~ year, data=FD[ii,], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
                } else {
                    points( landings ~ year, data=FD[ii,], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
                }
            }
        }
        axis(1, at=xlabels, labels=TRUE)   
        axis( 2 )
        legend("topleft", c("N-ENS X 5", "S-ENS", "4X X 5"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2)
    }

  
    if (toget=="timeseries_cpue") {
        plot( cpue ~ year, data=FD, type="n", ylab="Catch rate (kg/trap)", xlab="Year", xaxt="n", xlim=xrange, ylim=c(0, max(FD$cpue, na.rm=T)) ) 
        for (m in 1:length(regions)) {
            ii = which( FD$region == regions[m] )
            if (length(ii) > 0){
                points( cpue ~ year, data=FD[ii,], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
            }
        }
        axis(1, at=xlabels, labels=TRUE)   
        axis( 2 )
        legend("topleft", c("N-ENS", "S-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2)
    }
  
    if (toget=="timeseries_effort") {
        plot( effort ~ year, data=FD, type="n", ylab="Effort (1000 trap hauls)", xlab="Year", xaxt="n", xlim=xrange, ylim=c(0, max(FD$effort, na.rm=T)) ) 
        for (m in 1:length(regions)) {
            ii = which( FD$region == regions[m] )
            if (length(ii) > 0){
                points( effort ~ year, data=FD[ii,], type="b", col=cols[m], lwd=3, lty=lns[m], pch=pts[m])
            }
        }
        axis(1, at=xlabels, labels=TRUE)   
        axis( 2 )
        legend(x="topleft", c("N-ENS", "S-ENS", "4X"), bty="n", lty=lns, lwd=2, pch=pts, col=cols, cex=1.2)
    }

}