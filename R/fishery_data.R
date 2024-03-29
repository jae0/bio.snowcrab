

fishery_data = function( 
    toget = c("summary_annual", "summary_monthly", "summary_weekly", "shell_condition", "fraction_observed"), 
    regions = c("cfanorth", "cfasouth", "cfa4x"), y=NULL ) {
 
    require(data.table)
    
    out = list()

    fstat = setDT( snowcrab_landings_db() )

    if (regions[1]=="cfaall") {
        fstat$region = "cfaall"
    } else {
        fstat$region = NA
        for (rg in regions) {
            if ( rg=="cfanorth" ) cfa = c("cfa20", "cfa21", "cfa22", "cfanorth", "north")
            if ( rg=="cfasouth" ) cfa = c("cfa23", "cfa24", "cfasouth", "cfaslope")
            if ( rg=="cfa4x") {
                cfa = "cfa4x"
                # closed in 2018/2019
                fstat0 = fstat[1,] 
                fstat0$yr = fstat0$year= 2018
                fstat0$landings = 1e-9
                fstat0$cpue = 1e-9
                fstat0$effort = 1e-12
                fstat0$region= rg
                fstat = rbind(fstat, fstat0)
            }
            ii = which( fstat$cfa %in% cfa )
            if (length(ii) > 0) fstat$region[ii] = rg 
        }
    }

    if (is.null(y)) y = sort(unique(fstat$yr) )


    # note: "yr" is fishing year, in 4x: 1999-2000 is yr=1999
    fstat = fstat[ which(fstat$yr %in% y) , ]
    
    
    # same rule as in snowcrab_landings_db -- 650lbs/trap is a reasonable upper limit 
    ii = which( fstat$cpue > (650*0.454)) 
    if (length(ii) > 0 ) fstat$cpue[ ii] = NA  
 
    if (any( grepl("summary_annual", toget)) |  any( grepl("fraction_observed", toget)) ) {
        # default is annual summary

        summary_annual = fstat[, .(landings=sum(landings, na.rm=TRUE), cpue=mean(cpue, na.rm=TRUE)), by=.(region, yr) ]
        summary_annual = summary_annual[ which(!is.na(summary_annual$region)) , ]

        summary_annual$effort =  summary_annual$landings / summary_annual$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
        summary_annual = summary_annual[ order(region, yr),]
 #       setDF(summary_annual)

        summary_annual$landings = summary_annual$landings / 1000  # t
        summary_annual$effort = summary_annual$effort / 1000  # 1000 th

        # due to landings occurring in GULF region, some tweaks here (probably better in landings db): 
        # ii = which( summary_annual$region=="cfanorth" & summary_annual$yr==2013 ) 
        # if (length(ii) ==1) summary_annual$cpue[ii ] = 106 
        # ii = which( summary_annual$region=="cfanorth" & summary_annual$yr==2014 ) 
        # if (length(ii) ==1) summary_annual$cpue[ii ] = 104.5 
        # 
        # please update this (it is embeded in a function right now but needs to be moved to a CSV or DB lookup)
        tacs = snowcrab_tacs() #TODO

        summary_annual = merge(summary_annual, tacs, by=c("region", "yr"), all.x=TRUE, all.y=FALSE)
        summary_annual$landings = round(summary_annual$landings )
        summary_annual$effort = round(summary_annual$effort, 1 )
        summary_annual$cpue = round(summary_annual$cpue )
        summary_annual$TAC = round(summary_annual$TAC)
        summary_annual$Licenses = as.numeric(summary_annual$Licenses)
        summary_annual = summary_annual[order(summary_annual$region, summary_annual$yr),]
 
        out$summary_annual = summary_annual
    }
    

    if (any( grepl("summary_monthly", toget)) ) {
        # monthly stats
        smon = NA
        
        fstat$month = month( as.POSIXct(fstat$date.landed) )

        smon = fstat[, .(landings=sum(landings, na.rm=TRUE), cpue=mean(cpue, na.rm=TRUE)), by=.(region, yr, month) ]
            
        smon = smon[ which(!is.na(smon$region)) , ]

        smon$effort =  smon$landings / smon$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    #       setDF(smon)

        smon$landings = smon$landings / 1000  # t
        smon$effort = smon$effort / 1000  # 1000 th

        smon$landings = round(smon$landings )
        smon$effort = round(smon$effort, 1 )
        smon$cpue = round(smon$cpue )
 
        # add breaks in aug for 4x to plot correctly
        smon = rbind( smon, CJ( region="cfa4x", month=8, yr=y, cpue=NA, landings=NA, effort=NA))
        
        # add breaks  for north to plot correctly: month=6 (separate spring-fall)
        smon = rbind( smon, CJ( region="cfanorth", month=6.5, yr=y, cpue=NA, landings=NA, effort=NA))

        smon = smon[ order(region, yr, month),]
        out$summary_monthly = smon

        if (0) {
            
            region = "cfanorth"
            region = "cfasouth"
            region = "cfa4x"
        
            i = which(smon$region==region & is.finite(smon$month) )
            plot( landings ~ month, smon[ i,], col=smon$yr[i] )
            for (j in y) {
                k = which(smon$region==region & smon$yr==j)
                lines( landings ~ month, smon[k,], col=j )
            }

            i = which(smon$region==region & is.finite(smon$month) )
            plot( effort ~ month, smon[ i,], col=smon$yr[i] )
            for (j in y) {
                k = which(smon$region==region & smon$yr==j)
                lines( effort ~ month, smon[k,], col=j )
            }

            i = which(smon$region==region & is.finite(smon$month) )
            plot( cpue ~ month, smon[ i,], col=smon$yr[i] )
            for (j in y) {
                k = which(smon$region==region & smon$yr==j)
                lines( cpue ~ month, smon[k,], col=j )
            }

        }
    }

    if (any( grepl("summary_weekly", toget)) ) {
        # weekly stats
        swk = NA
        
        fstat$week = week( as.POSIXct(fstat$date.landed) )

        swk = fstat[, .(landings=sum(landings, na.rm=TRUE), cpue=mean(cpue, na.rm=TRUE)), by=.(region, yr, week) ]
            
        swk = swk[ which(!is.na(swk$region)) , ]

        swk$effort =  swk$landings / swk$cpue  ## estimate effort level as direct estimates are underestimates (due to improper logbook records)
    #       setDF(swk)

        swk$landings = swk$landings / 1000  # t
        swk$effort = swk$effort / 1000  # 1000 th

        swk$landings = round(swk$landings )
        swk$effort = round(swk$effort, 1 )
        swk$cpue = round(swk$cpue )

        # add breaks in aug for 4x to plot correctly
        swk = rbind( swk, CJ( region="cfa4x", week=floor(8/12*52)+0.5, yr=y, cpue=NA, landings=NA, effort=NA))

       # add breaks in aug for north to plot correctly
        swk = rbind( swk, CJ( region="cfanorth", week=26.5, yr=y, cpue=NA, landings=NA, effort=NA))

        swk = swk[ order(region, yr, week),]
        out$summary_weekly = swk

        if (0) {
            
            region = "cfanorth"
            region = "cfasouth"
            region = "cfa4x"

            i = which(swk$region==region & is.finite(swk$week) )
            plot( landings ~ week, swk[ i,], col=swk$yr[i] )
            for (j in y) {
                k = which(swk$region==region & swk$yr==j)
                lines( landings ~ week, swk[k,], col=j )
            }

            i = which(swk$region==region & is.finite(swk$week) )
            plot( effort ~ week, swk[ i,], col=swk$yr[i] )
            for (j in y) {
                k = which(swk$region==region & swk$yr==j)
                lines( effort ~ week, swk[k,], col=j )
            }

            i = which(swk$region==region & is.finite(swk$week) )
            plot( cpue ~ week, swk[ i,], col=swk$yr[i] )
            for (j in y) {
                k = which(swk$region==region & swk$yr==j)
                lines( cpue ~ week, swk[k,], col=j )
            }


        }
    }


    if ( any( grepl("shell_condition", toget)) | any( grepl("fraction_observed", toget)) )  {

        male = 0
        odb = observer.db("odb")
        setDT(odb)
        odb = odb[ which( odb$sex==male & odb$cw >= 95 & odb$cw < 170 & odb$prodcd_id=="0" & is.finite(odb$shell) ) ,]  # commerical sized crab only
        odb$region = NA
        for ( reg in regions ) {
            r = polygon_inside(x=odb, region=aegis.polygons::polygon_internal_code(reg), planar=FALSE)
            if (length(r)> 0)  odb$region[r] = reg
        }
    }

    if ( any( grepl("shell_condition", toget)) ) { 
        shell_condition = odb[ !is.na(odb$region), .N, by=.(region, fishyr, shell) ]
        shell_condition[, total:=sum(N, na.rm=TRUE), by=.(region, fishyr)]
        shell_condition$percent = round(shell_condition$N / shell_condition$total, 3) * 100
        out$shell_condition=shell_condition
    }

    if (any( grepl("fraction_observed", toget)) ) {

      if (!exists("summary_annual", out)) stop( "summary_annual not found")

      if (!exists("totmass_kept", odb)) odb$totmass_kept = odb$totmass  # temporary fix until observer db raw data gets refreshed (JC 2022)
      message("Make sure all observer local raw tables have been refreshed")
      
      fraction_observed = odb[ !is.na(odb$region), .(no_traps=.N, kept=sum(totmass_kept, na.rm=TRUE), caught=sum(totmass, na.rm=TRUE) ), by=.(region, fishyr) ]
      setnames(fraction_observed, "fishyr", "yr")

      fraction_observed = fraction_observed[ out$summary_annual, on=.(region, yr)]
      # in percentages:

      fraction_observed$observed_landings_pct = round( fraction_observed$kept / fraction_observed$landings *100, 2)  # values in metric tonnes
      fraction_observed$observed_effort_pct = round( (fraction_observed$no_traps/1000) / fraction_observed$effort * 100,2)  #  effort is in 1000 th

      fraction_observed$observed_landings_pct[ which(!is.finite(fraction_observed$observed_landings_pct)) ] = 0
      fraction_observed$observed_effort_pct[ which(!is.finite(fraction_observed$observed_effort_pct)) ] = 0

      out$fraction_observed = fraction_observed 

    }


    return(out)

}