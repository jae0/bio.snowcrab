 

size_distributions = function(
    p=p, 
    outdir=file.path(p$project.outputdir, "size_structure"), 
    toget="base_data",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    xrange=c(10, 160),
    np = 512,
    dx=NULL,
    bw=2,
    kernel="gaussian",
    density_offset = NULL,
    redo=FALSE,
    add_zeros=FALSE,
    pg=NULL,
    moving_average=FALSE,
    ti_window=c(-4,4),
    sigdigits=2,
    lowpassfilter=0,
    lowpassfilter2=0.001,
    plot_solutions = FALSE,
    tlevels=c(-2, 6),
    zlevels=c(0, 100),
    Y=NULL ) { 

    if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 
     
    if (0) {
        regions=c("cfanorth", "cfasouth", "cfa4x")
        xrange=c(0, 160)
        dx = 2
    }
    
    # sex codes
    # male = 0
    # female = 1
    # sex.unknown = 2

    # # maturity codes
    # immature = 0
    # mature = 1
    # mat.unknown = 2
 

    if (toget=="base_data") {

        fn = file.path( outdir, "size_distributions_base_data.RDS" )
        Z = NULL
        if (!redo) {
            if (file.exists(fn)) {
                Z = readRDS(fn)
                if (is.null(xrange)) xrange = attr(Z, "xrange") # leave alone
                if (is.null(dx)) dx = attr(Z, "dx") # leave alone
                breaks = seq(xrange[1], xrange[2], by=dx)
                mids = breaks[-length(breaks)] + dx/2
                Z$cwd = discretize_data( Z$cw, brks=breaks, labels=mids, resolution=dx )  
                Z = Z[ is.finite(cwd) ,]
                attr(Z, "xrange") = xrange
                attr(Z, "dx") = dx
            }
            return(Z)
        }
        set = snowcrab.db( DS="set.clean")
        det = snowcrab.db( DS="det.initial")
        setDT(set)
        setDT(det)
        set$sid = paste(set$trip, set$set, sep="~")
        det$sid = paste(det$trip, det$set, sep="~")
        set$year = set$yr 
        set$region = NA
        for ( region in regions ) {
            r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
            if (length(r) > 0) set$region[r] = region
        }
        set$space_id = NA
        Z = sf::st_as_sf( set[,.(lon, lat)], coords=c("lon", "lat") )
        st_crs(Z) = st_crs( projection_proj4string("lonlat_wgs84") )
        for (aoi in 1:nrow(pg)) {
            ks = which(!is.na( st_points_in_polygons(pts=Z, polys=pg[aoi, "AUID"], varname= "AUID" ) ))
            if (length(ks) > 0 ) set$space_id[ks] = pg$AUID[aoi]
        }
        set= set[!is.na(region), ]
        set = set[, .(sid, region, space_id, year, sa, t, z, timestamp, julian, lon, lat)]
        det = det[, .(sid, shell, cw, sex, mass, mat, gonad, durometer)]
        Z = det[ set, on=.(sid)]
        
        # trim a few strange data points
        o = lm( log(mass) ~ log(cw), Z)
        todrop = which(abs(o$residuals) > 0.5)
        Z = Z[-todrop,]

        Z$shell = factor( Z$shell )

        attr(Z, "xrange") = xrange
        attr(Z, "dx") = dx

        print(fn)
        saveRDS(Z, file=fn, compress=TRUE)
        return(Z)
    }


    if (toget=="tabulated_data") {
        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        fn = file.path( outdir, "size_distributions_tabulated_data.RDS" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) {
                M =readRDS(fn)
                if (!is.null(Y)) M = M[ year %in% Y, ]
                if (add_zeros) {
                    # merge zeros here so we do not have to store the massive intermediary file
                    # CJ required to get zero counts dim(N) # 171624960    
                    M = M[ CJ( region, year, sex, mat, cwd, sid, unique=TRUE ), 
                        on=.( region, year, sex, mat, cwd, sid ) ]
                }
                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
                return(M)
            }
        }
        Z = size_distributions(p=p, toget="base_data", xrange=xrange, dx=dx )
        # aggregate by cwd 
        M = Z[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, sex, mat, cwd, sid) ]
        Z = NULL
        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$cwd = as.factor(M$cwd)
        M$region = as.factor(M$region)
        saveRDS(M, file=fn, compress=TRUE)
        # return this way to add zeros, if required
        return( size_distributions(p=p, toget="tabulated_data", add_zeros=add_zeros, redo=FALSE ) )
    }


    if (toget == "simple_direct" ) {
        # no linking across time, beak by year to reduces ram use
        savedir = file.path(outdir, "simple_direct")
        if (!dir.exists(savedir)) dir.create(savedir, recursive=TRUE, showWarnings =FALSE) 
        if (!redo) {
            M = NULL
            if (!is.vector(Y)) stop("Y should be a year vector")
            for (yr in as.character(Y)) {
                fn = file.path( savedir, paste("size_distributions_simple_direct_", yr, ".RDS" ))
                if (file.exists(fn)) {
                    m = NULL
                    m = readRDS(fn)
                    m$year = yr
                    M = rbind( M, m)
                }
            }
            return(M)
        }

        Z = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, add_zeros=FALSE  )
        # NOTE without offset, this implicitly drops the zeros
        if (is.null(density_offset)) density_offset = min( Z$density[ which(Z$density>0) ] )
        message( "Density offset: ", density_offset )
        Z = NULL; gc()
        for (yr in as.character(Y)) {
            M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, Y=yr, add_zeros=TRUE  )
            M$log_den = log(M$density + density_offset) 
            M = M[ ,         
                .(  nsamples = .N,
                    number_mean = mean( N, na.rm=TRUE ),
                    number_sd = sd( N, na.rm=TRUE ),
                    sa_mean = mean( sa, na.rm=TRUE ),
                    sa_sd = sd( sa, na.rm=TRUE ),
                    mass = mean( mass, na.rm=TRUE),
                    mass_sd = sd( mass, na.rm=TRUE),
                    den   = mean(density, na.rm=TRUE), 
                    den_sd =sd(density, na.rm=TRUE), 
                    den_lb=mean(density, na.rm=TRUE) - 1.96*sd(density, na.rm=TRUE),
                    den_ub=mean(density, na.rm=TRUE) + 1.96*sd(density, na.rm=TRUE),
                    denl     = exp(mean(log_den, na.rm=TRUE))-density_offset, 
                    denl_log = log(exp(mean(log_den, na.rm=TRUE))-density_offset), 
                    den_sd_log = sd(log_den, na.rm=TRUE), 
                    den_lb_log=exp(mean(log_den, na.rm=TRUE)-density_offset - 1.96*sd(log_den, na.rm=TRUE)),
                    den_ub_log=exp(mean(log_den, na.rm=TRUE)-density_offset + 1.96*sd(log_den, na.rm=TRUE))
                ), 
                by= .( region, sex, mat, cwd)
            ]
            attr(M, "density_offset") = density_offset
            fn = file.path( savedir, paste("size_distributions_simple_direct_", yr, ".RDS" ))
            saveRDS(M, file=fn, compress=TRUE)
        }

        return(size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=Y, redo=FALSE))     
    }



    if (toget=="linear_model") {
        require(biglm)
        fn = file.path( outdir, "size_distributions_lm.RDS" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, Y=Y, add_zeros=TRUE )
        setDT(M)
        fit = biglm::bigglm( density ~ region:year:mat:cwd:sex - 1, data=M, family=gaussian(link="identity") )
        O = summary(fit)$coefficients
        res = tstrsplit(rownames(O), ":")
        setDT(res)
        setnames(res, new=c("region", "year", "mat", "cwd", "sex"))
        res[,region:=as.numeric(gsub( "region", "", region))]
        res[,year:=as.factor(gsub( "year", "", year))]
        res[,mat:=as.factor(gsub( "mat", "", mat))]
        res[,cwd:=as.numeric(gsub( "cwd", "", cwd))]
        res[,sex:=as.numeric(gsub( "sex", "", sex))]
        O = data.table(O)
        colnames(O) = c("density",  "density_se", "t", "p")
        O = cbind( res, O )
        saveRDS(O, file=fn, compress=TRUE )
        return(O)
    }


    if (toget=="poisson_glm") {
        stop("This takes far too long to use ... ")
        fn = file.path( outdir, "size_distributions_poisson_glm.RDS" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, Y=Y, add_zeros=TRUE )
        M$ID = as.factor( paste( M$region, M$year, M$sex, M$mat, M$cwd, sep="_") )
        # subset 
        ss = M[ region=="cfanorth" & sex=="0" & year %in% as.character(2015:2022), which=TRUE]
        # NOTE this is too large of a problem for glm
        fit = biglm::bigglm( N ~ ID - 1 +offset(log_sa), data=M[ss,], 
            family=poisson(link="log"), na.action="na.omit" )
        P = data.table( ID = names(coef(fit)), mean=coef(fit) )
        nm = matrix( unlist( strsplit(P$ID, "_")), ncol=5, byrow=TRUE)
        P = cbind( P, nm)
        names(P) = c("ID", "N", "region", "year", "sex", "mat", "cwd")
        P$region = gsub("^ID", "", P$region)
        O = list( fit=fit, P=P )
        saveRDS(O, file=fn, compress=TRUE )
        return(O)
    }


    if (toget=="poisson_inla") {
        stop("This takes far too long to use ... ")
        require(INLA)
        fn = file.path( outdir, "size_distributions_inla_poisson.RDS" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, add_zeros=TRUE )
        M$tag ="o"
        P = CJ( 
            N = NA,
            log_sa = 0,   # log(1) ... 1km^2
            region = regions,
            year = p$yrs,  
            sex = c("0", "1"), 
            mat = c("0", "1"),
            cwd = levels(M$cwd),
            tag= "p"
        )
        pn = names(P)
        M = copy( rbind( M[, ..pn ], P ) ) # a deep copy and not a reference
        # M$ID = paste( M$region, M$year, M$sex, M$mat, M$cwd, sep="_")
        # subset 
        ss = M[ region=="cfanorth" & sex=="0" & year %in% as.character(2015:2022), which=TRUE]
        fit = inla( N ~ year + mat + cwd + year:mat:cwd - 1 + offset(log_sa), data=M[ss,], 
            family="poisson", verbose=TRUE)
        iP = M[tag=="p", which=TRUE]
        P$N = fit$summary.fitted.values$mean[iP]
        P$Nsd = fit$summary.fitted.values$sd[iP]
        P$Nlb = fit$summary.fitted.values$"0.025quant"[iP]
        P$Nub = fit$summary.fitted.values$"0.975quant"[iP]
        O = list( fit=fit, P=P )
        saveRDS(O, file=fn, compress=TRUE )
        return(O)
    }


    if (toget=="kernel_density_weighted") {
  
        file_prefix = "kernel_densities"
        if (moving_average) file_prefix = "kernel_densities_moving_average"

        outdir = file.path( p$project.outputdir, "size_structure", paste(file_prefix, "_", round(bw,3), "_", np, sep="") )
        dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 
 
        xr = round( log(xrange), digits=2 ) 
        if(is.null(dx)) dx = diff(xr)/(np-1)
        xvals = seq( xr[1], xr[2], by=dx )

        if (!redo) {
            if (is.null(Y)) stop("year Y must be provided")
            M = NULL
            for (yr in as.character(Y) ) {
                fn = file.path( outdir, paste( file_prefix, "_", yr, ".csv", sep="" ) )
                if (file.exists(fn)) {
                    varnames <- try( data.table::fread(fn, nrows = 1, header = FALSE), silent = FALSE)
                    if (inherits(varnames, "try-error")) return(NULL)
                    if (ncol(varnames) <= 1) return(NULL)  # no data
                    m = data.table::fread( fn, header=TRUE )
                    if (colnames(m)[1] == "V1") m[[1]] = NULL 
                    setDT(m)
                    if (exists("X", m)) m$X = NULL  # should not be necessary .. but just in case
                }
                if (!is.null(m)) M = rbind(M, m)
            }

            M$sex = as.character(M$sex)
            M$mat = as.character(M$mat)
            M$au = as.character(M$au)
            xvals = attributes(M)$xvals

            if (exists("sid", M)) {
                set = snowcrab.db( DS="set.clean")
                set$sid = paste(set$trip, set$set, sep="~")
                setDT(set)
                set = set[ , .(sid, t, z, lon, lat ) ]
                set$region = NA
                for ( region in regions ) {
                    r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
                    if (length(r) > 0) set$region[r] = region
                }
                set$zi = cut( set$z, breaks=c(zlevels, Inf ), labels=zlevels )   
                set$ti = cut( set$t, breaks=c(tlevels, Inf ), labels=tlevels )   
                M = set[M, on="sid"]
            }

            attr(M, "xrange") = xrange
            attr(M, "xvals") = xvals
            attr(M, "xr") = xr
            attr(M, "bw") = bw
            attr(M, "dx") = dx

            return(M)            
        }

        sigdigits=3
        
        M = size_distributions( p=p, toget="base_data") # , xrange=xrange, dx=dx  not sent due to not being relevant
        M$sex = as.character(M$sex)
        M$mat = as.character(M$mat)
        M$logcw = log(M$cw)
        M$wt = 1.0 / M$sa

        yrs = sort( unique( M$year ) ) # 1996:present
        if (!is.null(Y)) yrs = as.numeric(Y)

        sexes = c("0", "1")  # 0 ,1, 2 male, female, unknown
        mats = c("0", "1")   # 0 ,1, 2  imm, mat, unknown

        nbs = attributes(pg)$nb$nbs

        if (moving_average) {
            # weekly basis
            M$ti = M$year + round(trunc(M$julian / 365 * 52 ) / 52, digits=sigdigits)         # discretize time quarterly
            # if (0) { yr=2022; wk=40; auid="360"; sx="1"; mt="1" }
            
            for (yr in yrs) {
                out1 = NULL
                out2 = NULL
                # print(yr)
                for (wk in 1:52) {
                    mti1 = yr + (wk + ti_window[1])/52
                    mti2 = yr + (wk + ti_window[2])/52
                    kt = M[ ti >= mti1 & ti <= mti2, which=TRUE ]
                    if ( length(kt) < 1)  next() 
                    auids = na.omit(unique(M$space_id[kt]))
                    for (auid in auids) {
                        ii = which(pg$AUID == auid)
                        if (length(ii) != 1 ) next()  # NA's  should not happen
                        aus = unique( c( pg$AUID[nbs[[ ii ]]], auid ) )
                        ka = intersect( kt, M[ space_id %in% aus, which=TRUE ] )
                        if (length(ka) < 1) next()
                        for (sx in sexes) {
                            ks = intersect( ka, M[ sex==sx, which=TRUE ] )
                            if (length(ks) < 1) next()
                            for (mt in mats) {
                                n =  intersect( ks, M[ mat==mt, which=TRUE] )
                                N =  length(n) 
                                if (N < 1) next() 
                                tout = paste("| sex: ", sx, "| mat: ", mt, "| au: ", auid, "|year: ", yr, "| week: ", wk, "| N: ", N ) 
                                message(tout )
                                uu = try( density( M$logcw[n], bw=bw, kernel=kernel, from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE ))
                                if (inherits(uu, "class-error")) next()
                                uu$y = uu$y / sum(uu$y) / dx  # density
                                
                                out1 = rbind( out1, data.table( sex=sx, mat=mt, au=auid, year=yr, wk=wk, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                                out2 = rbind( out2, data.table( t(uu$y))  )
                                res = NULL
                            } # mat
                        } # sex
                    }  # au
                }
                if (is.null(out1) | is.null(out2)) next()
                out = cbind(out1, out2)
                fnout  = file.path( outdir, paste( file_prefix, "_", yr, ".csv", sep="" ) )
                data.table::fwrite( cbind(out1, out2), file=fnout )
                print(fnout ) 
            }

        } else {
            # quarterly basis if (0) { yr=2022; season=40; au="360"; sex="1"; mat="1" }
            M$ti = M$year + round(trunc(M$julian / 365 * 4 ) / 4, digits=sigdigits)         # discretize time quarterly
            seasons = seq(0, 0.75, by=0.25)  #            
            for (yr in yrs) {
                out1 = NULL
                out2 = NULL
                for (season in seasons) {
                    mti = yr + season
                    kt = M[ ti == mti, which=TRUE ]
                    if ( length(kt) < 1)  next() 
                    sids = unique(M$sid[kt])
                    for (si in sids) {
                        ka = intersect( kt, M[ sid == si, which=TRUE ] )
                        if (length(ka) < 1) next()
                        au = unique(M$space_id[ka])[1]
                        for (sx in sexes) {
                            ks = intersect( ka, M[ sex==sx, which=TRUE ] )
                            if (length(ks) < 1) next()
                            for (mt in mats) {
                                n =  intersect( ks, M[ mat==mt, which=TRUE] )
                                N =  length(n) 
                                if (N < 1) next() 

                                tout = paste("|sid: ", si, "| sex: ", sx, "| mat: ", mt, "| au: ", au, "|year: ", yr, "| season: ", season, "| N: ", N ) 
                                message(tout )
                                uu = try( density( M$logcw[n], bw=bw, kernel=kernel, from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE ))
                                if (inherits(uu, "class-error")) next()
                                uu$y = uu$y / sum(uu$y) / dx  # density
                                
                                out1 = rbind( out1, data.table( sid=si, sex=sx, mat=mt, au=au, year=yr, season=season, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                                out2 = rbind( out2, data.table( t(uu$y))  )
                                res = NULL
                            } # mat
                        } # sex
                    }  # au
                }   # seasons  
                if (is.null(out1) | is.null(out2)) next()
                out = cbind(out1, out2)
                fnout  = file.path( outdir, paste( file_prefix, "_", yr, ".csv", sep="" ) )
                data.table::fwrite( cbind(out1, out2), file=fnout )
                print(fnout ) 
            }
        }
        return ( size_distributions(p=p, toget="kernel_density_weighted", 
            moving_average=moving_average,
            pg=pg, ti_window=ti_window, 
            bw=bw, np=np, xrange =xrange, Y=Y, redo=FALSE ))
    }

  
 
}

 

 
