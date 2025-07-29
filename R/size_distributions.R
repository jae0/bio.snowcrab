 

size_distributions = function(
    p=p, 
    outdir=file.path(p$project.outputdir, "size_structure"), 
    toget="",
    M=NULL,
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    # xrange=c(10, 150),
    np = 512,
    # dx=NULL,
    span = NULL,
    ldx=NULL,
    bw=2,
    kernel="gaussian",   
    redo=FALSE,
    add_zeros=FALSE,
    pg=NULL,
    pg_subset=NULL, 
    ti_window=c(-4,4),
    n_min=1,
    sigdigits=2,
    strata="yasm",
    lowpassfilter=0,
    lowpassfilter2=0.001,
    plot_solutions = FALSE,
    tlevels=c(-2, 6),
    zlevels=c(0, 100),
    Y=NULL ) { 

    if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 
      

    # note ranges in CW will be log transformed later
    if (is.null(span)) {
        span = function( sexid) {
            switch(sexid,
                male   = c( 5, 155, 40),
                female = c( 5, 95,  40)
            )
        }
    }
        # lspan[1]  = log(span[1]) 
        # lspan[2]  = log(span[2]) 

    # c( xrange[1], xrange[2], nx) 
    #   nx = floor( abs(diff(xrange) / dx) )
        
        

    # sex codes
    # male = 0
    # female = 1
    # sex.unknown = 2

    # # maturity codes
    # immature = 0
    # mature = 1
    # mat.unknown = 2
 


    if (toget == "rawdata" ) {
 
        if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 

        if (!redo) {
            M = NULL 
            fn = file.path( outdir, paste("size_distributions_rawdata", ".rdz", sep="" ))
            if (file.exists(fn)) {
                M = aegis::read_write_fast(fn) 
            }
            return(M)
        }
        
        set = snowcrab.db( DS="set.clean")
        setDT(set)
        set$sid = paste(set$trip, set$set, sep="~")
        set$year = set$yr 
        set$region = NA
        for ( region in regions ) {
            r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
            if (length(r) > 0) set$region[r] = region
        }
        set= set[!is.na(region), ]
        set = set[, .(sid, region, year, sa, t, z, timestamp, julian, lon, lat)]


        det = snowcrab.db( DS="det.initial")
        setDT(det)
        det$sid = paste(det$trip, det$set, sep="~")
        det = det[, .(sid, shell, cw, sex, mass, mat, gonad, durometer, chela, abdomen)]
        det = det[ mat %in% c("0", "1") & sex %in% c("0", "1"),]


        M = set[ det, on=.(sid)]
        
        det = set = NULL; gc()
         
        # trim a few strange data points
        o = lm( log(mass) ~ log(cw) + factor(sex) + factor(mat), M)
        todrop = which(abs(o$residuals) > 0.5)
        if (length(todrop)>0) {
            M[todrop, "mass"] = NA  # 17 in 2024
        }

            # these are global parameters
            # # sex codes
            # male = 0
            # female = 1
            # sex.unknown = 2

            # # maturity codes
            # immature = 0
            # mature = 1
            # mat.unknown = 2

        # female large
        todrop = M[ sex=="1" & cw > 90, which=TRUE ]  
        if (length(todrop)>0) {
            sex_female = which(is.finite(M$abdomen[todrop]))  # must be female if there is abdomen
            sex_male = which(!is.finite(M$abdomen[todrop]))  # must be male ? assumed
            
            M$sex[todrop][sex_male] = "0"  # recode to male 
            M = M[-todrop[sex_female],]
        }

        # female large imm
        # todrop = M[ sex=="1" & cw > 80 & mat=="0", which=TRUE ]  
        # if (length(todrop)>0) M = M[-todrop,]

        # female large mat
        # todrop = M[ sex=="1" & cw > 90 & mat=="1", which=TRUE ]  
        # if (length(todrop)>0) M = M[-todrop,]

        # female small mat
        todrop = M[ sex=="1" & cw <35 & mat=="1", which=TRUE ]  
        if (length(todrop)>0) M = M[-todrop,]


        # male large imm
        # todrop = M[ sex=="0" & cw > 135 & mat=="0", which=TRUE ]  
        # if (length(todrop)>0) M = M[-todrop,]

        # male large mat
        # todrop = M[ sex=="0" & cw > 150 & mat=="1", which=TRUE ]  
        # if (length(todrop)>0) M = M[-todrop,]

        # male small mat
        todrop = M[ sex=="0" & cw <49 & mat=="1", which=TRUE ]  
        if (length(todrop)>0) M = M[-todrop,]

        M$shell = factor( M$shell )
          
        fn = file.path( outdir, paste("size_distributions_rawdata", ".rdz", sep="" ))
        read_write_fast( data=M, fn=fn )
        return(M)
    }



    if (toget == "crude" ) {
 
        savedir = file.path(outdir, "crude")
        if (!dir.exists(savedir)) dir.create(savedir, recursive=TRUE, showWarnings =FALSE) 

        if (!redo) {
            M = NULL
            if (!is.vector(Y)) stop("Y should be a year vector")
            for (yr in as.character(Y)) {
                fn = file.path( savedir, paste("size_distributions_crude_", yr, ".rdz", sep="" ))
                if (file.exists(fn)) {
                    m = NULL
                    m = aegis::read_write_fast(fn)
                    m$year = yr
                    M = rbind( M, m)
                }
            }
            return(M)
        }
       
        P = size_distributions(p=p, toget="rawdata", outdir=outdir)
        
        sexid = list(male = "0", female = "1" )
        P$cwd = NA

        for (j in c("male", "female")) {
            k = which( P$sex==sexid[[j]] )
            P$cwd[k] = discretize_data( x=P$cw[k], span=p$span(j) )  
        }

        P = P[ is.finite(cwd) ,]

        if (!is.null(Y)) P = P[ year %in% Y, ]
        
        # aggregate by cwd 
        Z = P[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, sex, mat, cwd, sid) ]
        
        P = NULL;gc()

        Z$year = as.factor(Z$year)
        Z$region = as.factor(Z$region)
        Z$cwd = as.factor(Z$cwd)

        # merge zeros here so we do not have to store the massive intermediary file
        # CJ required to get zero counts dim(N) # 171624960    
        Z = Z[ CJ( region, year, sex, mat, cwd, sid, unique=TRUE ), 
                on=.( region, year, sex, mat, cwd, sid ) ]
        Z[ !is.finite(N),   "N"] = 0
        Z[ !is.finite(mass), "mass"] = 0 
        Z[ !is.finite(sa), "sa"] = 1 #dummy value
        Z$density = Z$N / Z$sa
        Z[ !is.finite(density), "density"] = 0  

        for (yr in as.character(Y)) {
            M = Z[ year==yr, ]
            M$log_density = log(M$density)  
            M = M[ ,         
                .(  nsamples = length(which(N>0)),
                    number_mean = mean( N, na.rm=TRUE ),
                    number_sd = sd( N, na.rm=TRUE ),
                    sa_mean = mean( sa, na.rm=TRUE ),
                    sa_sd = sd( sa, na.rm=TRUE ),
                    mass = mean( mass, na.rm=TRUE),
                    mass_sd = sd( mass, na.rm=TRUE),
                    den     = mean( density, na.rm=TRUE), 
                    den_sd  = sd( density, na.rm=TRUE), 
                    den_log_geo    = geometric_mean_sd( log_density, "mean" )  ,
                    den_log_geo_sd = geometric_mean_sd( log_density, "sd" ) 
                ), 
                by= .( region, sex, mat, cwd)
            ]

            M$den_lb = M$den - 1.96*M$deb_sd
            M$den_ub = M$den + 1.96*M$deb_sd
            
            M$denl = exp(M$den_log_geo)
            M$denl_lb = exp(M$den_log_geo - (1.96*M$den_log_geo_sd))
            M$denl_ub = exp(M$den_log_geo + (1.96*M$den_log_geo_sd))

            # capture NaN's due to no counts
            i = which(!is.finite(M$den_log_geo))
            if (length(i)>0) M$den_log_geo[i] = 0

            i = which(!is.finite(M$den_log_geo_sd))
            if (length(i)>0) M$den_log_geo_sd[i] = 0


            i = which(!is.finite(M$denl))
            if (length(i)>0) M$denl[i] = 0

            i = which(!is.finite(M$denl_lb))
            if (length(i)>0) M$denl_lb[i] = 0


            i = which(!is.finite(M$denl_ub))
            if (length(i)>0) M$denl_ub[i] = 0


            fn = file.path( savedir, paste("size_distributions_crude_", yr, ".rdz", sep="" ))
            read_write_fast( data=M, fn=fn )
        }

        return(size_distributions(p=p, toget="crude", span=span, Y=Y, outdir=outdir, redo=FALSE))     
    }



    if (toget=="base_data") {
        # georeferenced "crude" data
        savedir = file.path(outdir, "base_data")
        if (!dir.exists(savedir)) dir.create(savedir, recursive=TRUE, showWarnings =FALSE) 

        if (!redo) {
            M = NULL
            if (!is.vector(Y)) stop("Y should be a year vector")
            for (yr in as.character(Y)) {
                fn = file.path( savedir, paste("size_distributions_base_data_", yr, ".rdz", sep="" ))
                if (file.exists(fn)) {
                    m = NULL
                    m = aegis::read_write_fast(fn)
                    m$year = yr
                    M = rbind( M, m)
                }
            }
            return(M)
        }
 
        # add georeference to pg

        set = snowcrab.db( DS="set.clean")
        setDT(set)
        set$sid = paste(set$trip, set$set, sep="~")
        set = set[, .(sid, lon, lat)]
        set$space_id = NA
        Z = sf::st_as_sf( set[,.(lon, lat)], coords=c("lon", "lat") )
        st_crs(Z) = st_crs( projection_proj4string("lonlat_wgs84") )
        for (aoi in 1:nrow(pg)) {
            ks = which(!is.na( st_points_in_polygons(pts=Z, polys=pg[aoi, "AUID"], varname= "AUID" ) ))
            if (length(ks) > 0 ) set$space_id[ks] = pg$AUID[aoi]
        }
        set = set[, .(sid, space_id)]

        for (yr in Y) {
            S = copy(set)
            basedata = size_distributions(p=p, toget="crude", outdir=outdir, Y=yr)
            basedata = S[ basedata, on=.(sid)]
            fn = file.path( savedir, paste("size_distributions_base_data_", yr, ".rdz", sep="" ))
            read_write_fast( data=basedata, fn=fn )
            print(fn)
        }

        return(size_distributions(p=p, toget="base_data", span=span, Y=Y, outdir=outdir, redo=FALSE))
    }


    if (toget=="tabulated_data") {
        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        fn = file.path( outdir, "size_distributions_tabulated_data.rdz" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) {
                M =aegis::read_write_fast(fn)
                if (!is.null(Y)) M = M[ year %in% Y, ]
                if (add_zeros) {
                    # merge zeros here so we do not have to store the massive intermediary file
                    # CJ required to get zero counts dim(N) # 171624960    
                    M = M[ CJ( sex, mat, cwd, sid, unique=TRUE ), 
                         on=.( sex, mat, cwd, sid ) ]
                }
                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
                return(M)
            }
        }

        M = size_distributions(p=p, toget="base_data", span=span, outdir=outdir )

        # aggregate by cwd 
        M = M[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, sex, mat, cwd, sid) ]
        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$cwd = as.factor(M$cwd)
        read_write_fast( data=M, fn=fn )
        # return this way to add zeros, if required
        return( size_distributions(p=p, toget="tabulated_data", outdir=outdir, add_zeros=add_zeros, redo=FALSE ) )
    }



    if (toget=="tabulated_data_by_stage") {
        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        fn = file.path( outdir, "size_distributions_tabulated_data_by_stage.rdz" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) {
                M =aegis::read_write_fast(fn)
                if (!is.null(Y)) M = M[ year %in% Y, ]
                if (add_zeros) {
                    # merge zeros here so we do not have to store the massive intermediary file
                    # CJ required to get zero counts dim(N) # 171624960    
                    M = M[ CJ( stage, sid, unique=TRUE ), 
                        on=.(  stage, sid ) ]
                }
                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
                return(M)
            }
        }
        M = size_distributions(p=p, toget="base_data", span=span, outdir=outdir )
        # aggregate by cwd 

        mds = size_distributions(p=p, toget="modal_groups", outdir=outdir, redo=FALSE )
        
        M$stage = filter.stage( M, mds ) 
        M = M[!is.na(stage),]

        M = M[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, stage, sid) ]
        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$stage = as.factor(M$stage)
        read_write_fast( data=M, fn=fn )
        # return this way to add zeros, if required
        return( size_distributions(p=p, toget="tabulated_data_by_stage", outdir=outdir, add_zeros=add_zeros, redo=FALSE ) )
    }

 


    if (toget=="linear_model") {
        require(biglm)
        fn = file.path( outdir, "size_distributions_lm.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", span=span, Y=Y, outdir=outdir, add_zeros=TRUE )
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
        read_write_fast( data=O, fn=fn )
        return(O)
    }


    if (toget=="poisson_glm") {
        stop("This takes far too long to use ... ")
        fn = file.path( outdir, "size_distributions_poisson_glm.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", span=span, Y=Y, outdir=outdir, add_zeros=TRUE )
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
        read_write_fast( data=O, fn=fn )
        return(O)
    }


    if (toget=="poisson_inla") {
        stop("This takes far too long to use ... ")
        require(INLA)
        fn = file.path( outdir, "size_distributions_inla_poisson.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_distributions(p=p, toget="tabulated_data", span=span, outdir=outdir, add_zeros=TRUE )
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
        read_write_fast( data=O, fn=fn )
        return(O)
    }


    if (toget=="kernel_density_weighted") {
  
        kdtype = paste( "kernel_densities", strata, np, sep="_" )

        outdir = file.path( p$project.outputdir, "size_structure", kdtype )
        dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 
 
        fn0 = file.path( outdir, "kernel_density_weighted.rdz" )

        xr = round( log(xrange), digits=2 ) 
        if (is.null(ldx)) ldx = diff(xr)/(np-1)
        xvals = seq( xr[1], xr[2], by=ldx )
    
        if (!redo) {

            if (is.null(Y)) {
                # if no years, then current saved version:
                O = NULL
                if (!redo) {
                    if (file.exists(fn)) O =aegis::read_write_fast(fn0)
                    return(O)
                }
            }

            M = NULL
            
            for (yr in as.character(Y) ) {
                fn   = file.path( outdir, paste( "kd_", yr, ".rdz", sep="" ) )
                if (file.exists(fn)) {
                    m = NULL
                    m = read_write_fast( fn=fn )
                    # setDT(m)
                    # if (exists("X", m)) m$X = NULL  # should not be necessary .. but just in case
                    if (!is.null(m)) M = rbind(M, m)
                }
            }

            M$sex = as.character(M$sex)
            M$mat = as.character(M$mat)
            M$au = as.character(M$au)

            attr(M, "xrange") = xrange
            attr(M, "xvals") = xvals
            attr(M, "xr") = xr
            attr(M, "bw") = bw
            attr(M, "ldx") = ldx
            attr(M, "np") = np
            attr(M, "ti_window") = ti_window
            attr(M, "pg") = pg
            attr(M, "sigdigits") = sigdigits
            
            read_write_fast( data=M, fn=fn0 )

            return(M)            
        }

      
        M = size_distributions( p=p, toget="base_data", outdir=outdir) # , span=span  not sent due to not being relevant
        M$sex = as.character(M$sex)
        M$mat = as.character(M$mat)
        M$logcw = log(M$cw)
        M$wt = 1.0 / M$sa

        yrs = sort( unique( M$year ) ) # 1996:present
        if (!is.null(Y)) yrs = as.numeric(Y)

        sexes = c("0", "1")  # 0 ,1, 2 male, female, unknown
        mats = c("0", "1")   # 0 ,1, 2  imm, mat, unknown

        nbs = attributes(pg)$nb$nbs  # immediate neighbours only
     
        if (strata=="yasm") {
            # weekly basis
            M$ti = M$year + round(trunc(M$julian / 365 * 52 ) / 52, digits=sigdigits)         # discretize time quarterly
            # if (0) { yr=2022; wk=40; auid="360"; sx="1"; ma="1" }
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
                            for (ma in mats) {
                                n =  intersect( ks, M[ mat==ma, which=TRUE] )
                                N =  length(n) 
                                if (N < 1) next() 
                                tout = paste("| sex: ", sx, "| mat: ", ma, "| au: ", auid, "|year: ", yr, "| week: ", wk, "| N: ", N ) 
                                message(tout )
                                uu = try( density( M$logcw[n], bw=bw[[sx]][[ma]], kernel=kernel, from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE ))
                                if (inherits(uu, "class-error")) next()
                                uu$y = uu$y / sum(uu$y) / ldx  # density
                                out1 = rbind( out1, data.table( sex=sx, mat=ma, au=auid, year=yr, wk=wk, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                                out2 = rbind( out2, data.table( t(uu$y))  )
                                res = NULL
                            } # mat
                        } # sex
                    }  # au
                }
                if (is.null(out1) | is.null(out2)) next()
                out = cbind(out1, out2)
                fnout  = file.path( outdir, paste( "kd_", yr, ".rdz", sep="" ) )
                read_write_fast( data=out, fn=fnout )
                print(fnout ) 
            }

        } else if (strata=="smryzt")  {

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
                            for (ma in mats) {
                                n =  intersect( ks, M[ mat==ma, which=TRUE] )
                                N =  length(n) 
                                if (N < 1) next() 

                                tout = paste("|sid: ", si, "| sex: ", sx, "| mat: ", ma, "| au: ", au, "|year: ", yr, "| season: ", season, "| N: ", N ) 
                                message(tout )
                                uu = try( density( M$logcw[n], bw=bw[[sx]][[ma]], kernel=kernel, from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE ))
                                if (inherits(uu, "class-error")) next()
                                uu$y = uu$y / sum(uu$y) / ldx  # density
                                
                                out1 = rbind( out1, data.table( sid=si, sex=sx, mat=ma, au=au, year=yr, season=season, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                                out2 = rbind( out2, data.table( t(uu$y))  )
                                res = NULL
                            } # mat
                        } # sex
                    }  # au
                }   # seasons  
                if (is.null(out1) | is.null(out2)) next()
                out = cbind(out1, out2)
                fnout  = file.path( outdir, paste( "kd_", yr, ".rdz", sep="" ) )
                read_write_fast( data=out, fn=fnout )
                print(fnout ) 
            }
        }
        return ( size_distributions(p=p, toget="kernel_density_weighted", strata=strata, outdir=outdir, 
            pg=pg, ti_window=ti_window,  sigdigits=sigdigits,  
            bw=bw, np=np, xrange =xrange, Y=Y, redo=FALSE ))
    }
 
    # -----------------

    if (toget=="kernel_density_modes") { 

        fn = file.path( outdir, paste("size_distributions_summary_", strata, ".rdz", sep="") )
        
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }

        # spatial window is nearest-neighbours in spatial graph
 
        M = size_distributions(p=p, toget="kernel_density_weighted", 
          bw=bw, np=np, ldx=ldx, xrange=xrange, Y=years, strata=strata, pg=pg, 
          sigdigits=sigdigits, ti_window=ti_window,  outdir=outdir )   
   
        xrange =   attr(M, "xrange")
        xvals =   attr(M, "xvals")
        xr =   attr(M, "xr")
        bw =   attr(M, "bw")
        ldx =   attr(M, "ldx")
        np =   attr(M, "np")
        ti_window = attr(M, "ti_window")
        pg = attr(M, "pg")
        sigdigits = attr(M, "sigdigits")  

        zlevels=c(0, 100)
        tlevels= c(-2, 6)
        sexes=c("0", "1")
        mats=c("0", "1") 
        
        peaks = data.table()
        troughs = data.table()
        peak_values = data.table()
        trough_values = data.table()

        if (strata=="smryzt") {
            aus=c("cfanorth", "cfasouth", "cfa4x")
            
            if (exists("sid", M)) {
                set = snowcrab.db( DS="set.clean")
                set$sid = paste(set$trip, set$set, sep="~")
                setDT(set)
                set = set[ , .(sid, t, z, lon, lat ) ]
                set$region = NA
                for ( region in aus ) {
                    r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
                    if (length(r) > 0) set$region[r] = region
                }
                set$zi = cut( set$z, breaks=c(zlevels, Inf ), labels=zlevels )   
                set$ti = cut( set$t, breaks=c(tlevels, Inf ), labels=tlevels )   
                M = set[M, on="sid"]
            }

            K = aggregate_by( M, 
                agg_by = c("year", "sex", "mat", "region", "zi", "ti" ),  # strata
                xvals= xvals,
                recale_density_to_numerical_density=TRUE,  ### keep normalized to reduce scale issues
                agg_function = function(x) {exp(mean( log(x), na.rm=TRUE) ) }, # geometric_mean 
                add_offset=TRUE 
            )

            for ( s in sexes ) {
            for ( m in mats ) {
            for ( r in aus) {
            for ( y in years ) {
            for ( z in zlevels ) {
            for ( t in tlevels ) {
                vn = paste(y,s,m,r,z,t, sep="_" )
                if (!exists(vn, K)) next()

                mds = identify_modes( Z=as.vector(t(K[, ..vn])),
                    sigdigits=sigdigits, 
                    lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2, 
                    dx=ldx, X=xvals,
                    n_min=n_min, plot_solutions=TRUE )   
                if (is.null(mds)) next()
                if (inherits(mds, "try-error")) next()
                # out[[s]][[m]][[r]][[y]][[z]][[t]] = mds
                peaks = rbind(peaks, cbind(s, m, r, y, z, t, t(t(mds[["peaks"]])) ))
                troughs = rbind(peaks, cbind(s, m, r, y, z, t, t(t(mds[["troughs"]])) ))
                peak_values = rbind(peak_values, cbind(s, m, r, y, z, t, t(t(mds[["peak_values"]])) ))
                trough_values = rbind(trough_values, cbind(s, m, r, y, z, t, t(t(mds[["trough_values"]])) ))
            }}} }}}
            
            setnames(peaks, "V7", "peaks")
            setnames(troughs, "V7", "troughs")
            setnames(peak_values, "V7", "peak_values")
            setnames(trough_values, "V7", "trough_values")

        } else if (strata=="yasm" ) {
       
            aus=pg$AUID
            K = aggregate_by( M, 
                agg_by = c( "year", "au", "sex", "mat" ),  # strata
                xvals= xvals,
                recale_density_to_numerical_density=TRUE,  ### keep normalized to reduce scale issues
                agg_function = function(x) {exp(mean( log(x), na.rm=TRUE) ) }, # geometric_mean 
                add_offset=TRUE 
            )
            for ( s in sexes ) {
            for ( m in mats ) {
            for ( a in aus) {
            for ( y in years ) {
                vn = paste(y,a,s,m, sep="_" )
                if (!exists(vn, K)) next()
                mds = identify_modes( Z=as.vector(t(K[, ..vn])),
                  sigdigits=sigdigits, 
                    lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2, 
                    dx=ldx, X=xvals,
                    n_min=n_min, plot_solutions=TRUE)   
                if (is.null(mds)) next()
                if (inherits(mds, "try-error")) next()
                peaks = rbind(peaks, cbind(s, m, a, y, t(t(mds[["peaks"]])) ))
                troughs = rbind(peaks, cbind(s, m, a, y, t(t(mds[["troughs"]])) ))
                peak_values = rbind(peak_values, cbind(s, m, a, y, t(t(mds[["peak_values"]])) ))
                trough_values = rbind(trough_values, cbind(s, m, a, y, t(t(mds[["trough_values"]])) ))
            }}} }
    
            setnames(peaks, "V5", "peaks")
            setnames(troughs, "V5", "troughs")
            setnames(peak_values, "V5", "peak_values")
            setnames(trough_values, "V5", "trough_values")
            
        }
  
        peaks$peaks = as.numeric(peaks$peaks)
        peak_values$peak_values = as.numeric(peak_values$peak_values)
        troughs$troughs = as.numeric(troughs$troughs)
        trough_values$trough_values = as.numeric(trough_values$trough_values)
 
        O = list()
        O$peaks=peaks
        O$peak_values=peak_values
        O$troughs=troughs
        O$trough_values=trough_values
   
        vn = "peaks"
        wn = "peak_values"
        out = NULL
        dists = NULL
 
        # no aus ( agg across all space) .. mode of modes
        for (yr in years) {
        for (sx in c("0", "1")) {
        for (ma in c("0", "1")) {
            Z = unlist(O[[vn]][ s==sx & m==ma & y==yr & a %in% aus, ..vn])
            W = unlist(O[[wn]][ s==sx & m==ma & y==yr & a %in% aus, ..wn])
            if (length(Z) < 1) next()
            mds = identify_modes( 
                Z=Z, # W=W,
                n_min=n_min, 
                lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2,
                xvals=xvals, dx=ldx, bw=bw[[sx]][[ma]], sigdigits=sigdigits, plot=TRUE) 
            if (is.na(mds$N)) next()

            out = rbind( out, data.table( cw=mds$peaks, mat=ma, sex=sx, year=yr) )
            dists = rbind( dists, data.table( 
                cw=mds$u$x, 
                density=mds$u$y, 
                N=mds$N, 
                mat=ma, sex=sx, year=yr) )
            
        }}}
        
        O[["ysm"]] = list(peaks=out, densities=dists)
        
        out = NULL
        dists = NULL
        
        if (strata=="yasm") {
            
            for (yr in years) {
            for (sx in c("0", "1")) {
            for (ma in c("0", "1")) {
            for (au in aus) {

            Z = unlist(O[[vn]][ s==sx & m==ma & y==yr & a %in% au, ..vn])
            W = unlist(O[[wn]][ s==sx & m==ma & y==yr & a %in% au, ..wn])

            if (length(Z) < 1) next()
            mds = identify_modes( 
                Z=Z, # W=W,
                n_min=n_min, 
                lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2,
                xvals=xvals, dx=ldx, bw=bw[[sx]][[ma]], sigdigits=sigdigits, plot=TRUE) 
            if (is.na(mds$N)) next()

            out = rbind( out, data.table( cw=mds$peaks, mat=ma, sex=sx, year=yr, auid=au) )
            dists = rbind( dists, data.table( 
                cw=mds$u$x, 
                density=mds$u$y, 
                N=mds$N, 
                mat=ma, sex=sx, year=yr, auid=au) )
            
            }}}}
         
        } else if (strata=="smryzt") {
            
            regions = unique(O[[vn]]$r)
            tis = unique(O[[vn]]$t)
            zis = unique(O[[vn]]$z)

            for (yr in years) {
            for (sx in c("0", "1")) {
            for (ma in c("0", "1")) {
            for (re in regions ) {
            for (tt in tis) {
            for (zz in zis) {
            Z = unlist(O[[vn]][ s==sx & m==ma & y==yr & r==re & t==tt & z==zz, ..vn])
            if (length(Z) < 1) next()

            mds = identify_modes( 
                Z=Z,  n_min=n_min, 
                lowpassfilter=0.0, lowpassfilter2=0,
                xvals=xvals, dx=ldx, bw=bw[[sx]][[ma]], sigdigits=sigdigits, plot=TRUE) 
            if (is.na(mds$N)) next()

            out = rbind( out, data.table( cw=mds$peaks, mat=ma, sex=sx, year=yr) )
            dists = rbind( dists, data.table( 
                cw=mds$u$x, 
                density=mds$u$y, 
                N=mds$N, 
                mat=ma, sex=sx, year=yr, region=re, temp=tt, depth=zz) )
            
            }}}}}}
 
        }

        O[[strata]] = list(peaks=out, densities=dists)
 
        read_write_fast( data=O, fn=fn )
        return(O)
    } 

    # ---------------

    if ( toget %in% c("modal_groups", "modal_groups_models") ) {
        
        if (!redo) {
            if (toget =="modal_groups") {
                mds = NULL
                fn = file.path(outdir, "modes_summary.rdz")
                if (file.exists(fn)) mds = read_write_fast(fn)
                if (!is.null(mds)) return (mds)
            }
            if (toget =="modal_groups_models") {
                mds_models = NULL
                fn = file.path(outdir, "modes_models.rdz")
                if (file.exists(fn)) mds_models = aegis::read_write_fast(fn)
                if (!is.null(mds_models)) return (mds_models)
            }        
        }
        
        message("These solutions need to be checked carefully and cleaned appropriately ...")

        survey_size_freq_dir = file.path( p$annual.results, "figures", "size.freq", "survey")

        if (is.null(M)) M = size_distributions(p=p, toget="kernel_density_modes", strata=strata, outdir=outdir )
        
        MI = M[["ysm"]][["densities"]]
        MO = M[["ysm"]][["peaks"]]

        fn = file.path(survey_size_freq_dir, "modes_male_imm_allsolutions.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot(density~cw, MI[sex=="0" & mat=="0" , ], pch="." )
            abline(v=MO[ sex=="0" & mat=="0", cw ], col="gray", lwd=0.5 )
        dev.off()
        print(fn)

        fn = file.path(survey_size_freq_dir, "modes_male_mat_allsolutions.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot(density~cw, MI[sex=="0" & mat=="1" , ], pch=".")  # NOTE misses the largest size group
            abline(v=MO[ sex=="0" & mat=="1", cw ], col="gray" )
        dev.off()
        print(fn)

        fn = file.path(survey_size_freq_dir, "modes_female_imm_all_solutions.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot(density~cw, MI[sex=="1" & mat=="0" , ], pch=".")
            abline(v=MO[ sex=="1" & mat=="0", cw ], col="gray" )
        dev.off()
        print(fn)

        fn = file.path(survey_size_freq_dir, "modes_female_mat_all_solutions.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot(density~cw, MI[sex=="1" & mat=="1" , ], pch=".")
            abline(v=MO[ sex=="1" & mat=="1", cw ], col="gray" )
        dev.off()
        print(fn)

        # collect point estimates 
 
        fn = file.path(survey_size_freq_dir, "modes_male_imm.png" )
        png(filename=fn, width=1000,height=600, res=144)
        mi = identify_modes( Z = unlist(MO[ sex=="0" & mat=="0" , cw]),  
            lowpassfilter2=lowpassfilter2, xvals=xvals, dx=ldx, bw=bw[["0"]][["0"]], sigdigits=3, plot=TRUE) 
        abline(v=4, col="orange", lwd=2, lty="dashed") # likely a nonmode
        dev.off()
        print(fn)

        fn = file.path(survey_size_freq_dir, "modes_male_mat.png" )
        png(filename=fn, width=1000,height=600, res=144)
        mm = identify_modes( Z = unlist(MO[ sex=="0" & mat=="1" , cw]),  
            lowpassfilter2=lowpassfilter2, xvals=xvals, dx=ldx, bw=bw[["0"]][["1"]], sigdigits=3, plot=TRUE) 
        dev.off()
        print(fn)

        fn = file.path(survey_size_freq_dir, "modes_female_imm.png" )
        png(filename=fn, width=1000,height=600, res=144)
        fi = identify_modes( Z = unlist(MO[ sex=="1" & mat=="0" , cw]),  
            lowpassfilter2=lowpassfilter2, xvals=xvals, dx=ldx, bw=bw[["1"]][["0"]], sigdigits=3, plot=TRUE) 
        dev.off()
        print(fn)
    
        fn = file.path(survey_size_freq_dir, "modes_female_mat.png" )
        png(filename=fn, width=1000,height=600, res=144)
        fm = identify_modes( Z = unlist(MO[ sex=="1" & mat=="1" , cw]),  
            lowpassfilter2=lowpassfilter2, xvals=xvals, dx=ldx, bw=bw[["1"]][["1"]], sigdigits=3, plot=TRUE) 
        dev.off()
        print(fn)
    
        mds = rbind( 
            data.table( logcw = fi$peaks,   sex="f", mat= "i" ),
            data.table( logcw = fm$peaks,   sex="f", mat= "m" ),
            data.table( logcw = mi$peaks,   sex="m", mat= "i" ),
            data.table( logcw = mm$peaks,   sex="m", mat= "m" )
        )
        mds$cw = exp(mds$logcw)

        if (0) {
            # check if there are any strangemodes and remove
            plot( mds$logcw[ mds$sex=="f"])
            plot( mds$logcw[ mds$sex=="m"])  
        }

        # bad = mds[  logcw > 3.95 & logcw < 4.05 & sex=="m" & mat=="i", which=TRUE ]
        # if (length(bad)>0) mds = mds[-bad,]

        f = mds[ sex=="f", ][order(mat, cw),]
        f$seq = 1:nrow(f)

        fn = file.path(survey_size_freq_dir, "modes_female_growth_trajectory_empirical.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot( cw ~ seq, f)
            i =4:6  # hyp: imm just under corresponding mature size
            arrows(f$seq[i], f$cw[i], f$seq[i+3], f$cw[i+3], length=0.2, col= 1:3)
            i = f[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(f$seq[i], f$cw[i], f$seq[i+1], f$cw[i+1], length=0.2 )
        dev.off()
        print(fn)
        

        m = mds[ sex=="m", ][order(mat, cw),]
        m$seq = 1:nrow(m)


        fn = file.path(survey_size_freq_dir, "modes_male_growth_trajectory_empirical.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot( cw ~ seq, m)
            i = 5:7 # hyp: imm just under corresponding mature size
            arrows(m$seq[i], m$cw[i], m$seq[i+3], m$cw[i+3], length=0.2, col= 1:3)
            i = m[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(m$seq[i], m$cw[i], m$seq[i+1], m$cw[i+1], length=0.2 )
            # last immature group -> maturity is missing .. add it below
        dev.off()
        print(fn)
        

        # assign instar: imm patterns seems simple
        mds$instar = NA
        
        # female
        ii = mds[sex=="f" & mat=="i", which=TRUE]
        mds$instar[ii] = cw_to_instar( mds$logcw[ii], "f" ) 

        jj = mds[sex=="f" & mat=="m", which=TRUE]
        for (j in jj) {
            k = mds[sex=="f" & mat=="i" & logcw < mds$logcw[j] , which=TRUE]
            mds$instar[j] = max( mds$instar[ k] ) + 1
        }

     
        fn = file.path(survey_size_freq_dir, "modes_female_growth_trajectory_empirical_tweaked.png" )
        png(filename=fn, width=1000,height=600, res=144)
           # verify:
            f = mds[ sex=="f", ][order(mat, cw),]
            plot( cw ~ instar, f)
            j = f[ mat=="m", which=TRUE]
            for (i in j ){
                k = f[mat=="i" & instar==(f$instar[i] -1), which=TRUE ] 
                arrows(f$instar[i], f$cw[i], f$instar[k], f$cw[k], length=0.2, code=1, col="red")
            }
            i = f[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(f$instar[i], f$cw[i], f$instar[i+1], f$cw[i+1], length=0.2 )
        dev.off()
        print(fn)


        # male
        ii = mds[sex=="m" & mat=="i", which=TRUE]
        mds$instar[ii] = cw_to_instar( mds$logcw[ii], "m" ) 

        jj = mds[sex=="m" & mat=="m", which=TRUE]
        for (j in jj) {
            k = mds[sex=="m" & mat=="i" & logcw < mds$logcw[j] , which=TRUE]
            mds$instar[j] = max( mds$instar[ k] ) + 1
        }

 
     
        fn = file.path(survey_size_freq_dir, "modes_male_growth_trajectory_empirical_tweaked.png" )
        png(filename=fn, width=1000,height=600, res=144)
            # verify:
            m = mds[ sex=="m", ][order(mat, cw),]
            plot( cw ~ instar, m)
            j = m[ mat=="m", which=TRUE]
            for (i in j ){
                k = m[mat=="i" & instar==(m$instar[i] -1), which=TRUE ] 
                arrows(m$instar[i], m$cw[i], m$instar[k], m$cw[k], length=0.2, code=1, col="blue")
            }
            i = m[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(m$instar[i], m$cw[i], m$instar[i+1], m$cw[i+1], length=0.2 )
        dev.off()
        print(fn)

 
        mds$pred = NA  # pred bade on logcw0 and maturity
        mds$logcw0 = NA

        i = mds[sex=="f" & mat=="i", which=TRUE ] 
        for (j in i) {
            k = mds[sex=="f" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
            if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
        }
        
        i = mds[sex=="f" & mat=="m", which=TRUE ] 
        for (j in i) {
            k = mds[sex=="f" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
            if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
        }

        i = mds[sex=="m" & mat=="i", which=TRUE ] 
        for (j in i) {
            k = mds[sex=="m" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
            if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
        }
        
        i = mds[sex=="m" & mat=="m", which=TRUE ] 
        for (j in i) {
            k = mds[sex=="m" & mat=="i" & instar==(mds[j, instar]-1), which=TRUE ]
            if (length(k) >0) mds$logcw0[j] = mds$logcw[k]
        }


        of = lm( logcw ~ logcw0 * mat,  mds[ sex=="f",], na.action="na.omit")
        mds$pred[ which(mds$sex=="f")] = predict( of, mds[ sex=="f",] )
        summary(of)
        
        if (0) {
            plot(logcw~logcw0,  mds[ sex=="f",])
            points(logcw~logcw0,  mds[ sex=="f" & mat=="i",], col="red")
            points(pred~logcw0,  mds[ sex=="f" & mat=="i",], col="green")
            points(pred~logcw0,  mds[ sex=="f" & mat=="m",], col="purple")
        }

        om = lm( logcw ~ logcw0 * mat,  mds[ sex=="m",], na.action="na.omit")
        

        mds$pred[ which(mds$sex=="m")] = predict( om, mds[ sex=="m",] )
        summary(om)
        
        if (0) {
            plot(pred~logcw0,  mds[ sex=="m",])
            points(logcw~logcw0,  mds[ sex=="m" & mat=="m",], col="red", pch="+")
            points(logcw~logcw0,  mds[ sex=="m" & mat=="i",], col="red")
            points(pred~logcw0,  mds[ sex=="m" & mat=="i",], col="green")
            points(pred~logcw0,  mds[ sex=="m" & mat=="m",], col="purple", pch="+")
        }    

        # add unobserved instars: 1:4 and 13 Male
        oif = lm( logcw ~ instar, mds[sex=="f" & mat=="i", ], na.action="na.omit")
        omf = lm( logcw~ instar, mds[sex=="f" & mat=="m", ], na.action="na.omit")
        summary(oif) # Adjusted R-squared:  0.999
        summary(omf) # Adjusted R-squared:  0.977

        oim = lm( logcw~ instar, mds[sex=="m" & mat=="i", ], na.action="na.omit")
        omm = lm( logcw~ instar, mds[sex=="m" & mat=="m", ], na.action="na.omit")
        summary(oim) # Adjusted R-squared:  0.999
        summary(omm) # Adjusted R-squared:  0.999

        unobs= CJ(logcw=NA, sex=c("m", "f"), mat="i", cw=NA, instar=1:3, pred=NA, logcw0=NA)

        mds = rbind(mds, unobs)
        # mature male instar not represented (due to rarity vs instar 12)
        # add instar 13
        logcw12 = mds[sex=="m" & mat=="i" & instar==12, logcw]
        instar13 = data.table( NA, "m", "m", NA, 13, NA, logcw12)
        names(instar13) = names(mds)
        mds = rbind(mds, instar13 )
        
        mds$predicted = NA
        mds$predicted_se = NA

        i = mds[sex=="f" & mat=="i", which=TRUE]
        ip = predict(oif, mds[i,], se.fit=TRUE )
        mds$predicted[i] =ip$fit
        mds$predicted_se[i] =ip$se.fit

        i = mds[sex=="f" & mat=="m", which=TRUE]
        ip = predict(omf, mds[i,], se.fit=TRUE )
        mds$predicted[i] =ip$fit
        mds$predicted_se[i] =ip$se.fit
    
        i = mds[sex=="m" & mat=="i", which=TRUE]
        ip = predict(oim, mds[i,], se.fit=TRUE )
        mds$predicted[i] =ip$fit
        mds$predicted_se[i] =ip$se.fit

        i = mds[sex=="m" & mat=="m", which=TRUE]
        ip = predict(omm, mds[i,], se.fit=TRUE )
        mds$predicted[i] =ip$fit
        mds$predicted_se[i] =ip$se.fit

       # full predicted pattern:
        mds$cwmean = exp(mds$predicted)
        mds$cwlb = exp(mds$predicted - 1.96*mds$predicted_se )
        mds$cwub = exp(mds$predicted + 1.96*mds$predicted_se )
    
        fn = file.path(survey_size_freq_dir, "modes_growth_female.png" )
        png(filename=fn, width=1000, height=600, res=144)
            f = mds[ sex=="f", ][order(mat, instar),]
            plot( cwmean ~ instar, f, type="p" )
            j = f[ mat=="m", which=TRUE]
            for (i in j ){
                k = f[mat=="i" & instar==(f$instar[i] -1), which=TRUE ] 
                arrows(f$instar[i], f$cwmean[i], f$instar[k], f$cwmean[k], length=0.2, code=1, col="red")
            }
            i = f[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(f$instar[i], f$cwmean[i], f$instar[i+1], f$cwmean[i+1], length=0.2 )
        dev.off()
    
        fn = file.path(survey_size_freq_dir, "modes_growth_male.png" )
        png(filename=fn, width=1000,height=600, res=144)
            m = mds[ sex=="m", ][order(mat, instar),]
            plot( cwmean ~ instar, m, type="p" )
            j = m[ mat=="m", which=TRUE]
            for (i in j ){
                k = m[mat=="i" & instar==(m$instar[i] -1), which=TRUE ] 
                arrows(m$instar[i], m$cwmean[i], m$instar[k], m$cwmean[k], length=0.2, code=1, col="blue")
            }
            i = m[ mat=="i", which=TRUE]
            i = i[-length(i)]
            arrows(m$instar[i], m$cwmean[i], m$instar[i+1], m$cwmean[i+1], length=0.2 )
        dev.off()


        mds$diff = exp(mds$logcw) - exp(mds$logcw0) 
        mds$diff_prop = mds$diff / exp(mds$logcw0)  # fractional increase
        mds$diffp = exp(mds$pred) - exp(mds$logcw0) 

        fn = file.path(survey_size_freq_dir, "modes_growth_increment.png" )
        png(filename=fn, width=1000,height=600, res=144)
            plot(diff_prop ~ instar, mds, type="n")
            points(diff_prop ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=19 )
            points(diff_prop ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=21 )
            points(diff_prop ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=19 )
            points(diff_prop ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=21 )
        dev.off()

        if (0) {
            mean(mds$diff_prop[ mds$sex=="f"], na.rm=TRUE)  # 30 % increase each moult
            mean(mds$diff_prop[ mds$sex=="m"], na.rm=TRUE)  # 30 % increase each moult

            mean(mds$diff_prop[ mds$sex=="f" & mds$mat=="i" ], na.rm=TRUE)  # 36.7 % increase each moult
            mean(mds$diff_prop[ mds$sex=="m" & mds$mat=="i" ], na.rm=TRUE)  # 34.9 % increase each moult

            mean(mds$diff_prop[ mds$sex=="f" & mds$mat=="m" ], na.rm=TRUE)  # 17 % increase each moult
            mean(mds$diff_prop[ mds$sex=="m" & mds$mat=="m" ], na.rm=TRUE)  # 16 % increase each moult

            plot(diff ~ instar, mds, type="n")
            points(diff ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=19 )
            points(diff ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=21 )
            points(diff ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=19 )
            points(diff ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=21 )

            points(diffp ~ instar, mds[sex=="f" & mat=="i",], col="orange", pch=23 )
            points(diffp ~ instar, mds[sex=="f" & mat=="m",], col="darkred", pch=24 )
            points(diffp ~ instar, mds[sex=="m" & mat=="i",], col="green", pch=23 )
            points(diffp ~ instar, mds[sex=="m" & mat=="m",], col="darkblue", pch=24 )

        }

        # this is to distinguish between same maturity groups, unlike instar dtermination above
        mds$logcw_predicted0 = NA
        instar0 = mds$instar-1

        i = mds[sex=="f" & mat=="i", which=TRUE ] 
        mds$logcw_predicted0[i] = predict(oif, data.table(instar=instar0[i]) ) 
        
        i = mds[sex=="f" & mat=="m", which=TRUE ] 
        mds$logcw_predicted0[i] = predict(omf, data.table(instar=instar0[i]) ) 

        i = mds[sex=="m" & mat=="i", which=TRUE ] 
        mds$logcw_predicted0[i] = predict(oim, data.table(instar=instar0[i]) ) 
       
        i = mds[sex=="m" & mat=="m", which=TRUE ] 
        mds$logcw_predicted0[i] = predict(omm, data.table(instar=instar0[i]) ) 

        mds$diff_predicted =   mds$predicted - mds$logcw_predicted0
        mds$predicted_lb =  mds$predicted - mds$diff_predicted /2
        mds$predicted_ub =  mds$predicted + mds$diff_predicted /2
    
        instar = paste("0", mds$instar, sep="")
        instar = substring( instar, nchar(instar)-1, nchar(instar))
        mds$stage = paste(mds$sex, mds$mat, instar, sep="|")

        # save as rdata for use in julia
        fn = file.path(outdir, "modes_summary.rdz")
        read_write_fast(mds, fn=fn)

        mds_models  = list(oif, oim, omf, omm, of, om)
        fn = file.path(outdir, "modes_models.rdz")
        read_write_fast(mds_models, fn=fn)
    
        return(mds)
    }

 

    # ---------------

    if ( toget %in% c("modal_groups_carstm_inputs") ) {


      outputdir = file.path( p$modeldir, p$carstm_model_label )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )
    
      fn = file.path( outputdir, "carstm_inputs.rdz" )

      M = NULL
   
      if (!redo) {
          if (file.exists(fn)) M = aegis::read_write_fast(fn)
          if (!is.null(M)) return (M)
      }
    
      ## Base data
      mds = size_distributions(p=p, toget="modal_groups", outdir=outdir )

      M = size_distributions(p=p, toget="base_data", pg=pg, span=span, outdir=outdir )
      M$id = gsub("~", ".", M$sid)
      M = M[ year %in% p$yrs, ]

      M$cwd = discretize_data( x=M$cw, span=span )  

      M$mat[ M$mat=="2" & M$shell != "1" ] = "1"  # override
   
      M$stage = filter.stage( M, mds ) 
      M = M[!is.na(stage),]
 
      # aggregate by cwd 
      M = M[, .( N=.N ), by=.( id, sex, mat, cwd ) ]
      M = M[ CJ( id, sex, mat, cwd, unique=TRUE ), on=.( id, sex, mat, cwd ) ]
      M$N[ which(is.na(M$N))] = 0

      set = snowcrab.db( DS="set.clean")
      set$id = paste(set$trip, set$set, sep=".")
      setDT(set)
      set = set[ , .(id, t, z, lon, lat, plon, plat, yr, timestamp, julian, sa ) ]
      
    #   set$region = NA
    #   for ( region in c("cfanorth", "cfasouth", "cfa4x") ) {
    #       r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=FALSE)
    #       if (length(r) > 0) set$region[r] = region
    #   }

      set$space_id = NA
      Z = sf::st_as_sf( set[,.(lon, lat)], coords=c("lon", "lat") )
      st_crs(Z) = st_crs( projection_proj4string("lonlat_wgs84") )
      for (aoi in 1:nrow(pg)) {
          ks = which(!is.na( st_points_in_polygons(pts=Z, polys=pg[aoi, "AUID"], varname= "AUID" ) ))
          if (length(ks) > 0 ) set$space_id[ks] = pg$AUID[aoi]
      }

      M = set[M, on="id"]
      set=NULL

      M = M[ !is.na(space_id), ]
      M = M[ !is.na(yr), ]
 
      M$year = factor(M$yr)
     
      M$sa[ which(!is.finite(M$sa))] = 1 # dummy value
 
      # some survey timestamps extend into January (e.g., 2020) force them to be part of the correct "survey year", i.e., "yr"
      i = which(lubridate::month(M$timestamp)==1)
      if (length(i) > 0) M$timestamp[i] = M$timestamp[i] - lubridate::duration(month=1)

      M$tiyr = lubridate::decimal_date(M$timestamp)
   
      M$dyear = lubridate::decimal_date( M$timestamp ) - M$yr
      M$dyear[ M$dyear > 1] = 0.99  # a survey year can run into the next year, cap the seasonal compenent at the calendar year for modellng 

      # reduce size
      M = M[ which( M$lon > p$corners$lon[1] & M$lon < p$corners$lon[2]  & M$lat > p$corners$lat[1] & M$lat < p$corners$lat[2] ), ]
      # levelplot(z.mean~plon+plat, data=M, aspect="iso")

      # should already have this information ... just in case 
      require(aegis.survey)
      pSU = survey_parameters( yrs=1999:p$year.assessment )
      SU = survey_db( DS="set", p=pSU ) 
      
      oo = match( M$id, SU$id )
      
      mz = which(!is.finite(M$z))
      if (length(mz) > 0 ) M$z[mz] = SU$z[oo[mz]]

      mt = which(!is.finite(M$t))
      if (length(mt) > 0 ) M$t[mt] = SU$t[oo[mt]]
  
      pSU = SU = NULL
      gc()
      
      # data_offset is SA in km^2
      M$data_offset = M$sa 
      M$data_offset[which(!is.finite(M$data_offset))] = median(M$data_offset, na.rm=TRUE )  # just in case missing data
      M = M[ which(  is.finite(M$data_offset)   ),  ]

#  dim(M) # [1] 409472     20
      
      # basic space-time expansion
      M = carstm_prepare_inputdata( 
        p=p, M=M, sppoly=pg, 
        APS_data_offset=1, 
        retain_positions_outside_of_boundary = 25,  # centroid-point unit of p$aegis_proj4string_planar_km
        vars_to_retain=c("id", "N", "sa",  "cwd", "mat", "sex", "data.offset", "t", "z") 
      )

# dim(M) # 416696     19


      # expand also for a few other items:
      obs = M[tag=="observations",]
      aps = M[tag=="predictions",]
# dim(obs) # 404936     19
# dim(aps) #  11760    19

      aps$sex = NULL
      aps$mat = NULL
      aps$cwd = NULL

      n_aps = nrow(aps)
      mats = c("0", "1")
      n_mat = length(mats)
      aps = cbind( aps[ rep.int(1:n_aps, n_mat), ], rep.int( mats, rep(n_aps, n_mat )) )
      names(aps)[ncol(aps)] = "mat"

      n_aps = nrow(aps)
      sexes = c("0", "1")
      n_sex = length(sexes)
      aps = cbind( aps[ rep.int(1:n_aps, n_sex), ], rep.int( sexes, rep(n_aps, n_sex )) )
      names(aps)[ncol(aps)] = "sex"

      n_aps = nrow(aps)
      cwds = mids
      n_cwd = length(cwds)
      aps = cbind( aps[ rep.int(1:n_aps, n_cwd), ], rep.int( cwds, rep(n_aps, n_cwd )) )
      names(aps)[ncol(aps)] = "cwd"
 
      nms = intersect( names(obs), names(aps) )
      M = rbind(obs[,..nms], aps[, ..nms])


#      setDF(M)
      #   these vars being missing means zero-valued
      #   vars_to_zero = c( "density" )
      #   for ( vn in vars_to_zero ) {
      #     if (exists( vn, M)) {
      #       i = which( is.na(M[, vn] ) )
      #       if (length(i) > 0) M[i, vn] = 0 
      #     }
      #   }

      if ( exists("substrate.grainsize", M)) M$log.substrate.grainsize = log( M$substrate.grainsize )

      if (!exists("yr", M)) M$yr = M$year  # req for meanweights

      # IMPERATIVE: 
      i =  which(!is.finite(M$z))
      j =  which(!is.finite(M$t)) 

      if (length(j)>0 | length(i)>0) {
        warning( "Some areal units that have no information on key covariates ... you will need to drop these and do a sppoly/nb reduction with areal_units_neighbourhood_reset() :")
            print( "Missing depths:")
        print(unique(M$AUID[i]) )
        print( "Missing temperatures:")
        print(unique(M$AUID[j] ) )
      }


      if (0) {
        # Note used right now but if addtional survey data from groundfish used ...
        # predictions to: westeren 2a and NED
        gears_ref = "Nephrops"
        i = which(is.na(M$gear)) 
        M$gear[ i ] = gears_ref
        gears = unique(M$gear[-i])
        gears = c( gears_ref, setdiff( gears, gears_ref ) ) # reorder
        M$gear = as.numeric( factor( M$gear, levels=gears ) )
        attr( M$gear, "levels" ) = gears

        M$vessel = substring(M$id,1,3)
        M$id = NULL 

        vessels_ref = "xxxx"
        i = which(is.na(M$vessel) )
        M$vessel[ i ] = vessels_ref
        vessels = unique(M$vessel[-i])
        vessels = c( vessels_ref, setdiff( vessels, vessels_ref ) ) # reorder
        M$vessel= as.numeric( factor( M$vessel, levels=vessels ) )
        attr( M$vessel, "levels" ) = vessels
      }

  
      if (0) {
        # drop data without covariates 
        i = which(!is.finite( rowSums(M[, .(z, t, pca1, pca2 ) ] )) )
        if (length(i) > 0 ) {
          au = unique( M$AUID[i] )
          j = which( M$AUID %in% au )
          if (length(j) > 0 ) {

            plot( pg["npts"] , reset=FALSE, col=NA )
            plot( pg[j, "npts"] , add=TRUE, col="red" )
          
            M = M[ -j, ]
            pg = pg[ which(! pg$AUID %in% au ), ] 
            pg = areal_units_neighbourhood_reset( pg, snap=2 )
          }
        }
      }

    
    M = M[ is.finite(M$cwd), ]
    M = M[ M$sex %in% c("0", "1"), ]
    M = M[ M$mat %in% c("0", "1"), ]

      # imperative covariate(s)
      M = M[ which(is.finite(M$z)), ]  
      M = M[ which(is.finite(M$t)), ]  
  
      M$space = match( M$AUID, pg$AUID) # for bym/car .. must be numeric index matching neighbourhood graphs
      M$space_time = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 
      M$space_cyclic = M$space  # copy for space_time component (INLA does not like to re-use the same variable in a model formula) 

      M$time = match( M$year, p$yrs ) # copy for space_time component .. for groups, must be numeric index
      M$time_space = M$time    
       
            
        M$pa = NA
        M$pa[which(M$N>0 & M$tag=="observations")] = 1
        M$pa[which(M$N==0 & M$tag=="observations")] = 0

        M$N[ !is.finite(M$N) ] = 1  # prediction surface
        M$N[ M$N==0 ] = 1  # where observations are zero, assume effort was at least 1
        M$sa[ !is.finite(M$sa) ] = 1  # prediction surface
        M$Ntrials = round( M$N / M$sa )  # weight = density

      read_write_fast( data=M, fn=fn)
      return(M) 
    }



}

 

 
