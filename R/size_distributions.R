 

size_distributions = function(
    p=p, 
    outdir=file.path(p$project.outputdir, "size_structure"), 
    toget="base_data",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    xrange=c(10, 160),
    np = 512,
    dx=NULL,
    cwbr = 2, 
    bw=2,
    kernel="gaussian",
    density_offset = NULL,
    redo=FALSE,
    add_zeros=FALSE,
    pg=NULL,
    sigdigits=2,
    lowpassfilter=0,
    lowpassfilter2=0.001,
    plot_solutions = FALSE,
    tlevels=c(-2, 6),
    zlevels=c(0, 100),
    Y=NULL ) { 

    if (!dir.exists(outdir)) dir.create(outdir)
     
    if (0) {
        regions=c("cfanorth", "cfasouth", "cfa4x")
        xrange=c(0, 160)
        cwbr = 2
    }
    
    
    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2
 


    if (toget=="base_data") {

        fn = file.path( outdir, "size_distributions_base_data.RDS" )
        if (!redo) {
            Y = NULL
            if (file.exists(fn)) Y =readRDS(fn)
            return(Y)
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
        

        Y = det[ set, on=.(sid)]
        
        # trim a few strange data points
        o = lm( log(mass) ~ log(cw), Y)
        todrop = which(abs(o$residuals) > 0.5)
        Y = Y[-todrop,]

        breaks = seq(xrange[1], xrange[2], by=cwbr)
        mids = breaks[-length(breaks)] + cwbr/2

        Y$cwd = discretize_data( Y$cw, brks=breaks, labels=mids, resolution=cwbr )  
        
        Y = Y[ is.finite(cwd) ,]
        
        Y$shell = factor( Y$shell )

        print(fn)
        saveRDS(Y, file=fn, compress=TRUE)
        return(Y)
    }


    if (toget=="tabulated_data") {

        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        
        fn = file.path( outdir, "size_distributions_tabulated_data.RDS" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) M =readRDS(fn)
            if (add_zeros) {
                # do it here so we dot not have to store the massive intermediary file
                
                # CJ required to get zero counts dim(N) # 171624960    
                M = M[ CJ( region, year, sex, mat, cwd, sid, unique=TRUE ), 
                    on=.( region, year, sex, mat, cwd, sid ) ]

                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                        
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
            }
            return(M)
        }
        
        if (is.null(Y)) Y = size_distributions(p=p, toget="base_data" )
    
        # aggregate by cwd 
        M = Y[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, sex, mat, cwd, sid) ]
        Y = NULL

        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$cwd = as.factor(M$cwd)
        M$region = as.factor(M$region)
        saveRDS(M, file=fn, compress=TRUE)
        
        # return this way to add zeros, if required
        return( size_distributions(p=p, toget="tabulated_data", add_zeros=add_zeros, redo=FALSE ) )
    }


    if (toget == "simple_direct" ) {
 
        fn = file.path( outdir, "size_distributions_simple_direct.RDS" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) M =readRDS(fn)
            return(M)
        }
        
        if (is.null(Y)) Y = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE  )

        message("Make this an incremental update ...")
     
        # NOTE without offset, this implicitly drops the zeros
        if (is.null(density_offset)) density_offset = min( M$density[ which(M$density>0) ] )
        message( "Density offset: ", density_offset )

        M$log_den = log(M$density + density_offset) 
   
        # aggregate across sid
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
            by= .( region, year, sex, mat, cwd)
        ]
          
        attr(M, "density_offset") = density_offset
        saveRDS(M, file=fn, compress=TRUE)
        return(M)     
    }



    if (toget=="linear_model") {
        require(biglm)

        fn = file.path( outdir, "size_distributions_lm.RDS" )
        if (!redo) {
            O = NULL
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }

        M = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE )

        gc()
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

        stop("This takes far too long to run ... ")
 
        fn = file.path( outdir, "size_distributions_poisson_glm.RDS" )
        if (!redo) {
            O = NULL
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }

        M = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE )
 
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
       

       stop("This takes far too long to run ... ")
  
        require(INLA)

        fn = file.path( outdir, "size_distributions_inla_poisson.RDS" )
        if (!redo) {
            O = NULL
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }

        M = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE )
        M$tag ="o"

        gc()
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
 
        outdir = file.path( p$project.outputdir, "size_structure", paste("kernel_densities", "_", round(bw,3), "_", np, sep="") )
        dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 

        xr = round( log(xrange), digits=2 ) 
        if(is.null(dx)) dx = diff(xr)/(np-1)
        xvals = seq( xr[1], xr[2], by=dx )

        if (!redo) {
            if (is.null(Y)) stop("year Y must be provided")
            M = NULL
            for (yr in Y) {
                fn = file.path( outdir, paste( "kernel_densities_", yr, ".csv", sep="" ) )
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


        set = snowcrab.db( DS="set.clean")
        set$sid = paste(set$trip, set$set, sep="~")
        setDT(set)
        set = set[ , .(sid, yr, t, z, lon, lat ) ]
        set$region = NA
        for ( region in regions ) {
            r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
            if (length(r) > 0) set$region[r] = region
        }
        set$zi = cut( set$z, breaks=c(zlevels, 5000 ), labels=zlevels)   
        set$ti = cut( set$t, breaks=c(tlevels, 100 ), labels=tlevels )   

        # next aggregate distributions across time /space
        cols = colnames(M)
        cols = cols[grep("^V[[:digit:]]+",cols)]
        xvals = attributes(M)$xvals
 
        M = set[M, on="sid"]
 
            attr(M, "xrange") = xrange
            attr(M, "xvals") = xvals
            attr(M, "xr") = xr
            attr(M, "bw") = bw
            attr(M, "dx") = dx

            return(M)            
        }

        sigdigits=3

        nb = attributes(pg)$nb$nbs 
        aus = pg$AUID
        
        M = size_distributions( p=p, toget="base_data" )
        M$sex = as.character(M$sex)
        M$mat = as.character(M$mat)
        M$logcw = log(M$cw)
        M$wt = 1.0 / M$sa

        M$ti = M$year + round(trunc(M$julian / 365 * 4 ) / 4, digits=sigdigits)         # discretize time quarterly
        seasons = seq(0, 0.75, by=0.25)  
        yrs = sort( unique( M$year ) ) # 1996:present
        if (!is.null(Y)) yrs=Y

        sexes = c("0", "1")  # 0 ,1, 2 male, female, unknown
        mats = c("0", "1")   # 0 ,1, 2  imm, mat, unknown

        if (0) { yr=2022; season=40; au="360"; sex="1"; mat="1" }

        for (yr in yrs) {
            fnout  = file.path( outdir, paste( "kernel_densities_", yr, ".csv", sep="" ) )
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
                            uu = density( M$logcw[n], bw=bw, kernel=kernel, from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE )
                            uu$y = uu$y / sum(uu$y) / dx  # density
                            out1 = rbind( out1, data.table( sid=si, sex=sx, mat=mt, au=au, year=yr, season=season, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                            out2 = rbind( out2, data.table( t(uu$y))  )
                            res = NULL
                        } # mat
                    } # sex
                }  # au
            }   # time  
     
            data.table::fwrite( cbind(out1, out2), file=fnout )
            print(fnout ) 

        }
 
    }


 
}

 

 
