 

size_distributions = function(
    p=p, 
    outdir=file.path(p$project.outputdir, "size_structure"), 
    toget="base_data",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    xrange=c(0, 160),
    cwbr = 2, 
    bw=2,
    kernel="gaussian",
    kexp=1,
    density_offset = NULL,
    grad_method="simple",
    redo=FALSE,
    add_zeros=FALSE,
    pg=NULL,
    sigdigits=2,
    lowpassfilter=2,
    lowpassfilter2=0.001,
    n_neighbours=0,
    n_min=10,
    n_cutoff=0,
    plot_solutions = FALSE,
    ti_window = c(-3, 3),
    group=NULL,
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
                den_lb=exp(mean(log_den, na.rm=TRUE)-density_offset - 1.96*sd(log_den, na.rm=TRUE)),
                den_ub=exp(mean(log_den, na.rm=TRUE)-density_offset + 1.96*sd(log_den, na.rm=TRUE))
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
 
        outdir = file.path( p$project.outputdir, "size_structure", "posteriors_summaries" )
        dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 

        if (!redo) {
            fn = file.path( outdir, paste( "posterior_summaries_", Y, ".csv", sep="" ) )
            M = NULL
            if (file.exists(fn)) {
                M =read.csv(fn)
                setDT(M)
                M$X = NULL
                M$sex = as.character(M$sex)
                M$mat = as.character(M$mat)
                M$au = as.character(M$au)
            }

            return(M)
        }
        
        M = size_distributions( p=p, toget="base_data" )
        M$sex = as.character(M$sex)
        M$mat = as.character(M$mat)

        M$logcw = log(M$cw)

        M$wt = 1.0 / M$sa

        sigdigits=3
        
        # data_resolution is ~1 to 2 mm (observations error)
        bw =  round(median(log(22:142)-log(20:140)), digits=2 )   # approx SD for each interval

        # discretize time quarterly
        M$ti = M$year + round(trunc(M$julian / 365 * 4 ) / 4, digits=sigdigits) # weekly 

        # dimensionality of problem
        nb = attributes(pg)$nb$nbs 
        aus = pg$AUID
       
        # max(M$logcw) # 5.14166355650266
        # min(M$logcw) # 2.138889000323256
        xr = round( range(M$logcw, na.rm=TRUE), digits=2 ) 
        
        np = 256
        dx = diff(xr)/(np-1)
    
        seasons = seq(0, 0.75, by=0.25)  
      
        yrs = sort( unique( M$year ) ) # 1996:present
        if (!is.null(Y)) yrs=Y

        sexes = c("0", "1")  # 0 ,1, 2 male, female, unknown
        mats = c("0", "1")   # 0 ,1, 2  imm, mat, unknown

        if (0) { yr=2022; season=40; au="360"; sex="1"; mat="1" }

        for (yr in yrs) {
            
            fnout  = file.path( outdir, paste( "posterior_summaries_", yr, ".csv", sep="" ) )

            out1 = NULL
            out2 = NULL

            for (season in seasons) {
                mti = yr + season
                kt = M[ ti == mti, which=TRUE ]

                if ( length(kt) < 1)  next() 
                
                #for (au in aus) {
                sids = unique(M$sid[kt])
                for (si in sids) {
                    ka = intersect( kt, M[ sid == si, which=TRUE ] )
                    if (length(ka) < 1) next()
                 
                    au = unique(M$space_id[ka])[1]
                    for (sex in sexes) {
                        ks = intersect( ka, M[ sex==sex, which=TRUE ] )
                        if (length(ks) < 1) next()
                        for (mat in mats) {
                            n =  intersect( ks, M[ mat==mat, which=TRUE] )
                            N =  length(n) 
                            if (N < 1) next() 
                            tout = paste("|sid: ", si, "| sex: ", sex, "| mat: ", mat, "| au: ", au, "|year: ", yr, "| season: ", season, "| N: ", N ) 
                            message(tout )
                            uu = density( M$logcw[n], bw=bw, kernel="gaussian", from=xr[1], to=xr[2], n=np, weights=M$wt[n], na.rm=TRUE )
                            uu$y = uu$y / sum(uu$y) / dx  # density
                            out1 = rbind( out1, data.table( sid=si, sex=sex, mat=mat, au=au, year=yr, season=season, Nsample=N, Neffective=round( sum( M$wt[n]) ) )) 
                            out2 = rbind( out2, data.table( t(uu$y))  )
                            res = NULL
                        } # mat
                    } # sex
                }  # au
            }   # time  
     
            write.csv( file=fnout, cbind(out1, out2) )
            print(fnout ) 

        }
 
    }

    if (toget=="modal_analysis") {

        fn = file.path( outdir, "size_distributions_modes.RDS" )
        if (!redo) {
            Y = NULL
            if (file.exists(fn)) Y =readRDS(fn)
            if (is.null(group)) {
                return(Y)
            } else {
                return( Y[[group]] )
            }
        }

        if(0){
            kernel="gaussian"
            grad_method="Richardson"
            bw=0.01
            sigdigits=3
            ti_window = c(-6, 6)
            kexp=1
            n_min=10
            n_cutoff=3
            lowpassfilter=0.0001
            lowpassfilter2=0.0001
            n_neighbours=2
            plot_solutions=TRUE
            xrange=c(0, 160) 
            cwbr = 2
            kernel="gaussian"
            density_offset = NULL
            add_zeros=FALSE
            pg=NULL
            n_neighbours=0
        }
    
        # similar to kernel_density_weighted in approach but with aggregation across time and space windows
        
        auid = pg$AUID
     
        # max(M$logcw) # 5.14166355650266
        # min(M$logcw) # 2.138889000323256                          
   
        xr =  c(2.138889000323256 , 5.14166355650266)
        np = 256
        dx = diff(xr)/(np-1)
            
        for ( yr in yrs ) {
            Y = size_distributions(p=p, toget="kernel_density_weighted", Y=yr )
            # Y0 = Y[,1:8]
            for (season in c( 0, 0.25, 0.5, 0.75 )) {
                kt = Y[ season == season, which=TRUE ]
                if (length(kt)< nmin) next()

                for (sex in c("0", "1"))  {
                for (mat in c("0", "1"))  {
                    out = NULL
                    ksm = intersect( kt, which(  Y$mat == mat &  Y$sex == sex  ) )
                    if (length(ksm) < nmin) next()
                    aus = unique( Y$au[ksm] )

                    for (au in aus) {
                        aoi = au
                        if ( n_neighbours > 0 ) {
                            aui = which(auid %in% aoi )
                            aoi = identify_neigbours( nb=attributes(pg)$nb, index=aui, n_neighbours=n_neighbours  )
                        }
                        ka = which( Y$au %in% auid[aoi] )
                        k = intersect( ksm, ka ) 
                        if (length(k > n_min)) {

                        kd = Y[k, .(.SD * Neffective), .SDcols=patterns("^V[[:digit:]]+")  ] 
                        kd = colSums(kd, na.rm=TRUE) / sum(kd, na.rm=TRUE) / dx
                   
                        Ns = sum(Y[k,Nsample], na.rm=TRUE)
                        Ne = sum(Y[k,Neffective], na.rm=TRUE)
  
        }}}}}}

        out = list(male = mmat, female=fmat, imale=mimm, ifemale=fimm ) # immature = ximm 
        saveRDS( out, file=fn, compress=TRUE)

        return (out)
    }

 
}

 

 
