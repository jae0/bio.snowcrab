 

size_distributions = function(
    p=p, 
    outdir=file.path(p$project.outputdir, "size_structure"), 
    toget="base_data",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    xrange=c(0, 160),
    cwbr = 2, 
    density_offset = NULL   ,
    redo=FALSE,
    Y=NULL ) { 

    if (!dir.exists(outdir)) dir.create(outdir)
     
    # sex codes
    male = 0
    female = 1
    sex.unknown = 2

    # maturity codes
    immature = 0
    mature = 1
    mat.unknown = 2
 
    breaks = seq(xrange[1], xrange[2], by=cwbr)
    mids = breaks[-length(breaks)] + cwbr/2


    if (toget=="base_data") {

        fn = file.path( outdir, "size_distributions_base_data.RDS" )
        if (!redo) {
            Y = NULL
            if (file.exists(fn)) Y =readRDS(fn)
            return(Y)
        }

        set = snowcrab.db( DS="set.clean")
        set$sid = paste(set$trip, set$set, sep="~")

        det = snowcrab.db( DS="det.initial")
        det$sid = paste(det$trip, det$set, sep="~")

        set$year = set$yr 
        set$region = NA

        for ( region in regions ) {
            r = polygon_inside(x=set, region=aegis.polygons::polygon_internal_code(region), planar=F)
            if (length(r) > 0) set[r, "region"] = region
        }

        set= set[!is.na(set$region), ]
        set = set[, c("sid", "region", "year", "sa", "t", "z", "timestamp", "julian")]

        det = det[, c("sid", "shell", "cw", "sex", "mass", "mat", "gonad", "durometer")]
        
        setDT(set)
        setDT(det)

        Y = det[ set, on=.(sid)]
        
        # trim a few strange data points
        o = lm( log(mass) ~ log(cw), Y)
        todrop = which(abs(o$residuals) > 0.5)
        Y = Y[-todrop,]

        Y$cwd = as.numeric( as.character( cut( Y$cw, breaks=breaks, labels =mids ) ) )
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
        return(M)
    }

    if (toget == "simple_direct" ) {
 
        fn = file.path( outdir, "size_distributions_simple_direct.RDS" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) M =readRDS(fn)
            return(M)
        }
        
        if (is.null(Y)) Y = size_distributions(p=p, toget="tabulated_data"  )

        # do it here so we dot not have to store the massive intermediary file
        
        # CJ required to get zero counts dim(N) # 171624960    
        M = Y[ CJ( region, year, sex, mat, cwd, sid, unique=TRUE ), 
               on=.( region, year, sex, mat, cwd, sid ) ]

        M[ !is.finite(N),   "N"] = 0
        M[ !is.finite(mass), "mass"] = 0 
        M[ !is.finite(sa), "sa"] = 1 #dummy value
                
        M$density = M$N / M$sa
        M[ !is.finite(density), "density"] = 0  
     
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


    if (toget=="data_cubes") {
          
        require(data.cube)

        M = size_distributions(p=p, toget=Y, redo=redo)

        regs = unique(M$region)    # "cfanorth" "cfasouth" "cfa4x"
        yrs = sort( unique( M$year ) ) # 1996:present
        sexes = unique(M$sex)  # 0 ,1, 2 male, female, unknown
        mats = unique(M$mat)   # 0 ,1, 2  imm, mat, unknown
        # mids = 1.5, 2.5 ... 184.5
        shells = levels(M$shell) 
        
        data.cube.installed = try( require(data.cube), silent=TRUE )
        if (!data.cube.installed) {        
            # also stored on github
            install.packages("data.cube", repos = paste0("https://", c(
                "jangorecki.gitlab.io/data.cube",
                "cloud.r-project.org"
            )))
        }
 
        R = as.array( M[, .(region, year, sex, mat, cwd, den_log)],  measure="den_log",
            dimnames=list(region=regs, year=yrs, sex=sexes, mat=mats, cwd=mids ) ) 

        U = as.array( M[, .(region, year, sex, mat, cwd, den_sd_log)],  measure="den_sd_log",
            dimnames=list(region=regs, year=yrs, sex=sexes, mat=mats, cwd=mids ) ) 

        M = list(log_density=R, log_density_sd=U)      
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

        M = size_distributions(p=p, toget="tabulated_data" )

        gc()
        setDT(M)

        fit = biglm::bigglm( density ~ region:year:mat:cwd:sex - 1, data=M, family=gaussian(link="identity") )

        O = summary(fit)$coefficients

        orn = tstrsplit(rownames(O), ":")
        setDT(orn)
        setnames(orn, new=c("region", "year", "mat", "cwd", "sex"))

        orn[,region:=as.numeric(gsub( "region", "", region))]
        orn[,year:=as.factor(gsub( "year", "", year))]
        orn[,mat:=as.factor(gsub( "mat", "", mat))]
        orn[,cwd:=as.numeric(gsub( "cwd", "", cwd))]
        orn[,sex:=as.numeric(gsub( "sex", "", sex))]

        O = data.table(O)
        colnames(O) = c("density",  "density_se", "t", "p")
        O = cbind( orn, O )
        
        saveRDS(O, file=fn, compress=TRUE )

        return(O)
    }


    if (toget=="poisson_model") {

        stop("This takes far too long ... ")

        require(biglm)

        fn = file.path( outdir, "size_distributions_glm_poisson.RDS" )
        if (!redo) {
            O = NULL
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }

        M = size_distributions(p=p, toget="tabulated_data" )
        M$log_sa = log(M$sa)
        M$sa =NULL

        gc()
        setDF(M)
     
        fit = biglm::bigglm( N ~ region:year:mat:cwd:sex - 1 +offset(log_sa), data=M, 
            family=poisson(link="log") )
        O = coef(fit)

        saveRDS(O, file=fn, compress=TRUE )

        return(O)
    }



    if (toget=="poisson_inla") {
       
        require(INLA)

        fn = file.path( outdir, "size_distributions_inla_poisson.RDS" )
        if (!redo) {
            O = NULL
            if (file.exists(fn)) O =readRDS(fn)
            return(O)
        }

        M = size_distributions(p=p, toget="tabulated_data" )
        M$log_sa = log(M$sa)
        M$sa =NULL

        gc()
        setDF(M)
        M = M[which(M$year %in% as.character(c(2010:2015)) & M$region=="cfanorth" ),]

        fit = inla( N ~ year:mat:cwd:sex - 1 +offset(log_sa), data=M, 
            family="poisson" )

        saveRDS(fit, file=fn, compress=TRUE )

        return(fit)
    }



}

 

 
