size_distributions = function(
    p=p, 
    outdir=NULL,
    toget="",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    span = NULL,
    redo=FALSE,
    pg=NULL,
    Y=NULL ) { 

    if (is.null(outdir)) outdir = project.datadirectory( "bio.snowcrab", "output", "size_structure" ) 
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

    if (is.null(Y)) {
        if (exists("yrs", p)) Y = p$yrs
    }

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
            r = polygon_inside(x=set, region=region, planar=F)
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
            P$cwd[k] = discretize_data( x=P$cw[k], span=span(j) )  
        }

        P = P[ is.finite(cwd) ,]

        P = P[ year %in% Y, ]
        
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

            M$den_lb = M$den - 1.96*M$den_sd
            M$den_ub = M$den + 1.96*M$den_sd
            
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



}

 

 
