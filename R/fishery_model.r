
fishery_model = function(  p=NULL, DS="plot",  
  plotresults=TRUE, tag="default", areas=c("cfanorth", "cfasouth", "cfa4x"),   
  vname="", type="density", res=NULL, fit=NULL, fn=NULL, aulabels=c("N-ENS","S-ENS","4X"), save.plot=TRUE, ... ) {

#     sb = bio.snowcrab::fishery_model( DS="data_aggregated_timeseries", p=p  )


  if (tag=="default") {
    if (!is.null(p)) if (exists("tag", p)) tag = p$tag
  }

 


  if (DS=="logistic_parameters") {
    
    p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

    out = list()

    if (!exists("method", out)) out$method = "stan"  # "jags", etc.

    if (!exists("carstm_model_label", p)) p$carstm_model_label = "default"   

    if (!exists("outdir", out)) out$outdir = file.path( p$modeldir, p$carstm_model_label, "fishery_model_results", tag )

    if (!exists("fnres", out)) out$fnres  = file.path( out$outdir, paste( "logistics_model_results", p$year.assessment, out$method, tag, "RDS", sep=".") )

    if (!exists("fnfit", out)) out$fnfit  = file.path( out$outdir, paste( "logistics_model_fit", p$year.assessment, out$method, tag, "RDS", sep=".") )
    
    dir.create( out$outdir, showWarnings = FALSE, recursive = TRUE  )

    message( "Results will be saved to:", out$outdir)

    # observations
    if (!exists("fmdata", out)) {
      out$fmdata = fishery_model( DS="data_aggregated_timeseries", p=p  )
      oo = apply( out$fmdata$IOA, 2, range, na.rm=TRUE )
      for (i in 1:ncol(oo)) {
        out$fmdata$IOA[,i] = (out$fmdata$IOA[,i] - oo[1,i] )/ diff(oo[,i])  # force median 0.5 with most data inside 0,1
      }
      # out$fmdata$IOA_min = apply( out$fmdata$IOA, 2, min, na.rm=TRUE ) 
    }

    if (!exists("er", out$fmdata)) out$fmdata$er = 0.2  # target exploitation rate
    if (!exists("U", out$fmdata))  out$fmdata$U = ncol( out$fmdata$IOA)  # number of regions
    if (!exists("N", out$fmdata))  out$fmdata$N = nrow( out$fmdata$IOA)  # no years with data
    if (!exists("M", out$fmdata))  out$fmdata$M = 3 # no years for projections
    if (!exists("ty", out$fmdata)) out$fmdata$ty = which(p$yrs == 2004)  # index of the transition year (2004) between spring and fall surveys
    if (!exists("cfa4x", out$fmdata))  out$fmdata$cfa4x = 3 # column index of cfa4x
    if (!exists("eps",   out$fmdata))  out$fmdata$eps = 1e-9  # small non-zero number

    out$fmdata$missing = ifelse( is.finite(out$fmdata$IOA), 0, 1)
    out$fmdata$missing_n = colSums(out$fmdata$missing)
    out$fmdata$missing_ntot = sum(out$fmdata$missing_n)

    # this must be done last
    out$fmdata$IOA[ which(!is.finite(out$fmdata$IOA)) ] = 0 # reset NAs to 0 as stan does not take NAs
    out$fmdata$CAT[ which(!is.finite(out$fmdata$CAT)) ] = out$fmdata$eps  # remove NA's

    # priors
    if (!exists("Kmu", out$fmdata)) out$fmdata$Kmu =  c( 5.5, 65.0, 2.0 )   ## based upon prior historical analyses (when stmv and kriging were attempted)
    if (!exists("rmu", out$fmdata)) out$fmdata$rmu =  c( 1.0, 1.0, 1.0 )    ## biological constraint 
    if (!exists("qmu", out$fmdata)) out$fmdata$qmu =  c( 1.0, 1.0, 1.0 )    ## based upon video observations q is close to 1 .. but sampling locations can of course cause bias (avoiding rocks and bedrock)

    if (!exists("Ksd", out$fmdata)) out$fmdata$Ksd =  c( 0.25, 0.25, 0.25 ) * out$fmdata$Kmu   
    if (!exists("rsd", out$fmdata)) out$fmdata$rsd =  c( 0.1, 0.1, 0.1 ) * out$fmdata$rmu  # smaller SD's to encourage solutions closer to prior means
    if (!exists("qsd", out$fmdata)) out$fmdata$qsd =  c( 0.1, 0.1, 0.1 ) * out$fmdata$qmu   
 
    return(out)
  }


  if (DS=="stan_surplus_production_2019_model") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       
      }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.0> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=eps, upper=1.0> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps, upper=1.0> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps, upper=1.25> [M+N,U] bm;  // permit overshoot bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
              Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
              if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }
      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ beta( 8, 2 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ beta( 1, 100000 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] / K[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] /K[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3]/K[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j]/K[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3]/ K[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] / K[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( fmax(fmin(b0,0.99),eps), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

 
        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }
        }

 
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }


    "
    )
  }



  if (DS=="stan_surplus_production_2022_model_variation1_wider_qc_uniform") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=-1, upper=1> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ uniform( -0.5, 0.5 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 2 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }



  if (DS=="stan_surplus_production_2022_model_qc_uniform") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=0, upper=1> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.25> [U] bosd;  // observation error
        vector <lower=eps, upper=0.25> [U] bpsd;  // process error
        // vector <lower=eps, upper=0.25> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coefficient
        qc ~ uniform( 0, 1 ) ; // i.e., Y:b offset constant  ;  

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }

  if (DS=="stan_surplus_production_2022_model_variation1_wider_qc_normal") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=-1, upper=1> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ normal( 0, 0.25 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }



  if (DS=="stan_surplus_production_2022_model_variation1_wider_qc_normal_0") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=-1, upper=1> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ normal( 0, 0.000001 ) ; // i.e.,force to be 0
        
        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }


  if (DS=="stan_surplus_production_2022_model_qc_beta") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=0, upper=0.5> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ beta( 1, 2 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }




  if (DS=="stan_surplus_production_2022_model_qc_cauchy") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=0, upper=1.0> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ cauchy( 0, 0.1 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }


  if (DS=="stan_surplus_production_2022_model_qc_cauchy_wider") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        int NU;
        int Ndata;

        MN = M+N ;
        N1 = N+1;
        NU = N*U;
        Ndata = NU -  missing_ntot; 
       }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=-1.0, upper=1.0> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] rem_sd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
      }

      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ cauchy( 0, 0.1 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - CAT[1,j]/K[j]  , eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - CAT[i-1,j]/K[j] , eps) )), bosd[j] ) ;
          }
        }

        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - CAT[i,3]/K[3]  , eps) )) , bosd[3] ) ;
        }


        //  transition year (ty)
        for (j in 1:2) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0 , eps) ) ) , bosd[j] ) ; //NENS and SENS
        }
        ( Y[ty,3] / q[3] ) +  qc[3]  ~ normal( (  ( fmax( bm[ty,3]  - CAT[ty,3]/K[3] , eps) ) ) , bosd[3] ) ; //SENS

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - CAT[i,j]/K[j]   , eps) ) ) , bosd[j] ) ; //   fall surveys
          }
        }

        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        vector[Ndata] log_lik;
        int n;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
        for (j in 1:U) {
          for(i in 1:N) {
            B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
            C[i,j] = CAT[i,j] ;
          }

          for (i in N1:MN) {
            B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
            C[i,j] = er*bm[(i-1),j] * K[j] ;
          }

        }
        
        // for loo calcs (pointwise loglikelihoods)
        n=0;
        for (j in 1:U) {
          for (i in 1:N) {
            if ( missing[i,j] == 0 ) {
              n += 1;
              log_lik[n] = normal_lpdf( Y[i,j] | B[i,j], bosd[j] );
            } 
          }
        }

      }
    "
    )
  }





  if (DS=="stan_surplus_production_catch_observation") {
    return( "
      data {

        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        vector[U] qsd ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }

      transformed data {
        int MN;
        int N1;
        MN = M+N ;
        N1 = N+1;
      }

      parameters {
        vector <lower=eps> [U] K;
        vector <lower=0.25, upper=2.0> [U] r;  // biologically should be ~ 1 
        vector <lower=eps, upper=2.5> [U] q;  // multiplicative factor, unlikely to be >200%
        vector <lower=eps, upper=0.5> [U] qc;  //  offset .. unlikely to be off by > 50%
        vector <lower=eps, upper=0.5> [U] bosd;  // observation error
        vector <lower=eps, upper=0.5> [U] bpsd;  // process error
        vector <lower=eps, upper=0.5> [U] catsd;  // catch error
        vector <lower=eps> [U] b0;
        vector <lower=eps> [missing_ntot] IOAmissing;
        matrix <lower=eps> [M+N,U] bm;  // force bm max 1 
        matrix <lower=eps> [N,U] cat;  // estimated catch
        vector <lower=eps, upper=0.5> [U] catQ;  // multiplicative factor .. fraction misreported
      }

      transformed parameters {
        
        matrix[N,U] Y;  // index of abundance
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
                Y[i,j] = IOA[i,j]   ;  // translation of a zscore to a positive internal scale 
                if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }

      }

      model {

        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        b0 ~ normal( 0.5, 0.1 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;
        catsd ~ cauchy( 0, 0.1 ) ;
        catQ ~  beta( 1, 5 ) ;

        q ~ normal( qmu, qsd ) ; // i.e., Y:b scaling coeeficient
        qc ~ beta( 1, 5 ) ; // i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )

        // -------------------
        // biomass observation model
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018


        // catch observation model -- cat = CAT/K *catQ
        // spring surveys
        for (j in 1:2) {
          ( CAT[1,j]/K[j] ) * (1.0 + catQ[j]) ~ normal( cat[1,j], catsd[j] );  
          for (i in 2:(ty-1) ){
            ( CAT[i-1,j]/K[j] ) * (1.0 + catQ[j]) ~ normal( cat[i-1,j], catsd[j] );  
          }
        }
        for (i in 1:(ty-1) ){
          ( CAT[i,3]/K[3] ) * (1.0 + catQ[3]) ~ normal( cat[i,3], catsd[3] );  
        }
        //  transition year (ty)
        for (j in 1:2) {
          ( CAT[ty-1,j] + CAT[ty,j] ) / ( 2.0* K[j] )  * (1.0 + catQ[j]) ~ normal( cat[ty,j], catsd[j] );  
        }
        ( CAT[ty,3] / K[3] )  * (1.0 + catQ[3]) ~ normal( cat[ty,3], catsd[3] );  

        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( CAT[i,j]/K[j] ) * (1.0 + catQ[j]) ~ normal( cat[i,j], catsd[j] );  
          }
        }

        // biomass observation model 
        // spring surveys
        for (j in 1:2) {
          ( Y[1, j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[1,j] - cat[1,j], eps) )) , bosd[j] ) ;
          for (i in 2:(ty-1) ){
            ( Y[i, j] / q[j] ) +  qc[j]   ~ normal( ( (  fmax( bm[i,j] - cat[i-1,j], eps) )), bosd[j] ) ;
          }
        }
        for (i in 1:(ty-1) ){
          ( Y[i, 3] / q[3] ) +  qc[3]  ~ normal( ( ( fmax( bm[i,3] - cat[i,3], eps) )) , bosd[3] ) ;
        }

        //  transition year (ty)
        for (j in 1:3) {
          ( Y[ty,j] / q[j] ) +  qc[j]  ~ normal( ( ( fmax( bm[ty,j]  - cat[ty,j], eps) ) ) , bosd[j] ) ; 
        }
        // fall surveys
        for (j in 1:3) {
          for (i in (ty+1):N) {
            ( Y[i,j] / q[j] ) +  qc[j] ~ normal( ( ( fmax( bm[i,j] - cat[i,j], eps) ) ), bosd[j] ) ; //   fall surveys
          }
        }


        // -------------------
        // biomass process model
        // fmax .. force positive value .. initial conditions
        bm[1,] ~ normal( (b0), bpsd ) ;

        for (j in 1:U) {
          real o;
          for (i in 2:N) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - cat[i-1,j] )), bpsd[j] ) ;
          }
          for (i in N1:MN) {
            o = r[j] * fmax( 1.0 - bm[i-1,j], eps ) ; 
            bm[i,j] ~ normal( ( ( bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j]  )), bpsd[j] ) ;
          }
        }
      }

      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;

        // -------------------
        // fishing mortality
        // fall fisheries

         for (j in 1:U) {
           for (i in 1:N) {
             F[i,j] =  -log( fmax( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
           }
           for (i in N1:MN) {
             F[i,j] =  -log( fmax( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
           }
         }

 
        // -------------------
        // parameter estimates for output
        for (j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }

        // recaled estimates
         for (j in 1:U) {
           for(i in 1:N) {
             B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
             C[i,j] = CAT[i,j] ;
           }

           for (i in N1:MN) {
             B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
             C[i,j] = er*bm[(i-1),j] * K[j] ;
           }

         }

      }
    "
    )
  }



  if (DS=="stan_surplus_production_2020") {
    return( "
      data {
        int<lower=0> N; // no. years
        int<lower=0> U; // no. regions
        int<lower=0> M; // no. years to project
        int ty;
        real er ;
        real eps ;
        vector[U] Ksd;
        vector[U] rsd;
        vector[U] qsd;
        vector[U] Kmu ;
        vector[U] rmu ;
        vector[U] qmu ;
        matrix[N,U] CAT;
        matrix[N,U] IOA;
        matrix[N,U] missing;
        array[U] int missing_n;
        int missing_ntot;
      }
      transformed data {
        int MN;
        int N1;
        MN = M+N ;
        N1 = N+1;
      }
      parameters {
        vector <lower=eps>[U] K;
        vector <lower=eps>[U] r;
        vector <lower=eps>[U] q;
        vector <lower=eps>[U] qs;
        vector <lower=eps,upper=(1-eps)>[U] bosd;  // observation error
        vector <lower=eps,upper=(1-eps)>[U] bpsd;  // process error
        vector <lower=eps,upper=(1-eps)>[U] b0;
        vector <lower=eps>[missing_ntot] IOAmissing;
        matrix <lower=eps>[M+N,U] bm;
      }
      transformed parameters {
        matrix[N,U] Y;  // index of abundance
        matrix[N,U] Ymu;  // collator used to force positive values for lognormal
        matrix[MN,U] bmmu; // collator used to force positive values for lognormal
        matrix[MN,U] rem;  // observed catch
        // copy parameters to a new variable (Y) with imputed missing values
        {
          int ii;
          ii = 0;
          for (j in 1:U) {
            for (i in 1:N) {
              Y[i,j] = IOA[i,j];
              if ( missing[i,j] == 1 ) {
                ii = ii+1;
                Y[i,j] = IOAmissing[ii];
              }
            }
          }
        }
        // -------------------
        // removals (catch) observation model, standardized to K (assuming no errors in observation of catch!)
        for (j in 1:U) {
          rem[1:N,j] =  CAT[1:N,j]/K[j] ;
          rem[(N+1):MN,j] =  er*bm[ N:(MN-1),j] ;  // forecasts
        }
        // -------------------
        // observation model calcs and contraints:
        // Ymu = 'surveyed/observed' residual biomass at time of survey (Bsurveyed)
        // cfanorth(1) and cfasouth(2)
        //   This is slightly complicated because a fall / spring survey correction is required:
        //   B represents the total fishable biomass available in fishing year y
        //     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
        //     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
        // spring surveys from 1998 to 2003
        //   this is conceptualized in the following time line:
        //     '|' == start/end of each new fishing year
        //     Sf = Survey in fall
        //     Ss = Survey in spring
        //     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
        // Cfa 4X -- fall/winter fishery
        //    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018
        for (j in 1:2) {
          Ymu[1,j]        = qs[j] * bm[1,j] - rem[1,j] ; // starting year approximation
          Ymu[2:(ty-1),j] = qs[j] * bm[2:(ty-1),j] - rem[1:(ty-2),j] ; //spring surveys
          Ymu[ty,j]       = q[j]  * bm[ty,j] - (rem[(ty-1),j] + rem[ty,j] )/2.0  ; // transition year .. approximation
          Ymu[(ty+1):N,j] = q[j]  * bm[(ty+1):N,j] - rem[(ty+1):N,j] ;   // fall surveys
        }
        {
          int k;
          k=3;
          Ymu[1,k]        = q[k] * bm[1,k]   - rem[1,k] ; // starting year approximation ymu[1991] = bm[1991]-rem[1991]
          Ymu[2:(ty-1),k] = q[k] * bm[2:(ty-1),k] - rem[2:(ty-1),k];
          Ymu[ty:N,k]     = q[k] * bm[ty:N,k] - rem[ty:N,k];
        }
        for (j in 1:U) {
          for (i in 1:N) {
            Ymu[i,j] = K[j] * fmax( Ymu[i,j], eps); // force positive value
          }
        }
        // -------------------
        // process model calcs and constraints
        for (j in 1:U) {
          bmmu[1,j] = b0[j] ; // biomass at first year
          for (i in 2:MN) {
            bmmu[i,j] = bm[i-1,j] * ( 1.0 + r[j]*(1-bm[i-1,j]) ) - rem[i-1,j] ;
          }
        }
        for (j in 1:U) {
          for (i in 1:MN) {
            bmmu[i,j] = fmax(bmmu[i,j], eps);  // force positive value
          }
        }
      }
      model {
        // -------------------
        // priors for parameters
        K ~ normal( Kmu, Ksd )  ;
        r ~ normal( rmu, rsd )  ;
        q ~ normal( qmu, qsd )  ;
        qs ~ normal( qmu, qsd )  ;
        b0 ~ beta( 8, 2 ) ; // starting b prior to first catch event
        bosd ~ cauchy( 0, 0.1 ) ;  // slightly informative .. center of mass between (0,1)
        bpsd ~ cauchy( 0, 0.1 ) ;
        // -------------------
        // biomass observation model
        for (j in 1:U) {
          log(Y[1:N,j]) ~ normal( log(Ymu[1:N,j]), bosd[j] ) ;  // stan thinks Y is being transformed due to attempt to impute missing values .. ignore
        }
        // -------------------
        // biomass process model
        for (j in 1:U) {
          log(bm[1:MN,j]) ~ normal( log(bmmu[1:MN,j]), bpsd[j] ) ;
        }
        // could have used lognormal but this parameterization is 10X faster and more stable
        target += - log(fabs(Y));  // required due to log transf above
        target += - log(fabs(bm));
      }
      generated quantities {
        vector[U] MSY;
        vector[U] BMSY;
        vector[U] FMSY;
        matrix[MN,U] B;
        matrix[MN,U] C;
        matrix[MN,U] F;
        matrix[M,U] TAC;
        // -------------------
        // fishing mortality
        // fall fisheries
         for (j in 1:3) {
           for (i in 1:N) {
             F[i,j] =  1.0 - rem[i,j] / bm[i,j]  ;
           }
         }
         for (j in 1:U) {
           for (i in N1:MN) {
             F[i,j] =  1.0 - er * bm[i-1,j] / bm[i,j]  ;
           }
           for (i in 1:MN) {
             F[i,j] =  -log( fmax( F[i,j], eps) )  ;
           }
         }
        // -------------------
        // parameter estimates for output
        for(j in 1:U) {
           MSY[j]    = r[j]* exp(K[j]) / 4 ; // maximum height of of the latent productivity (yield)
           BMSY[j]   = exp(K[j])/2 ; // biomass at MSY
           FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; // fishing mortality at MSY
        }
        // recaled estimates
         for(j in 1:U) {
           for(i in 1:MN) {
             B[i,j] = (bm[i,j] - rem[i,j]) * K[j] ;
             C[i,j] = rem[i,j]*K[j] ;
           }
           for(i in 1:M) {
             TAC[i,j] = rem[N+i,j]*K[j] ;
           }
         }
      }
    "
    )
  }




  if (DS=="data_aggregated_timeseries" ) {

    cfanorth =  1 # column index
    cfasouth =  2 # column index
    cfa4x =  3 # column index

    landings = bio.snowcrab::logbook.db(DS="logbook") 
      # NOTE:: message( "Fishing 'yr' for CFA 4X has been set to starting year:: 2001-2002 -> 2001, etc.")
      # year is year of capture
      # yr is "fishing year" relative to the assessment cycle
    landings = landings[ which (landings$cfa %in% c( "cfanorth", "cfasouth", "cfa4x" ) ) , ]
    L = tapply( landings$landings, INDEX=landings[,c("yr", "cfa")], FUN=sum, na.rm=T )
    nL = nrow(L)

    cfaall = tapply( landings$landings, INDEX=landings[,c("yr")], FUN=sum, na.rm=T )
    L = cbind( L, cfaall )
    L = L / 1000/1000  # convert to kt pN$fishery_model_label = "stan_surplus_production_2022_model_qc_cauchy"
    L[ !is.finite(L)] = 0

    L = as.data.frame( L[ match( p$yrs, rownames(L) ), areas ] )

    # biomass data: post-fishery biomass are determined by survey B)
    B = aggregate_simulations( fn=carstm_filenames( p, returnvalue="filename", fn="aggregated_timeseries" ) )$RES

    rownames(B) = B$yrs
    B = as.data.frame( B[ match( p$yrs, B$yrs ), areas ] )

    # cfa4x have had no estimates prior to 2004

    cfanorth.baddata = which( p$yrs <= 2004 )
#    B[ cfanorth.baddata, cfanorth ] = NA

    cfasouth.baddata = which( p$yrs <= 2004 )
#    B[ cfasouth.baddata, cfasouth ] = NA

    cfa.nodata =   which( p$yrs <= 2004 )
    B[ cfa.nodata , cfa4x ] = NA

    sb = list(
      IOA = as.matrix(B), # observed index of abundance
      CAT = as.matrix(L)  # catches  , assume 20% handling mortality and illegal landings
    )

    return(sb)
  }


 
  if (DS=="logistic_model" ) {

    message( "Output location is: ", p$fishery_model$outdir )
    dir.create( p$fishery_model$outdir, recursive=T, showWarnings=F )


    if (is.null(fit)) {
      p$fishery_model$stancode$compile()
      fit = p$fishery_model$stancode$sample( ... )
    }

    if (0) {
        # Posterior summary statistics .. the $sample() method creates R6 CmdStanMCMC objects,
        fit$summary()  #  summarise_draws() from posterior package
        fit$summary("K", "r", "q")

        # this is a draws_array object from the posterior package
        # str(fit$sampler_diagnostics())

        diagnostics_df = as_draws_df(fit$sampler_diagnostics())
        print(diagnostics_df)

        fit$cmdstan_diagnose()
        fit$cmdstan_summary()

        # Posterior draws .. $draws()  a 3-D array (iteration x chain x variable)
        draws_array = fit$draws()
        str(draws_array)

        draws_df = as_draws_df(draws_array) # as_draws_matrix() for matrix
        print(draws_df)

    }

    fit$save_object( file = p$fishery_model$fnfit )   #  save this way due to R-lazy loading; RDS file

    res = list( mcmc=stan_extract( as_draws_df(fit$draws() ) ), p=p )

    read_write_fast( data=res, file=p$fishery_model$fnres)

    return(res)
  }



  if (DS=="samples" ) {
    res = NULL
    if (file.exists(p$fishery_model$fnres)) res = aegis::read_write_fast(p$fishery_model$fnres)
    return(res)
  }


  if (DS=="fit" ) {
    fit = NULL
    if (file.exists(p$fishery_model$fnfit)) fit = aegis::read_write_fast(p$fishery_model$fnfit)
    return(fit)
  }


  if (DS=="plot") {

    y = res$mcmc
    sb= res$p$fishery_model$fmdata
    
    if (is.null(fn)) {
      outdir = pN$fishery_model$outdir 
    } else{
      outdir = dirname(fn)
    }

    ntacs = sb$nProj
    yrs0 = res$p$yrs
    yrs = c( yrs0, (max(yrs0)+c(1:sb$M) ) )
    yrs.last = max(yrs0) + 0.5
    ndata = length(yrs0)
    hdat = 1:ndata


    if (vname =="r.ts") {
      # catch this first as the layout is different
        if (is.null(fn)) fn=file.path(outdir, "r.ts.density.png" )

        br = 75

        plot.new()
        layout( matrix(c(1:(sb$N*3)), ncol=3, nrow=sb$N ))
        par(mar = c(1., 1., 0.65, 0.75))

        u = apply( y[["r"]], c(2, 3), median, na.rm=T  )

        for (i in 1:3) {
        for (yr in 1:sb$N) {
          qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
          qs = signif( qs, 3 )
          pdat = u[ yr,i ]
          xrange = range( pdat, na.rm=T )
          postdat = hist( pdat, breaks=br, plot=FALSE )
          yrange = range( 0, postdat$density, na.rm=T ) * 1.02
          hist( pdat, freq=FALSE, breaks=br, xlim=xrange, ylim=yrange, main="", xlab="", ylab="Density", col="lightgray", border="gray")
          YR = rownames(sb$IOA) [yr]
          legend( "topright", bty="n", legend=paste( aulabels[i], YR, "r", " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ), cex=0.9 )
        }}

      if (save.plot) savePlot( filename=fn, type="png" )
      return( fn)

    }


    plot.new()
    layout( matrix(c(1,2,3), 3, 1 ))
    par(mar = c(4.4, 4.4, 0.65, 0.75))

    prr = NULL
    prr$class ="none"  # no priors by default


    if ( type=="density" ) {  # default

      if ( vname=="K" ) {
        
        if (is.null(fn)) fn=file.path(outdir, "K.density.png" )
        u = y[[vname]]

        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
           pdat = u[,i]
           prr=NULL
           prr$class="normal"

          # E(X) = exp(mu + 1/2 sigma^2)
          # med(X) = exp(mu)
          # Var(X) = exp(2*mu + sigma^2)*(exp(sigma^2) - 1)
          # CV = sqrt( Var(X) ) / E(X) = sqrt(exp(sigma^2) - 1) ~ sigma; sigma < 1/2
          # SD(X) = sqrt( exp(2*mu + sigma^2)*(exp(sigma^2) - 1) )
          #  or   = CV * E(X)

           prr$mean= sb$Kmu[i]
           prr$sd =  sb$Ksd[i]
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ... )

          legend( "topright", bty="n", legend=paste( aulabels[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ))
      }
    }

      if ( vname=="r" ) {
        if (is.null(fn)) fn=file.path(outdir, "r.density.png" )

        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {    prr=NULL
          prr=NULL
          prr$class='normal'
          prr$mean=sb$rmu[i]
          prr$sd= sb$rsd[i]
          pdat = u[,i]
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n", legend=paste( aulabels[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" ) )
      }}



      if ( vname=="q" ) {
        if (is.null(fn)) fn=file.path(outdir, "q.density.png" )

        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class="normal"
          prr$mean=sb$qmu[i]
          prr$sd=sb$qsd[i]
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n", legend=paste( aulabels[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" )
      )}}


      if ( vname=="qc" ) {
        if (is.null(fn)) fn=file.path(outdir, "qc.density.png" )
        u = y[[vname]]
        qc = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class="cauchy"
          prr$a = 0
          prr$b = 0.1
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n", legend=paste( aulabels[i], "\n", vname, " = ", qc[2,i], " {", qc[1,i], ", ",  qc[3,i], "}  ", sep="" )
      )}}


      if ( vname=="qs" ) {
        if (is.null(fn)) fn=file.path(outdir, "qs.density.png" )
        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        for (i in 1:3) {
          pdat = u[,i]
          # prr=NULL
          # prr$class="normal"
          # prr$mean=sb$q0x[i]
          # prr$sd=sb$q0x[i]*sb$cv
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n", legend=paste( aulabels[i], "\n", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}  ", sep="" )
      )}}


      if ( vname=="BMSY" ) {
        if (is.null(fn)) fn=file.path(outdir, "BMSY.density.png" )
        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class="none"
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n",
            legend=paste( aulabels[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )
      )}}

      if ( vname=="FMSY" ) {
        if (is.null(fn)) fn=file.path(outdir, "FMSY.density.png" )
        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class="none"
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n",
            legend=paste( aulabels[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )
      )}}


      if ( vname=="bosd" ) {
        if (is.null(fn)) fn=file.path(outdir, "bosd.density.png" )
        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class='uniform'
          prr$max=3
          prr$min=0
          #prr$class="lognormal"
          #prr$meanlog=sb$bomup
          #prr$sdlog=sqrt(sb$bosdp)
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n",
            legend=paste( aulabels[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )
      )}}

      if ( vname=="bpsd" ) {
        if (is.null(fn)) fn=file.path(outdir, "bpsd.density.png" )
        u = y[[vname]]
        qs = apply( u, 2, quantile, probs=c(0.025, 0.5, 0.975) )
        qs = signif( qs, 3 )
        for (i in 1:3) {
          pdat = u[,i]
          prr=NULL
          prr$class='uniform'
          prr$max=3
          prr$min=0
          #prr$class="lognormal"
          #prr$meanlog=sb$bpmup
          #prr$sdlog=sqrt(sb$bpsdp)
          plot.freq.distribution.prior.posterior( prior=prr, posterior=pdat, ...  )
          legend( "topright", bty="n",
            legend=paste( aulabels[i], " ", vname, " = ", qs[2,i], " {", qs[1,i], ", ",  qs[3,i], "}", sep="" )
      )}}


    }

    # --------------

    if ( type=="timeseries" ) {

    plot.new()
    layout( matrix(c(1,2,3), 3, 1 ))
    par(mar = c(4.4, 4.4, 0.65, 0.75))

      if (vname=="biomass") {
        
        if (is.null(fn))  fn=file.path(outdir, "biomass.timeseries.png" )

        u = y[["B"]]

        for (i in 1:3) {
          meanval = apply( u[,1:sb$N,i], 2, mean, na.rm=T  )

          prs = seq( from=0.025, to=0.975, length.out=600)
          Bq =  apply( u[,1:sb$N,i], 2, quantile, probs=prs, na.rm=T )

          #yran = range(c(0, Bq, sb$IOA[,i] ), na.rm=T )*1.01
          yran = range(c(0, Bq ), na.rm=T )*1.01
          plot( yrs0, Bq[1,], type="n", ylim=yran, xlim=range(yrs0), xlab="", ylab=""  ) #change xlim to yrs0 to remove 3 yr projection
          cols = gray.colors( floor(length( prs)/2) )
          cols2 = c(cols[length(cols):1], cols )
          for ( j in 1:length(prs) ) {
            lines ( yrs0, Bq[j,], lwd=4, col=cols2[j] )
          }
          # lines( yrs0, B, lwd=3, col="darkgreen" )
          #abline (v=yrs.last , lwd=2, lty="dashed" ) #can comment out this line if not providing forward projection
          if (i==2) title( ylab="Fishable biomass (kt)" )
          if (i==3) title( xlab="Year" )
          #points( yrs0, qIOA, pch=20, col="darkgreen" )
          #lines ( yrs0, qIOA, lwd=3, col="darkgreen", lty="dashed" )
          lines ( yrs0, meanval, lwd=2, col="blue", lty="dotted" )
          #points( yrs0, IOA, pch=20, col="darkred" )
          #lines( yrs0, IOA, lwd=3, lty="dotdash", col="red" )
          # legend( "topright", bty="n", legend=aulabels[i])
        }
      }

      if (vname=="rem") {

        if (is.null(fn))  fn=file.path(outdir, "rem.timeseries.png" ) 
        
        u = y[[vname]]

        for (i in 1:3) {
          meanval = apply( u[,,i], 2, mean, na.rm=T  )

          prs = seq( from=0.025, to=0.975, length.out=600)
          Bq =  apply( u[,,i], 2, quantile, probs=prs, na.rm=T  )

          #yran = range(c(0, Bq, sb$IOA[,i] ), na.rm=T )*1.01
          yran = range(c(0, Bq ), na.rm=T )*1.01
          plot( yrs, Bq[1,], type="n", ylim=yran, xlim=range(yrs0), xlab="", ylab=""  ) #change xlim to yrs0 to remove 3 yr projection
          cols = gray.colors( floor(length( prs)/2) )
          cols2 = c(cols[length(cols):1], cols )
          for ( j in 1:length(prs) ) {
            lines ( yrs, Bq[j,], lwd=4, col=cols2[j] )
          }
          # lines( yrs, rem, lwd=3, col="darkgreen" )
          #abline (v=yrs.last , lwd=2, lty="dashed" ) #can comment out this line if not providing forward projection
          if (i==2) title( ylab="Fishable biomass (kt)" )
          if (i==3) title( xlab="Year" )
          #points( yrs0, qIOA, pch=20, col="darkgreen" )
          #lines ( yrs0, qIOA, lwd=3, col="darkgreen", lty="dashed" )
          lines ( yrs, meanval, lwd=2, col="blue", lty="dotted" )
          #points( yrs0, IOA, pch=20, col="darkred" )
          #lines( yrs0, IOA, lwd=3, lty="dotdash", col="red" )
          legend( "topright", bty="n", legend=aulabels[i])
        }
      }

      if (vname=="fishingmortality") {

        if (is.null(fn))  fn=file.path(outdir, "fishingmortality.timeseries.png" ) 

        Fmsy = apply( y[["FMSY"]], 2, mean, na.rm=T )
      
        prs = seq( from=0.025, to=0.975, length.out=600)

        F = y[["F"]]
   
        Fi = apply( F[, 1:sb$N, ] , c(2,3), quantile, probs=prs, na.rm=T )
        
        for (i in 1:3) {
          yran = range(c(0, max(c(Fi,Fmsy))), na.rm=T )*1.05
          yran = pmin( yran, 1.2 )
          plot( yrs0, Fi[1,,i], type="n", ylim=yran, xlab="", ylab="" )
          cols = gray.colors( floor(length( prs)/2) )
          cols2 = c(cols[length(cols):1], cols )
          if (i %in% c(1,2)){
          for ( j in 1:length(prs) ) {
            lines ( yrs0, Fi[j,,i], lwd=4, col=cols2[j] )
          }}
          if (i==3){
            for ( j in 1:length(prs) ) {
            lines ( yrs0-0.2, Fi[j,,i], lwd=4, col=cols2[j] )
          }}
          if (i==2) title( ylab="Fishing mortality" )
          if (i==3) title( xlab="Year" )
          # legend( "topright", bty="n", legend=aulabels[i])
          abline (h=-log(1-0.2), lwd=2, lty="dashed" )
          abline (h=Fmsy[i], lwd=2, lty="solid", col="red" )
        }
      }

    }

    if (type=="hcr") {
      if (vname=="default") {

        if (is.null(fn))  fn=file.path(outdir, "hcr.default.png" )
        
        B =  apply(  y[["B"]], c(2,3), median, na.rm=T  )
        F =  apply( y[["F"]], c(2,3), median, na.rm=T  )
        K =  apply( y[["K"]], c(2), median, na.rm=T  )
        FMSY = apply( y[["FMSY"]], c(2), median, na.rm=T  )
        BMSY = apply( y[["BMSY"]], c(2), median, na.rm=T  )

        for (i in 1:3 ) {
          ylims = c(0, min( 1, max( FMSY[i] * 1.25, F[hdat,i] ) ) )
          plot( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),
            ylim=ylims, col="darkorange", cex=0.8, lwd=2, xlab="", ylab="", pch=20 )


          # nn = as.matrix( cbind( Bx=as.vector( y$B[,ndata,i] ), Fx = as.vector( y$F[,ndata,i] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)

          if (i==3) title( xlab="Fishable biomass (kt)" )
          if (i==2) title( ylab="Fishing mortality" )

          F30 = -log(1-0.3)
          F10 = -log(1-0.1)

          Fref =  0.22
          Bmsy = K[i] * 0.5
          Bref = K[i] * 0.2
          BK = K[i]
          BK25 = K[i] * .25
          Fhistorical = mean( F[hdat,i], na.rm=T )
          Bhistorical = mean( B[hdat,i], na.rm=T )
          yl = 0.05

          polygon(x=c(Bmsy,Bmsy*2,Bmsy*2, Bmsy),y=c(-0.08,-0.1,FMSY[i],FMSY[i]),col='lightgreen',border=NA)
          polygon(x=c(Bmsy/2,Bmsy,Bmsy, Bmsy/2),y=c(-0.08,-0.1,FMSY[i],FMSY[i]),col='lightgoldenrod',border=NA)
          polygon(x=c(0,Bmsy/2,Bmsy/2, 0),y=c(-0.08,-0.1,FMSY[i],FMSY[i]),col='darksalmon',border=NA)

#might need adjustment below to offset F vs B. Need to plot F against PREfishery biomass
          lines( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),
            ylim=ylims, col='blue', cex=0.8, lwd=2, xlab="", ylab="", pch=20 )

          abline (h=Fref, lty="solid", col="gray", lwd=2 )

          abline (h=F10, lty="dotted", col="gray")
          # text( 0.05*K[i], F10, "10% HR", pos=1 )

          abline (h=F30, lty="dotted", col="gray")
          # text( 0.05*K[i], F30, "30% HR", pos=1 )


          abline (h=FMSY[i], lty="dashed", col="red" )

          # abline (h=Fhistorical, lty="dashed")
          # text( 0.05*K[i], Fhistorical, "Mean", pos=1, lwd=2 )

          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)

          abline (v=Bmsy, lty="dotted")

          abline (v=BK, lty="dotted")

          abline (v=BK25, lty="dotted")

          text( Bmsy-0.01*K[i], yl, "K/2" , srt=90, pos=3)
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)
          text( 0.05*K[i], Fref, "20% HR", pos=1 )
          text( 0.05*K[i], FMSY[i], "FMSY", pos=3, lwd=2, col="red" )
          if (i %in% c(1,2)){
            text( B[1:(ndata-1),i], F[1:(ndata-1),i],  labels=yrs0[-ndata], pos=3, cex= 0.8 )
            points( B[ndata,i], F[ndata,i],  pch=21, bg='darkorange', cex= 1.4 )
            text( B[ndata,i], F[ndata,i],  labels=yrs0[ndata], pos=3, cex= 1.4, font=2 )

            text( 0, ylims[2]*0.9,  labels=aulabels[i], pos=3, cex= 0.85 )
          }
          if (i==3){
          text( B[1:(ndata-1),i], F[1:(ndata-1),i],  labels=(yrs0[-ndata]), pos=3, cex= 0.8 )
          points( B[ndata,i], F[ndata,i],  pch=21, bg='darkorange', cex= 1.4 )
          text( B[ndata,i], F[ndata,i],  labels=yrs0[ndata], pos=3, cex= 1.4, font=2 )
          text( 0, ylims[2]*0.9,  labels=aulabels[i], pos=3, cex= 0.85 )
          }
          # abline (v=Bhistorical, lty="dashed")
          # text( Bhistorical-0.01*K[i], yl, "Mean" , srt=90, pos=3,  lwd=2)
        }
      #Enable below to see the annual F estimates for inclusion in the document
          print("F for N-ENS" )
          print(F[hdat,1] )
          print("F for S-ENS" )
          print(F[hdat,2] )
        print("F for 4X" )
        print(F[hdat,3] )



      }

      if (vname=="default.unmodelled") {

        if (is.null(fn))  fn=file.path(outdir, "hcr.default.unmodelled.png" )


        B =  sb$IOA
        F =  apply( y[["F"]], c(2,3), median, na.rm=T  )

        areas = c("cfa4x", "cfasouth", "cfanorth" )
        regions = c("4X", "S-ENS", "N-ENS")

        td = exploitationrates(p=p, areas=areas, labels=regions, CFA4X.exclude.year.assessment=FALSE )


        K =  apply( y[["K"]], c(2), median, na.rm=T  )
        FMSY = apply( y[["FMSY"]], c(2), median, na.rm=T  )
        BMSY = apply( y[["BMSY"]], c(2), median, na.rm=T  )

        for (i in 1:3 ) {
          ylims = c(0, FMSY[i] * 1.25)
          plot( B[hdat,i], F[hdat,i],  type="b", xlim=c(0, K[i] * 1.1 ),
            ylim=ylims, col="darkorange", cex=0.8, lwd=2, xlab="", ylab="", pch=20 )

          # nn = as.matrix( cbind( Bx=as.vector( y$B[,ndata,i] ), Fx = as.vector( y$F[,ndata,i] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)

          if (i==3) title( xlab="Fishable biomass (kt)" )
          if (i==2) title( ylab="Fishing mortality" )

          F30 = -log(1-0.3)
          F10 = -log(1-0.1)

          Fref =  0.22
          Bmsy = BMSY[i]
          Bref = K[i] * 0.2
          BK = K[i]
          BK25 = K[i] * .25
          Fhistorical = mean( F[hdat,i], na.rm=T )
          Bhistorical = mean( B[hdat,i], na.rm=T )
          yl = 0.05


          abline (h=Fref, lty="solid", col="gray", lwd=2 )

          abline (h=F10, lty="dotted", col="gray")
          # text( 0.05*K[i], F10, "10% HR", pos=1 )

          abline (h=F30, lty="dotted", col="gray")
          # text( 0.05*K[i], F30, "30% HR", pos=1 )


          abline (h=FMSY[i], lty="dashed", col="red" )

          # abline (h=Fhistorical, lty="dashed")
          # text( 0.05*K[i], Fhistorical, "Mean", pos=1, lwd=2 )

          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)

          abline (v=Bmsy, lty="dotted")

          abline (v=BK, lty="dotted")

          abline (v=BK25, lty="dotted")

          text( 0.05*K[i], Fref, "20% HR", pos=1 )
          text( 0.05*K[i], FMSY[i], "FMSY", pos=3, lwd=2, col="red" )
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)
          text( Bmsy-0.01*K[i], yl, "K/2" , srt=90, pos=3)
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)
          text( B[hdat,i], F[hdat,i],  labels=yrs0, pos=3, cex= 0.8 )

          text( 0, ylims[2]*0.9,  labels=aulabels[i], pos=3, cex= 0.85 )



          # abline (v=Bhistorical, lty="dashed")
          # text( Bhistorical-0.01*K[i], yl, "Mean" , srt=90, pos=3,  lwd=2)

        }
      }

      if (vname=="simple") {

        if (is.null(fn))  fn=file.path(outdir, "hcr.simple.png" )

        require(car)

        B =  apply(  y[["B"]], c(2,3), median, na.rm=T  )
        F =  apply( y[["F"]], c(2,3), median, na.rm=T  )
        C =  apply( y[["C"]], c(2,3), median, na.rm=T  )
        K =  apply( y[["K"]], c(2), median, na.rm=T  )
        FMSY = apply( y[["FMSY"]], c(2), median, na.rm=T  )
        BMSY = apply( y[["BMSY"]], c(2), median, na.rm=T  )

        aulabels = c("N-ENS", "S-ENS", "4X")

        for (i in 1:3 ) {
          ylims = max(C[,i] )* c(0, 1.1)
          plot( B[hdat,i], C[hdat,i],  type="l", xlim=c(0, K[i]*1.05  ),
            ylim=ylims, xlab="", ylab="",  lwd=2, col="darkorange" )

          abline(0,0.1, lty="dotted", lwd=2, col="gray" )
          abline(0,0.2, lwd=3, col="gray" )
          abline(0,0.3, lty="dotted", lwd=2, col="gray" )

          points( B[hdat,i], C[hdat,i], col="orange", cex=0.8, pch=20 )
          text( B[hdat,i], C[hdat,i],  labels=yrs0, pos=3, cex= 0.85 )

          text( 0, ylims[2]*0.9,  labels=aulabels[i], pos=3, cex= 0.85 )

          # nn = as.matrix( cbind( Bx=as.vector( y$B[,ndata,i] ), Fx = as.vector( y$F[,ndata,i] ) ))
          # ellipse.2d(nn[,1], nn[,2], pv=0.05, sc=30)

          if (i==3) title( xlab="Fishable biomass (kt)" )
          if (i==2) title( ylab="Catch (kt)" )

          Cmsy = ( exp( FMSY[i] ) - 1)
          Cref = ( exp( FMSY[i] * 0.2 ) - 1) * K[i]
          Bmsy = BMSY[i]
          Bref = K[i] * 0.2
          BK = K[i]
          BK25 = K[i] * .25
          # Chistorical = mean( C[hdat,i], na.rm=T )
          # Bhistorical = mean( B[hdat,i], na.rm=T )
          yl = 0.1 * max(C[hdat,i])

          # abline (h=Fref, lty="dotted")
          # text( 0.25, Fref, "Target\n (0.2 * FMSY) ", pos=1 )

          abline (0, Cmsy, lty="dotted", col="red")
          # text( 0.1*K[i], Cmsy, "FMSY", pos=1 )

          # abline (h=Chistorical, lty="dashed")
          # text( 0.1*K[i], Chistorical, "Mean", pos=3, lwd=2 )

          # abline (v=Bref, lty="dotted")
          # text( Bref-0.2, 0.25, "Lower biomass reference point\n (LBRP = 0.2 * BMSY)" , srt=90, pos=3)


          abline (v=BK, lty="dotted")
          text( BK-0.01*K[i], yl, "K" , srt=90, pos=3)

          abline (v=BK/2, lty="dotted")
          text( BK/2-0.01*K[i], yl, "K/2" , srt=90, pos=3)

          abline (v=BK25, lty="dotted")
          text( BK25-0.01*K[i], yl, "K/4" , srt=90, pos=3)

        }


      }
    }

    if (type=="diagnostic.catch") {

    }


    if (type=="diagnostic.phase") {
      
      if (is.null(fn))  fn=file.path(outdir, "diagnostic.phase.png" )

        B =  apply(  y[["B"]], c(2,3), median, na.rm=T  )
        K =  apply( y[["K"]], c(2), median, na.rm=T  )

      for (i in 1:3 ) {
        plot( B[1:ndata-1,i], B[2:ndata,i],  type="b", xlab="t", ylab="t+1",
          xlim=c(0, K[i] * 1.25 ), ylim=max(K[i] )* c(0, 1.25), lwd=2, col="darkorange" )

         # abline(0,0.1, lty="dotted", lwd=2, col="gray" )
         abline( coef=c(0,1) )
       #text( B[1:ndata-1,], B[2:ndata,],  labels=yrs4 , pos=4, cex=0.8 )
      }
    }


    if (type=="diagnostic.errors") {

       if (is.null(fn))  fn=file.path(outdir, "diagnostic.errors.png" )

      # observation vs process error
      graphics.off()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))
      require(car)

      eP = y[["bpsd"]]  
      eO = y[["bosd"]]
      for (i in 1:3 ) {
          plot( eP[,,i], eO[,,i],  type="p", pch=22 )
          if (i==2) title( ylab="Process error (SD)" )
          if (i==3) title( xlab="Observation error (SD)" )

      }


    }


    if (type=="diagnostic.production") {
     
      B =  apply(  y[["B"]], c(2,3), median, na.rm=T  )
      P =  apply( y[["P"]], c(2,3), median, na.rm=T  )
      C =  apply( y[["C"]], c(2,3), median, na.rm=T  )
      K =  apply( y[["K"]], c(2), median, na.rm=T  )
      FMSY = apply( y[["FMSY"]], c(2), median, na.rm=T  )
      BMSY = apply( y[["BMSY"]], c(2), median, na.rm=T  )
      MSY = apply( y[["MSY"]], c(2), median, na.rm=T  )


      # production vs biomass
      plot.new()
      layout( matrix(c(1,2,3), 3, 1 ))
      par(mar = c(5, 4, 0, 2))
      for (i in 1:3) {
        plot( B[,i], P[,i], type="n", pch=20, ylim=c(0, max( c(P[,i], MSY[i]))*1.1), xlim=c(0,K[i]*1.05), xlab="Biomass; kt", ylab="Yield; kt"  )
        a = MSY[i] / (BMSY[i])^2
        curve( -a*(x-BMSY[i])^2 + MSY[i], from=0, to=K[i], add=TRUE, lwd=3, col="gray" )
        abline(v=BMSY[i], lty="dotted", lwd=3, col="gray")
        points( B[,i], P[,i], type="p", pch=20  )
        text( B[,i], P[,i], yrs0, pos=3, cex=0.8 )
        # abline(h=0)
      }

    }

    if (is.null(fn)) fn=file.path(outdir, "diagnostic.production.png" )
  
    if(save.plot) savePlot( filename=fn, type="png" )
    return( fn)

  }


}



