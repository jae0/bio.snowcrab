

using Zygote, Turing

@model function fishery_model_turing( Y, Kmu, Ksd, CAT; ty=6, M=3, er=0.2, eps=1e-9 )
    
  # -------------------
  # parameters
  
  N, U = size(Y)
  
  
  # -------------------
  # priors
  
  K  ~ arraydist( [ TruncatedNormal( Kmu[i], Ksd[i], eps, Inf) for i in 1:U] )  ; # (mu, sd)
  r ~ filldist( TruncatedNormal( 1.0, 0.1, 0.25, 2.0), 3 )  # (mu, sd)
  b0 ~ filldist( TruncatedNormal( 0.5, 0.1, eps, 1.0), 3) ; # starting b prior to first catch event
  bosd ~ filldist( truncated( Cauchy( 0.0, 0.5 ), eps, 1.0), 3) ;  # slightly informative .. center of mass between (0,1)
  bpsd ~ filldist( truncated( Cauchy( 0.0, 0.5 ), eps, 1.0), 3) ;
  q ~ filldist( TruncatedNormal( 1.0, 0.1, eps, 2.5), 3) ; # i.e., Y:b scaling coeeficient
  qc ~ filldist( truncated( Cauchy( 0.0, 0.1 ), -1.0, 1.0), 3) ; # i.e., Y:b offset constant  ; plot(dbeta(seq(0,1,by=0.1), 1, 10 ) )
  
  # -------------------
  # biomass process model
  # db/dt = r K b (1-b) - c ; b, c are normalized by K  

  # max .. force positive value .. initial conditions
  
  bm = tzeros(Real, N+M, U)
  
  for j in 1:U 
    bm[1,j] ~ TruncatedNormal( b0[j], bpsd[j], eps, 1.0)  ;
    for i in 2:N 
      o = r[j] * ( 1.0 - bm[i-1,j] ) ; 
      bm[i,j] ~ TruncatedNormal(   bm[i-1,j] * ( 1.0 + o ) - CAT[i-1,j]/K[j], bpsd[j], eps, 1.0)  ;
    end
    for i in (N+1):(M+N) 
      o = r[j] * ( 1.0 - bm[i-1,j] ) ; 
      bm[i,j] ~ TruncatedNormal(   bm[i-1,j] * ( 1.0 + o ) - er*bm[(i-1),j], bpsd[j], eps, 1.0)  ;
    end
  end
  
  # -------------------
  # biomass observation model
  # cfanorth(1) and cfasouth(2)
  #   This is slightly complicated because a fall / spring survey correction is required:
  #   B represents the total fishable biomass available in fishing year y
  #     in fall surveys:    Btot(t) = Bsurveyed(t) + removals(t)
  #     in spring surveys:  Btot(t) = Bsurveyed(t) + removals(t-1)
  # spring surveys from 1998 to 2003
  #   this is conceptualized in the following time line:
  #     '|' == start/end of each new fishing year
  #     Sf = Survey in fall
  #     Ss = Survey in spring
  #     |...(t-2)...|.Ss..(t-1)...|...(t=2004)..Sf.|...(t+1).Sf..|...(t+2)..Sf.|...
  # Cfa 4X -- fall/winter fishery
  #    Btot(t) = Bsurveyed(t) + removals(t)  ## .. 2018-2019 -> 2018
  
  
  # -------------------
  # spring surveys
  for j in 1:2 
    ys = ( Y[1, j] / q[j] ) +  qc[j]
    ys ~ TruncatedNormal( bm[1,j] - CAT[1,j]/K[j] , bosd[j], eps, 1.0) ;
    for i in 2:(ty-1) 
      ys = ( Y[i, j] / q[j] ) +  qc[j]
      ys  ~ TruncatedNormal( bm[i,j] - CAT[i-1,j]/K[j] , bosd[j], eps, 1.0)  ;
    end
  end
  
  for i in 1:(ty-1)  
    ys = ( Y[i, 3] / q[3] ) +  qc[3]
    ys  ~ TruncatedNormal( bm[i,3] - CAT[i,3]/K[3], bosd[3], eps, 1.0)  ;
  end
  
  
  # -------------------
  #  transition year (ty)
  for j in 1:2 
    ys = ( Y[ty,j] / q[j] ) +  qc[j]
    ys  ~ TruncatedNormal(  bm[ty,j]  - (CAT[ty-1,j]/K[j]  + CAT[ty,j]/K[j] ) / 2.0  , bosd[j], eps, 1.0)  ; #NENS and SENS
  end
  ys = ( Y[ty,3] / q[3] ) +  qc[3]
  ys  ~ TruncatedNormal(  bm[ty,3]  - CAT[ty,3]/K[3] , bosd[3], eps, 1.0)  ; #SENS
  
  # fall surveys
  for j in 1:U 
    for i in (ty+1):N 
      ys = ( Y[i,j] / q[j] ) +  qc[j]
      ys ~ TruncatedNormal(  bm[i,j] - CAT[i,j]/K[j], bosd[j], eps, 1.0)  ; #   fall surveys
    end
  end
  
  
  
  # -------------------
  # fishing mortality
  # fall fisheries
  F = tzeros(Real, N+M, U)
  B = tzeros(Real, N+M, U)
  C = tzeros(Real, N+M, U)
  MSY = tzeros(Real,  U)
  BMSY = tzeros(Real,  U);
  FMSY = tzeros(Real,  U);
  
  for j in 1:U 
    for i in 1:N 
      F[i,j] =  -log( max( 1.0 - CAT[i,j]/K[j]  / bm[i,j], eps) )  ;
    end
    for i in (N+1):(M+N)  
      F[i,j] =  -log( max( 1.0 - er * bm[i-1,j] / bm[i,j], eps) )  ;
    end
  end
  
  # -------------------
  # parameter estimates for output
  for j in 1:U 
    MSY[j]    = r[j]* exp(K[j]) / 4 ; # maximum height of of the latent productivity (yield)
    BMSY[j]   = exp(K[j])/2 ; # biomass at MSY
    FMSY[j]   = 2.0 * MSY[j] / exp(K[j]) ; # fishing mortality at MSY
  end
  
  # recaled estimates
  for j in 1:U 
    for i in 1:N 
      B[i,j] = bm[i,j] * K[j] - CAT[i,j] ;
      C[i,j] = CAT[i,j] ;
    end
    
    for i in (N+1):(M+N) 
      B[i,j] = (bm[i,j] - er*bm[(i-1),j]) * K[j] ;
      C[i,j] = er*bm[(i-1),j] * K[j] ;
    end
    
  end
  
end 

# end @model fishery 
 