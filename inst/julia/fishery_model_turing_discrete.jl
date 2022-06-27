
# https://catalyst.sciml.ai/stable/tutorials/compositional_modeling/

# using Revise
# using Distributed

# using Turing
# Add four processes to use for sampling.
# addprocs(4)

using Turing, Tracker, Memoization
Turing.setadbackend(:tracker)

# Turing.setrdcache(true)
# Turing.emptyrdcache()

using Catalyst, DifferentialEquations, DiffEqBase, DiffEqJump, StochasticDiffEq
using MKL, Plots, StatsPlots, GraphViz, Latexify
 
# Turing.setadbackend(:forwarddiff)

using  RData  

fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn.RData"
o = load( fndat, convert=true)
Y = o["Y"]
Ksd = o["Ksd"]
Kmu = o["Kmu"]
removals = o["CAT"]
 

# To clarify, “reverse mode” AD is efficient when you have a functions f(x) with small number of outputs fi and many inputs xj (in computing ∂fi/∂xj ), i.e. for functions mapping x∈Rm to f∈Rn with n≪m . (For example, in neural-network training where you want the derivative of one loss function ( n=1 ) with respect to millions ( m ) of network parameters. (The “manual” application of such a technique is also known as an adjoint method 37, and in the neural-net case it is called backpropagation.)

# In contrast, forward-mode AD (as in ForwardDiff.jl) is better when there is a small number of inputs and a large number of outputs, i.e. when n≫m , i.e. when you are computing many functions of a few variables. (It essentially corresponds to “manual” application of the chain rule in the most obvious way.)

# Zygote and ReverseDiff are both reverse-mode AD, but while ReverseDiff pushes custom types through your code to compute the backward pass (hence your code must be written to accept generic types), Zygote effectively rewrites the source code of your functions and works through more arbitrary code. 
  

  # initial values for Logistic!
  # almost working but solution decay to negative numbers though it is supposed to be bounded ..  
  au = 2
  
  pseason = 0.8
  M = 3  # no years to project
  eps = 1e-9
  er = 0.2
  
  kmu = Kmu[au]
  ksd = Ksd[au]
   
  tspan = (0.0, 24.0 )
  si = Y[:,au]  # "survey index"
  N = length(si)

  removed =  Integer.( floor.( removals[:,au]  ))
  fish_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  survey_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  dt = 0.1
  tstops = survey_time

 
  
   
  @model function fishery_model_turing_incremental_discrete( si, kmu, ksd, removed, N=length(si) )
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K  ~  TruncatedNormal( kmu, ksd, 1e-9, Inf)   ; # (mu, sd)
    r ~  TruncatedNormal( 1.0, 0.25, 0.25, 2.0)   # (mu, sd)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    q ~  TruncatedNormal( 1.0, 0.2, 0.1, 2.5)  ; # i.e., Y:b scaling coefficient
    qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # ymean = fished (postfishery) abundance
    ymean =  tzeros(Real, N)
    ymean[1] ~ Beta( 10.0, 1.0)  ; # starting b prior to first catch event

    for i in 2:N
      ymean[i] ~ TruncatedNormal( r * ymean[i-1] * ( 1.0 - ymean[i-1] ) - removed[i-1]/K, bpsd, eps, 1.0)  ;
    end
      #  if msol.retcode != :Success
      #    Turing.@addlogprob! -Inf
      #    return nothing
      #  end

    ys = tzeros(Real, N)

    for i in 1:N
      si[i] ~ TruncatedNormal( (ymean[i] - qc)*q, bosd, 0.0, 1.0 )  ;
      ys[i] = ymean[i] * K
    end
    # @show (r, K, bosd, bpsd, q, qc)
    # @show ymean
  end
 
   
  

  fmod = fishery_model_turing_incremental_discrete( si, kmu, ksd, removed  )
   
  # Prior predictive check: 
  prior_res = sample(   fmod, Prior(), 10, nwarmup = 10, nchains = 3 );
  missing_data = Vector{Missing}(missing, N)
  mod_preds = fmod( missing_data, kmu, ksd, removed )
  prior_check = predict( fmod, prior_res )
  summarystats( prior_check )
  plot(prior_check)

  # Posterior sampling
  using Threads
  res = sample(   fmod, NUTS(0.95), MCMCThreads(), 2000, nwarmup = 1000, nchains = 3 );

  # res  =  sample( fmod,  Turing.SMC(), 10 )
  # res  =  sample( fmod,  Turing.NUTS(0.95), 10 )

  using MCMCChains, AdvancedHMC

  summarystats(res)

  density( res[:,:r,:])
  density( res[:,:K,:])
  density( res[:,:q,:])
  density( res[:,:qc,:])
  

  group( res, :ymean)  ##== res(:, 5:20, :) == res[[:ymean]]
  
  corner(res)

  plot( res[[:ymean]] )

  plot(
    traceplot(res),
    meanplot(res),
    density(res),
    histogram(res),
    mixeddensity(res),
    autocorplot(res),
    dpi=300, size=(840,600)
  )

  plot(res, seriestype=(:meanplot, :autocorplot), dpi=300)


  plot(res[:,[Symbol("ymean[$i]") for i in 1:10],:])

  using ArviZ
  using PyPlot



# ArviZ ships with style sheets!
ArviZ.use_style("arviz-darkgrid")
  

plot_autocorr(res; var_names=["r", "K"]);
gcf()

idata = from_mcmcchains( res; library="Turing" )

Plots.plot( survey_time , summarystats(idata.posterior; var_names=["ymean"]).mean )
Plots.plot!( survey_time , si )


Plots.plot( survey_time , summarystats(idata.posterior; var_names=["ys"]).mean )



oo = Array(res, (size(res)[1], size(res)[2]) )

  kmean = mean( res[[3]].value .*  res[[2]].value)

  pm = [mean( res[[1]].value ), mean( res[[2]].value ) ]

  odeprob = ODEProblem( Logistic!, [kmean], tspan, pm )
  bm = solve( odeprob, Rosenbrock23(); saveat=0.1, callback=cb, reltol=1e-15,abstol=1e-15 )
  plot!(bm, label="ode-fishing")

  prob = remake( prob, u0=[kmean], p=pm)
  bm = solve(prob, SSAStepper(),  saveat=0.1, callback=cb  ) #
  plot!(bm, lw=2, label="jump-fishing" )
 
  plot(; legend=false)
  posterior_samples = sample(res[[1, 2, 3]], 50; replace=false)
  for p in eachrow(Array(posterior_samples))
      bmpost = solve(prob, Tsit5(); p=pm, saveat=0.1)
      plot!(bmpost; alpha=0.1, color="#BBBBBB")
  end

  # Plot simulation and noisy observations.
  bm =  solve( prob, AutoTsit5(Rosenbrock23()); saveat=0.1, callback=cb, reltol=1e-16, abstol=1e-16 ) 
  plot!(bm; color=[1 2], linewidth=1)
  scatter!(bm.t, si'; color=[1 2])


end



  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(sN+M)
  B = zeros(N+M)
  C = zeros(N+M)

  C[1:N] = removals ./ K
  C[(N+1):(M+N)] = er .* bm[(N):(M+N-1)]
  C = 1.0 .- C / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = bm .* K  
 



bm = tzeros(Real, N+M, U)
  
bm[1] ~ TruncatedNormal( b0, bpsd, eps, 1.0)  ;
for i in 2:N 
  o = r * ( 1.0 - bm[i-1] ) ; 
  bm[i] ~ TruncatedNormal(   bm[i-1] * ( 1.0 + o ) - CAT[i-1]/K, bpsd, eps, 1.0)  ;
end
for i in (N+1):(M+N) 
  o = r * ( 1.0 - bm[i-1] ) ; 
  bm[i] ~ TruncatedNormal(   bm[i-1] * ( 1.0 + o ) - er*bm[(i-1)], bpsd, eps, 1.0)  ;
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


# north and south
if j in 1:2 {
  # spring surveys
  ys = ( Y[1, j] / q ) +  qc
  ys ~ TruncatedNormal( bm[1,j] - CAT[1,j]/K , bosd, eps, 1.0) ;
  for i in 2:(ty-1) 
    ys = ( Y[i, j] / q ) +  qc
    ys  ~ TruncatedNormal( bm[i,j] - CAT[i-1,j]/K , bosd, eps, 1.0)  ;
  end
  #  transition year (ty)
  ys = ( Y[ty,j] / q ) +  qc
  ys  ~ TruncatedNormal(  bm[ty,j]  - (CAT[ty-1,j]/K  + CAT[ty,j]/K ) / 2.0  , bosd, eps, 1.0)  ; #NENS and SENS
  # fall surveys
  for j in 1:U 
    for i in (ty+1):N 
      ys = ( Y[i,j] / q ) +  qc
      ys ~ TruncatedNormal(  bm[i,j] - CAT[i,j]/K, bosd, eps, 1.0)  ; #   fall surveys
    end
  end

end

# cfa4X
if j ==3
  # spring surveys
  for i in 1:(ty-1)  
    ys = ( Y[i, 3] / q[3] ) +  qc[3]
    ys  ~ TruncatedNormal( bm[i,3] - CAT[i,3]/K[3], bosd[3], eps, 1.0)  ;
  end
  #  transition year (ty)
  ys = ( Y[ty,3] / q[3] ) +  qc[3]
  ys  ~ TruncatedNormal(  bm[ty,3]  - CAT[ty,3]/K[3] , bosd[3], eps, 1.0)  ; #SENS
  # fall surveys
  for j in 1:U 
    for i in (ty+1):N 
      ys = ( Y[i,j] / q ) +  qc
      ys ~ TruncatedNormal(  bm[i,j] - CAT[i,j]/K, bosd, eps, 1.0)  ; #   fall surveys
    end
  end

end  

# deterministic computations: 
F = zeros(sN+M)
B = zeros(N+M)
C = zeros(N+M)

C[1:N] = CAT[1:N] ./ K
C[(N+1):(M+N)] = er .* bm[(N):(M+N-1)]
C = 1.0 -. C / bm

F =  -log( max.(C, eps) )  ;

   
# -------------------
# parameter estimates for output
MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
BMSY   = exp(K)/2 ; # biomass at MSY
FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY


# recaled estimates

B[1:N] = bm[1:N] *. K - CAT[1:N] ;
B[(N+1):(M+N)] = (bm[(N+1):(M+N)] - C[(N):(M+N-1)]) *. K ;


 