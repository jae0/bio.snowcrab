

dir = expanduser("~/julia/snowcrab/")  # The directory of your package, for you maybe "C:\something"  
push!(LOAD_PATH, dir)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(dir)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL", 
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  
  "Plots", "StatsPlots", "MultivariateStats"
]

# using ArviZ
# using PyPlot

# # ArviZ ships with style sheets!
# ArviZ.use_style("arviz-darkgrid")

#  Pkg.add( pkgs ) # add required packages

for pk in pkgs; @eval using $(Symbol(pk)); end


# https://catalyst.sciml.ai/stable/tutorials/compositional_modeling/


Turing.setprogress!(false);
# Turing.setadbackend(:zygote)
Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)

# Turing.setrdcache(true)
# Turing.emptyrdcache()


fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn.RData"
o = load( fndat, convert=true)
Y = o["Y"]
Ksd = o["Ksd"]
Kmu = o["Kmu"]
removals = o["CAT"]
 
 
rs = @reaction_network begin
  r, X--> 2X
  r*X/K, X --> 0
end r K 



  if false

    @show rs_ode = convert(ODESystem, rs)

    # as an SPDE:
    sprob = SDEProblem(rs, u0, tspan, p)
    ssol  = solve(sprob, LambaEM(), reltol=1e-3)
    plot(ssol,lw=2,title="Adaptive SDE: Birth-Death Process")

    Graph(rs)

    # as ode
    u0 = [:X => 1.0, :Y => 1.0, :XY => 1.0, :Z1 => 1.0, :Z2 => 1.0]
    oprob = ODEProblem(rn, u0, tspan, [])

    species(rs)
    parameters(rs)
    reactions(rs)
    latexify(rs, starred=true)
    latexify( )
    g = Graph(rs)
    odesys = convert(ODESystem, repressilator)

    using Plots, GraphRecipes
    plot(TreePlot(rs), method=:tree, fontsize=12, nodeshape=:ellipse)
 
  end

  # initial values for Logistic!
  # almost working but solution decay to negative numbers though it is supposed to be bounded ..


  au = 2
  
  pseason = 0.8
  M = 3  # no years to project
  eps = 1e-9
  er = 0.2
  
  kmu = Kmu[au] * 1000
  ksd = Ksd[au] * 1000
   
  tspan = (0.0, 24.0 )
  si = Y[:,au]  # "survey index"
  N = length(si)

  removed =  Integer.( floor.( removals[:,au] .* 1000 ))
  fish_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  survey_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  dt = 0.1
  tstops = survey_time


  # specify u0 and tspan via functions
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 24.0 )  #p[4] tspan[0], p[5] is dtspan
   
  u0_func(p) = Integer(floor(p[3] ))
  tspan_func(p) = (p[4], p[5])
  
  tspan = (p[4], p[5])

  function affect_fishing!(integrator)
    i = findall(t -> t==integrator.t, fish_time)
    integrator.u[1] -=  removed[i[1]] 
    reset_aggregated_jumps!(integrator)
  end
 
  cb = PresetTimeCallback( fish_time, affect_fishing! )
  
  
  
  p = ( 0.99, 61560.4181, 0.94*61560, 0.0, 24.0 )  #p[4] tspan[0], p[5] is dtspan
   
  dprob = DiscreteProblem(rs, [u0_func(p)], tspan_func(p), p[1:2] )
  
  # basc dynamics with no fishing
  prob = JumpProblem(rs, dprob, Direct(),  saveat=dt )





  # simple run
  bm = solve( prob, SSAStepper(), saveat=0.1  )
  plot!(bm, lw=2, title=("Gillespie solution for $au"), label="simple"  )
  
  # simple run with fishing
  bm = solve( prob, SSAStepper(),  saveat=0.1, callback=cb  ) #
  plot!(bm, lw=2, label="fishing" )
  

  
  # testing use of remake
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 30.0 )  #p[4] tspan[0], p[5] is dtspan
  prob2 = remake( prob, u0=[u0_func(p)], tspan=tspan_func(p), p=p[1:2])
  bm2 = solve(prob2, SSAStepper(),  saveat=0.1, callback=cb    ) 
  plot!(bm2, lw=2, label="test")
  
 
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 30.0 )  #p[4] tspan[0], p[5] is dtspan
  prob2 = remake( prob, u0=[u0_func(p)], tspan=tspan_func(p), p=p[1:2])
  bm = solve( prob2, SSAStepper(),  saveat=0.1, callback=cb  ) #
   
  plot!(bm, lw=2, label="fishing - lower b0" )

   
   
  @model function fishery_model_turing_incremental_ssa( si, kmu, ksd, removed, prob, N=length(si), ::Type{T} = Float64) where {T}
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
     # priors
     K  ~  TruncatedNormal( kmu, ksd, kmu/10.0, kmu*10.0)   ; # (mu, sd)
     r ~  TruncatedNormal( 1.0, 0.25, 0.25, 2.0)   # (mu, sd)
     bosd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
     bpsd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
     q ~  TruncatedNormal( 1.0, 0.25, 0.1, 2.5)  ; # i.e., Y:b scaling coeeficient
     qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
     # initial conditions
     ymean =  Vector{T}(undef, N)
     ymean[1] ~ Beta( 10.0, 1.0)  ; # starting b prior to first catch event
     # process model
     for i in 2:N
       #p = T[r, K, ymean[i-1]*K, i-1.0, i*1.0 ]
        # @show (typeof(p[3]), length(p[3]))
        #  if typeof( p[3] ) !== Float64
        #   Turing.@addlogprob! -Inf
        #   return nothing
        #  end
        jp = remake( prob, u0=T[ymean[i-1]*K], tspan=(i-1.0, i*1.0 ), p=T[r, K] )
        msol = solve( jp, SSAStepper(), callback=cb  ) 
        j = findall(t -> t==survey_time[i], msol.t)
          if length(j) > 0 
            ymean[i] ~ TruncatedNormal( msol.u[j[1]][1]/K, bpsd, 1e-9, 1.2)  ; 
          else
            Turing.@addlogprob! -Inf
            return nothing
          end

      end
     # @show ymean
     # observation model
    #  for i in 1:N
    #    si[i] ~ TruncatedNormal( (ymean[i] + qc) * q, bosd, 1e-9, 1.2 )
    #  end 

    @. si ~ TruncatedNormal( (ymean .+ qc) .* q, bosd, 1e-9, 1.2 ) 
  end
 
    
   
  # @model function fishery_model_turing_ssa( si, kmu, ksd, removed, prob, N=length(si)  ) 
  #   # single global model is still to slow to use for parameter estimation
  #   # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  #    # priors
  #    K  ~  TruncatedNormal( kmu, ksd, kmu/10.0, kmu*10.0)   ; # (mu, sd)
  #    r ~  TruncatedNormal( 1.0, 0.25, 0.25, 2.0)   # (mu, sd)
  #    bosd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
  #    bpsd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
  #    q ~  TruncatedNormal( 1.0, 0.25, 0.1, 2.5)  ; # i.e., Y:b scaling coeeficient
  #    qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
  #    # initial conditions
  #    ymean =  tzeros(N)
  #    ymean[1] ~ Beta( 10.0, 1.0)  ; # starting b prior to first catch event
  #    p = (r, K, ymean[1]*K, 0.0, 25.0 )
  #    msol = solve( remake( prob, u0=[ymean[1]*K], tspan=tspan_func(p), p=[r, K] ), SSAStepper(), callback=cb, tstops=survey_time  ) 
  #    # process model
  #    for i in 2:N
  #      j = findall(t -> t==survey_time[i], msol.t)
  #      if length(j) > 0
  #       ymean[i] ~ TruncatedNormal( msol.u[j[1]][1]/K, bpsd, 1e-9, 1.2)  ; 
  #      end
  #    end
  #    # observation model
  #    @. si ~ TruncatedNormal( (ymean + qc) * q, bosd, 1e-9, 1.2 ) 
  # end
 
  # too slow to use
  # fmod = fishery_model_turing_ssa( si, kmu, ksd, removed, prob  )

  fmod = fishery_model_turing_incremental_ssa( si, kmu, ksd, removed, prob  )

  # test run and spin up compilation
  res  =  sample( fmod,  Turing.MH(), 10 )
  res  =  sample( fmod,  Turing.NUTS(0.65), 10 , max_depth=12 ) 

  res  =  sample( fmod,  Turing.PG(), 10 )

  res  =  sample( fmod,  Turing.MH(), 1000,  thinning=10 , discard_initial=500  )

  res  =  sample( fmod,  Turing.MH(), MCMCThreads(), 1000, 4, thinning=100, discard_initial=1000 )

  res  =  sample( fmod,  Turing.NUTS(0.65), 1000 , max_depth=12 ) 
 
  
  
  
  # prob = JumpProblem( rs, dprob, Direct(),  callback=cb, saveat=dt  ) # tstops=tstops, 
  # res[:,[Symbol("ymean[1]") ]]

  for u in sample(1:100, 10) 
    u0 = res[u,:K,1] * res[u,:"ymean[1]",1]
    prob2 = remake( prob, u0=[u0], tspan=(0.0, 24.0), p=[res[u,:r,1], res[u,:K,1]] )
    msol = solve( prob2, SSAStepper(), saveat=0.1, callback=cb  ) #
    plot!(msol; alpha=0.2, color="#CCCCCC")
  end
  plot!(; legend=false)

  k0 = Integer(floor(mean( res[[:"ymean[1]"]].value ) *  mean( res[[:K]].value ) ))

  pm = [mean( res[[:r]].value ), mean( res[[:K]].value ) ]

  msol = solve( remake( prob, u0=[k0], tspan=(0.0, 24.0), p=pm ), SSAStepper(), saveat=0.1, callback=cb  ) #
  plot!(msol, label="ssa-fishing")


  msol = solve( remake( prob, u0=[k0], tspan=(0.0, 24.0), p=pm ), SSAStepper(), saveat=0.1   ) #
  plot!(msol, label="ssa-nofishing")

  # back transform si to normal scale 
  yhat = ( si ./ mean(res[[:"q"]].value) .- mean(res[[:"qc"]].value)) .* mean(res[[:"K"]].value) 
  scatter!(1:N, yhat   ; color=[1 2])
  plot!(1:N, yhat  ; color=[1 2])


  
  # look at predictions:
  si_pred = Vector{Union{Missing, Float64}}(undef, length(si))
  prob2 = JumpProblem(rs, dprob, Direct(),  saveat=dt , callback=cb )
  fmod_pred = fmod( si_pred, kmu, ksd, removed, prob  ) 
 
  predictions = predict(fmod_pred, res)
  y_pred = vec(mean(Array(group(predictions, :si)); dims = 1));
  
  plot( si, y_pred )
  sum(abs2, si - y_pred) ≤ 0.1



# Generate a MLE estimate.
mle_estimate = optimize(fmod, MLE())

# Generate a MAP estimate.
map_estimate = optimize(fmod, MAP())









plot_autocorr(res; var_names=["r", "K"]);
 

idata = from_mcmcchains( res; library="Turing" )

Plots.plot( survey_time , summarystats(idata.posterior; var_names=["ymean"]).mean )
Plots.plot!( survey_time , si )


Plots.plot!( survey_time , summarystats(idata.posterior; var_names=["ymean"]).mean .* mean( res[[:K]].value ), legend=:bottomright, ylim=(0,55000) )




  # Plot observations.
  bm = is .* K . ./ q .+ qc 

  plot!( fish_time, bm; color=[1 2], linewidth=1)
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


 