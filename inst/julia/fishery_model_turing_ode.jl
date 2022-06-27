

dir = expanduser("~/julia/snowcrab/")  # The directory of your package, for you maybe "C:\something"  
push!(LOAD_PATH, dir)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(dir)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL", "LazyArrays",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "ArviZ",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",   
  "Plots", "StatsPlots", "MultivariateStats"
]

  # using LazyArrays

# using Turing
# using DifferentialEquations, LinearAlgebra
# using MKL

# # Load StatsPlots for visualizations and diagnostics.
# using StatsPlots, Plots

# Turing.setadbackend(:forwarddiff)

# using ArviZ
# using PyPlot

# # ArviZ ships with style sheets!
# ArviZ.use_style("arviz-darkgrid")

#  Pkg.add( pkgs ) # add required packages

for pk in pkgs; @eval using $(Symbol(pk)); end


# https://catalyst.sciml.ai/stable/tutorials/compositional_modeling/


Turing.setprogress!(false);
Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)

# Turing.setrdcache(true)
# Turing.emptyrdcache()
 

# To clarify, “reverse mode” AD is efficient when you have a functions f(x) with small number of outputs fi and many inputs xj (in computing ∂fi/∂xj ), i.e. for functions mapping x∈Rm to f∈Rn with n≪m . (For example, in neural-network training where you want the derivative of one loss function ( n=1 ) with respect to millions ( m ) of network parameters. (The “manual” application of such a technique is also known as an adjoint method 37, and in the neural-net case it is called backpropagation.)

# In contrast, forward-mode AD (as in ForwardDiff.jl) is better when there is a small number of inputs and a large number of outputs, i.e. when n≫m , i.e. when you are computing many functions of a few variables. (It essentially corresponds to “manual” application of the chain rule in the most obvious way.)

# Zygote and ReverseDiff are both reverse-mode AD, but while ReverseDiff pushes custom types through your code to compute the backward pass (hence your code must be written to accept generic types), Zygote effectively rewrites the source code of your functions and works through more arbitrary code. 
 
# Turing.setadbackend(:forwarddiff)
 
fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn.RData"
o = load( fndat, convert=true)
Y = o["Y"]
Ksd = o["Ksd"]
Kmu = o["Kmu"]
removals = o["CAT"]
 

function Logistic!( du, u, p, t)
  du[1] = p[1] * u[1]  * (1.0 - u[1]/p[2])  # specific rate
end


if false
  # initial values for Logistic! 
  u0 = [0.1 ]
  p = (1.0, 1.0)
  tspan = (0.0, 5.0)
  prob = ODEProblem( Logistic!, u0,  tspan, p )
  res =  solve( prob, Tsit5(), saveat=0.1 )
  Plots.plot( res )
  # plot( res.t, reduce(hcat, res.u)' )  #res.u is an Vector{Vector{}} 
end




if false
  # initial values for Logistic!
  # almost working but solution decay to negative numbers though it is supposed to be bounded ..  
  au = 2
  
  pseason=0.8
  M = 3  # no years to project
  eps = 1e-9
  er = 0.2
  
  kmu = Kmu[au] * 1000
  ksd = Ksd[au] * 1000
   
  #p = ( 1.0, kmu, 0.75*kmu, 0.0, 24.0 )  # r,K,u0,t0,tinterval
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 10.0 )  #p[4]  , p[5]  

  si = Y[:,au]  # "survey index"
  N = length(si)
  dt = 0.1
  # specify u0 and tspan via functions
  
  u0_func(p) = p[3] 
  tspan_func(p) = (p[4], p[5])
  
  tspan = (0.0, 24.0)


  # musrandom start 2t add a time element for the b0 estimate .. cannot send directly a new initial condition 
  # so this simply adds an element for the first time step (0.1)
  removed = removals[:,au] .* 1000
  survey_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  fish_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  tstops = survey_time
 
  function affect_fishing!(integrator)
    i = findall(t -> t==integrator.t, fish_time)
    integrator.u[1] -=  removed[i[1] ] 
  end

  cb =  PresetTimeCallback( fish_time, affect_fishing! )
  
 

  prob = ODEProblem( Logistic!, [0.1], tspan, p, save_positions=(false, false) )
  msol =  solve( prob, Tsit5(), saveat=dt )
  plot( msol, label="ode, no fishing" )
  
  
  prob = ODEProblem( Logistic!, [p[3]], tspan_func, p, save_positions=(false, false) )
  msol =  solve( prob, Tsit5(), saveat=dt )
  plot!( msol, label="ode, no fishing random start 1" )

  
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 24.0 )  #p[4] tspan[0], p[5] is dtspan
  prob = remake(prob; f=Logistic!, u0=[u0_func(p)], tspan=tspan_func(p), p=p )
  msol =  solve( prob, Tsit5(), saveat=dt )
  plot!( msol, label="ode, no fishing random start 2" )
 
  
  p = ( 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 24.0 )  #p[4] tspan[0], p[5] is dtspan
  prob = remake(prob; f=Logistic!, u0=[u0_func(p)], tspan=tspan_func(p), p=p )
  msol =  solve( prob, Tsit5(), saveat=dt, callback=cb  )
  plot!( msol, label="ode, fishing", legend = :bottomright) 
 
  
  @model function fishery_model_turing_incremental_ode( si, kmu, ksd, removed, prob, N=length(si), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K  ~  TruncatedNormal( kmu, ksd, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    r ~  TruncatedNormal( 1.0, 0.1, 0.25, 2.0)   # (mu, sd)
    bpsd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
    q ~  TruncatedNormal( 1.0, 0.1, 1.0e-9, 2.5)  ; # i.e., Y:b scaling coeeficient
    qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # initial conditions
    ymean =  Vector{T}(undef, N)
    ymean[1] ~ Beta( 10.0, 1.0)  ; # starting b prior to first catch event

    # process model
    for i in 2:N
      p = T[ r, K, ymean[i-1]*K, i-1.0, i*1.0 ]
      msol = solve( 
        remake( prob, f=Logistic!, u0=[u0_func(p)], tspan=tspan_func(p), p=p[1:2] ), 
        Rosenbrock23(), 
        callback=cb,  
        saveat=dt 
      ) 
      if msol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
      end
      ym = max( last(msol.u)[1]/K, 1e-9 )
      ymean[i] ~ TruncatedNormal( ym, bpsd, 1e-9, 1.25)  ; 
      # @show ymean
    end

    # observation model
    @. si ~ TruncatedNormal( (ymean *q) + qc, bosd, 0.0, 1.25 ) 

  end

  
  @model function fishery_model_turing_ode( si, kmu, ksd, removed, prob, N=length(si), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K  ~  TruncatedNormal( kmu, ksd, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    r ~  TruncatedNormal( 1.0, 0.1, 0.25, 2.0)   # (mu, sd)
    bpsd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 10.0 )  ;  # slightly informative .. center of mass between (0,1)
    q ~  TruncatedNormal( 1.0, 0.1, 1.0e-9, 2.5)  ; # i.e., Y:b scaling coeeficient
    qc ~  TruncatedNormal( 0.0, 0.25, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # initial conditions
    ymean =  Vector{T}(undef, N)
    ymean[1] ~ Beta( 10.0, 1.0)  ; # starting b prior to first catch event

    p = T[ r, K, ymean[1]*K, 0.0, 24.0 ]

    # process model
    msol = solve( 
      remake( prob, f=Logistic!, u0=[u0_func(p)], tspan=tspan_func(p), p=p[1:2] ), 
      Rosenbrock23(), 
      callback=cb,  
      saveat=dt 
    ) 
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

    for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        ym = msol.u[j[1]][1]
        if typeof( ym ) !== Float64
          ym = ymean[1]*K
        end
        ym = max( ym, 1e-9 )
        ymean[i] ~ TruncatedNormal( ym/K, bpsd, 1e-9, 1.2)  ; 
      end
    end
 
    # observation model
    @. si ~ TruncatedNormal( (ymean *q) + qc, bosd, 0.0, 1.25 ) 

  end

  p = [ 1.0, kmu, rand(Beta(10,1))*kmu, 0.0, 23.0 ]  #p[4] tspan[0], p[5] is dtspan
  prob = ODEProblem( Logistic!, [0.1], tspan, p[1:2], saveat=dt, callback=cb    )


  prob = remake(prob; f=Logistic!, u0=[u0_func(p)], tspan=tspan_func(p), p=p[1:2] )
 
  fmod = fishery_model_turing_incremental_ode( si, kmu, ksd, removed, prob )
 
  res  =  sample( fishery_model_turing_incremental_ode( si, kmu, ksd, removed, prob ),  Turing.MH(), 10 )

  res  =  sample( fishery_model_turing_incremental_ode( si, kmu, ksd, removed, prob ),  Turing.NUTS(0.9), 1000 )
 


  fmod = fishery_model_turing_ode( si, kmu, ksd, removed, prob )

  res  =  sample( fishery_model_turing_ode( si, kmu, ksd, removed, prob ),  Turing.MH(), 10 )

  res  =  sample( fishery_model_turing_ode( si, kmu, ksd, removed, prob ),  Turing.NUTS(0.9), 1000 )
 


  plot!(; legend=false, ylim=(0,80000))
  for u in 500:1000  
    u0 = res[u,:K,1] * res[u,:"ymean[1]",1]
    prob2 = remake( prob, u0=[u0], tspan=(0.0, 24.0), p=[res[u,:r,1], res[u,:K,1]] )
    msol = solve(prob2, Tsit5(); saveat=0.1, callback=cb )
    plot!(msol; alpha=0.2, color="#DDDDDD")
  end
  
  u0 = Integer(floor(mean( res[[:"ymean[1]"]].value ) *  mean( res[[:K]].value ) ))

  pm = [mean( res[[:r]].value ), mean( res[[:K]].value ) ]

  msol = solve( remake( prob, u0=[u0], tspan=(0.0, 24.0), p=pm ), Tsit5(), saveat=0.1, callback=cb  ) #
  plot!(msol, label="ssa-fishing")


  prob = ODEProblem( Logistic!, [0.1], tspan, p[1:2], saveat=dt  )
  msol = solve( remake( prob, u0=[u0], tspan=(0.0, 24.0), p=pm ), Tsit5(); saveat=0.1  ) #
  plot!(msol, label="ssa-nofishing")

   
  # back transform si to normal scale 
  yhat = ( si  .- mean(res[[:"qc"]].value)) ./ mean(res[[:"q"]].value).* mean(res[[:"K"]].value) 
  scatter!(1:N, yhat   ; color=[1 2])
  plot!(1:N, yhat  ; color=[4 2])


  u = zeros(N)
  for i in 1:N
    u[i] = mean( res[:,Symbol("ymean[$i]"),:] .* res[:,:K,:] ) 
  end
  scatter!(1:N, u  ; color=[5 2])
 
  
  idata = from_mcmcchains( uu; library="Turing" )
 

# ArviZ ships with style sheets!
ArviZ.use_style("arviz-darkgrid")
  

plot_autocorr(res; var_names=["r", "K"]);
 

idata = from_mcmcchains( res; library="Turing" )

Plots.plot( survey_time , summarystats(idata.posterior; var_names=["ymean"]).mean )
Plots.plot!( survey_time , si )


Plots.plot!( survey_time , summarystats(idata.posterior; var_names=["ymean"]).mean .* mean( res[[:K]].value ) )


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


 
using Flux, DiffEqFlux
params = Flux.params(p)

msol =  solve( prob, Tsit5(), callback=cb, saveat=dt)  

function predict_rd() # Our 1-layer "neural network"
  solve(prob,Tsit5(),p=p,saveat=dt)[1] # override with new parameters
end

loss_rd() = sum(abs2,x-1 for x in predict_rd()) # loss function (squared absolute) x-1

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cbflux = function () #callback

  # must add a time element for the b0 estimate .. cannot send directly a new initial condition 
  # so this simply adds an element for the first time step (0.1)
  removed = removals[:,au] .* 1000
  survey_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  fish_time =  collect( 1.0:1.0:N ) .+ pseason  # time of observations for survey
  tstops = survey_time

  function Logistic!( du, u, p, t)
    du[1] = p[1] * u[1]  * (1.0 - u[1]/p[2])  # specific rate
  end
 
  function affect_init!(integrator)
    integrator.u[1] = u0[1] #b0
    reset_aggregated_jumps!(integrator)
  end

  function affect_fishing!(integrator)
    i = findall(t -> t==integrator.t, fish_time)
    integrator.u[1] -=  removed[i[1] ] 
    # reset_aggregated_jumps!(integrator)
  end

  cb = CallbackSet(  
    PresetTimeCallback( 0.1, affect_init! ),
    PresetTimeCallback( fish_time, affect_fishing! )
  )
 
  prob = ODEProblem( Logistic!, u0, tspan, p, save_positions=(false, false) )
  msol =  solve( prob, Tsit5(), saveat=dt )
  plot( msol, label="ode, no fishing" )
  
 function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=p),Tsit5(),saveat=dt), ylim=(0,kmu*2)))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_rd, params, data, opt, cb = cbflux)

m = Chain(
  Conv((2,2), 1=>16, relu),
  x -> maxpool(x, (2,2)),
  Conv((2,2), 16=>8, relu),
  x -> maxpool(x, (2,2)),
  x -> reshape(x, :, size(x, 4)),
  x -> solve(prob,Tsit5(),u0=x,saveat=0.1)[1,:],
  Dense(288, 10), softmax) |> gpu


