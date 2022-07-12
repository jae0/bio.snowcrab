
# ------------------------------
# Delay SSA size structured

# https://github.com/palmtree2013/DelaySSAToolkit.jl
# https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/
# NOTE::: require 03.snowcrab_carstm.r to be completed 


# from https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/:

# "A major assumption behind the majority of stochastic models of biochemical kinetics is the memoryless hypothesis, i.e., the stochastic dynamics of the reactants is only influenced by the current state of the system, which implies that the waiting times for reaction events obey exponential distributions. Gillespie developed a stochastic simulation algorithm (SSA) to simulate stochastic dynamics for such systems [1]. While this Markovian assumption considerably simplifies model analysis, it is dubious for modelling certain non-elementary reaction events that encapsulate multiple intermediate reaction steps [2]."

# For a few number of jumps, DelayRejection and DelayDirect will often perform better than other aggregators.

# For large numbers of jumps with sparse chain like structures and similar jump rates, for example continuous time random walks, DelayDirectCR and DelayMNRM often have the best performance.


dir = expanduser("~/julia/snowcrab/")  # The directory of your package, for you maybe "C:\something"  
push!(LOAD_PATH, dir)  # add the directory to the load path, so it can be found

import Pkg  # or using Pkg
Pkg.activate(dir)  # so now you activate the package
# Pkg.activate(@__DIR__()) #  same folder as the file itself.

Base.active_project()  # to make sure it's the package you meant to activate, print the path to console so you get a visual confirmation it's the package you meant to use

pkgs = [ 
  "Revise", "RData", "MKL",  "LazyArrays", "Flux", "StatsBase", "StaticArrays", "ForwardDiff", "DiffResults",
  "Turing", "Zygote", "Memoization", "ModelingToolkit", "Distributions", "DynamicPPL",
  "Catalyst", "DifferentialEquations", "LinearAlgebra",  "Interpolations", "DelaySSAToolkit",
  "Plots", "StatsPlots", "MultivariateStats"
]
 
for pk in pkgs; @eval using $(Symbol(pk)); end

#  Pkg.add( pkgs ) # add required packages

 


# ------------------------------
# Part 1 -- construct basic data and parameter list defining the main characteristics of the study
 
get_data_with_RCall = false

if get_data_with_RCall

    using RCall

    # typing <$> in Julia's  command prompt starts an R session.  
    
    $
      
    {
        # this is R-code that creates local RData file with required data
        # type <backspace> to escape back to julia
        source( file.path( code_root, "bio_startup.R" )  )
        require(bio.snowcrab)   # loadfunctions("bio.snowcrab")
        fishery_model_data_inputs( year.assessment=2021, type="size_structured_numerical_dynamics" )
    }

    # now back in Julia, fetch data into julia's workspace (replace fndat with the  filenane printed above )
    @rget Y  
    @rget Kmu 
    @rget removals 
    @rget ty
    @rget habitat_surface_areas
    
    # mechanism to run the rest if self contained
    include("/home/jae/bio/bio.snowcrab/inst/julia/fishery_model_turing_ode.jl")

else
    using  RData
    
    fndat = "/home/jae/bio.data/bio.snowcrab/modelled/1999_present_fb/fishery_model_results/turing1/biodyn_number_size_struct.RData"
    o = load( fndat, convert=true)
    Y = o["Y"]
    Kmu = o["Kmu"]
    removals = o["L"]
    habitat_surface_areas = o["hsa"]
end

plot( Y[:,:yr], Y[:,:cfasouth_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfasouth_M1] )
plot!( Y[:,:yr] .+2, Y[:,:cfasouth_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfasouth_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfasouth_M4] )


plot( Y[:,:yr], Y[:,:cfanorth_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfanorth_M1] )
plot!( Y[:,:yr] .+2 , Y[:,:cfanorth_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfanorth_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfanorth_M4] )


plot( Y[:,:yr], Y[:,:cfa4x_M0] )
plot!( Y[:,:yr] .+1 , Y[:,:cfa4x_M1] )
plot!( Y[:,:yr] .+2 , Y[:,:cfa4x_M2] )
plot!( Y[:,:yr] .+3, Y[:,:cfa4x_M3] )
plot!( Y[:,:yr] .+4, Y[:,:cfa4x_M4] )



# ------------------------------

Turing.setprogress!(false);
Turing.setadbackend(:zygote)
# Turing.setadbackend(:forwarddiff)
# Turing.setadbackend(:reversediff)
# Turing.setadbackend(:tracker)
 
 
 
function size_structured!( du, u, h, p, t)
  # M0, M1, M2, M3, M4, F  = u
  bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau, hsa  = p
  tr10 = vM1 * u[2]   # transition 1 -> 0; h index same as u index (offset by 1)   
  tr21 = vM2 * u[3]   # transitiom 2 -> 1
  tr32 = vM3 * u[4]   # transitiom 3 -> 2
  tr43 = vM4 * u[5]   # transitiom 4 -> 3
  f8  = h(p, t-8)[6]        # no fem 8 yrs ago
  du[1] = tr10             s- (mM0 * u[1]) * (u[1]/ K0)        
  du[2] = tr21      - tr10 - (mM1 * u[2]) * (u[2]/ K1)  
  du[3] = tr32      - tr21 - (mM2 * u[3]) * (u[3]/ K2) 
  du[4] = tr43      - tr32 - (mM3 * u[4]) * (u[4]/ K3) 
  du[5] = bM4 * f8 - tr43  - (mM4 * u[5]) * (u[5]/ K4)  
  du[6] = bF  * f8         - (mF0 * u[6]) * (u[6]/ KF)   # fem mat simple logistic with lag tau and density dep on present numbers
end

 
function size_structured_hsa!( du, u, h, p, t)
  # M0, M1, M2, M3, M4, F  = u
  bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau, hsa  = p
  tr10 = vM1 * u[2]   # transition 1 -> 0; h index same as u index (offset by 1)   
  tr21 = vM2 * u[3]   # transitiom 2 -> 1
  tr32 = vM3 * u[4]   # transitiom 3 -> 2
  tr43 = vM4 * u[5]   # transitiom 4 -> 3
  f8  = h(p, t-8)[6]        # no fem 8 yrs ago
  du[1] = tr10             - (mM0 * u[1]) * (u[1]/ K0) * hsa(t)       
  du[2] = tr21      - tr10 - (mM1 * u[2]) * (u[2]/ K1) * hsa(t) 
  du[3] = tr32      - tr21 - (mM2 * u[3]) * (u[3]/ K2) * hsa(t)
  du[4] = tr43      - tr32 - (mM3 * u[4]) * (u[4]/ K3) * hsa(t)
  du[5] = bM4 * f8 - tr43  - (mM4 * u[5]) * (u[5]/ K4) * hsa(t) 
  du[6] = bF  * f8         - (mF0 * u[6]) * (u[6]/ KF) * hsa(t)  # fem mat simple logistic with lag tau and density dep on present numbers
end

nS = 6 # n components

# history function 0.5 default
# h(p,t) = ones( nS ) .* 0.5  #values of u before t0
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
 
tau = 1  # delay

lags = [tau]

solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()


if false
  # testing DDE version -- initial values for size_structured! 
  
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* 100

  bM4=3.0; bF=1.2;
  K0=100; K1=100; K2=100;  K3=100; K4=100; KF=100;
  mM0=0.2; mM1=0.2; mM2=0.2; mM3=0.2; mM4=0.2; mF0=0.2;
  vM1=0.8; vM2=0.8; vM3=0.8; vM4=0.8;  
  tau=1
  
  p = ( bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau )

  tspan = (0.0, 100.0)
  
  prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
  res =  solve( prob,  solver, saveat=0.1 )
  plot( res )
  # plot( res.t, reduce(hcat, res.u)' )  #res.u is an Vector{Vector{}} 
 
end




# -------------------------
# other parameters

au = 2  # cfasouth
eps = 1e-9

tspan = (1990.0, 2030.0)

# convert to number .. 0.56 is ave mean weight of fb
kmu = Kmu[au] * 1000 *1000 / 0.56


M0 = Y[:,:cfasouth_M0]  # "survey index"
M1 = Y[:,:cfasouth_M1]  # "survey index"
M2 = Y[:,:cfasouth_M2]  # "survey index"
M3 = Y[:,:cfasouth_M3]  # "survey index"
M4 = Y[:,:cfasouth_M4]  # "survey index"
F = Y[:,:cfasouth_f_mat]  # "survey index"
 
N = length(M0)
dt = 0.1

survey_time = Y[:,:yrs]   # time of observations for survey

# this only adds habitat space  ... predation is also a useful one .. 
# speed is the issue 
forcing_time = survey
external_forcing = 1.0 - habitat_surface_areas[:,cfasouth_hsa] / max( habitat_surface_areas[:,cfasouth_hsa] ) # invert
hsa = scale(interpolate(external_forcing, BSpline(Linear())), forcing_time)


fish_time = removals[:,:ts]
removed = removals[:,:cfasouth]

function affect_fishing!(integrator)
  i = findall(t -> t==integrator.t, fish_time)
  integrator.u[1] -=  removed[ i[1] ] 
end

cb =  PresetTimeCallback( fish_time, affect_fishing! )
 
nS = 6 # n components

# history function 0.5 default
# h(p,t) = ones( nS ) .* 0.5  #values of u before t0
h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5 *kmu
 
tau = 1  # delay

lags = [tau]


solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()






# ---------------
if false
    #test run
    u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* kmu
  
    bM4=2.0; bF=1.2;
    K0=kmu; K1=kmu; K2=kmu;  K3=kmu; K4=kmu; KF=kmu;
    mM0=0.1; mM1=0.3; mM2=0.3; mM3=0.4; mM4=0.4; mF0=0.2;
    vM1=0.75; vM2=0.9; vM3=0.9; vM4=0.9;  
    tau=1

    forcing_time = survey_time
    external_forcing = ones(length(forcing_time))  # turns it off
    hsa = scale(interpolate(external_forcing, BSpline(Linear())), forcing_time)
 
    
    K = kmu  
      
    p = ( bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau, hsa )

    prob = DDEProblem( size_structured!, u0, h, tspan, p; constant_lags=lags )
    msol =  solve( prob,  solver, saveat=dt )
    plot( msol, label="dde, no fishing", legend=:left )
     
    prob = remake(prob; u0=u0, h=h, tspan=tspan, p=p, constant_lags=lags  )
    msol =  solve( prob,  solver, saveat=dt, callback=cb  )
    plot( msol, label="dde, fishing", legend = :bottomright) # , ylim=(0,kmu) 

    
    scatter(time_forcing, data_forcing,   label = "discrete forcing")
    plot!(time_forcing, g1_cst(time_forcing), label = "forcing1",  line = (:dot,   :red))

    
    prob = DDEProblem( size_structured_hsa!, u0, h, tspan, p; constant_lags=lags )
    msol =  solve( prob,  solver, saveat=dt )
    plot( msol, label="dde, hsa, no fishing", legend=:left )
     
end

 
# ---------------

@model function fishery_model_turing_dde( M0, M1, M2, M3, M4, F, kmu, tspan, prob, N=length(M0), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K0  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K1  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K2  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K3  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K4  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    KF  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    
    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

    q0 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q1 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q2 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q3 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q4 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    qF ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient

    qc0 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc1 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc2 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc3 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc4 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qcF ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # initial conditions
    m0 =  Vector{T}(undef, N)
    m1 =  Vector{T}(undef, N)
    m2 =  Vector{T}(undef, N)
    m3 =  Vector{T}(undef, N)
    m4 =  Vector{T}(undef, N)
    f0 =  Vector{T}(undef, N)

    m0[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m1[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m2[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m3[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m4[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    f0[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event

    # birth rate of M4 from F_8
    bM4 ~  TruncatedNormal( 1.0, 0.25, 0.25, 3.0)   # (mu, sd)
    bF  ~  TruncatedNormal( 1.0, 0.25, 0.25, 3.0)   # (mu, sd)
    
    # mortality
    mM0 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM1 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM2 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM3 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM4 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mF0 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    
    # transition rates
    vM1 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM2 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM3 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM4 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)

    # process model
    u0 = T[ m0[1]*K0, m1[1]*K1, m2[1]*K2, m3[1]*K3, m4[1]*K4, f0[1]*KF ] 

    p = ( bM4, bF, 
      K0, K1, K2, K3, K4, KF, 
      mM0, mM1, mM2, mM3, mM4, mF0, 
      vM1, vM2, vM3, vM4, 
      tau, hsa )
    
    msol = solve( remake( prob; u0=u0, h=h, tspan=tspan, p=p ), solver, callback=cb, saveat=dt )
    if msol.retcode != :Success
      Turing.@addlogprob! -Inf
      return nothing
    end

   for i in 2:N
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        m0[i] ~ TruncatedNormal( max( msol.u[j[1]][1] / K0, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m1[i] ~ TruncatedNormal( max( msol.u[j[1]][2] / K1, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m2[i] ~ TruncatedNormal( max( msol.u[j[1]][3] / K2, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m3[i] ~ TruncatedNormal( max( msol.u[j[1]][4] / K3, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m4[i] ~ TruncatedNormal( max( msol.u[j[1]][5] / K4, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        f0[i] ~ TruncatedNormal( max( msol.u[j[1]][6] / KF, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
      end
    end
  
    # observation model
    # @. M0 ~ TruncatedNormal( (f0 *q) + qc, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M0 ~ TruncatedNormal( (m0 + qc0) * q0, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M1 ~ TruncatedNormal( (m1 + qc1) * q1, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M2 ~ TruncatedNormal( (m2 + qc2) * q2, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M3 ~ TruncatedNormal( (m3 + qc3) * q3, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M4 ~ TruncatedNormal( (m4 + qc4) * q4, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. F  ~ TruncatedNormal( (f0 + qcF) * qF, bosd, 0.0, 1.25 )  # M0 (0,1)  
    
  
end


@model function fishery_model_turing_incremental_dde( M0, M1, M2, M3, M4, F, kmu, tspan, prob, N=length(M0), ::Type{T}=Float64 ) where {T}  
    # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
    # priors
    K0  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K1  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K2  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K3  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    K4  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)
    KF  ~  TruncatedNormal( kmu, kmu*0.25, kmu/10.0, kmu*10.0)   ; # (mu, sd)

    bpsd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)
    bosd ~  Beta( 1.0, 5.0 )  ;  # slightly informative .. center of mass between (0,1)

    q0 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q1 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q2 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q3 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    q4 ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient
    qF ~  TruncatedNormal( 1.0, 0.1, 0.1, 10.0)  ; # i.e., Y:b scaling coeeficient

    qc0 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc1 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc2 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc3 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qc4 ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    qcF ~  TruncatedNormal( 0.0, 0.1, -1.0, 1.0)  ; # i.e., Y:b offset constant   
    
    # initial conditions
    m0 =  Vector{T}(undef, N)
    m1 =  Vector{T}(undef, N)
    m2 =  Vector{T}(undef, N)
    m3 =  Vector{T}(undef, N)
    m4 =  Vector{T}(undef, N)
    f0 =  Vector{T}(undef, N)

    m0[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m1[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m2[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m3[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    m4[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event
    f0[1] ~  truncated( Cauchy( 0.5, 1.0), 0.1, 1.0 )  ; # starting b prior to first catch event

    # birth rate of M4 from F_8
    bM4 ~  TruncatedNormal( 1.0, 0.25, 0.25, 3.0)   # (mu, sd)
    bF  ~  TruncatedNormal( 1.0, 0.25, 0.25, 3.0)   # (mu, sd)
    
    # mortality
    mM0 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM1 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM2 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM3 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mM4 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 
    mF0 ~  TruncatedNormal( 0.2, 0.1, 1e-9, 0.9) 

    # transition rates
    vM1 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM2 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM3 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)
    vM4 ~  TruncatedNormal( 0.5, 0.1, 1e-9, 1.0)

    # process model
    t0 = floor(survey_time[1])

    for i in 2:N
      tsp = (t0+i-1.1, t0+i+0.1 )
   
      u0 = T[ m0[1]*K0, m1[1]*K1, m2[1]*K2, m3[1]*K3, m4[1]*K4, f0[1]*KF ] 
  
      p = ( bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau, hsa )

      msol = solve( 
        remake( prob; u0=u0, h=h, tspan=tsp, p=p, constant_lags=lags   ), 
        solver, 
        callback=cb,  
        saveat=dt 
      )
      if msol.retcode != :Success
        Turing.@addlogprob! -Inf
        return nothing
      end
      j = findall(t -> t==survey_time[i], msol.t)
      if length(j) > 0
        m0[i] ~ TruncatedNormal( max( msol.u[j[1]][1] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m1[i] ~ TruncatedNormal( max( msol.u[j[1]][2] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m2[i] ~ TruncatedNormal( max( msol.u[j[1]][3] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m3[i] ~ TruncatedNormal( max( msol.u[j[1]][4] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        m4[i] ~ TruncatedNormal( max( msol.u[j[1]][5] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
        f0[i] ~ TruncatedNormal( max( msol.u[j[1]][6] / K, 1e-9 ), bpsd, 1e-9, 1.25)  ; 
      end
    end
  
    # observation model
    # @. M0 ~ TruncatedNormal( (f0 *q) + qc, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M0 ~ TruncatedNormal( (m0 + qc0) * q0, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M1 ~ TruncatedNormal( (m1 + qc1) * q1, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M2 ~ TruncatedNormal( (m2 + qc2) * q2, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M3 ~ TruncatedNormal( (m3 + qc3) * q3, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. M4 ~ TruncatedNormal( (m4 + qc4) * q4, bosd, 0.0, 1.25 )  # M0 (0,1)  
    @. F  ~ TruncatedNormal( (f0 + qcF) * qF, bosd, 0.0, 1.25 )  # M0 (0,1)  

end
 
  

#  ----
#  run
  u0 = [ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ] .* kmu
 
  K = kmu  
  
  p = ( bM4, bF, K0, K1, K2, K3, K4, KF, mM0, mM1, mM2, mM3, mM4, mF0, vM1, vM2, vM3, vM4, tau, hsa )
  
  prob = DDEProblem( size_structured!, u0, h, tspan, p, constant_lags=lags )
  
  fmod = fishery_model_turing_incremental_dde( M0, M1, M2, M3, M4, F, kmu,  tspan, prob )

  fmod = fishery_model_turing_dde( M0, M1, M2, M3, M4, F, kmu, tspan, prob )


# testing
res  =  sample( fmod,  Turing.MH(), 3 )
res  =  sample( fmod,  Turing.NUTS( 3, 0.65), 3 )



n_samples = 1000
n_adapts = 500

res  =  sample( fmod,  Turing.NUTS(n_adapts, 0.65), n_samples )


t0 = floor(survey_time[1])

v = 1 # female
v = 6 # M0

for u in 1:100
  for  i in 1:N
    u0 = res[u,:K,1] * [ 
      res[u,Symbol("m0[$i]"),1],
      res[u,Symbol("m1[$i]"),1],
      res[u,Symbol("m2[$i]"),1],
      res[u,Symbol("m3[$i]"),1],
      res[u,Symbol("m4[$i]"),1],
      res[u,Symbol("f0[$i]"),1]
    ]
    
    tspn = ( t0+i-1.1, t0+i+0.1 )
    
    p = ( 
      res[u,:r,1], res[u,:K,1], res[u,:bM4,1], 
      res[u,:mM0,1], 
      res[u,:mM1,1], 
      res[u,:mM2,1], 
      res[u,:mM3,1], 
      res[u,:mM4,1], 
      res[u,:vM1,1], 
      res[u,:vM2,1], 
      res[u,:vM3,1], 
      res[u,:vM4,1], 
      tau, hsa 
    )
  
    msol = solve( 
      remake( prob, u0=u0, h=h, tspan=tspn, p=p  ), 
      solver, 
      callback=cb,  
      saveat=dt   ) #
    # plot!( msol; alpha=0.05, color=[3 4] )
    plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.05 ) 

  end
end


plot!(; legend=false, xlim=(1997,2023) )

u0 = mean( res[[:K]].value ) * [ 
  mean( res[[:Symbol("m0[$i]")]].value ),
  mean( res[[:Symbol("m1[$i]")]].value ),
  mean( res[[:Symbol("m2[$i]")]].value ),
  mean( res[[:Symbol("m3[$i]")]].value ),
  mean( res[[:Symbol("m4[$i]")]].value ),
  mean( res[[:Symbol("f0[$i]")]].value )
]

pm = ( 
  mean( res[[:r]].value ), mean( res[[:K]].value ), mean( res[[:bM4]].value ), 
  mean( res[[:mM0]].value ), 
  mean( res[[:mM1]].value ), 
  mean( res[[:mM2]].value ), 
  mean( res[[:mM3]].value ), 
  mean( res[[:mM4]].value ), 
  mean( res[[:vM1]].value ),
  mean( res[[:vM2]].value ),
  mean( res[[:vM3]].value ),
  mean( res[[:vM4]].value ),
  tau, hsa 
)

msol = solve( remake( prob, u0=u0, h=h, tspan=tspan, p=pm ), solver, callback=cb, saveat=dt )  
plot!(msol, label="ode-mean-fishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 
plot!(; legend=false, xlim=(1997,2023) )

prob2 = DDEProblem( size_structured!, u0, h, tspan, pm, saveat=dt )
msol = solve( prob2,  solver, saveat=dt ) #  effective nullify callbacks
plot!(msol, label="ode-mean-nofishing")
v = 1
plot!( msol.t, reduce(hcat, msol.u)'[:,v], color=[v v] , alpha=0.5 ) 

# back transform M0 to normal scale 
yhat = ( M0 .* mean(res[[:q0]].value) .- mean(res[[:qc0]].value )) .* mean(res[[:K]].value) 
scatter!(survey_time, yhat   ; color=[1 2])
plot!(survey_time, yhat  ; color=[1 2])
plot!(; legend=false, xlim=(1997,2023) )


# sample and plot means from model
w = zeros(N)
  for u in 1:length(res)  
    for i in 1:N
      w[i] = res[u,:K0,1] * res[u,Symbol("m0[$i]"),1]
    end
    plot!(survey_time, w  ;  alpha=0.1, color=[5 2])
  end

  plot!(; legend=false, xlim=(1997,2023) )


  u = zeros(N)
  v = zeros(N)
  for  i in 1:N
    u[i] = mean( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
    v[i] = std( res[:,Symbol("m0[$i]"),:] .* res[:,:K,:] ) 
  end
  scatter!(survey_time, u  ; color=[3 2])
  
  
  
# misc computed quantities

# params need to be named  .. return only that which is specified by "return()", below
pm = ( 
  r=r, K=K, bM4=bM4, 
  mM0=mM0, mM1=mM1, mM2=mM2, mM3=mM3, mM4=mM4, 
  vM1=vM1, vM2=vM2, vM3=vM3, vM4=vM4  
)

@model function fm_test( M0, M1, M2, M3, M4, F, kmu, tspan, prob, N=length(M0), ::Type{T}=Float64 ) where {T}  
  
  # deterministic computations: do from similations:
  M=3
  er=0.2

  F = zeros(N+M)
  B = zeros(N+M)
  C = zeros(N+M)

  C[1:N] = removed ./ K
  C[(N+1):(M+N)] = er .* bm[(N):(M+N-1)]
  C = 1.0 .- C / bm

  F =  -log( max.(C, eps) )  ;
  
  # parameter estimates for output
  MSY    = r* exp(K) / 4 ; # maximum height of of the latent productivity (yield)
  BMSY   = exp(K)/2 ; # biomass at MSY
  FMSY   = 2.0 * MSY / exp(K) ; # fishing mortality at MSY
 
  # recaled estimates
  B = bm .* K  
 
  
  # recaled estimates
  
  B[1:N] = bm[1:N] *. K - L[1:N] ;
  B[(N+1):(M+N)] = (bm[(N+1):(M+N)] - C[(N):(M+N-1)]) *. K ;

  
  return( test=r+1, )
end
fmod2 = fm_test(M0, M1, M2, M3, M4, F, kmu, tspan, prob )

gq = generated_quantities(fmod2, pm)

# gq = generated_quantities(res, values(p), keys(p))
# m = vec(getindex.( gq, 1))
# density!(m, lab="generated quantity (VI)")
# vline!([0], lab="true value")

 
  
  # look at predictions:
  M0_pred = Vector{Union{Missing, Float64}}(undef, length(M0))
  fmod_pred = fmod( M0_pred, kmu,  tspan, prob  ) 
 
  predictions = predict(fmod_pred, res)
  y_pred = vec(mean(Array(group(predictions, :M0)); dims = 1));
  
  plot( M0, y_pred )
  sum(abs2, M0 - y_pred) ≤ 0.1


 
 
 
do_variational_inference = false
if do_variational_inference
  # to do Variational Inference (an sd term goes less than 0 .. not sure of the cause ):
   
  res_vi =  vi(fmod, Turing.ADVI( 10, 1000));
 
     # Run sampler, collect results. @doc(Variational.ADVI) : 
     # samples_per_step::Int64
     # Number of samples used to estimate the ELBO in each optimization step.
     # max_iters::Int64
     # Maximum number of gradient steps.
  
  res_vi_samples = rand( res_vi, 1000)  # sample via simulation
  
  p1 = histogram(res_vi_samples[1, :]; bins=100, normed=true, alpha=0.2, color=:blue, label="")
  density!(res_vi_samples[1, :]; label="s (ADVI)", color=:blue, linewidth=2)
  density!(res, :s; label="s (NUTS)", color=:green, linewidth=2)
  vline!([var(x)]; label="s (data)", color=:black)
  vline!([mean(res_vi_samples[1, :])]; color=:blue, label="")
  
  p2 = histogram(res_vi_samples[2, :]; bins=100, normed=true, alpha=0.2, color=:blue, label="")
  density!(res_vi_samples[2, :]; label="m (ADVI)", color=:blue, linewidth=2)
  density!(res, :m; label="m (NUTS)", color=:green, linewidth=2)
  vline!([mean(x)]; color=:black, label="m (data)")
  vline!([mean(res_vi_samples[2, :])]; color=:blue, label="")
  
  plot(p1, p2; layout=(2, 1), size=(900, 500))
  
  
  do_maximum_likelihood = false
  if do_maximum_likelihood
    res_mle = Turing.optimize(fmod, MLE())
    res_mle = Turing.optimize(fmod, MLE())
    res_mle = Turing.optimize(fmod, MLE(), NelderMead())
    res_mle = Turing.optimize(fmod, MLE(), SimulatedAnnealing())
    res_mle = Turing.optimize(fmod, MLE(), ParticleSwarm())
    res_mle = Turing.optimize(fmod, MLE(), Newton())
    res_mle = Turing.optimize(fmod, MLE(), AcceleratedGradientDescent())
    res_mle = Turing.optimize(fmod, MLE(), Newton(), Optim.Options(iterations=10_000, allow_f_increases=true))

    using StatsBase
    coeftable(res_mle)
end


do_maximum_aposteriori = false
if do_maximum_aposteriori
  res_map = Turing.optimize(fmod, MAP())
end


   



 
using Flux, DiffEqFlux
params = Flux.params(p)

msol =  solve( prob,  solver, callback=cb, saveat=dt)  

function predict_rd() # Our 1-layer "neural network"
  solve(prob, solver,p=p,saveat=dt)[1] # override with new parameters
end

loss_rd() = sum(abs2,x-1 for x in predict_rd()) # loss function (squared absolute) x-1

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cbflux = function () #callback
  # function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=p), solver,saveat=dt), ylim=(0,kmu*2)))
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
  x -> solve(prob, solver,u0=x,saveat=0.1)[1,:],
  Dense(288, 10), softmax) |> gpu


