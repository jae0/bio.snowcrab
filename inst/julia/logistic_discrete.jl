
theme(:default)  # defaults for graphics
#theme(:vibrant)
#theme(:bright)


# to start a graphics window
# gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=4)
gr()

 
Y = o["Y"][∈(yrs).(o["Y"].yrs), :]
removals = o["L"][∈(yrs).(o["L"].yrs), :]

Kmu = [5.0, 60.0, 1.25]
 
nT = length(yrs)
nP = 5  # number of predictions into future (with no fishing)
nM = nP + nT  # total number of prediction years
dt = 1  # time resolution of solutions .. discrete annual so 1 year
nS = 1 # no state variables, not used 
no_digits = 3  # time floating point rounding 
smallnumber = 1.0e-9 # floating point value sufficient to assume 0 valued
    
# "survey index"
S = Y[:,Symbol("$aulab"  )]

# scale index where required
Smean = mean(skipmissing(S))
Sstd = std( skipmissing(S))
Smin = minimum(skipmissing(S))
Smax = maximum( skipmissing(S))
Srange = Smax - Smin 

SminFraction = Smin ./ Srange  # used as informative prior mean in some runs

# id index
ki = aulab=="cfanorth" ? 1 :
     aulab=="cfasouth" ? 2 :
     aulab=="cfa4x"    ? 3 :
     0  # default

kmu = Kmu[ki]  

survey_time = Y[:,:yrs]   # time of observations for survey
# prediction_time = floor.(vcat( survey_time, collect(1:nP) .+ maximum(survey_time) ) )     
fish_time = round.( round.( removals[:,:ts] ./ dt; digits=0 ) .* dt; digits=no_digits)   # time of observations for landings

removed = removals[:,Symbol("$aulab")]
  
smallnumber = 1.0 / kmu / 10.0  # floating point value of sufficient to assume 0 valued

no_digits = 3  # time floating point rounding

survey_time =  round.( round.( Y[:,:yrs] ./ dt; digits=0 ) .* dt ; digits=no_digits)    # time of observations for survey
 
predtime =  4.0/12.0  # predictions ("m") are prefishery .. arbitrarily setting to 4/12
prediction_time =
  floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  #yrs
  round( predtime/dt ; digits=no_digits)   # april (m== prefishery)

iok = findall( !ismissing, S )
 
# basic params for "logistic_discrete"
PM = (
  yrs=yrs,
  nS = nS, 
  nT = length(yrs),
  nP = nP,  # number of predictions into future (with no fishing)
  nM = nM,  # total number of prediction years
  K = (kmu, 0.25*kmu, kmu/5.0, kmu*5.0 ),
  r = (1.0, 0.1, 0.5, 1.5),
  bpsd = ( 0.1, 0.05, 0.01, 0.5 ),
  bosd = ( 0.1, 0.05, 0.01, 0.5 ),
  q1 = (  1.0, 0.1,  0.01, 10.0 ),
  q0 = ( SminFraction, 0.1, -1.0, 1.0 ),
  m0 = ( 0.9, 0.2, 0.1, 1.0), 
  mlim =(0.0, 1.0),
  removed=removed,
  S = S,
  iok = iok,
  yeartransition = 0
) 

# run model estimations / overrides
Turing.setprogress!(false);

if (aulab=="cfanorth") | (aulab=="cfasouth")
  PM = @set PM.yeartransition = 6
end

PM = @set PM.K = (kmu, 0.25*kmu, kmu/10.0, kmu*10.0 )
PM = @set PM.r = (1.0, 0.1, 0.25, 2.0)
PM = @set PM.bpsd = (0, 0.1, 1.0e-9, 0.5) 
PM = @set PM.bosd = (0, kmu*0.5, 1.0e-9, kmu/2.0)
PM = @set PM.q1 = (1.0, 0.5,  1.0e-1, 10.0)
PM = @set PM.m0 = ( 5, 5)
PM = @set PM.mlim = ( 0.0, 1.25)

fmod = logistic_discrete_turing_historical( PM )  # q1 only

# Turing NUTS-specific default options  ..  see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
n_adapts, n_samples, n_chains = 10000, 10000, 4  # kind of high .. but mixing is poor in this model parameterization
target_acceptance_rate = 0.99
max_depth=14  ## too high and it become impossibly slow
init_ϵ = 0.01

# by default use NUTS sampler ... SMC is another good option if NUTS is too slow
turing_sampler = Turing.NUTS(n_samples, target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )
print( "\n\nSampling for: ", aulab, year_assessment, "\n\n" )


# Julia-code to generate model solutions and posteriors. 
# This is an abbrevariated version of:
# https://github.com/jae0/dynamical_model/blob/master/snowcrab/04.snowcrab_fishery_model.md 
 
res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains ) # sample in parallel 

# if threading is not working (MSWindows?) try:
# res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

# mcmc save file name and location
res_fn = joinpath( model_outdir, string("results_turing", "_", aulab, ".hdf5" ) )  
@save res_fn res
print( "\n\n", "Results file:",  res_fn, "\n\n" )

