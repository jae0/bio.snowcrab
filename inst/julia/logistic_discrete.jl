 
# data

Y = o["Y"][∈(yrs).(o["Y"].yrs), :]  # index of abundance (kt)
removals = o["L"][∈(yrs).(o["L"].yrs), :]  # landings (kt)

nT = length(yrs) 
nP = 0  # number of predictions into future (with no fishing)
nM = nP + nT  # total number of prediction years
dt = 1  # time resolution of solutions .. discrete annual so 1 year
nS = 1 # no state variables, not used 
no_digits = 3  # time floating point rounding 
smallnumber = 1.0e-9 # floating point value sufficient to assume 0 valued
    
# "survey index"
S = Y[:,Symbol("$aulab"  )]  # select subarea inputs for modelling
 
# id index for each subarea
ki = aulab=="cfanorth" ? 1 :
     aulab=="cfasouth" ? 2 :
     aulab=="cfa4x"    ? 3 :
     0  # default

Kmu = [5.0, 55.0, 1.5] # vector of guess of magnitude of carrying capacity (kt)
kmu = Kmu[ki]  

smallnumber = 1.0 / kmu / 10.0  # floating point value of sufficient to assume 0 valued
no_digits = 3  # time floating point rounding

survey_time = Y[:,:yrs]    # time of observations for survey .. for plotting

# adjust survey time to be approximately in spring or autumn
spring = findall( x -> x < 2004, survey_time )
fall   = findall( x -> x >= 2004, survey_time )
survey_time[spring] = survey_time[spring] .+ 5.0/12.0  # arbitrary ..  "spring"
survey_time[fall]   = survey_time[fall] .+ 10.0/12.0   # time of survey in "fall" ; all of 4X sep-mar gets called "fall" for simplicity in plot
survey_time =  round.(  survey_time ; digits=no_digits)    # time of observations for survey
 
removed = removals[:,Symbol("$aulab")]

# predictions ("m") are "prefishery" .. arbitrarily setting to 4/12
predtime =  4.0/ 12.0

postfishery = 11.9 / 12.0   
 
prediction_time = floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  round( predtime/dt ; digits=no_digits)   # april (m== prefishery)
yrs_pred_report = findall( x -> (x >= 1999.0) & (x <= (Real(year_assessment)+nP +1.0)), prediction_time )
prediction_time_ss = prediction_time[yrs_pred_report]

postfishery_time = floor.( vcat( collect(minimum(yrs) : (maximum(yrs)+nP) ) )  ) .+  round( postfishery/dt ; digits=no_digits)   # april (m== prefishery)
yrs_pred_report = findall( x -> (x >= 1999.0) & (x <= (Real(year_assessment)+nP +1.0)), postfishery_time )
postfishery_time_ss = postfishery_time[yrs_pred_report]

iok = findall( !ismissing, S )  # index of data locations

# basic parameters for "logistic_discrete_historical"
PM = (
  yrs=yrs,
  nS = nS, 
  nT = length(yrs),
  nP = nP,  # number of predictions into future (with no fishing)
  nM = nM,  # total number of prediction years
  K = kmu .* (1.0, 0.25, 1.0/5.0, 5.0 ),
  r = (1.0, 0.1, 0.01, 2.0),
  bpsd = ( 0.25, 0.1, 0.01, 0.5 ),
  bosd = ( 0.25, 0.1, 0.01, 0.5 ),
  q1 = ( 1.0, 0.25, 0.01, 2.0 ),
  mlim =( 0.001, 0.999 ),
  m0 = (0.5, 0.1, 0.001, 1.25 ),
  removed=removed,
  S = S,
  iok = iok,
  yeartransition = 0
) 

if (aulab=="cfanorth") | (aulab=="cfasouth")
  PM = @set PM.yeartransition = 6
end
  
if (! @isdefined outformat )
  outformat ="png"
end

if (! @isdefined model_variation )
  model_variation = "logistic_discrete_historical"   
end


if (! @isdefined outputs_dir )
  outputs_dir = joinpath( homedir(), "bio.data", "bio.snowcrab", "fishery_model", string(year_assessment), model_variation )
end

mkpath( outputs_dir )


# set up Turing model
fmod = logistic_discrete_turing_historical( PM )  # q1 only

# Turing NUTS-specific default options  ..  see write up here: https://turing.ml/dev/docs/using-turing/sampler-viz
n_adapts, n_samples, n_chains = 10000, 10000, 4  # kind of high .. but mixing is poor in this model parameterization

target_acceptance_rate, max_depth, init_ϵ = 0.75, 9, 0.05

# by default use NUTS sampler ... SMC is another good option if NUTS is too slow
turing_sampler = Turing.NUTS( target_acceptance_rate; max_depth=max_depth, init_ϵ=init_ϵ )
print( "\n\nSampling for: ", aulab, ": ", year_assessment, "\n\n" )

Turing.setprogress!(true);

# generate model solutions and posteriors. This is an abbrevariated version of:
# https://github.com/jae0/dynamical_model/blob/master/snowcrab/04.snowcrab_fishery_model.md 
print( "\nPrior predictive sampling\n" )
prior = sample( fmod, Prior(), n_samples) # Prior predictive check
@show prior

print( "\nPosterior predictive sampling\n" )
res  =  sample( fmod, turing_sampler, MCMCThreads(), n_samples, n_chains; drop_warmup=true ) # sample in parallel 
@show res

# if threading is not working (MSWindows?) try:
# res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

# mcmc save file name and location
res_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, ".hdf5" ) )  
@save res_fn res
# @load res_fn res
print( "\n\n", "Results file:",  res_fn, "\n\n" )

summary_directory = joinpath( outputs_dir, string("results_turing", "_", aulab ) )
print( "\n\n", "Posteriors of transformed parameters saved at: ", summary_directory, "\n\n" )

# n scaled, n unscaled, biomass of fb with and without fishing, model_traces, model_times 
m, num, bio, trace, trace_bio, trace_time = fishery_model_predictions(res ) 

# fishing (kt), relative Fishing mortality, instantaneous fishing mortality:
Fkt, FR, FM = fishery_model_mortality() 

# save a few data files as semicolon-delimited CSV's for use outside Julia
summary_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, "_summary", ".csv" ) )  
CSV.write( summary_fn,  summarize( res ), delim=";" )  # use semicolon as , also used in parm names
  
bio_fn1 = joinpath( outputs_dir, string("results_turing", "_", aulab, "_bio_fishing", ".csv" ) )  
CSV.write( bio_fn1,  DataFrame(bio[:,:], :auto), delim=";" )  # use semicolon as , also used in parm names

fm_fn = joinpath( outputs_dir, string("results_turing", "_", aulab, "_fm", ".csv" ) )  
CSV.write( fm_fn,  DataFrame( FM, :auto), delim=";" )  # use semicolon as , also used in parm names

print( "\n\n", "Plots being created at: ",  outputs_dir, "\n\n" )

# annual snapshots of biomass (kt) prefishery
pl = fishery_model_plot( toplot=("survey", "fishing" ) )
savefig(pl, joinpath( outputs_dir, string("plot_predictions_", aulab, ".", outformat) )  )
 
# annual snapshots of biomass (kt) postfishery
pl = fishery_model_plot( toplot=("survey", "postfishery" ) )
savefig(pl, joinpath( outputs_dir, string("plot_predictions_postfishery_", aulab, ".", outformat) )  )
 
# plot fishing mortality, , 
pl = fishery_model_plot( toplot="fishing_mortality" )
savefig(pl, joinpath( outputs_dir, string("plot_fishing_mortality_", aulab, ".", outformat) )  )
 
# HCR plot
pl = fishery_model_plot( toplot="harvest_control_rule", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("plot_hcr_", aulab, ".", outformat) )  )

# HCR plot
pl = fishery_model_plot( toplot="harvest_control_rule_postfishery", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("plot_hcr_postfishery_", aulab, ".", outformat) )  )

# grey is prior, purple is posterior 
L = truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])
pl = plot_prior_posterior( "r", prior, res )
save(joinpath( outputs_dir, string("plot_prior_r_", aulab, ".", outformat) ), pl  )

L = truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4]) 
pl = plot_prior_posterior( "K", prior, res, bw=0.04)
save(joinpath( outputs_dir, string("plot_prior_K_", aulab, ".", outformat) ), pl  )
 
L = truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    
pl = plot_prior_posterior( "q1", prior, res)
save(joinpath( outputs_dir, string("plot_prior_q1_", aulab, ".", outformat) ), pl  )
 
L = truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] )
pl = plot_prior_posterior( "bpsd", prior, res)
save(joinpath( outputs_dir, string("plot_prior_bpsd_", aulab, ".", outformat) ), pl  )
 
L = truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] ) 
pl = plot_prior_posterior( "bosd", prior, res)
save(joinpath( outputs_dir, string("plot_prior_bosd_", aulab, ".", outformat) ), pl  )
 
pl = fishery_model_plot( toplot="surplus_production", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("plot_surplus_production_", aulab, ".", outformat) )  )
 
pl = fishery_model_plot( toplot="state_space", n_sample=1000 ) #, alphav=0.01 )  # hcr
savefig(pl, joinpath( outputs_dir, string("plot_state_space_", aulab, ".", outformat) )  )
 


print( "\n\n", string("Fishery model for ", aulab, " completed."), "\n\n" )

