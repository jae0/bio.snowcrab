Turing.@model function logistic_discrete_turing_historical( PM )
  # biomass process model: dn/dt = r n (1-n/K) - removed ; b, removed are not normalized by K  
  # priors 
  K ~ truncated(Normal( PM.K[1], PM.K[2]), PM.K[3], PM.K[4])  
  r ~ truncated(Normal( PM.r[1], PM.r[2]), PM.r[3], PM.r[4])   # (mu, sd)
  bpsd ~ truncated(Normal( PM.bpsd[1], PM.bpsd[2]), PM.bpsd[3], PM.bpsd[4] ) # slightly informative .. center of mass between (0,0.5)
  bosd ~ truncated(Normal( PM.bosd[1], PM.bosd[2]), PM.bosd[3], PM.bosd[4] ) # slightly informative .. center of mass between (0,0.5)
  q1 ~ truncated(Normal( PM.q1[1], PM.q1[2]), PM.q1[3], PM.q1[4] )    

  # m's are "total available for fishery" (latent truth)
  m = tzeros( PM.nM )
  m[1] ~ truncated(Normal(PM.m0[1], PM.m0[2]), PM.m0[3], PM.m0[4])   ; # starting b prior to first catch event

  for i in 2:PM.nT
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ) - PM.removed[i-1]/K, bpsd), PM.mlim[1], PM.mlim[2])  ;
  end

  for i in (PM.nT+1):PM.nM
    m[i] ~ truncated(Normal( m[i-1] + r * m[i-1] * ( 1.0 - m[i-1] ), bpsd), PM.mlim[1], PM.mlim[2])  ; # predict with no removals (prefishery)
  end
 
  if any( x -> x < 0.0, m)  # permit overshoot
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood
  # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale
  # observation model: S = (m - q0)/ q1   <=>   m = S * q1 + q0  
  # see function: abundance_from_index      
  if PM.yeartransition == 0
    # 4X
    for i in PM.iok
      PM.S[i] ~ Normal(  (K * m[i] - PM.removed[i]) / q1, K * bosd )  ; # fall survey
    end

  else
    # cfanorth and cfasouth:
    # spring to fall survey: transition year = 2004
    # spring = 1:5
    # fall = 6:last
    
    for i in PM.iok
      if  i < PM.yeartransition
        PM.S[i] ~ Normal(  K * m[i] / q1, K * bosd )  ;  # spring survey
      elseif i == PM.yeartransition
        PM.S[i] ~ Normal(  ( K * m[i] - (PM.removed[i-1] + PM.removed[i]) / 2.0) / q1, K * bosd )  ;  # transition year  .. averaging should be done before .. less computation 
      else
        PM.S[i] ~ Normal(  ( K * m[i] - PM.removed[i] ) / q1 , K * bosd )  ; # fall survey
      end
    end
  end

end



function expand_grid(; kws...)
  names, vals = keys(kws), values(kws)
  return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end
 

 
function discretize_decimal( x, delta=0.01 ) 
  num_digits = Int(ceil( log10(1.0 / delta)) )   # time floating point rounding
  out = round.( round.( x ./ delta; digits=0 ) .* delta; digits=num_digits)
  return out
end

   
function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=-1 )

  nchains = size(res)[3]
  nsims = size(res)[1]
 
  if n_sample == -1
    # do all
    n_sample = nchains * nsims
    oo = expand_grid( sims=1:nsims, chains=1:nchains)
  else
    oo = DataFrame( sims=rand(1:nsims, n_sample), chains=rand(1:nchains, n_sample) )
  end

  md = zeros(nM, n_sample) 
  mb = zeros(nM, n_sample)
   
  z = 0

  while z <= n_sample 
    z += 1
    z > n_sample && break
    j = oo[z, :sims]  # nsims
    l = oo[z, :chains] # nchains
    for i in 1:nM
      md[i,z] = res[j, Symbol("m[$i]"), l]
      mb[i,z] = md[i,z] * res[j, Symbol("K"), l]
    end
  end

  # additional nothings to keep same expectations as continuous models
  return (md, nothing, mb, nothing, nothing, nothing )  

end


function logistic_discrete_reference_points(r, K)
  expK = exp.(K) 
  msy   = r .* expK ./ 4.0 ; # maximum height of of the latent productivity (yield)
  bmsy  = expK ./ 2.0 ; # biomass at MSY
  fmsy  = 2.0 .* msy ./ expK ; # fishing mortality at MSY
  return (msy, bmsy, fmsy)
end


function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end

 

function showall( x )
    # print everything to console
    show(stdout, "text/plain", x) # display all estimates
end


function plots_diagnostic( res, vn="K" ) 
  gr()
  pl = plot()
  pl = density!(pl, res[ Symbol(vn) ])
  return pl
end


function fishery_model_mortality(; removed=removed, bio=bio, survey_time=survey_time )    
  fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete 
  Fkt = removed
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  # FM[ FM .< eps(0.0)] .= zero(eltype(FM))
  return ( Fkt, FR, FM  )
end


function abundance_from_index( Sai, res  )
  # map S <=> m  where S = observation index on unit scale; m = latent, scaled abundance on unit scale

  K =  vec( res[:,Symbol("K"),:] )'
  q1 =  vec( res[:,Symbol("q1"),:] )'
  S_m = Sai .* q1   # already on K scale
  
  return S_m
end



function fishery_model_plot(; toplot=("fishing", "survey"), n_sample=min(250, size(bio)[2]),
  res=res, bio=bio, FM=FM, 
  S=S,
  prediction_time=prediction_time, prediction_time_ss=prediction_time_ss, survey_time=survey_time, yrs=yrs, 
  alphav=0.075, pl= plot(), time_range=(floor(minimum(prediction_time_ss))-1.0, ceil(maximum(prediction_time_ss))+0.5  )
)
 
  nsims = size(bio)[2]
  ss = rand(1:nsims, n_sample)  # sample index

  if any(isequal.("trace", toplot))  
    @warn "trace is not valid for a discrete model"
    
  end 

  if any(isequal.("nofishing", toplot))  
    @warn "nofishing not implemented"
    
  end 

  # extract sims (with fishing) -- prefishery biomass
  # plot biomass
  if any(isequal.("fishing", toplot))   # this means there is fishing occuring ( and plot pre-fishery biomass )
    g = bio   # [ yr,  sim ]
    pl = plot!(pl, prediction_time_ss, g[:,ss] ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time_ss, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end
 

  # extract sims (with fishing)
  # plot biomass
  if any(isequal.("postfishery", toplot))   # this means there is fishing occuring ( and plot post-fishery biomass )
    g = bio .- PM.removed    # [ yr,  sim ] .- landings
    pl = plot!(pl, postfishery_time_ss, g[:,ss] ;  alpha=alphav, color=:orange)
    pl = plot!(pl, postfishery_time_ss, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("footprint", toplot))  
    @warn "footprint not implemented"
    
  end
   

  if any(isequal.("survey", toplot))    # scale index (at survey time)
    # map S -> m and then multiply by K
    # where S=observation on unit scale; m=latent, scaled abundance on unit scale
    S_m = abundance_from_index( S, res  )
    S_K = mean(S_m, dims=2)  # average by year
    pl = plot!(pl, survey_time, S_K, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, S_K, markersize=4, color=:darkgray)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )

  end
   

  if any(isequal.("fishing_mortality", toplot))  
 
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, prediction_time_ss, FM[:,ss] ;  alpha=0.02, color=:lightslateblue)
    pl = plot!(pl, prediction_time_ss, FMmean ;  alpha=0.8, color=:slateblue, lw=4)
    pl = plot!(pl, ylim=(0, ub ) )
    pl = plot!(pl ; legend=false )
    pl = plot!(pl; xlim=time_range )

  end


  if any(isequal.("fishing_mortality_vs_footprint", toplot))  
    @warn "footprint not implemented"
    

  end


  if any(isequal.("harvest_control_rule_footprint", toplot))  
    @warn "footprint not implemented"
    

  end
   

  if any(isequal.("harvest_control_rule", toplot))  

    r = vec( Array(res[:, Symbol("r"), :]) )
    K = vec( Array(res[:, Symbol("K"), :]) ) 
    (msy, bmsy, fmsy) = logistic_discrete_reference_points(r, K)

    pl = hline!(pl, fmsy[ss]; alpha=0.01, color=:lightgray )
    pl = hline!(pl, [mean(fmsy)];  alpha=0.6, color=:darkgray, lw=5 )
    pl = hline!(pl, [quantile(fmsy, 0.975)];  alpha=0.5, color=:gray, lw=2, line=:dash )
    pl = hline!(pl, [quantile(fmsy, 0.025)];  alpha=0.5, color=:gray, lw=2, line=:dash )
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(prediction_time_ss)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
    fb = bio[1:length(prediction_time_ss),:]
    fb_mean = mean(fb, dims=2)
    fm_mean = mean(FM, dims=2)
  
    fbbb = [quantile(fb[nt,:], 0.025), quantile(fb[nt,:], 0.975) ]

    FMbb = [quantile(FM[nt,:], 0.975), quantile(FM[nt,:], 0.025) ]
     
    pl = scatter!(pl, [fb[nt,:]], [FM[nt,:]] ;  alpha=0.01, color=:magenta, markersize=2.5, markerstrokewidth=0)
    pl = scatter!(pl, fbbb, FMbb;  alpha=0.5, color=:magenta, markershape=:star, markersize=6, markerstrokewidth=1)

    pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.9, color=:gold, markersize=8, markerstrokewidth=1)
    
    pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)
    pl = scatter!(pl,  fb_mean, fm_mean;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0  )
    pl = scatter!(pl,  fb_mean .+0.051, fm_mean .-0.0025;  alpha=0.8, color=colours,  markersize=0, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, prediction_time_ss), :top, :left, pointsize=8) )

    ub = max( quantile(K, 0.75), maximum( fb_mean ), maximum(fmsy) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, quantile(fmsy, 0.975)  ) )
  
  end
   
 

  if any(isequal.("harvest_control_rule_postfishery", toplot))  

    r = vec( Array(res[:, Symbol("r"), :]) )
    K = vec( Array(res[:, Symbol("K"), :]) ) 
    (msy, bmsy, fmsy) = logistic_discrete_reference_points(r, K)

    pl = hline!(pl, fmsy[ss]; alpha=0.01, color=:lightgray )
    pl = hline!(pl, [mean(fmsy)];  alpha=0.6, color=:darkgray, lw=5 )
    pl = hline!(pl, [quantile(fmsy, 0.975)];  alpha=0.5, color=:gray, lw=2, line=:dash )
    pl = hline!(pl, [quantile(fmsy, 0.025)];  alpha=0.5, color=:gray, lw=2, line=:dash )
  
    pl = vline!(pl, K[ss];  alpha=0.05, color=:limegreen )
    pl = vline!(pl, K[ss]./2;  alpha=0.05, color=:darkkhaki )
    pl = vline!(pl, K[ss]./4;  alpha=0.05, color=:darkred )
  
    pl = vline!(pl, [mean(K)];  alpha=0.6, color=:chartreuse4, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)];  alpha=0.5, color=:chartreuse4, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/2.0];  alpha=0.6, color=:darkkhaki, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/2.0;  alpha=0.5, color=:darkkhaki, lw=2, line=:dash )
  
    pl = vline!(pl, [mean(K)/4.0];  alpha=0.6, color=:darkred, lw=5 )
    pl = vline!(pl, [quantile(K, 0.975)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
    pl = vline!(pl, [quantile(K, 0.025)]/4.0;  alpha=0.5, color=:darkred, lw=2, line=:dash )
  
    nt = length(prediction_time_ss)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
    fb = bio[1:length(prediction_time_ss),:] .- PM.removed
    fb_mean = mean(fb, dims=2)
    fm_mean = mean(FM, dims=2)
  
    fbbb = [quantile(fb[nt,:], 0.025), quantile(fb[nt,:], 0.975) ]

    FMbb = [quantile(FM[nt,:], 0.975), quantile(FM[nt,:], 0.025) ]
     
    pl = scatter!(pl, [fb[nt,:]], [FM[nt,:]] ;  alpha=0.01, color=:magenta, markersize=2.5, markerstrokewidth=0)
    pl = scatter!(pl, fbbb, FMbb;  alpha=0.5, color=:magenta, markershape=:star, markersize=6, markerstrokewidth=1)

    pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.9, color=:gold, markersize=8, markerstrokewidth=1)
    
    pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)
    pl = scatter!(pl,  fb_mean, fm_mean;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0  )
    pl = scatter!(pl,  fb_mean .+0.051, fm_mean .-0.0025;  alpha=0.8, color=colours,  markersize=0, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, prediction_time_ss), :top, :left, pointsize=8) )

    ub = max( quantile(K, 0.75), maximum( fb_mean ), maximum(fmsy) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, quantile(fmsy, 0.975)  ) )
  
  end
   
 

  if any(isequal.("prediction_space", toplot))  
    # b(t+1) vs b(t) * K / q1
    
    K = vec( Array(res[:, Symbol("K"), :]) )[ss] 
    q1 = vec( Array(res[:, Symbol("q1"), :]) )[ss] 

    t1 = 2:length(survey_time)
    t0 = 1:(length(survey_time)-1)

    b1 = bio[t1,ss] ./ K' .* q1'  # b(t+1)
    b0 = S[t0] ./ K' .* q1'   # b(t)

    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:(nt-1), (nt-1) )]
    pl = scatter!(pl, b0, b1;  alpha=0.2, color=colours, markersize=2.5, markerstrokewidth=0)

  end


  if any(isequal.("state_space", toplot))  
    # b(t+1) vs b(t) * K / q1
    K = vec( Array(res[:, Symbol("K"), :]) )[ss]
    q1 = vec( Array(res[:, Symbol("q1"), :]) )[ss]

    t1 = 2:length(survey_time)
    t0 = 1:(length(survey_time)-1)

    b = bio[:,ss]  ./ K' .* q1' 
    b1 = b[ t1,:]  # b(t+1)
    b0 = b[ t0,:]  # b(t)

    nsims = size(bio)[2]
    ss = rand(1:nsims, 1000)  # sample index
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:(nt-1), (nt-1) )]
    pl = scatter!(pl, b0, b1; legend=false, alpha=0.2, color=colours, markersize=2.5, markerstrokewidth=0)

  end

 

  if any(isequal.("surplus_production", toplot))  
    # $S(B) = rB (1-B/K)$ (i.e., "Schaefer" 1954 form )

    K = vec( Array(res[:, Symbol("K"), :]) )[ss]
    r = vec( Array(res[:, Symbol("r"), :]) )[ss]
    q1 = vec( Array(res[:, Symbol("q1"), :]) )[ss]

    b = bio[:,ss]  ./ K' .* q1'
    sp = r' .* b  .* (1.0 .- b) # surplus prod
 
    nsims = size(bio)[2]
    ss = rand(1:nsims, 1000)  # sample index
  
    nt = length(survey_time)
    colours = get(ColorSchemes.tab20c, 1:nt, :extrema )[rand(1:nt, nt )]
    for i in 1:nt
      pl = scatter!(pl, b[i,:], sp[i,:];  alpha=0.3, color=colours[i], markersize=2.0, markerstrokewidth=0)
    end
    pl = plot(pl; legend=false, xlab="Biomass", ylab="Surplus Production")
  
  end

  return(pl)
 
end


function plot_prior_posterior( vn, prior, posterior; bw=0.02 )
  pri =  vec(collect( prior[:,Symbol(vn),:] ))
  pos =  vec(collect( posterior[:,Symbol(vn),:] ))
  pl = plot(ylabel="Density", xlabel=vn ) 
  pl = density!(pl, pri,  fill=true, color = :slateblue, fillalpha=0.25, bandwidth = bw, lw=0, label="Prior")
  pl = density!(pl, pos,  fill=true, color = :purple, fillalpha=0.5, bandwidth = bw, lw=0, label="Posterior")
  return(pl)
end

  