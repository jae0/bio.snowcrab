using Turing
  

function size_structured_dde!( du, u, h, p, t )
  # here u, du are actual numbers .. not normalized by K due to use of callbacks

  # b, K, d, d2, v, tau, hsa  = p
  b, K, d, d2, v  = p
 
  @inbounds begin

      # this break down seems to speed it up a bit ... not sure why
      br =  b .* h(p, t-8.0)[6]  
      tr =  v .* h(p, t-1.0)[2:5]
      dh = d2 .* u ./ hsa(t, 1:6)
      dr = d .* u  .+  dh .* u 
      
      du[1] = tr[1] * K[2] / K[1]            - dr[1]       # note:
      du[2] = tr[2] * K[3] / K[2]   - tr[1]  - dr[2]
      du[3] = tr[3] * K[4] / K[3]   - tr[2]  - dr[3]
      du[4] = tr[4] * K[5] / K[4]   - tr[3]  - dr[4]
      du[5] = br[1] * K[6] / K[5]   - tr[4]  - dr[5]
      du[6] = br[2]                          - dr[6]      # fem mat simple logistic with lag tau and density dep on present numbers
    
    end

end


function dde_parameters()
    # these are dummy initial values .. just to get things started
    b=[1.7, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .*kmu;
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    v=[0.65, 0.68, 0.61, 0.79];
    d2 = [0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    # params = ( b, K, d, d2, v, tau, hsa)
    params = ( b, K, d, d2, v)
    return params
end



function β( mode, conc )
  # alternate parameterization of beta distribution 
  # conc = α + β     https://en.wikipedia.org/wiki/Beta_distribution
  beta1 = mode *( conc - 2  ) + 1.0
  beta2 = (1.0 - mode) * ( conc - 2  ) + 1.0
  Beta( beta1, beta2 ) 
end 



@model function size_structured_dde_turing_generic( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
 
  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~ filldist(  Gamma(2.0, 0.05),  nS ) # #  working: β(0.1, 10.0);  plot(x->pdf(β(0.01, 8), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
  =#
  d ~   filldist( LogNormal( -1.508, 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -1.508, 0.5 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
  =#
  v ~   filldist( LogNormal( 0.3782, 0.5 ),  4 ) # transition rates # plot(x->pdf(β(0.99, 10), x), xlim=(0,1))  

  u0 ~  filldist( Beta(1, 1), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v )

  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      abstol=1.0e-10, 
      reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
      end
  end

end


@model function size_structured_dde_turing_reference( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
 
  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~ filldist(  Gamma(2.0, 0.1),  nS ) # #  working: β(0.1, 10.0);  plot(x->pdf(β(0.01, 8), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    log( exp(0.5)-1.0 ) = log(0.6487 ) = -0.4328 .. 50%
  =#
  d ~   filldist( LogNormal( -1.508, 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -0.7096, 0.5 ), nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.80)-1.0) = log(1.226)  = 0.2034  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.691) = 0.5255    .. ~95% (moult) transition rate per year (mode) 
    
  =#
  v ~   filldist( LogNormal( 0.2034, 0.5 ),  4 ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  

  u0 ~  filldist( Beta(2, 2), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v  )

  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      # abstol=1.0e-10, 
      # reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[k] )  # observation and process error combined
      end
  end

end


@model function size_structured_dde_turing_testing( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
 
  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~  arraydist(Gamma.(2.0, Scv) ) # #  working: β(0.1, 10.0);  plot(x->pdf(Gamma(2.0, .1), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  d ~   filldist( LogNormal( -1.508 , 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -0.7096, 0.5 ), nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  v ~   filldist( LogNormal( 0.3782, 0.5 ),  4 ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  

  u0 ~  filldist( Beta(1.0, 1.0), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v  )
 
  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      # abstol=1.0e-10, 
      # reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[Si[i],k] )  # observation and process error combined
      end
  end

end
 
# ------------------------------


@model function size_structured_dde_turing_north( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
 
  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~  arraydist(Gamma.(2.0, Scv) ) # #  working: β(0.1, 10.0);  plot(x->pdf(Gamma(2.0, .1), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  d ~   filldist( LogNormal( -1.508 , 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -0.7096, 0.5 ), nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  v ~   filldist( LogNormal( 0.3782, 0.5 ),  4 ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  

  u0 ~  filldist( β(0.8, 8.0), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v  )
  
  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      # abstol=1.0e-10, 
      # reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[Si[i],k] )  # observation and process error combined
      end
  end

end
 
 
# ------------------------------

@model function size_structured_dde_turing_south( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
 

  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~  arraydist(Gamma.(2.0, Scv) ) # #  working: β(0.1, 10.0);  plot(x->pdf(Gamma(2.0, .1), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  d ~   filldist( LogNormal( -1.050 , 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -0.7096  , 0.5 ), nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  v ~   filldist( LogNormal( 0.3782, 0.5 ),  4 ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  

  u0 ~  filldist( β(0.8, 8.0), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v  )
 
  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      # abstol=1.0e-10, 
      # reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[Si[i],k] )  # observation and process error combined
      end
  end

end
 
 
# ------------------------------


@model function size_structured_dde_turing_4x( S, kmu, tspan, prob, nS, 
  solver=MethodOfSteps(Tsit5()), dt = 0.01)
  
  # plot(x->pdf(LogNormal(log(kmu),0.2), x), xlim=(0,kmu*5))
  K ~ filldist( LogNormal(log(kmu), 0.25), nS )  # kmu is max of a multiyear group , serves as upper bound for all
  q ~ filldist( Normal( 1.0, 0.1 ), nS )
  qc ~ arraydist([Normal( -SminFraction[i], 0.1) for i in 1:nS])  # informative prior on relative height 
 
  model_sd ~  arraydist(Gamma.(2.0, Scv) ) # #  working: β(0.1, 10.0);  plot(x->pdf(Gamma(2.0, .1), x), xlim=(0,1)) # uniform 

  # lognormal (1,1) has a mode at 1, with a large variability 
  b ~   filldist( LogNormal(  1.0, 1.0 ),  2 )   # centered on 1; plot(x->pdf(LogNormal(1.0, 1.0), x), xlim=(0,10)) # mode of 5

  #= note: 
    log( exp(0.1)-1.0 ) = log(0.1052 ) = -2.252  .. 10% mortality (mode)
    log( exp(0.2)-1.0 ) = log(0.2214 ) = -1.508  .. 20% mortality (mode)
    log( exp(0.3)-1.0 ) = log(0.3499 ) = -1.050  .. 30% 
    log( exp(0.4)-1.0 ) = log(0.4918 ) = -0.7096 .. 40%
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  d ~   filldist( LogNormal( -1.508 , 0.25 ), nS ) # plot(x->pdf(LogNormal(0.2, 1.0), x), xlim=(0, 2)) 
  d2 ~  filldist( LogNormal( -0.7096, 0.5 ), nS ) # plot(x->pdf(LogNormal(-0.7096, 0.25 ), x), xlim=(0, 2)) 

  #= note: 
    log( exp(0.90)-1.0) = log(1.46)  = 0.3782  .. ~90% (moult) transition rate per year (mode)
    log( exp(0.95)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    log( exp(0.99)-1.0) = log(1.586) = 0.461   .. ~95% (moult) transition rate per year (mode) 
    # Can go higher than 100% .. as the size-based categories are imperfect and there is also input from other groups
  =#
  v ~   filldist( LogNormal( 0.3782, 0.5 ),  4 ) # transition rates # plot(x->pdf(LogNormal( 0.3782, 0.5 ), x), xlim=(0,1))  

  u0 ~  filldist( β(0.2, 12.0), nS )  # plot(x->pdf(Beta(1, 1), x), xlim=(0,1)) # uniform 

  # pm = ( b, K, d, d2, v, tau, hsa )
  pm = ( b, K, d, d2, v  )
 
  # process model
  
  prob_new = remake( prob; u0=u0, h=h, tspan=tspan, p=pm )

  msol = solve(
      prob_new,
      solver, 
      callback=cb,
      isoutofdomain=(y,p,t)->any(x -> x<0.0, y),  # permit exceeding K
      # abstol=1.0e-10, 
      # reltol=1.0e-10, 
    #  dt=dt,  
      saveat=dt
  )

  # @show msol.retcode
  if msol.retcode != :Success
    Turing.@addlogprob! -Inf
    return nothing
  end

  # likelihood of the data
  for i in 1:nSI
      ii = findall(x->x==survey_time[Si[i]], msol.t)[1]
      for k in 1:nS
          S[Si[i],k] ~ Normal( msol.u[ii][k] * q[k] + qc[k], model_sd[Si[i],k] )  # observation and process error combined
      end
  end

end
 
# ------------------------------


function fishery_model_test( test=("basic", "random_external_forcing", "fishing", "nofishing") )

  ## test DifferentialEquations DDE model 
  gr()
  theme(:default)
 
  if any( occursin.( r"basic", test )  )

    ##  h, hsa, cb, tau, etc. are defined in the *_environment.jl file
    b=[0.8, 0.5]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* kmu;
    v=[0.65, 0.68, 0.61, 0.79];
    d=[0.15, 0.11, 0.14, 0.17, 0.16, 0.19];
    d2 =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
    
    u0 = [ 0.65, 0.6, 0.52, 0.62, 0.58, 0.32 ]  ; 
    tau=[1.0] 
    # params = ( b, K, d, d2, v, tau, hsa )
    params = ( b, K, d, d2, v  )
 
    prob1 = DDEProblem( size_structured_dde!, u0, h, tspan, params, constant_lags=tau  )  # tau=[1]
    out = msol1 =  solve( prob1,  solver, callback=cb, saveat=dt )
    pl = plot()
    pl = plot( pl, msol1, ; legend=false, xlim=(1999,2021), title="basic" )

  elseif  any( occursin.( r"random_external_forcing", test )  )

    # testing DDE version -- initial values for size_structured_dde! 
    # basic run
    ks = 1000
    u0 = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 ] ;
    b=[ 0.6, 0.6 ]
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    d=[0.4, 0.4, 0.4, 0.4, 0.4, 0.4];
    v=[0.9, 0.9, 0.9, 0.9];  
    d2=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5] 
    tau=[1.0] 
  
    survey_time = 1999:2021
    
    external_forcing = ones(length(survey_time),6)  # turns it off
    # external_forcing = rand(length(survey_time),6) # random
    
    efc = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc, 1999:2021, 1:6 )
    
    # p = ( b, K, d, d2, v, tau, hsa )   
    p = ( b, K, d, d2, v )   
    
    tspan = (1990.0, 2050.0)
    nS = length(u0)  # n components
    
    # history function for time before t0:: 0.5 default
    # h(p,t) = ones( nS ) .* 0.5  #values of u before t0
    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS)  .*0.5
    tau = [1.0]  # delay
    
    solver = MethodOfSteps(Tsit5())  # solver; BS3() and Vern6() also RK4()
    prob2 = DDEProblem( size_structured_dde! , u0, h, tspan, p; constant_lags=tau )
    out = msol2 =  solve( prob2,  solver, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol2 ; title="basic" )
  
  elseif  any( occursin.( r"nofishing", test ) )

    ks = 1.0e9
    K=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    b=[0.6, 0.4]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2=[0.6, 0.6, 0.6, 0.6, 0.6, 0.6] 
    tau=[1.0] 

    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* 0.5
   
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    external_forcing = ones(length(survey_time),6)  # turn off external forcing  (wrt viable habitat)
    efc1 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc1, 1999:2021, 1:6 )
    
    # p = ( b, K, d, d2, v, tau, hsa )   
    p = ( b, K, d, d2, v )   
    
    prob3 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    out = msol3 =  solve( prob3,  solver, saveat=dt  ) #, isoutofdomain=(y,p,t)->any(x->(x<0)|(x>1), y) )
    pl = plot()

    pl = plot!( pl, msol3, title="dde, with random hsa, with fishing") 

    i = 1; # fb
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_1" ) 

    i = 6; # fem
    pl = plot!( pl, msol3.t, reduce(hcat, msol3.u)'[:,i], color=[1 1] , alpha=0.75, lw=5, title="fishing_6" ) 
  
   
  elseif  any( occursin.( r"fishing", test )  | occursin( r"nofishing", test ) )

    ks = 1.0e9
    K = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] .* ks;
    u0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
    b=[0.8, 0.5]
    v=[0.8066539263878958, 0.7165358852484025, 0.8341124383106499, 0.7857601054678678]
    d=[0.20223639702656727, 0.18864211980978104, 0.18063928527606177, 0.23030220100440996, 0.19713968752681676, 0.40151610614035915]
    d2=[0.4, 0.4, 0.4, 0.4, 0.4, 0.4] 
    tau=[1.0] 

    h(p, t; idxs=nothing) = typeof(idxs) <: Number ? 1.0 : ones(nS) .* ks
   
    tspan = (1990.0, 2025.0)
  
    dt = 0.001
    nS = length(u0)  # n components
   
    # add random viable habitat forcing
    external_forcing = rand(length(survey_time),6) # random
    efc2 = extrapolate( interpolate( external_forcing, (BSpline(Linear()), NoInterp()) ), Interpolations.Flat() )
    hsa = Interpolations.scale(efc2, 1999:2021, 1:6 )
   
    #p = ( b, K, d, d2, v, tau, hsa)   
    p = ( b, K, d, d2, v )   
   
    prob4 = DDEProblem( size_structured_dde!, u0, h, tspan, p; constant_lags=tau )
    # prob4 = remake( prob; u0=u0, h=h, tspan=tspan, p=p )
      
    out = msol4 =  solve( prob4,  solver, saveat=dt )
    pl = plot()

    pl = plot!( pl, msol4, title="dde, with random hsa, no fishing" )

    i = 1; 
    pl = plot!( pl, msol4.t, reduce(hcat, msol4.u)'[:,i], color=[3 3] , alpha=0.75, lw=5 ) 
    
    i = 6; 
    pl = plot!( pl, msol4.t, reduce(hcat, msol4.u)'[:,i], color=[1 1] , alpha=0.75, lw=5 ) 
    
  end

  pl = plot!( pl; legend=true, xlim=(1990,2021) )

  return (msol, pl) 

end
 
 

# -----------

function abundance_from_index( aindex, res, k )
  u = size(res)
  v = aindex .* ones(length(aindex), u[1]) # expand to matrix
  yhat = ( v' .- res[:,Symbol("qc[$k]"),:]  )  ./  res[:,Symbol("q[$k]"),:] .* res[:,Symbol("K[$k]"),:]  
  return yhat'
end


# -----------

  


function fishery_model_predictions_timeseries( num; prediction_time, plot_k )
  gk = num[:,plot_k,:,1]
  pl = plot()
  pl = plot!(pl, prediction_time, gk;  alpha=0.02, color=:lightslateblue)
  pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
  pl = plot!(pl; legend=false )
  pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

  # back transform S to normal scale .. do sims too (TODO)
  yhat = ( S[:,plot_k] .- mean(res[:,Symbol("qc[$plot_k]"),:]  ) ) ./ mean(res[:,Symbol("q[$plot_k]"),:]) .* mean(res[:,Symbol("K[$plot_k]"),:]  )

  pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
  pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
  pl = plot!(pl; legend=false )
  return (gk, pl)
end

# -----------





function fishing_mortality_instantaneous( removed, abundance )
  -log(  1.0 - (removed  / abundance)  )  ;
end


function removals_aggregate( removed, fish_year )
  landings_aggregated = DataFrame( yr=floor.(fish_year), rem = removed );
  out = combine(groupby(landings_aggregated,:yr),[:rem ] .=> sum )
  stimes = DataFrame( yr=floor.(survey_time) )
  out = leftjoin(stimes, out, on=:yr)
  sort!(out, :yr)
  oo = findall(x->ismissing(x), out[:,:rem_sum])
  if length(oo) > 0
    out[ oo[1], :rem_sum ] = 0.0
  end
  return(out)
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
  # vn = "b[1]"; density!(res[ Symbol(vn) ])
  # vn = "b[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[1]"; density!(res[ Symbol(vn) ])
  # vn = "d[2]"; density!(res[ Symbol(vn) ])
  # vn = "d[3]"; density!(res[ Symbol(vn) ])
  # vn = "d[4]"; density!(res[ Symbol(vn) ])
  # vn = "d[5]"; density!(res[ Symbol(vn) ])
  # vn = "d[6]"; density!(res[ Symbol(vn) ])
  # vn = "K[1]"; density!(res[ Symbol(vn) ])
  # vn = "K[2]"; density!(res[ Symbol(vn) ])
  # vn = "K[3]"; density!(res[ Symbol(vn) ])
  # vn = "K[4]"; density!(res[ Symbol(vn) ])
  # vn = "K[5]"; density!(res[ Symbol(vn) ])
  # vn = "K[6]"; density!(res[ Symbol(vn) ])
end


# ----------


function fishery_model_inference( fmod; rejection_rate=0.65, 
  n_adapts=1000, n_samples=1000, n_chains=1, max_depth=7,   
  turing_sampler=Turing.NUTS(n_adapts, rejection_rate; max_depth=max_depth ) 
)

  res  =  sample( fmod, turing_sampler, MCMCThreads(), 
    n_samples, n_chains #, thinning=thinning, discard_initial=discard_initial  
  )
   
  # other call options /approaches

  # res_means = FillArrays.Fill(summarize(res).nt[2], n_chains)
  # res  =  sample( fmod, Turing.NUTS(), 10, init_params=res_means )   # , thinning=5, discard_initial=30  )
  # res  =  sample( fmod, turing_sampler_test, n_sample_test, thinning=5, discard_initial=n_adapts_test ) # to see progress -- about 5 min
  # res = sample( fmod, turing_sampler_test, MCMCThreads(), n_adapts_test, n_chains_test, n_adapts_test ) # MCMCThreads(), n_sample_test, n_chains, n_adapts
  # res = fishery_model_inference( fmod, n_adapts=n_adapts_test, n_samples=n_sample_test, n_chains=n_chains_test, max_depth=7  ) # same thing

  # if on windows and threads are not working, use single processor mode:
  # res = mapreduce(c -> sample(fmod, turing_sampler, n_samples), chainscat, 1:n_chains)

  showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

  return res
end

# this repetition is not an error: this is Julia's type inference (in this case addition of "init_params")
function fishery_model_inference( fmod; rejection_rate=0.65, 
  n_adapts=1000, n_samples=1000, n_chains=1, max_depth=7,
  turing_sampler=Turing.NUTS(n_adapts, rejection_rate; max_depth=max_depth ),    
  init_params=NaN 
)
    res  =  sample( fmod, turing_sampler, MCMCThreads(), 
      n_samples, n_chains, init_params=init_params #, thinning=thinning, discard_initial=discard_initial  
    )
   showall(summarize(res ) )  # show(stdout, "text/plain", summarize(res)) # display all estimates

  return res
end




# -------------------


function fishery_model_predictions( res; prediction_time=prediction_time, n_sample=100 )

  nchains = size(res)[3]
  nsims = size(res)[1]
  
  md = zeros(nM, nS, n_sample, 2)  # number normalized
  mn = zeros(nM, nS, n_sample, 2)  # numbers
  mb = mn[:,1,:,:]  # biomass of first class

  trace_time = Vector{Vector{Float64}}()

  out10 = Vector{Vector{Float64}}()
  out11 = Vector{Vector{Float64}}()
  out12 = Vector{Vector{Float64}}()
  
  out1 = Vector{Vector{Float64}}()
  out2 = Vector{Vector{Float64}}()
  out3 = Vector{Vector{Float64}}()
  out4 = Vector{Vector{Float64}}()
  out5 = Vector{Vector{Float64}}()
  out6 = Vector{Vector{Float64}}()

  ntries = 0
  z = 0

  while z <= n_sample 
    ntries += 1
    ntries > n_sample*10 && break
    z >= n_sample && break

    j = rand(1:nsims)  # nsims
    l = rand(1:nchains) #nchains
    b = [ res[j, Symbol("b[$k]"), l] for k in 1:2]
    K = [ res[j, Symbol("K[$k]"), l] for k in 1:nS]
    v = [ res[j, Symbol("v[$k]"), l] for k in 1:4]
    d = [ res[j, Symbol("d[$k]"), l] for k in 1:nS]
    d2=[ res[j, Symbol("d2[$k]"), l] for k in 1:nS]
    u0 = [ res[j, Symbol("u0[$k]"), l] for k in 1:nS]

    # pm = ( b, K, d, d2, v, tau, hsa )
    pm = ( b, K, d, d2, v )

    prb = remake( prob; u0=u0 , h=h, tspan=tspan, p=pm )
    msol1 = solve( prb, solver, callback=cb, saveat=dt, dt=dt  )

    z==0 && push!(trace_time, msol1.t)

    msol0 = solve( prb, solver, saveat=trace_time[1] ) # no call backs
    
    if msol1.retcode == :Success && msol0.retcode == :Success
        z += 1
        
        # annual
        for i in 1:nM
            ii = findall(x->x==prediction_time[i], trace_time[1]) 
            if length(ii) > 0  
              ii = ii[1]
              sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(msol1.t[ii])  ./ 1000.0 ./ 1000.0  :  scale_factor   # n to kt
              md[i,:,z,1] = msol1.u[ii]  # with fishing
              md[i,:,z,2] = msol0.u[ii]  # no fishing
              mn[i,:,z,1] = msol1.u[ii]  .* K # with fishing scaled to K
              mn[i,:,z,2] = msol0.u[ii]  .* K # no fishing scaled to K
              mb[i,z,1] = mn[i,1,z,1]  .* sf  # biomass of state var 1 
              mb[i,z,2] = mn[i,1,z,2]  .* sf
            end
        end  # end for

        # traces for plotting, etc 
        sf  = nameof(typeof(mw)) == :ScaledInterpolation ? mw(trace_time[1])  ./ 1000.0 ./ 1000.0 :  scale_factor
        b10 = vec( reduce(hcat, msol0.u)'[:,1]) .* K[1] .* sf
        b11 = vec( reduce(hcat, msol1.u)'[:,1]) .* K[1] .* sf
        b12 = ( b10 .- b11 ) ./ b10 

        push!(out10, b10)
        push!(out11, b11)
        push!(out12, b12)

        push!(out1, vec( reduce(hcat, msol1.u)'[:,1]) .* K[1])
        push!(out2, vec( reduce(hcat, msol1.u)'[:,2]) .* K[2])
        push!(out3, vec( reduce(hcat, msol1.u)'[:,3]) .* K[3])
        push!(out4, vec( reduce(hcat, msol1.u)'[:,4]) .* K[4])
        push!(out5, vec( reduce(hcat, msol1.u)'[:,5]) .* K[5])
        push!(out6, vec( reduce(hcat, msol1.u)'[:,6]) .* K[6])

      end # if
  end  # while
 
  if z < n_sample 
    @warn  "Insufficient number of solutions" 
  end

  trace = (out1, out2, out3, out4, out5, out6 )
  trace_bio = (out10, out11, out12)
  return (md, mn, mb, trace, trace_bio, trace_time[1] )

end



# -----------


function fishery_model_mortality( ; removed=removed, bio=bio, survey_time=survey_time, fish_year=fish_year )   
  fb = bio[1:length(survey_time),:,1]  # the last 1 is for size struct; no effect in discrete 
  removed_annual_kt = removals_aggregate( removed, fish_year )
  Fkt = removed_annual_kt[:,:rem_sum] ./1000.0 ./ 1000.0  # removal in kg -> kt
  FR =  Fkt ./ ( Fkt .+  fb )  # relative F
  FM = -1 .* log.(  1.0 .- min.( FR, 0.99) )  # instantaneous F
  # FM[ FM .< eps(0.0)] .= zero(eltype(FM))
  return ( Fkt, FR, FM  )
end

# -----------


function fishery_model_plot(; toplot=("fishing", "nofishing", "survey"), n_sample=200,
  res=res, bio=bio, num=num, trace=trace, trace_bio=trace_bio, FM=FM, 
  S=S, si=1, scale_factor=scale_factor, 
  prediction_time=prediction_time, survey_time=survey_time, trace_time=trace_time, yrs=yrs, 
  alphav=0.05, pl= Plots.plot(), time_range=(floor(minimum(survey_time))-1.0, ceil(maximum(survey_time))+1.0 )
)
  
  nsims = size(bio)[2]
  ss = rand(1:nsims, n_sample)  # sample index

  # extract sims (with fishing)
  # plot biomass
  if any(isequal.("fishing", toplot))  
    g = bio[:,ss,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:orange)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkorange, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("nofishing", toplot))  
    g = bio[:,ss,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:lime)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:limegreen, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
  end

  if any(isequal.("footprint", toplot))  
    g1 = bio[:,ss,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,ss,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    pl = plot!(pl, prediction_time, g ;  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(g, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(g)*1.01 ) )
    pl = plot!(pl; xlim=time_range )
    
  end

  if any(isequal.("survey", toplot))  
    # back transform S to normal scale .. do sims too (TODO)
    k = 1
    qcm = mean(res[:,Symbol("qc[$k]"),:])
    qm  = mean(res[:,Symbol("q[$k]"),:])
    Km = mean(res[:,Symbol("K[$k]"),:]  )
    yhat = ( S[:,k] .- qcm  ) ./ qm .* Km   # abundance_from_index S[:,1]  
    if nameof(typeof(mw)) == :ScaledInterpolation
      yhat = yhat .* mw(yrs) ./ 1000.0  ./ 1000.0
    else
      yhat = yhat .* scale_factor
    end
    pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, yhat, markersize=4, color=:darkgray)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )

  end
  

  if any(isequal.("number", toplot))  

    gk = num[:,si,ss,1]
    pl = plot!(pl, prediction_time, gk;  alpha=alphav, color=:lightslateblue)
    pl = plot!(pl, prediction_time, mean(gk, dims=2);  alpha=0.8, color=:darkslateblue, lw=4)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; ylim=(0, maximum(gk)*1.01 ) )

    # back transform S to normal scale .. do sims too (TODO)
    yhat = ( S[:,si] .- mean(res[:,Symbol("qc[$si]"),:]  ) ) ./ mean(res[:,Symbol("q[$si]"),:]) .* mean(res[:,Symbol("K[$si]"),:]  )

    pl = plot!(pl, survey_time, yhat, color=:gray, lw=2 )
    pl = scatter!(pl, survey_time, yhat, markersize=4, color=:grey)
    pl = plot!(pl; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("trace", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2][ss], alpha=alphav, lw=1, color=:orange )
    pl = plot!( pl, trace_time, trace_bio[1][ss], alpha=alphav, lw=1, color=:lime )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end

  if any(isequal.("trace_fishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[2][ss], alpha=alphav, lw=1, color=:orange )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end
 
  if any(isequal.("trace_nofishing", toplot))  
    pl = plot!( pl, trace_time, trace_bio[1][ss], alpha=alphav, lw=1, color=:lime )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end


  if any(isequal.("trace_footprint", toplot))  
    pl = plot!( pl, trace_time, trace_bio[3][ss], alpha=alphav, lw=1, color=:lightslateblue )
    pl =  plot!(pl; legend=false )
    pl =  plot!(pl; xlim=time_range )
  end

  if any(isequal.("fishing_mortality", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    pl = plot!(pl, survey_time, FM[:,ss] ;  alpha=0.02, color=:lightslateblue)
    pl = plot!(pl, survey_time, FMmean ;  alpha=0.8, color=:slateblue, lw=4)
    pl = plot!(pl, ylim=(0, ub ) )
    pl = plot!(pl ; legend=false )
    pl = plot!(pl; xlim=time_range )
  end


  if any(isequal.("fishing_mortality_vs_footprint", toplot))  
    FMmean = mean( FM, dims=2)
    FMmean[isnan.(FMmean)] .= zero(eltype(FM))
    ub = maximum(FMmean) * 1.1
    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
    g = g[1:length(survey_time),:]
    pl = scatter!(pl, FM[:,ss], g[:,ss];  alpha=alphav, color=:lightslateblue)
    pl = scatter!(pl, FMmean, mean(g, dims=2);  
      alpha=0.8, color=:darkslateblue, lw=4, markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4))
    pl = plot!(pl ; legend=false )
  end


  if any(isequal.("harvest_control_rule_footprint", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
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
  
    nt = length(survey_time)
    colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
  
    fb_mean = mean(fb, dims=2)

    g1 = bio[:,:,1]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g2 = bio[:,:,2]   # [ yr,  sim, (with fishing=1; nofishing=2) ]
    g = ( g2 - g1 ) ./ g2
  
    g = g[1:length(survey_time),:]
  
    g_mean = mean(g, dims=2)
  
    # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
    pl = plot!(pl, fb_mean, g_mean ;  alpha=0.8, color=:slateblue, lw=3)
  
    pl = scatter!(pl,  fb_mean, g_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
    pl = scatter!(pl,  [fb_mean[nt]], [g_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
    
    ub = max( quantile(K, 0.95), maximum( fb_mean ) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(g_mean ) * 1.05  ) )
 
  end
   

  if any(isequal.("harvest_control_rule", toplot))  
    fb = bio[1:length(survey_time),:,1] 
 
    # mean weight by year
    sf = nameof(typeof(mw)) == :ScaledInterpolation ?  mw(yrs) ./ 1000.0  ./ 1000.0 : scale_factor
  
    # sample and plot posterior K
    K = vec( Array(res[:, Symbol("K[1]"), :]) ) .* mean(sf)  # convert to biomass 
  
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
  
    nt = length(survey_time)
    colours = get(colorschemes[:tab20c], 1:nt, :extrema )[rand(1:nt, nt)]
  
    # scatter!( fb, FM ;  alpha=0.3, color=colours, markersize=4, markerstrokewidth=0)
  
    fb_mean = mean(fb, dims=2)
    fm_mean = mean(FM, dims=2)
  
    # scatter!( [fb[nt,:]], [FM[nt,:]] ;  alpha=0.3, color=:yellow, markersize=6, markerstrokewidth=0)
    pl = plot!(pl, fb_mean, fm_mean ;  alpha=0.8, color=:slateblue, lw=3)
  
    pl = scatter!(pl,  fb_mean, fm_mean ;  alpha=0.8, color=colours,  markersize=4, markerstrokewidth=0,
      series_annotations = text.(trunc.(Int, survey_time), :top, :left, pointsize=4) )
    pl = scatter!(pl,  [fb_mean[nt]], [fm_mean[nt]] ;  alpha=0.8, color=:yellow, markersize=8, markerstrokewidth=1)
    
    ub = max( quantile(K, 0.95), maximum( fb_mean ) ) * 1.05
    pl = plot!(pl; legend=false, xlim=(0, ub ), ylim=(0, maximum(fm_mean ) * 1.05  ) )
 
  end
   
  return(pl)

end



