carstm_estimation_runtypes = function( snowcrab_filter_class=NULL, subtype=NULL ) {

  runtypes = list()
 

    
  runtypes[["R0"]] = list( 
    label= "default_R0",  
    thetaN=c( -0.798, 4.762, 0.620, 0.194, -0.927, 4.137, 2.724, 0.601, 1.521, 0.590, -0.219, 0.735  ) ,
    thetaW=c( 5.059, 8.307, 1.268, 12.584, 12.290, 11.263, 16.387, 13.426, 8.426, 0.735, 5.515, 4.292, 2.332 ) ,
    thetaH= c( -0.790, 3.663, 5.293, 1.211, -1.382, 3.752, 4.514, -1.619, 1.632, -0.959, 2.544, 2.728 ) 
  )


  runtypes[["fb"]] = list(  
    label= "default_fb",  
    thetaH=c( 0.187, 1.603, 2.131, 1.059, -0.015, 4.235, 4.792, -1.670, 2.050, -0.914, 3.126, 2.419 )
  )
  
 
  runtypes[["recruits"]] = list( 
    label= "default_male_recruits",  
    thetaH=c( 0.187, 1.603, 2.131, 1.059, -0.015, 4.235, 4.792, -1.670, 2.050, -0.914, 3.126, 2.419 )
  )

  # crab from instar 6 (~lower limit of detection of 20mm) 
  # to before they begin to mature sexually (instar 8) 
  runtypes[["pre.recruits.i6_i8"]] = list( 
    label= "default_MF_pre.recruits.i6_i8",  
    thetaH=c( 0.391, 1.382, 3.975, 1.217, -0.270, 3.304, 5.373, -1.659, 2.000, -0.963, 3.381, 2.453  )
  )

  runtypes[["imm"]] = list( 
    label= "default_imm",  
    thetaH=c( 1.070, 1.125, 2.441, 2.433, -1.111, 3.215, 4.590, -1.624, 2.283, -0.353, 0.089, 2.646 )
  )
 

  runtypes[["m.adolescent"]] = list( 
    label= "default_male_adolescent",  
    thetaH=c( 0.391, 1.382, 3.975, 1.217, -0.270, 3.304, 5.373, -1.659, 2.000, -0.963, 3.381, 2.453  )
  )
  
  runtypes[["m.mat"]] = list( 
    label= "default_male_mature",  
    thetaH=c( 0.187, 1.603, 2.131, 1.059, -0.015, 4.235, 4.792, -1.670, 2.050, -0.914, 3.126, 2.419  )
  )

  runtypes[["f.mat"]] = list( 
    label= "default_female_mature",  
    thetaH=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )

  runtypes[["f.adolescent"]] = list( 
    label= "default_female_adolescent",  
    thetaH=c( 1.038, 0.904, 1.795, 1.545, -1.054, 3.598, 3.079, -0.998, 1.770, -0.778, 2.761, 2.096  )
  )

  runtypes[["primiparous"]] = list( 
    label= "default_female_primiparous",  
    thetaH=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )

  runtypes[["multiparous"]] = list( 
    label= "default_female_multiparous",  
    thetaH=c( 1.822, 3.625, -0.537, 0.218, 2.792, 4.227, 4.881, -0.184, 3.138, -0.082, 2.903, 1.623 )
  )

  if ( !is.null(snowcrab_filter_class) ) {
    if (!exists( snowcrab_filter_class, runtypes )) {
      runtypes[[snowcrab_filter_class]] = list(label=snowcrab_filter_class)
    } 
    runtypes = runtypes[[snowcrab_filter_class]]
    if (!exists("label", runtypes)) runtypes$label = snowcrab_filter_class
    if (!is.null(subtype)) {
      runtypes$label = paste(runtypes$label, subtype, sep="_")
      if (grepl( "negative_binomial", subtype )) {
        runtypes$family$H = "nbinomial"
      }
    } 
  }
    

  return (runtypes)
  
}
