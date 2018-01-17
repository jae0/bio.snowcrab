
  interpolation.db = function( ip=NULL, DS=NULL, p=NULL, 
    varnames = c("snowcrab.large.males_abundance", "snowcrab.large.males_presence_absence", annot.cex=2) ) {

    if (DS %in% c( "biomass", "biomass.redo" )) {
    
      outdir = file.path( project.datadirectory("bio.snowcrab"), "modelled", "biomass" )
      dir.create(path=outdir, recursive=T, showWarnings=F)
      fn = file.path( outdir, "biomass.rdata" )

      if (DS=="biomass") {
        B = NULL
        if ( file.exists( fn) ) load(fn)
        return (B)
      }

      
      set = aegis::survey.db( p=p, DS="set.filter" ) # mature male > 95 mm 
 
      ii = which( set$totmass > 0 )
      qs = quantile( set$totmass[ii], probs=p$stmv_quantile_bounds, na.rm=TRUE )
      qs = log(qs) # 

      bm = snowcrab_stmv( p=p, DS="baseline", ret="mean", varnames=varnames )
      m = bm[[1]]  # biomass
      h = bm[[2]]  # habitat
      bm= NULL
      
      ll = which(h < p$habitat.threshold.quantile )
      if (length(ll) > 0 ) h[ll] = 0
      m = log( exp(m) * h ) # bm[[2]] is serving as weight/probabilities
      # #respect the bounds of input data (no extrapolation)
      
      qq = which( m < qs[1] )
      if (length(qq) > 0 ) m[qq] = NA
      rr = which( m > qs[2] )
      if (length(rr) > 0 ) m[rr] = qs[2]

      if(0) {
        bs = bathymetry.db(p=p, DS="baseline")
        levelplot( m[,16] ~ plon+plat, bs, aspect="iso")
        for (i in 1:16) print(levelplot( m[,i] ~ plon+plat, bs, aspect="iso"))
      }

      # more range checks
      s = snowcrab_stmv( p=p, DS="baseline", ret="sd", varnames=varnames )
      # range checks
      s = log( exp(s[[1]]) * s[[2]] ) # s[[2]] is serving as weight/probabilities
      
      sq = quantile(s, probs=p$stmv_quantile_bounds[2], na.rm=TRUE ) 
      # s[which(s > sq)] = sq  # cap upper bound of sd

      # mm = which(( m - 1.96*s ) < 0 )
      # if (length(mm)>0) {
      #   m[mm] = NA
      #   s[mm] = NA 
      # }

      # limit range of extrapolation to within a given distance from survey stations .. annual basis
      set = snowcrab.db( DS="set.clean")
      bs = bathymetry.db(p=p, DS="baseline")
      bb = array_map( "xy->1", bs, gridparams=p$gridparams )

      if (0) {
        # annual mask 
        for (iy in 1:p$ny ) {
          S = set[ which(set$yr==p$yrs[iy]), c("plon", "plat") ] 
          S = S[ !duplicated(S),]
          nn = array_map( "xy->1", S, gridparams=p$gridparams )
          overlap = match(  nn, bb )
          overlap = overlap[ which( is.finite( overlap ))] 
          o = bs[overlap,] 
          # add corners as buffer
          ot = t(o) # reshape to make addition simpler using R's cycling rules
          o = rbind( t( ot + p$threshold.distance*c( 1, 1) ), 
                     t( ot + p$threshold.distance*c(-1,-1) ),
                     t( ot + p$threshold.distance*c( 1,-1) ),
                     t( ot + p$threshold.distance*c(-1, 1) )
          )
          o = o[ !duplicated(o),]
          boundary= non_convex_hull( o, alpha=p$threshold.distance*4, plot=FALSE )
          outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
          m[outside.polygon,iy] = NA
          s[outside.polygon,iy] = NA
          h[outside.polygon,iy] = NA
        }
      }

      # a static mask
      S = set[ , c("plon", "plat") ] 
      S = S[ !duplicated(S),]
      nn = array_map( "xy->1", S, gridparams=p$gridparams )
      overlap = match(  nn, bb )
      overlap = overlap[ which( is.finite( overlap ))] 
      o = bs[overlap,] 
      # add corners as buffer
      ot = t(o) # reshape to make addition simpler using R's cycling rules
      o = rbind( t( ot + p$threshold.distance*c( 1, 1) ), 
                 t( ot + p$threshold.distance*c(-1,-1) ),
                 t( ot + p$threshold.distance*c( 1,-1) ),
                 t( ot + p$threshold.distance*c(-1, 1) )
      )
      o = o[ !duplicated(o),]
      boundary= non_convex_hull( o, alpha=25, plot=FALSE )
      outside.polygon = which( point.in.polygon( bs[,1], bs[,2], boundary[,1], boundary[,2] ) == 0 )
      m[outside.polygon,] = NA
      s[outside.polygon,] = NA
      h[outside.polygon,] = NA

      B = list( m=m, s=s, h=h )
      save( B, file=fn, compress=TRUE )
      B = NULL

      projectdir = file.path(p$data_root, "maps", "fishable.biomass", p$spatial.domain )
      datarange = seq( qs[1]-0.1, qs[2]+0.1, length.out=150)
      cols = color.code( "seis", datarange )
      m[which(!is.finite(m))] = qs[1]
      for (iy in 1:p$ny) {
        y = p$yrs[iy]
        outfn = paste( "prediction.abundance.mean", y, sep=".")
        xyz = cbind( bs[, c("plon", "plat")], m[,iy] )
        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=T, depthcontours=T, pts=NULL, annot=y,
          annot.cex=annot.cex, corners=p$planar.corners, at=datarange,
          col.regions=cols, rez=c(p$pres,p$pres) )
        print(lp)
        dev.off()
      }

      datarange = seq( 0, max(s,na.rm=TRUE), length.out=150)
      cols = color.code( "seis", datarange )
      s[which(!is.finite(s))] = 0.1
      s[which(s < 0.1)] = 0.1
      
      for (iy in 1:p$ny) {
        y = p$yrs[iy]
        outfn = paste( "prediction.abundance.sd", y, sep=".")
        xyz = cbind( bs[, c("plon", "plat")], s[,iy] )
        dir.create (projectdir, showWarnings=FALSE, recursive =TRUE)
        fn = file.path( projectdir, paste(outfn, "png", sep="." ) )
        png( filename=fn, width=3072, height=2304, pointsize=40, res=300 )
        lp = aegis::aegis_map( xyz=xyz, cfa.regions=T, depthcontours=T, pts=NULL, annot=y,
          annot.cex=annot.cex, corners=p$planar.corners, at=datarange,
          col.regions=cols, rez=c(p$pres,p$pres) )
        print(lp)
        dev.off()
      }

      return(fn)
    }

  
    # ------------------


    if (DS=="timeseries") {
      bm = interpolation.db(p=p, DS="biomass")
      bm$m = exp(bm$m) / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
      bm$s = exp(bm$s) / 10^3  # kg/km^2 to t/km^2  .. required for biomass.summary.db
      
      bs = bathymetry.db( p=p, DS="baseline")
      
      bq = min( bm$m[ bm$m > 0 ], na.rm=T )

      K = NULL
      nreg = length(p$regions)
      for (r in 1:nreg ){
        aoi = aegis::polygon_inside(x=bs[ , c("plon", "plat")], region=p$regions[r], planar=T)
        aoi = intersect( aoi, which( bs$plon > 250 ) )
        out = matrix( NA, nrow=p$ny, ncol=3) 
        
        for (y in 1:p$ny) {
          iHabitat = which( bm$h[,y] > p$habitat.threshold.quantile ) # any area with biomass > lowest threshold
          iHabitatRegion = intersect( aoi, iHabitat )
          out[ y, 1] = sum( bm$m[iHabitatRegion,y] , na.rm=TRUE ) # abundance weighted by Pr
          out[ y, 2] = sqrt( sum( (bm$s[iHabitatRegion,y])^2 , na.rm=TRUE ) )
          out[ y, 3] = sum( bm$h[iHabitatRegion,y] ) * (p$pres*p$pres)
        }
        
        ok = as.data.frame( out )
        names( ok) = c("total", "total.sd", "sa.region")
        
        ok$log.total = log(ok$total)
        ok$total.sd.ln = log(ok$total.sd) # as above
        ok$region = p$regions[r]
        ok$yr = p$yrs

        ok$lbound = ok$total - ok$total.sd*1.96 # normal assumption 
        ok$ubound = ok$total + ok$total.sd*1.96 
        K = rbind(K, ok)
      }

      return( K )      
    }


    if ( DS %in% c( "interpolation.simulation" ) ) {
      message(" simulation-based results are not ready at present")
      message(" defaulting to simple estimates based upon assymptotic assumptions" )
      message(" only R0.mass is supported for now" )
      
      out = NULL  
      if ( p$vars.to.model == "R0.mass" ) out = interpolation.db( p=p, DS="timeseries" )
      
      return(out)
    }
 

  # ---------


    if (DS =="habitat.temperatures") {

      bm = interpolation.db( p=p, DS="biomass" )
      ps = snowcrab_stmv(p=p, DS="output_data" ) # , voi=p$selection$name

      temp = ps$t * bm$h

      K = NULL
      nreg = length(p$regions)
      for (r in 1:nreg ){
        aoi = aegis::polygon_inside(x=bs[ , c("plon", "plat")], region=p$regions[r], planar=T)
        aoi = intersect( aoi, which( bs$plon > 250 ) )
        out = matrix( NA, nrow=p$ny, ncol=2) 
        
        for (y in 1:p$ny) {
          iHabitat = which( bm$h[,y] > p$habitat.threshold.quantile ) # any area with bm > lowest threshold
          iHabitatRegion = intersect( aoi, iHabitat )
          out[ y, 1] = mean( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
          out[ y, 2] = sd( temp[iHabitatRegion,y] , na.rm=TRUE ) # temperature weighted by Pr
        }
        
        ok = as.data.frame( out )
        names( ok) = c("temperature", "temperature.sd")
        ok$region = p$regions[r]
        ok$yr = p$yrs
        ok$lbound = ok$temperature - ok$temperature.sd*1.96 # normal assumption 
        ok$ubound = ok$temperature + ok$temperature.sd*1.96 
        K = rbind(K, ok)
      }

      return( K )      

    }


    return ("Completed" )
  }


