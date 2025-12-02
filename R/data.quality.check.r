
  data.quality.check = function( type, p ) {
    
    out = NULL
    out_message = NULL
    out_toprint = NULL

    # duplicated stations
    if (type=="stations") {
      set = snowcrab.db( DS="setInitial" )
      must.be.unique = paste( set$yr, set$station, sep="~" )
      i = which(must.be.unique %in% must.be.unique[which(duplicated(must.be.unique))]) # all i
      if (length(i) > 0) {      
        out = set[ i, c( "trip", "set", "station", "yr" ) ]
        j = which(is.finite(out$station))
        if (length(j) > 0) {
          out = out[ j, ]
          out = out[ order(out$yr, out$station) ,]
          out_toprint = paste( out$trip, out$set, out$station)
          out_message = "----- Duplicated trip, sets, stations ----- "
        }
      }
    }

    # counts of stations by area
    if(type=="count.stations") {
      set = snowcrab.db( DS="setInitial" )
      years = sort( unique( set$yr ) )
      nyears = length(years)
      nregions = length(p$regions)
      res = matrix( NA, nrow=nyears, ncol=nregions)
      for (r in 1:nregions) {
        nr =NULL
        nr = polygon_inside(x=set, region=p$regions[r], planar=F)
        for (y in 1:nyears) {
          ni = which( set$yr==years[y] )
          if (length(ni) > 0) {      
            res[y,r] = length( unique( intersect (nr, ni) ) )
          }
      }}

      out = as.data.frame(res)
      names(out) = c(p$regions)
      out$yr = years
      out = out[ order(out$yr), c("yr", p$regions)]
      out_message = "----- Number of stations: trip/sets ----- "
    }

    # positional information
    if(type=="position") {
      set = snowcrab.db( DS="setInitial" )
      inside = polygon_inside( set[, c("lon", "lat") ], "cfaall")
      if (length (inside) == nrow(set) ) {
        out_message = "All data are within positional bounds"
      } else {
        outside = setdiff( 1:nrow(set), inside )
        out = set[outside,]
        out_toprint = paste( out$trip, out$set, out$station)
        out_message = "----- The above are out of the cfa bounds ----- "  
      }
    }

    # positional information
    if(type=="position.difference") {
      set = snowcrab.db( DS="setInitial" )
      set.sub = split(set, set$station)
      
      for(i in 1:length(set.sub)){
        sta.match = set.sub[[i]]
        ave.lon.start = mean(sta.match$lon)
        ave.lon.end = mean(sta.match$lon1)
        ave.lat.start = mean(sta.match$lat)
        ave.lat.end = mean(sta.match$lat1)
        out = NULL
 
        for(j in 1:nrow(sta.match)){
          if(is.na(ave.lon.start) || is.na(ave.lat.start) || is.na(ave.lon.end) || is.na(ave.lat.end)){
            print("Contains NA positions: ")
            print( sta.match[ j, ] )
            out= rbind(out, sta.match[ j, ])
          }
          else{
          if(abs(sta.match$lon[j] - ave.lon.start) > .25){
            print(paste("Average start longitude: ", ave.lon.start, sep = ""))
            print( sta.match[ j, ] )
            out= rbind(out, sta.match[ j, ])
          }
          if(abs(sta.match$lon1[j] - ave.lon.end) > .25){
            print(paste("Average end longitude: ", ave.lon.end, sep = ""))
            print( sta.match[ j, ] )
            out= rbind(out, sta.match[ j, ])
          }
          if(abs(sta.match$lat[j] - ave.lat.start) > .25){
            print(paste("Average start latitude: ", ave.lat.start, sep = ""))
            print( sta.match[ j, ] )
            out= rbind(out, sta.match[ j, ])
          }
          if(abs(sta.match$lat1[j] - ave.lat.end) > .25){
            print(paste("Average end latitude: ", ave.lat.end, sep = ""))
            print( sta.match[ j, ] )
            out= rbind(out, sta.match[ j, ])
          }
          }
        }
      }
      
      out = out[ order(out$yr, out$trip, out$station) ,]

      out_toprint = paste( out$yr, out$trip, out$set, out$station)

      out_message = "----- Stations that are outside of historical positions -----"
      
    }
      

    if (type=="biologicals_fishno") {
       det = snowcrab.db( DS="det.rawdata" )
       names( det ) = rename.bio.snowcrab.variables(names(det) )
       i = which( !is.finite(det$crabno) )
       setDT(det)
       out = det[i, ][order(sdate),]
       out_toprint = paste( out$trip, out$set,  out$sdate, out$crabno)
       out_message = "----- Trip/sets that have no fishno -----"
    }


    if (type=="biologicals_duplicates") {
      det = snowcrab.db( DS="det.rawdata" )
      names( det ) = rename.bio.snowcrab.variables(names(det) )
      uid = paste(det$trip, det$set, det$sdate, det$crabno, sep="~")
      ii = which(duplicated( uid))
      setDT(det)
      ndups = length(ii)
      if (ndups > 0) {
        out = det[ uid %in% uid[ii] , ][order(sdate),]
        out_toprint = paste( out$sdate, out$trip, out$set, out$crabno)
        out_message = "----- Duplicated data found in det raw data... fix these above in ISDB .. -----"
      }
    }



    if (type=="biologicals_morphology") {

        set = snowcrab.db( DS="set.clean")
        det = snowcrab.db( DS="det.initial")
        setDT(set)
        setDT(det)
        set$sid = paste(set$trip, set$set, sep="~")
        det$sid = paste(det$trip, det$set, sep="~")
        set$year = set$yr 
        set$region = NA
        for ( region in regions ) {
            r = polygon_inside(x=set, region=region, planar=F)
            if (length(r) > 0) set$region[r] = region
        }
        set= set[!is.na(region), ]
        set = set[, .(sid, region, year, sa, t, z, timestamp, julian, lon, lat)]
        det = det[, .(sid, shell, cw, sex, mass, mat, gonad, durometer, abdomen, chela )]
        basedata = det[ set, on=.(sid)]
        
        out =list()

        # trim a few strange data points
        o = lm( log(mass) ~ log(cw) + as.factor(sex), basedata)
        out$allometry = basedata[which(abs(o$residuals) > 0.5),]

        large_females = basedata[ sex=="1" & cw > 90, which=TRUE ]  
        if (length(large_females)>0) {
            sex_male = which(!is.finite(basedata$abdomen[large_females]))  # must be male
            sex_female = which(is.finite(basedata$abdomen[large_females]))  # must be female
            basedata$sex[large_females][sex_male] = "0"  # recode to male 
            out$large_females = basedata[large_females[sex_female],]
        }
   
        #Cw.errors: Carapace Width below 5 or greater than 185
        out$cw_lt5_gt185 = det[ which(det$cw<5 | det$cw>185 ),] 

        #Chela.errors: Chela less than 1 or greater than 50
        out$chela_lt1_gt50 = det[which(det$chela < 1 | det$chela > 50  ),] 

        #Abdomen.errors:Abdomen less than 1 and greater than 66
        out$abdomen_lt1_gt66 = det[which(det$abdomen < 1 | det$abdomen > 66 ),]
  
        #Mass.errors: Mass less than 1 or greater than 1500
        out$mass.out_lt1_gt1500 = det[which( det$mass < 1 | det$mass > 1500  ),]

        #Sex.a: Indeterminate sex based on measurements taken (abdomen values where sex=male)
        out$sex.abdomen.male = det[which(is.finite( det$abdomen ) & det$sex==male),] 

        #Sex.c: Indeterminate sex based on measurements taken (chela values where sex=female
        out$sex.chela.female = det[which(is.finite( det$chela ) & det$sex==female),] 

        #Mat.errors: Unknown Maturity
        out$mat.unknown_withdata = det[which(det$mat ==2 & (is.finite(det$chela) | is.finite(det$abdomen))),]
    

 
          # these are global parameters
          # # sex codes
          # male = 0
          # female = 1
          # sex.unknown = 2

          # # maturity codes
          # immature = 0
          # mature = 1
          # mat.unknown = 2

        # out$large_females_immature = basedata[ sex=="1" & cw > 80 & mat=="0", ] 
        # out$large_females_mature = basedata[ sex=="1" & cw > 90 & mat=="1", ] 
        out$small_female_mature = basedata[ sex=="1" & cw <35 & mat=="1", ] 

        # out$large_males_immature = basedata[ sex=="0" & cw > 135 & mat=="0", ] 
        # out$large_male_mature = basedata[ sex=="0" & cw > 150 & mat=="1", ] 
        out$small_male_mature = basedata[ sex=="0" & cw <49 & mat=="1", ] 
         
        out_toprint = paste( names(out) )
        out_message = "----- Strange size data det raw data... fix these above in ISDB .. -----"
     
    }


    if (type=="seabird.load") {
        SS = snowcrab.db( DS="set.clean")
        i = which( is.na( SS$seabird_uid) & SS$yr %in% p$seabird.yToload & SS$yr >= 2012 )
        if (length(i) > 0) {
          out_message = "----- Missing seabird matches ----- "
          out = SS[ i,] 
          out = out[ order(out$yr, out$trip, out$station) , ] 
          out_toprint = paste( out$trip, out$set, out$station)
        }
    }


    if (type=="minilog.load") {
        SS = snowcrab.db( DS="set.clean")
        i = which( is.na( SS$minilog_uid) & SS$yr %in% p$minilog.yToload & SS$yr >= 2004 )
        if (length(i) > 0) {
          out = SS[ i,] 
          out = out[ order(out$yr, out$trip, out$station) , ] 
          out_toprint = paste( out$trip, out$set, out$station)
          out_message = "----- Missing minilog matches -----"
      }
    }

    if (type=="netmind.load") {
        SS = snowcrab.db( DS="set.clean")
        i = which( is.na( SS$netmind_uid) & SS$yr %in% p$netmind.yToload & SS$yr >= 2004 )
        if (length(i) > 0) {
          out =  SS[ i,] 
          out = out[ order(out$yr, out$trip, out$station) , ] 
          out_toprint = paste( out$trip, out$set, out$station)
          out_message = "----- Missing netmind matches ----- "
        }
    }


    if (type=="tow.duration") {
      set = snowcrab.db( DS="set.clean" )
      i = which( ( set$dt > 9  | set$dt < 3.5 )  & set$yr >=2004 )
      if  (length(i)>0 ) {
        out = set[i, c("trip", "set", "station", "dt", "timestamp", "yr")] 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_toprint = paste( out$trip, out$set, out$station)
        out_message = "----- The above have rather short/long tow times (dt) -----" 
      }
    }


    if (type=="tow.distance") {
      # expected = 2 knots * 5 min = 2 * 1.852 * 5/60 = 0.309 km ( so a good range is {-25%, +75%} = (0.232, 0.5408)
      set = snowcrab.db( DS="set.clean" )
      i = which( ( set$distance > 0.541  | set$distance < 0.232 )  & set$yr >=2004 )
      if  (length(i)>0 ) {
        out = set[i, c("trip", "set", "station", "distance", "timestamp", "yr")] 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_toprint = paste( out$trip, out$set, out$station)
        out_message = "----- The above have rather short/long tow distances -----" 
      }
    }


    if (type=="netmind.timestamp") {
      #  check times/data and merge remaining data using datestamps and {station, set}
      set = snowcrab.db( DS="set.clean" )
      nm = netmind.db("stats")
      nm$netmind.timestamp = nm$t0
      set = merge(
        set[,c(
          "trip", "set", "lon", "lat", "timestamp", "seabird_uid", "minilog_uid", "netmind_uid", "yr")],
        nm[ ,c(
          "netmind_uid","netmind.timestamp", "slon", "slat" )], 
        by="netmind_uid", all.x=TRUE, all.y=FALSE )
      
      time.diff = difftime( set$netmind.timestamp, set$timestamp )
      time.thresh = lubridate::minutes(30)
      i = which( abs( time.diff ) > time.thresh )
      if (length(i)>0) {
        out = set[i, ] 
        out = out[ order(out$timestamp, out$trip, out$set) , ] 
        out_toprint = paste( out$trip, out$set )
        out_message = "----- Potential date/time mismatches above -----"
      }
    }

    # netmind mismatches
    if(type=="netmind.mismatches") {
      set = snowcrab.db( DS="set.clean" )
      i = which( set$yr > 2004 & (set$netmind_uid==""| is.na(set$netmind_uid) )  )
      if ( length (i) > 0 ) {
        out = set[i,c("trip", "set", "station", "t0", "timestamp", "yr") ] 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_message =  "----- No netmind matches for the above sets -----"
      }
    }


    # poor minilog matches
    if (type=="minilog") {
      set = snowcrab.db( DS="set.clean" )
      must.be.unique = set$t0
      i = which(must.be.unique %in% must.be.unique[which(duplicated(must.be.unique))]) # all i
      if ( length (i) > 0 ) {
        x = set[ i, c( "trip", "set", "station", "t0", "yr" ) ]
        dup.t0 = x[ is.finite(x$t0), ]
        out =  dup.t0 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_message =  "----- Duplicated minilog times above -----" 
      }
    }

    # minilog mismatches
    if (type=="minilog.mismatches") {
      set = snowcrab.db( DS="set.clean" )
      i = which( set$yr > 2004 & (set$minilog_uid==""| is.na(set$minilog_uid) )  )
      if ( length (i) > 0 ) {
        out = set[i,c("trip", "set", "station", "t0", "timestamp", "yr") ] 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_message = "----- No minilog matches for the above sets -----"
      }
    }

    if (type=="minilog.dateproblems") {
      set = snowcrab.db( DS="set.clean" )
      time.thresh = lubridate::hours(1)
      i = which( abs( difftime( set$t0, set$timestamp )) > time.thresh )
      if ( length (i) > 0 ) {
        out = set[i, c("trip", "set", "station", "t0", "timestamp", "yr" )] 
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_message = "----- Date mismatches with Trip ID -----" 
      }
    }

    if (type=='na.spread') {
      set = snowcrab.db( DS="set.clean" )
      i = which(is.na(set$spread))
      if ( length (i) > 0 ) {
        out = set[i,c('trip','set','station', 'spread', "yr" )]
        out = out[ order(out$yr, out$trip, out$station) , ] 
        out_message = "----- NA in tow spread ----- "  
      }
    }

    if (type=='na.distance') {
      set = snowcrab.db( DS="set.clean" )
    	i = which(is.na(set$distance))
    	out = set[i,c('trip','set','station','distance' , "yr")]
      out = out[ order(out$yr, out$trip, out$station) , ] 
    	out_message = "----- NA in tow distance ----- " 
    }

    if (!is.null(out)) {
      
      message("\n")

      if (!is.null( out_toprint )) {
        print(  out_toprint ) 
      } else {
        print(   out )
      }

      if (!is.null(out_message)) message( out_message )

    } else {
      
      message("\n----- no issues detected -----\n" )
    
    }

    return(out)
  }

