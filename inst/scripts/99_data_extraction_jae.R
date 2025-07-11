#  these are convenience scripts to extract data from a MSWindows virtualbox session into linux (for Jae) with no   depenencies except DBI, Roracle


require(DBI)
require(ROracle)


fn_root = "C:/Users/choij/Desktop/datadump" 

year.assessment = 2024

yrs = 2022:(year.assessment + 1)  # add one to capture new year's data (for 4X)




fn.loc =  file.path( fn_root, "data", "trawl", "SNCRABSETS" )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )

con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
            # believeNRows=F required for oracle db's
for ( YR in yrs ) {
    fny = file.path( fn.loc, paste( YR,"rdz", sep="."))
    SNCRABSETS = NULL
    SNCRABSETS = ROracle::dbGetQuery(con, paste("select * from SNCRABSETS
                                                where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
                                                OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
    #Remove stations from previous years assesment
    ind = which(year(SNCRABSETS$BOARD_DATE)==YR & month(SNCRABSETS$BOARD_DATE) == 1) 
    if(length(ind)>0){
        SNCRABSETS = SNCRABSETS[-ind,]
    }
    if(nrow(SNCRABSETS) == 0){
        print(paste("No sets for ", YR)) 
    }
    else{
    read_write_fast( SNCRABSETS, file=fny, compress=TRUE)
    gc()  # garbage collection
    print(YR)
    }
}

ROracle::dbDisconnect(con)



fn.loc =  file.path( fn_root, "data", "trawl", "SNCRABDETAILS" )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE  )

con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

for ( YR in yrs ) {
    fny = file.path( fn.loc, paste( YR,"rdz", sep="."))
    SNCRABDETAILS = NULL
    #in following line replaced sqlQuery (Rrawdata) with  ROracle::dbGetQuery (ROracle)
    SNCRABDETAILS = ROracle::dbGetQuery(con,
        paste("select * from SNCRABDETAILS
                where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
                                                OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
    read_write_fast( SNCRABDETAILS, file=fny, compress=TRUE)
    gc()  # garbage collection
    print(YR)
}

ROracle::dbDisconnect(con)



fn.loc =  file.path( fn_root, "data", "trawl", "SNTRAWLBYCATCH" )

dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE  )

con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

for ( YR in yrs ) {
    fny = file.path( fn.loc, paste( YR,"rdz", sep="."))
    SNTRAWLBYCATCH = NULL
    #in following line replaced sqlQuery (Rrawdata) with  ROracle::dbGetQuery (ROracle)
    SNTRAWLBYCATCH = ROracle::dbGetQuery(con,
        paste("select * from SNTRAWLBYCATCH
                where EXTRACT(YEAR from BOARD_DATE) = ", YR , "
                                                OR (EXTRACT(YEAR from BOARD_DATE) = ", YR+1 , " AND EXTRACT(MONTH FROM Board_DATE)=1)") )
    read_write_fast( SNTRAWLBYCATCH, file=fny, compress=TRUE)
    gc()  # garbage collection
    print(YR)
}

ROracle::dbDisconnect(con)



fn.loc =  file.path( fn_root, "data", "logbook", "datadump" )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )
con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

for ( YR in yrs ) {
    fny = file.path( fn.loc, paste( YR,"rdz", sep="."))
    query = paste(
        "SELECT * from marfissci.marfis_crab ",
        "where target_spc=705",
        "AND EXTRACT(YEAR from DATE_LANDED) = ", YR )
    logbook = NULL
    #in following line replaced sqlQuery (Rrawdata) with  dbGetQuery (ROracle)
    logbook = ROracle::dbGetQuery(con, query )
    read_write_fast( logbook, file=fny, compress=T)
    gc()  # garbage collection
    print(YR)
}
ROracle::dbDisconnect(con)



fn.loc =  file.path( fn_root, "data", "logbook"  )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )
filename.licence = file.path( fn.loc, "lic.datadump.rdz" )
con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
lic = ROracle::dbGetQuery(con, "select * from marfissci.licence_areas")
read_write_fast(lic, file=filename.licence, compress=T)
ROracle::dbDisconnect(con)


fn.loc =  file.path( fn_root, "data", "logbook"  )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )
filename.areas = file.path( fn.loc, "areas.datadump.rdz" )
con=ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)
areas = ROracle::dbGetQuery(con, "select * from marfissci.areas")
read_write_fast(areas, file=filename.areas, compress=T)
ROracle::dbDisconnect(con)


fn.loc =  file.path( fn_root, "data", "observer", "datadump" )
dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )
con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

for ( YR in yrs ) {
fny = file.path( fn.loc, paste( YR, "rdz", sep="."))
odbq = paste(
    "SELECT s.LATITUDE, s.LONGITUDE, s.LANDING_DATE, s.SET_NO, s.PRODCD_ID, s.EST_CATCH, s.EST_KEPT_WT," ,
    "s.NUM_HOOK_HAUL, d.BOARD_DATE, d.FISH_NO, d.SEXCD_ID, d.FISH_LENGTH, " ,
    "d.FEMALE_ABDOMEN, d.CHELA_HEIGHT, d.SHELLCOND_CD, d.DUROMETRE, d.TRIP_ID, d.TRIP  " ,
    "FROM SNOWCRAB.SNCRABDETAILS_OBS d, SNOWCRAB.SNCRABSETS_OBS s " ,
    "WHERE d.TRIP_ID = s.TRIP_ID  " ,
    "AND d.SET_NO = s.SET_NO  " ,
    "AND d.FISH_NO Is Not Null" ,
    "AND EXTRACT(YEAR from d.BOARD_DATE) = ", YR )
odb = NULL
odb = ROracle::dbGetQuery(con, odbq )
read_write_fast( odb, file=fny, compress=T)
gc()  # garbage collection
print(YR)
}
ROracle::dbDisconnect(con)
   

# bycatch
fn.loc =  file.path( fn_root, "data", "observer", "bycatch" )
    dir.create( fn.loc, recursive = TRUE, showWarnings = FALSE )
 
con = ROracle::dbConnect(DBI::dbDriver("Oracle"),dbname=oracle.snowcrab.server , username=oracle.snowcrab.user, password=oracle.snowcrab.password, believeNRows=F)

for ( YR in yrs ) {
    fny = file.path( fn.loc, paste( YR, "rdz", sep="."))
    odbq = paste(
        "SELECT trip.trip_id, trip.trip, trip.board_date, trip.landing_date, st.set_no, vess.vessel_name, vess.license_no, vess.cfv,",
        "  isp.LATITUDE, isp.LONGITUDE, isp.DEPTH, isp.WATER_TEMPERATURE, ", 
        "  sc.common, st.comarea_id, st.nafarea_id,  ", 
        "  ca.speccd_id, ca.est_num_caught,  ca.est_kept_wt, ca.est_discard_wt, st.NUM_HOOK_HAUL,",
        "  fish.fish_no , fish.fish_length, fish.fish_weight ", 
        "FROM istrips trip, isgears gr, isfishsets st, iscatches ca, isfish fish, issetprofile isp, isvessels vess,  isobservercodes o, isspeciescodes sc ", 
        "WHERE trip.tripcd_id = 2509",
        "AND vess.vess_id = trip.vess_id",
        "AND vess.license_no = trip.license_no",
        "AND o.obscd_id = trip.obscd_id",
        "AND trip.trip_id = gr.trip_Id",
        "AND st.fishset_id = isp.fishset_id ",
        "AND isp.fishset_id = ca.fishset_id ",
        "AND ca.speccd_id = sc.speccd_id ",
        "AND trip.tripcd_id=2509 ",
        "AND st.specscd_id=2526 ",
        "AND isp.pntcd_id=4 ",
        "AND (trip.trip_id = st.trip_id AND gr.gear_id = st.gear_id)",
        "AND ca.catch_id = fish.catch_id(+)",
        "AND EXTRACT(YEAR from trip.board_date) = ", YR )

    odb = NULL
    odb = ROracle::dbGetQuery(con, odbq )
    read_write_fast( odb, file=fny, compress=T)
    gc()  # garbage collection
    print(YR)
}
ROracle::dbDisconnect(con)



# temperature
basedir = fn_root
loc.bottomdatabase = file.path( basedir, "archive", "bottomdatabase"  )
dir.create( loc.bottomdatabase, recursive=T, showWarnings=F )

con = ROracle::dbConnect( DBI::dbDriver("Oracle"),
    username = oracle.snowcrab.user,
    password = oracle.snowcrab.password,
    dbname = oracle.snowcrab.server
)

res = ROracle::dbSendQuery( con, "ALTER SESSION SET NLS_DATE_FORMAT = 'YYYY-MM-DD'")
res = ROracle::dbSendQuery( con, "ALTER SESSION SET NLS_TIMESTAMP_FORMAT = 'YYYY-MM-DD HH24:MI:SSXFF'")
res = ROracle::dbSendQuery( con, "ALTER SESSION SET NLS_TIMESTAMP_TZ_FORMAT = 'YYYY-MM-DD HH24:MI:SSXFF TZR'")

ythreshold =1970
    
    yrs = c(2020:2024)
    
      for ( yt in yrs ) {
        # yt = 2023
        if (yt < ythreshold) {
          message( "Warning: ",  yt, "is not in this database ... skipping" )
          next()
        }
        Z = NULL
        Z = ROracle::dbGetQuery( con,  paste(
          " select * " ,
          " from SC_TEMP_MERGE " ,
          " where EXTRACT(YEAR from SC_TEMP_MERGE.T_DATE) =", yt
        ) )

        if (nrow(Z) > 0 ) {
          if (ncol(Z) == 6) names(Z) = c("project", "date", "lat", "lon", "t_uid", "t" )
          if (ncol(Z) == 7) names(Z) = c("project", "date", "lat", "lon", "t_uid", "t", "z" )
          Z$yr = yt
          Z$dyear = lubridate::decimal_date( Z$date ) - Z$yr
        } else {
          Z = NULL
        }
        if (!is.null(Z)) {
          if ( nrow(Z) > 0  ) {
            fn = file.path(  loc.bottomdatabase, paste("bottom", yt, "rdz", sep="."))
            print (fn)
            read_write_fast( Z, file=fn, compress=T)
          }
        }
      }


# groundfish data
outdir = fn_root

fn.root =  file.path( outdir,  "gscat" )
dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    for ( YR in yrs ) {
    fny = file.path( fn.root, paste( YR,"rdz", sep="."))
    gscat = ROracle::dbGetQuery( connect,  paste(
        "select i.*, substr(mission,4,4) year " ,
        " from groundfish.gscat i " ,
        " where substr(MISSION,4,4)=", YR)
    )

    names(gscat) =  tolower( names(gscat) )
    print(fny)
    read_write_fast(gscat, file=fny, compress=T)
    gc()  # garbage collection
    print(YR)
    }

ROracle::dbDisconnect(connect)

    

connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fn.root =  file.path( outdir,  "gsdet" )
    dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

    for ( YR in yrs ) {
    fny = file.path( fn.root, paste( YR,"rdz", sep="."))
    gsdet = ROracle::dbGetQuery( connect,  paste(
        "select i.*, substr(mission,4,4) year" ,
        " from groundfish.gsdet i " ,
        " where substr(mission,4,4)=", YR)
    )
    names(gsdet) =  tolower( names(gsdet) )
    gsdet$mission = as.character( gsdet$mission )
    read_write_fast(gsdet, file=fny, compress=T)
    print(fny)
    gc()  # garbage collection
    print(YR)
    }

ROracle::dbDisconnect(connect)

      

connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fn.root =  file.path(outdir,  "gsinf" )
    dir.create( fn.root, recursive = TRUE, showWarnings = FALSE  )

    for ( YR in yrs ) {
    fny = file.path( fn.root, paste( YR,"rdz", sep="."))
    gsinf = ROracle::dbGetQuery( connect,  paste(
        "select * from groundfish.gsinf where EXTRACT(YEAR from SDATE) = ", YR )
    )
    names(gsinf) =  tolower( names(gsinf) )
    read_write_fast(gsinf, file=fny, compress=T)
    print(fny)
    gc()  # garbage collection
    print(YR)
    }
ROracle::dbDisconnect(connect)


connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fn = file.path( outdir,   "gsgear.rdz")
    gsgear =  ROracle::dbGetQuery(connect, "select * from groundfish.gsgear", as.is=T)
    names(gsgear) =  tolower( names(gsgear) )
    read_write_fast(gsgear, file=fn, compress=T)
    print(fn)

ROracle::dbDisconnect(connect)


    
connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fn = file.path( outdir,   "gslist.rdz")
    gslist = ROracle::dbGetQuery(connect, "select * from groundfish.gs_survey_list")
    names(gslist) =  tolower( names(gslist) )
    read_write_fast(gslist, file=fn, compress=T)
    print(fn)

ROracle::dbDisconnect(connect)


connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fnmiss = file.path( outdir,   "gsmissions.rdz")
    gsmissions = ROracle::dbGetQuery(connect, "select MISSION, VESEL, CRUNO from groundfish.gsmissions")
    names(gsmissions) =  tolower( names(gsmissions) )
    read_write_fast(gsmissions, file=fnmiss, compress=T)
    print(fnmiss)

ROracle::dbDisconnect(connect)


connect = ROracle::dbConnect( DBI::dbDriver("Oracle"), dbname=oracle.groundfish.server, username=oracle.personal.user, password=oracle.personal.password, believeNRows=F)

    fnspc = file.path( outdir,   "spcodes.rdz" )
    spcodes = ROracle::dbGetQuery(connect, "select * from groundfish.gsspecies", as.is=T)
    names(spcodes) =  tolower( names(spcodes) )
    read_write_fast(spcodes, file=fnspc)
    print( fnspc )
    print("Should follow up with a refresh of the taxonomy.db " )

ROracle::dbDisconnect(connect)


require(aegis)
require(stringr)

# copy paste the functions in the file or: 
source( file.path( "bio", "bio.snowcrab", "R", "ca_database_functions.r") ) 

tac.db =  CA.getTable("TAC")
lic.hist = CA.getTable("NUMLICENSE")
lic.db = CA.getTable("LIC_SUMMARY")

dir.create( file.path( fn_root, "data", "CA", "CA_db" ), recursive=TRUE,  showWarnings=FALSE  )
read_write_fast( 
    list( tac.db=tac.db, lic.hist=lic.hist, lic.db=lic.db), 
    fn=file.path( fn_root, "data", "CA", "CA_db", "CA.rdz" ) 
)  

# copy to ~/tmp in host system and then in command terminal

rsync -av ~/tmp/datadump/data ~/bio.data/bio.snowcrab/ 
rsync -av ~/tmp/datadump/archive/bottomdatabase  ~/bio.data/aegis/temperature/archive/
rsync -av ~/tmp/datadump/gs* ~/bio.data/aegis/groundfish/data/
rsync -av ~/tmp/datadump/spcodes* ~/bio.data/aegis/groundfish/data/


# save on tethys

rsync -av ~/tmp/datadump/data jae@142.2.22.74:/archive/bio.data/bio.snowcrab/ 
rsync -av ~/tmp/datadump/archive/bottomdatabase  jae@142.2.22.74:/archive/bio.data/aegis/temperature/archive/
rsync -av ~/tmp/datadump/gs* jae@142.2.22.74:/archive/bio.data/aegis/groundfish/data/
rsync -av ~/tmp/datadump/spcodes* jae@142.2.22.74:/archive/bio.data/aegis/groundfish/data/

