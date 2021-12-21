Stock assessment of Canada's Maritimes Region snow crab (Chionoectes oplio) leveraging aegis*, bio*, and stm* packages.


Installation:


1. To install:

```
  install.packages( "remotes", ask=FALSE, dependencies=TRUE ) # to inter-operate with github
  remotes::install_github( "jae0/aegis" ) # to bootstrap by installing directly from github
  remotes::install_github( "jae0/bio.snowcrab") # install bio.snowcrab and other required packages
```



2. Then, you need to have an Rprofile set up properly. Use the following, being careful to define the required R-global variables (see also: https://github.com/jae0/aegis/src/master/R/project.Rprofile.example.r):


```.

libPaths("~/R")
homedir = path.expand("~")
tmpdir = file.path( homedir, "tmp" )
work_root = file.path( homedir, "work" )    ### replace with correct path to work directory (local temporary storage)
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

# store your passwords and login here and make sure they are secure
try ( source( file.path( homedir, ".passwords" ) ) )

require( aegis )
```
 
For usage, examples can be found in aegis.*.

Here is a more expanded .Rprofile version that I use:


```
## options

ncols = Sys.getenv("COLUMNS")
if(nzchar(ncols)) options(width = as.integer(ncols))


repositories = getOption("repos")
repositories["CRAN"] = "https://cloud.r-project.org/"
repositories["INLA"] = "https://inla.r-inla-download.org/R/testing"
#repositories["CRAN"] = "http://mirror.its.dal.ca/cran/"
#repositories["INLA"] = "https://inla.r-inla-download.org/R/stable"

options(
  repos = repositories,
  papersize = "letter",
  width = 120,
  browser = "/usr/bin/chromium",
  mc.cores = parallel::detectCores(),
  pdfviewer="zathura",
  deparse.max.lines = 2,
  prompt="> ",
  digits=9,
  showWarnCalls=TRUE,
  showErrorCalls=TRUE,
  show.signif.stars=FALSE
)


# directory/file paths
homedir = path.expand("~")
tmpdir = file.path( homedir, "tmp" )
work_root = file.path( homedir, "tmp" )		 ### replace with correct path
code_root = file.path( homedir, "bio" )   ### replace with correct path
data_root = file.path( homedir, "bio.data" )   ### replace with correct path

setwd( code_root )  # this is to make TAB-completion easier inside of R ...
try( source( file.path( homedir, ".passwords" ) ) )  # could also just use REnviron for this ..

set.seed(12345)



# # alter quit/save default options
quit_without_saving = function() .Internal(quit("no", 0, FALSE))
close_all_graphics_windows = function(...) { graphics.off(...); invisible(gc()) }
load_bio_enviroment = function() { try( source( file.path( code_root, "bio_startup.R" ) ) ) }  # various user-specific local options


pkgs_reqd = unique( c(
  "akima", 
  "alphahull", 
  "BayesX", 
  "bigmemory",  
  "biglm", 
  "Cairo", 
  "colorspace",  
  "DBI", 
  "deSolve", 
  "emdbook", 
  "fasterize", 
  "grid",
  "GADMTools" , 
  "Hmisc", 
  "interp", 
  "raster",
  "FactoMineR",
  "fields", 
  "fftwtools", 
  "fftw", 
  "gstat",
  "geosphere", 
  "GillespieSSA", 
  "googlesheets4", 
  "lubridate",  
  "lattice",
  "maps", 
  "mapdata", 
  "maptools", 
  "maptools", 
  "maps", 
  "mapdata", 
  "nimble",  
  "numDeriv", 
  "pacman",
  "parallel",  
  "parallel", 
  "RColorBrewer",
  "RandomFields", 
  "RandomFieldsUtils",
  "R2BayesX",  
  "reshape2", 
  "rgeos",  
  "rgdal", 
  "splancs", 
  "SimInf", 
  "sf", 
  "sp", 
  "spdep",  
  "term", 
  "tmap",
  "truncnorm",
  "vegan" 
) )


pkgs_local_git = unique( c(
  "adapt",  
  "aegis",  
  "aegis.odemod",  
  "aegis.bathymetry",    
  "aegis.coastline",    
  "aegis.condition",    
  "aegis.metabolism",    
  "aegis.mpa",    
  "aegis.polygons",    
  "aegis.sizespectrum",    
  "aegis.substrate",    
  "aegis.survey",    
  "aegis.speciesarea",    
  "aegis.speciescomposition",    
  "aegis.temperature",    
  "bio.taxonomy",
  "bio.snowcrab",
  "carstm", 
  "ecomod",
  "netmensuration",  
  "stmv"
) )


pkg_install_git = function(source="github", ...) {
  #\\ add install_github flags e.g. force=TRUE to call if desired
  pkgsInstalled = .packages(all.available = TRUE)
  if ( ! "remotes" %in% pkgsInstalled ) install.packages( "remotes", dependencies=TRUE )
  bio = data.frame( libname=pkgs_local_git  )
  bio$local = file.path(code_root, bio$libname) 
  bio$github = paste( "jae0", bio$libname, sep="/") 
  bio$bitbucket = paste( "ecomod", bio$libname, sep="/") 
  if (source=="local") {
    for ( pkg in unique(bio$local) ) try( remotes::install_git( pkg, dependencies=FALSE, upgrade=TRUE, ... ) )
  } else if (source=="github"){
    for ( pkg in unique(bio$github) )  try( remotes::install_github( pkg, ... ) )
  } else if (source=="bitbucket"){
    for ( pkg in unique(bio$bitbucket) )  try( remotes::install_github( pkg, ... ) )
  } 
}

 
pkg_update = function(...)  update.packages( ask=FALSE, checkBuilt=TRUE )

pkg_reinstall = function(lib = .libPaths()[1], pkgs = as.data.frame(installed.packages(lib), stringsAsFactors=FALSE)$Package) {
  install.packages( lib=lib, pkgs=pkgs, type = 'source', dependencies=TRUE )
}


makeActiveBinding(".cc", close_all_graphics_windows, .GlobalEnv)
makeActiveBinding(".qq", quit_without_saving, .GlobalEnv)
makeActiveBinding(".bio", load_bio_enviroment, .GlobalEnv)
makeActiveBinding(".gi", pkg_install_git, .GlobalEnv)
makeActiveBinding(".pk", pkg_update, .GlobalEnv)
makeActiveBinding(".pkr", pkg_reinstall, .GlobalEnv)
 
libs0 = search()  # initial libs in memory
garbage_collect_libraries = function() {
  newobjs = setdiff( search(), libs0 )
  libs_new = newobjs[ grepl( "^package[:]", newobjs ) ] 
  # libs_new = gsub("package:", "", libs_new)
  if (length(libs_new)>0)  for ( pkg in libs_new) detach(pkg, character.only=TRUE )
  invisible( gc() )
}
makeActiveBinding(".gcl", garbage_collect_libraries, .GlobalEnv)

vars0 = ls()
garbage_collect_memory = function() {
  newvars = setdiff( ls(), vars0 )
  if (length(newvars)>0)  for ( nv in newvars) rm(list=nv)
  invisible( gc() )
}
makeActiveBinding(".gcm", garbage_collect_memory, .GlobalEnv)

environment_reinitialize = function() {
  close_all_graphics_windows() 
  garbage_collect_memory()
  garbage_collect_libraries()
  load_bio_enviroment()
}
makeActiveBinding(".init", environment_reinitialize, .GlobalEnv)


# ---- main init complete



.First = function(){
  # make some functions operate without parentheses:
  # store your passwords and login here and make sure they are with secure permissisions
  cat("---\nR session started at", date(), "with the following directories: \n")
  cat("\n  homedir = ", homedir )
  cat("\n  work_root = ", work_root )
  cat("\n  code_root = ", code_root )
  cat("\n  data_root = ", data_root )
  cat("\n")
  cat("\nShortcuts (enter command without parenthesese):\n")
  cat("\n  To (re)load environment: .bio" )
  cat("\n  To (re)initialize environment: .init" )
  cat("\n  To detach libraries: .gcl" )
  cat("\n  To clean memory: .gcm" )
  cat("\n  To quit: .qq" )
  cat("\n---\n\n")
}

.Last = function(){
  cat("\nSession ended at ", date(), "\n\n")
}



# ---- load bio specific libs:

# uncomment this out if you do want to load the environment automatically  
# -- this loads aegis, and any other libraries ..
# -- separation helps prevent race conditions when installing libraries
# .bio  


if (0) {

  # .bio starts this .. saved in file: file.path( code_root, "bio_startup.R" )  

    pkgsInstalled = .packages(all.available = TRUE)

    # specialized packages
    if (!"remotes" %in% pkgsInstalled ) utils::install.packages("remotes", dependencies=TRUE, ask=FALSE)
    if (!"INLA" %in% pkgsInstalled ) utils::install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
    # if (!"devtools" %in% pkgsInstalled ) utils::install.packages("devtools", dependencies=TRUE, ask=FALSE)
    # if (!"posterior" %in% pkgsInstalled ) remotes::install_github("stan-dev/posterior")
    # if (!"cmdstanr" %in% pkgsInstalled ) remotes::install_github("stan-dev/cmdstanr")
    # if (!"ROracle" %in% pkgsInstalled ) remotes::install_github("cran/ROracle", args=" --configure-args='--with-oci-lib=/usr/lib --with-oci-inc=/usr/include'", force=TRUE)
    # if (!"sf" %in% pkgsInstalled ) remotes::install_github("r-spatial/sf")


    if (!"aegis" %in% pkgsInstalled ) {
      message( "The package, aegis is missing. Install right now? (y/n):")
      ordl = readline()
      if (ordl=="y") remotes::install_github( "jae0/aegis")
    }


    # library
    # library("rstan") # observe startup messages
    # .libPaths("~/R/x86_64-pc-linux-gnu-library/4.0/")
    # "posterior", "bayesplot",  # not working right now 8 jun 2021
    # for (pkg in pkgs_reqd) if (! pkg %in% pkgsInstalled ) utils::install.packages( pkg, dependencies=TRUE, ask=FALSE)


    require(INLA)
    # inla.setOption(mkl=FALSE)
    # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
    inla.setOption(pardiso.license="~/paradiso.license")
    inla.pardiso.check()
    inla_num.threads =  floor(parallel::detectCores()/2 )
    inla.setOption(num.threads=paste(inla_num.threads, ":2", sep="") )


    # to force use, run with: control.compute=list(openmp.strategy="pardiso")

    suppressMessages( require( aegis ) )

}

```
