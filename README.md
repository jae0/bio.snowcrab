Stock assessment of Canada's Maritimes Region snow crab (*Chionoectes oplio*) leveraging aegis*, bio*, and stmv packages.

This project, bio.snowcrab is used by the Maritimes snow crab group to:

  - Assimilate data from various sources: 
    - Minilog, Seabird, Marport, GPS data streams
    - Snow crab survey data, entered into the At-Sea-Observer data base system
    - Logbook data of landings and effort ( 100% monitored )
    - At-sea-observed fisheries data ( 5% of catch monitored )
  - QA/QC
  - Model spatiotemporal variations of number (Poisson process) and mean body weight (Gaussian process) using Bayesian Conditional Autoregressive Models. Heavy lifting is done by INLA[https://www.r-inla.org/].
  - Model aggregate timeseries by region as Bayesian biomass dynamics autoregressive process (Logistic form) using STAN and obtain biological parameters that guide Harvest control rules. 
  - Generate routine figures, tables, reports and presentations.

Much of this project is generic and can be easily adapted for other species. Usage is shown in the scripts found in inst/scripts/0*.R. They represent the backbone of the assessment.  

Quarto and Rmarkdown documents can be found in inst/markdown/. They are meant
to be copied to a work directory such as bio.data/bio.snowcrab/assessments/ where
they can be run to generate reports on demand.


There is heavy reliance upon aegis.bathymetry, aegis.polygons, aegis.surveys and aegis.temperature. Though not necessary, they help inform the broader ecosystem-based approach that has been used with snow crab assessments since 2004 (when we received the mandate in Maritimes Region).  Examples of their usage are found in the individual aegis.* projects inst/scripts/0*.R files.


---

Installation:

1. To install:

```
  install.packages( "remotes", ask=FALSE, dependencies=TRUE ) # to inter-operate with github
  remotes::install_github( "jae0/aegis" ) # to bootstrap by installing directly from github
  remotes::install_github( "jae0/bio.snowcrab") # install bio.snowcrab and other required packages
```

2. Then, you need to have an Rprofile set up properly. Use the following, being careful to define the required R-global variables (see also: https://github.com/jae0/aegis/src/master/R/project.Rprofile.example.r):

```.

homedir = path.expand("~")
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

require( aegis )

# store your passwords and login here and make sure they are secure
# try ( source( file.path( homedir, ".passwords" ) ) )


```

A more expanded version, similar to what I use, can be found below:

https://github.com/jae0/aegis/blob/master/inst/scripts/example_Rprofile.R


 
