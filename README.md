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
