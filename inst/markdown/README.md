
The files in this directory are different markdown/Rmarkdown/Quarto reports and script sequences to be used for the analysis of the snow crab data. 

If you have GNU Make installed then it is as simple as running from the installation directory:

make snowcrab_working_paper YR=2023 

Add other options as required (eg, SOURCE (source directory), WK (work directory), etc..)

Change target to options in the Makefile. 


If you do not have GNU Make installed then you can run the commands in the Makfile directly from a command prompt:

But they will need to be copied to a working directory as there is usually write access is not permitted in the library directory. Copy or link also the markdown/media to the same work directory and then run the required commands. 

Use the Makefile commands as hints. E.g. for snowcrab_presentation_methods.html :

	ln -sf ~/bio/bio.snowcrab/inst/markdown/media ~/bio.data/bio.snowcrab/reports//media
	cp ~/bio/bio.snowcrab/inst/markdown/snowcrab_presentation_methods.rmd ~/bio.data/bio.snowcrab/reports/

	Rscript -e "rmarkdown::render('~/bio.data/bio.snowcrab/reports//snowcrab_presentation_methods.rmd', params=list( year.assessment=2023), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='~/bio.data/bio.snowcrab/reports/' ) " 
    



