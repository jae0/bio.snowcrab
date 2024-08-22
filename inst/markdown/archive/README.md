# Notes -- archived reports using rmarkdown (mostly)

The files in this directory are markdown (Rmarkdown) reports and script sequences to be used for the analysis of the snow crab data. 

## Viewing:

The can be viewed with a web-browser (some might require Markdown viewing extensions), or VScode's markdown preview or Rstudio's viewer. 

## Rendering

If you have GNU Make installed then it is as simple as running from the installation directory:

```
	make snowcrab_working_paper YR=2023 
```

This will copy relevant files to a temporary directory and run the rendering programs there. This is because modifying the source directory is unwise as it is a git repsitory or user write-permissions are usually constrained. 

Add other options as required (eg, SOURCE (source directory), WK (work directory), etc... See inside the header of the Makefile for options.)

-- OR --

If you do not have GNU Make installed then you can run the commands inside the Makefile directly from a command prompt. But again the relevant file need to be copied or sym-linked to a working directory as there is usually write access is not permitted in the library directory. Copy or link also the markdown/media to the same work directory and then run the required commands. 

Use the Makefile commands as hints. E.g. for snowcrab_presentation_methods.html :

```
	ln -sf ~/bio/bio.snowcrab/inst/markdown/media ~/bio.data/bio.snowcrab/assessments/media
	cp ~/bio/bio.snowcrab/inst/markdown/snowcrab_presentation_methods.rmd ~/bio.data/bio.snowcrab/assessments/

	Rscript -e "rmarkdown::render('~/bio.data/bio.snowcrab/assessments/snowcrab_presentation_methods.rmd', params=list( year.assessment=2023), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)', output_dir='~/bio.data/bio.snowcrab/assessments/' ) " 
```

-- OR -- 

In Rstudio or MS VScode there looks to be funcationality to render these documents. For Rmarkdown copy the report file and rename extension to *.Rmd or if Quarto, then *.qmd. (In testing, the file name extension seems to required at least for Quarto). Configuration of the code editors are platform-specific and a moving target so check current online resources. 

Do your tests with quick_test.md : 

```
	make rmarkdown quick_test YR=2023 
```



