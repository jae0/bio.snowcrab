# Notes

The files in this directory are markdown that can be rendered into other formats via Quarto (mostly), Rmarkdown (some old presentations), and pandoc. 

## Viewing:

They can be viewed with a web-browser (some might require Markdown viewing extensions), directly from Github, or VScode's markdown preview or Rstudio's viewer. 

## Rendering

If you have GNU Make installed then it is as simple as running from the  directory where you have downloaded this project: 

```shell
# location of the markdown files
cd ~/bio/bio.snowcrab/inst/markdown   

# {via Quarto}
make quarto FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=html 

```

This will copy relevant files to a temporary directory (WK) and run the rendering programs there. This is because modifying the source directory is unwise as it is a git repsitory or user write-permissions are usually constrained. 

Add other options as required (eg, SOURCE (source directory), WK (work directory), etc... See inside the header of the Makefile for options.)

It is posible to process dierctly with Rmarkdown or pandoc, but formattng markups are not the same and so additional tweaking would be requierd:

```shell
# {via Rmarkdown -- but might need additional tweaks as Quarto formats are sufficiently different}
make rmarkdown FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=bookdown::pdf_document2  DOCEXTENSION=pdf 

# {via pandoc directly ... your mileage will vary}
make pdf FN=02_fishery_summary  
```

-- OR --

If you do not have GNU Make installed then you can run the commands inside the Makefile directly from a command prompt. But again the relevant file need to be copied or symbolically-linked to a working directory as write access usually is not permitted in the library directory. Copy or link also the markdown/media to the same work directory and then run the required commands. 

Use the Makefile commands as hints. E.g. for snowcrab_presentation_methods.html:

```shell
# this is an Rmarkdown example, where the embeddings are still Rmarkdown code

ln -sf ~/bio/bio.snowcrab/inst/markdown/media ~/bio.data/bio.snowcrab/assessments/media

cp ~/bio/bio.snowcrab/inst/markdown/snowcrab_presentation_methods.rmd ~/bio.data/bio.snowcrab/assessments/

Rscript -e "rmarkdown::render('~/bio.data/bio.snowcrab/assessments/snowcrab_presentation_methods.rmd', params=list( year.assessment=2023), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)', output_dir='~/bio.data/bio.snowcrab/assessments/' ) " 
```

-- OR -- 

In Rstudio or VScode there is functionality to render these documents in place. For Rmarkdown copy the report file and rename extension to *.Rmd; or if Quarto, then *.qmd. Using Rmarkdown will require more reformatting and tweaks as Quarto and Rmarkdown are not exactly the same.

In testing, the file name extension seems to required at least for Quarto. Configuration of the code editors are platform-specific and a moving target so check current online resources. 

Do your tests with quick_test.md : 

```
	make quarto quick_test YR=2023 
```



