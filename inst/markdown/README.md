# Notes

The files in this directory are markdown that can be rendered into other formats via Quarto (via knittr and pandoc). 

The main data assimilation and modelling steps are odd numbered:

- ![01_snowcrab_data.md](01_snowcrab_data.md)
- ![03_biomass_indx_carstm.md](03_biomass_indx_carstm.md)
- ![05_snowcrab_fishery_model_turing.md](05_snowcrab_fishery_model_turing.md)


The main summary documents are even numbered:

- ![02_fishery_summary.md](02_fishery_summary.md)
- ![02_survey_summary.md](02_survey_summary.md)
- ![04_ecosystem_summary.md](04_ecosystem_summary.md)
- ![06_assessment_summary.md](06_assessment_summary.md)

Un-numbered files are presentations and operational documents, also Quarto/markdown formatted. They are basic templates that are intended to get you a working starting pointb that can be rapidly adpated for a given purposey. Rmarkdown/beamer seems to work better for presentations at present. Quarto formatting tweaks are needed for presentations (for now).


## Viewing:

Markdown documents can be viewed with:

- a web-browser (some might require Markdown viewing extensions)
- directly from Github (again through a web browser)
- VScode's markdown preview 
- Rstudio's viewer
- any text editor

## Rendering

Note that VSCode has a Makefile plugin with which you can run the makefile directly. But you will need quarto and gnu-make installed on your system (work fine on Linux but it has not been tested in MSWindows yet). 

### Quarto

If you have GNU Make installed then it is as simple as running from the  directory where you have downloaded this project: 

```shell
# location of the markdown files
cd ~/bio/bio.snowcrab/inst/markdown   

# {via Quarto}
make quarto FN=02_fishery_summary YR=2024 SOURCE=~/bio/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments DOCEXTENSION=html PARAMS="-P year_assessment:2024 -P sens:1" 

```

This will copy relevant files to a temporary directory (WK) and run the rendering programs there. This is because modifying the source directory is unwise as it is a git repsitory or user write-permissions are usually constrained. 

Add other options as required (eg, SOURCE (source directory), WK (work directory), etc... See inside the header of the Makefile for options.) 

Note: also the awkward passing of parameters to Quarto ... it is also duplcated to ensure Rmarkdown will also work.


### Rmarkdown

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

Rscript -e "rmarkdown::render('~/bio.data/bio.snowcrab/assessments/snowcrab_presentation_methods.rmd', params=list( year_assessment=2023), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)', output_dir='~/bio.data/bio.snowcrab/assessments/' ) " 


```

### VScode or Rstudio

In Rstudio or VScode there is functionality to render these documents in place. For Rmarkdown copy the report file and rename extension to *.Rmd; or if Quarto, then *.qmd. Using Rmarkdown will require more reformatting and tweaks as Quarto and Rmarkdown are not exactly the same.

In testing, the file name extension seems to required at least for Quarto. Configuration of the code editors are platform-specific and a moving target so check current online resources. 


### Pandoc

You can do everything directly through pandoc. That is covered extensively on the web (but still a moving target; the hints in the Makefile should be enough...), so do your own research.


## Testing

Do your tests with quick_test.md : 

```
	make quarto quick_test YR=2023 
```



