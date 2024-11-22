---
title: "quick_test to try a page before adding to main document "
subtitle: "this makes debugging a page easier and faster"
author: "Jae S. Choi"
# footnote: "jae.choi@dfo-mpo.gc.ca"
institute: "Bedford Institute of Oceanography, DFO Science"
# date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  beamer_presentation:
    theme: "metropolis"
    colortheme: "seagull"
    fonttheme: "professionalfonts"
    fig_caption: yes
    # latex_engine: pdflatex
    latex_engine: lualatex 
    keep_tex: true
classoption: 
  - aspectratio=169 #16:9 wide
  - t  # top align
header-includes: 
  - \usepackage{graphicx}
  - \usepackage[font={scriptsize}, labelfont={bf}]{caption}
  - \usepackage{float}
  - \usepackage{multicol}
  # - \usepackage{subfig}
  # - \newcommand{\btiny}{\begin{tiny}}
  # - \newcommand{\etiny}{\end{tiny}}
params:
  year.assessment: 2023
  media_loc: "media"
  debugging: FALSE
  loc_dde: ""
--- 


<!-- Preamble

This is a Markdown document ... To create HTML or PDF, etc, run: 


# for presentations to PDF (via beamer):
# note: section separation with '#' can confuse rmarkdown
  
  make rmarkdown FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=beamer_presentation DOCEXTENSION=pdf # {via Rmarkdown}
 
  make rmarkdown FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=html_document DOCEXTENSION=html # {via Rmarkdown}

# for html documents including presentations:
  make quarto FN=quick_test YR=2023 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=pdf # {via Quarto}


  make pdf FN=quick_test  # {via pandoc}

Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
   
Huge > huge > LARGE > Large > large > normalsize > small > footnotesize > scriptsize > tiny

-->



<!-- Set up R-environment -->

```{r setup, include=FALSE}
  require(knitr)
  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
#     dev.args = list(type = "cairo"),
    fig.retina = 2,
    dpi=192
  )

require(spsUtil)

quietly = spsUtil::quiet

require(ggplot2)
require(aegis)  # basic helper tools

year.assessment = 2024  # change this as appropriate
year_previous = year.assessment - 1

p = bio.snowcrab::load.environment( year.assessment=year.assessment )  

SCD = project.datadirectory("bio.snowcrab")
media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year.assessment)

#regions = c("cfanorth", "cfasouth", "cfa4x")
regions = c("cfanorth", "cfa23",  "cfa24", "cfa4x")
nregions = length(regions)


FD = fishery_data()  # mass in tonnes
fda = FD$summary_annual

library(tidyverse)
region = c(replicate(100,"A"),replicate(100,"B"),replicate(100,"C"))
variable = c(replicate(50,"var_a"),replicate(50,"var_b"),replicate(50,"var_a"),replicate(50,"var_b"),replicate(50,"var_a"),replicate(50,"var_b"))
value = runif(300, min = 0, max = 100)
df = data.frame(region = region, variable = variable, value = value)


```
  



# testing header creation ...  working


```{r}
#| echo: false
#| results: asis

regions_values <- unique(df$region)
for (region_value in regions_values){
  cat("## ", region_value, "\n")
  dfplot <- df |> filter(region == region_value)
  p = ggplot(dfplot, aes(x=variable, y=value)) + geom_boxplot() 
  print(p)
  cat("\n\n")
} 


```

 