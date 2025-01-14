---
title: "quick_test to try a page before adding to main document "
subtitle: "this makes debugging a page easier and faster"
author:
  - name: 
      given: Snow Crab Unit
      family: DFO Science
    # orcid: 0000-0003-3632-5723 
    # email: jae.choi@dfo-mpo.gc.ca
    # email: choi.jae.seok@gmail.com
    # corresponding: true
    affiliation: 
      - name: Bedford Institute of Oceanography, Fisheries and Oceans Canada
        city: Dartmouth
        state: NS
        # url: www.bio.gc.ca
# date: "`r format(Sys.time(), '%d %B, %Y')`"
date: last-modified
date-format: "YYYY-MM-D"
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
  year_assessment: 2024
  year_start: 1999
  media_loc: "media"
  debugging: FALSE
  loc_dde: ""
--- 


<!-- Preamble

This is a Markdown document ... To create HTML or PDF, etc, run: 


# for presentations to PDF (via beamer):
# note: section separation with '#' can confuse rmarkdown
  
  make rmarkdown FN=quick_test YR=2024 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=beamer_presentation DOCEXTENSION=pdf # {via Rmarkdown}
 
  make rmarkdown FN=quick_test YR=2024 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCTYPE=html_document DOCEXTENSION=html # {via Rmarkdown}

# for html documents including presentations:
  make quarto FN=quick_test YR=2024 SOURCE=~/projects/bio.snowcrab/inst/markdown WK=~/bio.data/bio.snowcrab/assessments  DOCEXTENSION=pdf  PARAMS="-P year_assessment:2024"  # {via Quarto}


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

year_assessment = params$year_assessment
year_previous = year_assessment - 1

p = bio.snowcrab::load.environment( year.assessment=year_assessment )  

SCD = project.datadirectory("bio.snowcrab")
media_loc = project.codedirectory("bio.snowcrab", "inst", "markdown", "media")

require(gt)  # table formatting

outtabledir = file.path( p$annual.results, "tables" )

years = as.character(1996: year_assessment)

lregions = list( region=c("cfanorth", "cfasouth", "cfa4x") )
nregions = length(regions[[1]] )
reg_labels = c("N-ENS", "S-ENS", "CFA 4X")  # formatted for label

regions = unlist(lregions)
nregions = length(regions)

FD = fishery_data( regions=regions)  # mass in tonnes
fda = FD$summary_annual

library(tidyverse)
region = c(replicate(100,"cfanorth"), replicate(100,"cfasouth"), replicate(100,"cfa4x"))
variable = c(replicate(50,"var_a"),replicate(50,"var_b"),replicate(50,"var_a"),replicate(50,"var_b"),replicate(50,"var_a"),replicate(50,"var_b"))
value = runif(300, min = 0, max = 100)
value2 = runif(300, min = -1, max = 1)
value3 = rnorm(300, 0, 1)
df = data.frame(region = region, variable = variable, value = value)


```
  



# testing header creation with *asis* ...  working
```{r}
#| echo: false
#| results: asis

for (reg in regions){
  cat("## ", reg, "\n")
  dfplot <- df |> filter(region == reg)
  p = ggplot(dfplot, aes(x=variable, y=value)) + geom_boxplot() 
  print(p)
  cat("\n\n")
} 



#### multi-panel figure example
```{r}
#| eval: true
#| echo: false 
#| output: true
#| label: map-XXX 
#| fig-cap: "XXXXX"
#| fig-dpi: 144
#| fig-height: 4
#| fig.show: hold
#| layout-ncol: 2
  
loc = file.path( SCD, "output", "maps", "survey.locations" )
years = year_assessment + c(0:-3)
fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```

Figure caption does notwork well for multiple figures ... XXX .

$~$



#### multi-panel table example
```{r}
#| eval: true
#| echo: false
#| output: true
#| results: asis
#| label: table-XXX
#| tbl-cap: "Table XXX"
#| layout-ncol: 2

for (r in 1:nregions) {
  reg = regions[r]
  REG = reg_labels[r]
  cat("#### ", REG, "\n")
  oo = df[ region==reg, .(
    Nstations = .N, 
    Nmale = sum(value, na.rm=TRUE) ,
    Nfemale = sum(value2, na.rm=TRUE) 
    ), by=.(yr)]
  oo$Total = oo$Nmale + oo$Nfemale 
  names(oo) = c("Year", "No. stations", "No. male", "No. female", "No. male mature", "No. female mature", "No. total"  )
  oo = oo[order(Year), ]

  out = gt::gt(oo) |> gt::tab_options(table.font.size = 10, data_row.padding = gt::px(1), 
    summary_row.padding = gt::px(1), grand_summary_row.padding = gt::px(1), 
    footnotes.padding = gt::px(1), source_notes.padding = gt::px(1), 
    row_group.padding = gt::px(1))
  print(out)
  cat("\n\n")
} 
 
```



$~$
   


```

 