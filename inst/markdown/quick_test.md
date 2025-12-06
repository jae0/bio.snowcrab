---
title: "quick_test to try a page before adding to main document "
subtitle: "this makes debugging a page easier and faster"
metadata-files:
  - _metadata.yml
params:
  year_assessment: 2024
  year_start: 1999
  data_loc:  "~/bio.data/bio.snowcrab"
  sens: 1
  debugging: FALSE
  model_variation: logistic_discrete_historical
profile: nens

--- 


<!-- 

make quarto FN=quick_test.md DOCTYPE=revealjs  PARAMS="-P year_assessment:2024"  --directory=~/bio/bio.snowcrab/inst/markdown

-->


```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup

  require(knitr)

  knitr::opts_chunk$set(
    root.dir = data_root,
    echo = FALSE,
    out.width="6.2in",
    fig.retina = 2,
    dpi=192
  )

  # things to load into memory (in next step) via _load_results.qmd
  toget = c( "fishery_results", "fishery_model" )  

```


::: {.content-visible when-profile="nens"}

```{r}
#| label: setup-observer-data-nens
#| eval: true
#| output: false
  region = "cfanorth"
  REGION = "N-ENS"

```

:::

::: {.content-visible when-profile="sens"}

```{r}
#| label: setup-observer-data-sens
#| eval: true
#| output: false

  region = "cfasouth"
  REGION = "S-ENS"
```
:::

::: {.content-visible when-profile="4x"}

```{r}
#| label: setup-observer-data-4x
#| eval: true
#| output: false

  region = "cfa4x"
  REGION = "4X"
```

:::


{{< include _load_results.qmd >}}  


## testing profile

```{r}
#| echo: false 
#| eval: true 
#| output: true

print(REGION)

print(regions)

```

<!--

## testing header creation with *asis* ...  working

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
```


## multi-panel figure example

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
SCD = project.datadirectory("bio.snowcrab")
 
loc = file.path( SCD, "output", "maps", "survey.locations" )
years = year_assessment + c(0:-3)
fn = check_file_exists( file.path( loc, paste( "survey.locations", years, "png", sep=".") ))
include_graphics( fn )
```

Figure caption does notwork well for multiple figures ... XXX .
 


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
   
-->

 