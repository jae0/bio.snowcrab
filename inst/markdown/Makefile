YR?=2024
DATADIR?=~/bio.data/bio.snowcrab
WK?=$(DATADIR)/assessments
FN?=FN_none
METHOD?=quarto
DOCTYPE?=html
PARAMS?="-P year_assessment:$(YR) -P model_variation:logistic_discrete_historical" # for quarto
# PARAMS?="year_assessment=$(YR), media:w_loc=$(DATADIR)/media, debugging=FALSE"   # for Rmarkdown

# warning: file, dir, and year variables names can include terminal or starting spaces so strip them  ..
$DATADIR=(strip $(DATADIR))
$FN=(strip $(FN))
$WK=(strip $(WK))
$YR=(strip $(YR))

# above are params that can be sent with make command
SOURCEDIR := $(realpath $(dir $(abspath $(MAKEFILE_LIST)) ))
FNM := $(basename $(notdir $(FN)))

WKY := $(WK)/$(YR)

# If the first argument is "run"...
ifeq (run, $(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

 
test:
	@echo "$(SOURCEDIR)"
	@echo "$(DATADIR)"
	@echo "$(WKY)"

# declare non-file targets:
.PHONY: dummy run all install clean cleanall distclean TAGS info check watch

.ONESHELL: copy # Applies to every target in the file!  .. to permit changing directory

dummy :
	@echo need to pass parameters to make
 
run : 
	@echo prog $(RUN_ARGS)  # just an example

clean :
	cd $(WKY)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl  

cleanall :
	cd $(WKY)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl *.epub *.html *.pdf  *.bib *.csl

preparefiles : 
	mkdir -p $(WKY)
	cp -f $(SOURCEDIR)/Makefile $(WKY)/
	cp -f $(SOURCEDIR)/*.yml $(WKY)/
	cp -f $(SOURCEDIR)/*.bib $(WKY)/
	cp -f $(SOURCEDIR)/*.csl $(WKY)/
	cp -f $(SOURCEDIR)/$(FN) $(WKY)/

quarto : preparefiles
	cd $(WKY)
	cp -f $(WKY)/$(FN) $(WKY)/$(FNM).qmd 	
	quarto render $(FNM).qmd --to $(DOCTYPE) $(PARAMS)
	ls $(WKY)/$(FNM).*
	
rmarkdown : preparefiles
	cd $(WKY)
	cp -f $(WKY)/$(FN) $(WKY)/$(FNM).rmd
	Rscript -e "rmarkdown::render('$(WKY)/$(FNM).rmd', params=list($(PARAMS)), output_format=$(DOCTYPE), output_dir=$(WKY))" 	
	ls $(WKY)/$(FNM).*
	
latex : preparefiles
	cd $(WKY)
	pandoc -s $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FNM).tex 
	ls $(WKY)/$(FNM).*
	
	
latex2pdf : latex
	cd $(WKY)
	pdflatex $(FNM).tex
	biber $(FNM) 
	pdflatex $(FNM).tex
	pdflatex $(FNM).tex
	ls $(WKY)/$(FNM).*
	
	
pdf : preparefiles
	cd $(WKY)
	pandoc -s $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FNM).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FNM).tex
	ls $(WKY)/$(FNM).*
	
	
html : preparefiles
	cd $(WKY)
	pandoc $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FNM).html \
#	--mathjax \
	--embed-resources --standalone
	ls $(WKY)/$(FNM).*
	
	
epub : preparefiles
	cd $(WKY)
	pandoc $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FNM).epub \
	--epub-cover-image=<cover-image> \
	--toc
	ls $(WKY)/$(FNM).*
	
	
docx : preparefiles
	cd $(WKY)
	pandoc $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FNM).docx \
	ls $(WKY)/$(FNM).*
	
	
odt : preparefiles
	cd $(WKY)
	pandoc $(FN) \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FNM).odt \
	ls $(WKY)/$(FNM).*
	
	
git:
	git commit -m"update"
	git checkout master
	git merge develop
	git push
	git pull
	git checkout develop
	git merge master
	git status

 
watch :
	while true; do make -q || make; sleep 0.5; done


 