YR?=2023
SOURCE:=${PWD}
# SOURCE?=~/projects/bio.snowcrab/inst/markdown
# SOURCE := $(abspath $(lastword $(MAKEFILE_LIST)) )
MEDIA?=$(SOURCE)/media
DATADIR?=~/bio.data/bio.snowcrab
WK?=~/bio.data/bio.snowcrab/assessments
FN?=FN_none
METHOD?=quarto
DOCTYPE?=beamer_presentation
# above are params that can be sent with make command
# warning: variables names include terminal spaces so match sure to keep it tidy ..
FNM := $(basename $(FN))
WKY := $(WK)/$(YR)

# If the first argument is "run"...
ifeq (run, $(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif
 
# export BSTINPUTS=.:media//:
# rmarkdown and rmarkdown docs needs to be copied to elsewhere, and media for csl's etc

test:
	@echo "$(SOURCE)"

# declare non-file targets:
.PHONY: dummy run all install clean cleanall distclean TAGS info check watch

.ONESHELL: copy # Applies to every target in the file!  .. to permit changing directory

dummy :
	@echo need to pass parameters to make
 
run : 
	@echo prog $(RUN_ARGS)  # just an example

clean :
	cd $(WKY)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl *.rmd *.qmd *.tex 

cleanall :
	cd $(WKY)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl *.docx *.odt *.epub *.html *.pdf *.rmd *.qmd *.tex *.md
 
quarto :
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/$(FNM).qmd 	
	ln -sf $(SOURCE)/media $(WKY)/media
	quarto render $(FNM).qmd --to html  
	@echo $(WKY)/$(FNM).html
	firefox $(WKY)/$(FNM).html  &
		
rmarkdown :
	$(info *****************************************************************************)
	$(info *** NOTE: for Beamer presentations, you want DOCTYPE=beamer_presentation  ***)
	$(info *** Rmarkdown document types: https://rmarkdown.rstudio.com/lesson-9.html ***)
	$(info *****************************************************************************)
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/$(FNM).rmd
	ln -sf $(SOURCE)/media $(WKY)/media
	Rscript -e "rmarkdown::render('$(WKY)/$(FNM).rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE  ), output_format='$(DOCTYPE)', output_dir='$(WKY)' ) " 
	@echo $(WKY)/$(FNM).pdf
	zathura $(WKY)/$(FNM).pdf &

latex2pdf:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).tex $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pdflatex $(FNM).tex
	biber $(FNM) 
	pdflatex $(FNM).tex
	pdflatex $(FNM).tex
	zathura $(FNM).pdf &

latex:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc -s $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FNM).tex 
	code $(FNM).tex &

pdf: 
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc -s $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FNM).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FNM).tex
	zathura $(FNM).pdf &

html:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FNM).html \
#	--mathjax \
	--embed-resources --standalone
	firefox $(FNM).html &

epub:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FNM).epub \
	--epub-cover-image=<cover-image> \
	--toc
	libreoffice $(FNM).epub &

docx:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FNM).docx \
	libreoffice $(FNM).docx &


odt:
	mkdir -p $(WKY)
	cd $(WKY)
	cp $(SOURCE)/Makefile $(WKY)/Makefile 	
	cp $(SOURCE)/$(FNM).md $(WKY)/	
	ln -sf $(SOURCE)/media $(WKY)/media
	pandoc $(FNM).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FNM).odt \
	libreoffice $(FNM).odt &

  
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


 