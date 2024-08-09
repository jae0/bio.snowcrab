YR?=2023
SOURCE := ${PWD}
# SOURCE?=~/projects/bio.snowcrab/inst/markdown
# SOURCE := $(abspath $(lastword $(MAKEFILE_LIST)) )
MEDIA?=$(SOURCE)/media
DATADIR?=~/bio.data/bio.snowcrab
WK?=$(DATADIR)/reports
DDE?=$(DATADIR)/fishery_model/$(YR)/size_structured_dde_normalized
FN?=FN_none
METHOD=quarto
# above are params that can be sent with make command
# warning: variables names include terminal spaces so match sure to keep it tidy ..

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
	cd $(WK)
	rm -rf *.log *.docx *.odt *.epub *.html 

cleanall :
	cd $(WK)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl  *.docx *.odt *.epub *.html *.pdf

quick_test :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.md $(WK)/	
	@echo $(WK)/$@.pdf
	Rscript -e "rmarkdown::render('$(WK)/$@.md', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) " && zathura $(WK)/$@.pdf &
 
quarto :
	ln -sf $(SOURCE)/media $(WK)/media
	cd $(WK)
	@echo $(WK)/$(FN).html
	cp $(SOURCE)/$(FN).md $(WK)/$(FN).qmd 	
	quarto render $(FN).qmd --to html  
	firefox $(WK)/$(FN).html  &
		
rmarkdown :
	ln -sf $(SOURCE)/media $(WK)/media
	cd $(WK)
	@echo $(WK)/$(FN).pdf
	cp $(SOURCE)/$(FN).md $(WK)/$(FN).rmd
	Rscript -e "rmarkdown::render('$(WK)/$(FN).rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) " 
	zathura $(WK)/$(FN).pdf &

latex2pdf:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).tex $(WK)/	
	cd $(WK)
	pdflatex $(FN).tex
	biber $(FN) 
	pdflatex $(FN).tex
	pdflatex $(FN).tex
	zathura $(FN).pdf &

latex:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc -s $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FN).tex 
	code $(FN).tex &

pdf: 
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc -s $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FN).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FN).tex
	zathura $(FN).pdf &

html:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FN).html \
#	--mathjax \
	--embed-resources --standalone
	firefox $(FN).html &

epub:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FN).epub \
	--epub-cover-image=<cover-image> \
	--toc
	libreoffice $(FN).epub &

docx:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FN).docx \
	libreoffice $(FN).docx &


odt:
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$(FN).md $(WK)/	
	cd $(WK)
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FN).odt \
	libreoffice $(FN).odt &

  
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


 