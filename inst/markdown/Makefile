YR?=2023
SOURCE?=~/bio/bio.snowcrab/inst/markdown
MEDIA?=$(SOURCE)/media
DATADIR?=~/bio.data/bio.snowcrab
WK?=$(DATADIR)/reports
DDE?=$(DATADIR)/fishery_model/2023/size_structured_dde_normalized
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)) )
FN?=FN_none
# warning: variables names include terminal spaces ..

# If the first argument is "run"...
ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif
 
# export BSTINPUTS=.:media//:
# rmarkdown and rmarkdown docs needs to be copied to elsewhere, and media for csl's etc


# declare non-file targets:
.PHONY: dummy run all install clean cleanall distclean TAGS info check watch

.ONESHELL: copy # Applies to every targets in the file!  .. to permit chaning directory

dummy :
	@echo need to pass parameters to make

 
run : 
	@echo prog $(RUN_ARGS)

clean :
	cd $(WK)
	rm -rf *.log *.docx *.odt *.epub *.html 

cleanall :
	cd $(WK)
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl  *.docx *.odt *.epub *.html *.pdf

all : snowcrab_presentation_general_summary snowcrab_presentation_assessment snowcrab_working_paper snowcrab_sar feeding_networks_report 02.tables 51_bycatch

snowcrab_presentation_general_summary :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.rmd $(WK)/	
	@echo $(WK)/$@.pdf
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)' ), output_dir='$(WK)' ) " && zathura $(WK)/$@.pdf &

snowcrab_presentation_assessment :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.rmd $(WK)/	
	@echo $(WK)/$@.pdf
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) " && zathura $(WK)/$@.pdf &

snowcrab_working_paper :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.rmd $(WK)/	
	@echo $(WK)/$@.pdf
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) " && zathura $(WK)/$@.pdf &

snowcrab_sar :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.rmd $(WK)/	
	@echo $(WK)/$@.docx
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) " &&  libreoffice $(WK)/$@.docx  &

feeding_networks_report :
	ln -sf $(SOURCE)/media $(WK)/media
	cp $(SOURCE)/$@.qmd $(WK)/	
	cd $(WK)
	@echo $(WK)/$@.html
	quarto render $@.qmd --to html && firefox $(WK)/$@.html  &


51_bycatch :
	ln -sf $(SOURCE)/media $(WK)/media
	cd $(WK)
	@echo $(WK)/$@.html
	quarto render $@.qmd --to html && firefox $(WK)/$@.html  &


02.tables :
	ln -sf $(SOURCE)/media $(WK)/media
	cd $(WK)
	@echo $(WK)/$@.html
	quarto render $@.qmd --to html  && firefox $(WK)/$@.html  &


pdf:
	pdflatex $(FN).tex
	biber $(FN) 
	pdflatex $(FN).tex
	pdflatex $(FN).tex
	zathura $(FN).pdf &

md2latex:
	pandoc -s $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FN).tex 
	code $(FN).tex &

md2pdf: 
	pandoc -s $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FN).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FN).tex
	zathura $(FN).pdf &

md2html:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FN).html \
	--mathjax \
	--self-contained
	firefox $(FN).html &

md2epub:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FN).epub \
	--epub-cover-image=<cover-image> \
	--toc
	libreoffice $(FN).epub &

md2docx:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FN).docx \
	libreoffice $(FN).docx &


md2odt:
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

