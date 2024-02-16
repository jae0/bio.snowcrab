YR?=2023
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))  # loc of this Makefile
WK?=~/bio.data/bio.snowcrab/reports
MEDIA?=$(WK)/media
FN?=FN_none
SOURCE?=~/bio/bio.snowcrab/inst/markdown
DDE?=/home/jae/bio.data/bio.snowcrab/fishery_model/2023/size_structured_dde_normalized

# If the first argument is "run"...
ifeq (run,$(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

# export BSTINPUTS=.:media//:

# declare non-file targets:
.PHONY: dummy run all install clean cleanall distclean TAGS info check watch

.ONESHELL: copy # Applies to every targets in the file!  .. to permit chaning directory

dummy :
	@echo need to pass parameters to make

copy:
	cp $(SOURCE)/Makefile $(WK)/
	cp -rf $(SOURCE)/*.md $(WK)/	
	cp -rf $(SOURCE)/*.qmd $(WK)/	
	cp -rf $(SOURCE)/*.rmd $(WK)/	
	@echo file have been copied to $(WK). Change to that directory and continue with reports;

run : 
	@echo prog $(RUN_ARGS)

clean :
	rm -rf *.log *.docx *.odt *.epub *.html 

cleanall :
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *.spl  *.docx *.odt *.epub *.html *.pdf

all : snowcrab_presentation_general_summary snowcrab_working_paper snowcrab_sar feeding_networks_report 02.tables 51_bycatch

snowcrab_presentation_general_summary : 
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) "
	zathura $@.pdf

snowcrab_presentation_assessment : 
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) "
	zathura $@.pdf

snowcrab_working_paper : 
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) "
	zathura $@.pdf

snowcrab_sar :
	Rscript -e "rmarkdown::render('$(WK)/$@.rmd', params=list( year.assessment=$(YR), media_loc='$(MEDIA)', debugging=FALSE , loc_dde='$(DDE)'  ), output_dir='$(WK)' ) "
	# zathura $@.pdf
	libreoffice $@.docx

feeding_networks_report :
	quarto render feeding_networks_report.qmd --to html
	firefox feeding_networks_report.html


51_bycatch :
	quarto render 51_bycatch.qmd --to html
	firefox 51_bycatch.html


02.tables :
	quarto render 02.tables.qmd --to html
	firefox 02.tables.html


pdf: latex  
	pdflatex $(FN).tex
	biber $(FN) 
	pdflatex $(FN).tex
	pdflatex $(FN).tex

view: 
	zathura $(FN).pdf 

md2pdf: 
	pandoc -s $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FN).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FN).tex

md2html:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FN).html \
	--mathjax \
	--self-contained

md2epub:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FN).epub \
	--epub-cover-image=<cover-image> \
	--toc

md2docx:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FN).docx \


md2odt:
	pandoc $(FN).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FN).odt \

  
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

