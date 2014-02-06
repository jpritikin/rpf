PDFS=flexmirt-plots.pdf

all: $(PDFS) 

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary
	-rm -r figure cache

%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit2pdf('$*.Rnw')"

%.R: %.Rnw
	$(R_HOME)/bin/Rscript -e "library(knitr); knit('$*.Rnw', tangle=TRUE)"
