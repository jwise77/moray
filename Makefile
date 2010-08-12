
USE_BIBTEX=yes

all: ms.pdf

# General rules on how to make postscript or PDF versions of a document

%.pdf: %.ps
	ps2pdf $*.ps $@

%.ps: %.dvi
	dvips -t letter -o $@ $*.dvi

%.dvi: %.tex
	latex $*.tex
ifeq ($(USE_BIBTEX),yes)
	-bibtex $*
endif
	latex $*.tex
	latex $*.tex
	latex $*.tex


.PHONY : clean veryclean

clean:
	rm -f *.bbl *.log *.aux *.out *.dvi *.blg *.toc *.lof *.lot *.tbx *.fgx *~

veryclean: clean
	rm -f *.ps *.pdf
