all: pdf

pdf: article.pdf

article.pdf: article.tex
	pdflatex article
	bibtex article
	pdflatex article
	pdflatex article

clean:
	rm -rf article.tex article.aux article.bbl article.blg article.log article.out article-*.pdf
