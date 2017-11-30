#include ../../setenv.mk
SHELL := cmd

all: pdf

chm:
	taskkill /FI "windowtitle eq nmblas"
	doxygen
	copy html\*.chm .


pdf: chm
	cd latex &&	call make.bat
	copy .\latex\refman.pdf .\nmblas.pdf

clean: 
	-rm -f -r -v html latex nmblas.chm nmblas.pdf *.html
	
open:
	www.bat
	
	#http://rc-module.github.io/nmblas
	
	#pdf
	#cp ./latex/refman.pdf ./html/nmblas.pdf
	#copy nmblas.html html\nmblas.html
	#nmblas.html 
github: clean chm 
	cd html &&\
	git init . && git checkout --orphan gh-pages &&\
	git add * &&\
	git commit -m "documentation-for-gh-pages" &&\
	git remote add origin https://github.com/RC-MODULE/nmblas.git  &&\
	git push -u -f origin gh-pages	


nmblas.html:	
	lua summary.lua
	xsltproc -o nmblas.html nmblas.xsl nmblas.xml
	
git-summary: nmblas.html 
	copy nmblas.html html\nmblas.html
	cd html &&\
	git add * &&\
	git commit -m "documentation-for-gh-pages" &&\
	git push -u -f origin gh-pages	
	
	