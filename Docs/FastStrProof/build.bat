@echo off

REM Build script for making paper from .tex files.

set TEX_FILE=FastStrProof


xelatex -synctex=1 -halt-on-error -file-line-error %TEX_FILE% &&	^
bibtex8 -W -c cp1256fa %TEX_FILE% &&								^
xelatex -synctex=1 -halt-on-error -file-line-error %TEX_FILE% &&	^
xelatex -synctex=1 -halt-on-error -file-line-error %TEX_FILE% &&	^
echo Done! &&														^
start %TEX_FILE%.pdf
