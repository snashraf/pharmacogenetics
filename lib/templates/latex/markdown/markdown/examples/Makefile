.PHONY: all clean
AUXFILES=*.tmp *.tui *.tuo *.mp *.tuc *.markdown.in *.markdown.out \
	*.markdown.lua *.log *.aux context-mkii.tex context-mkiv.tex \
	latex-pdftex.tex latex-xetex.tex latex-luatex.tex
AUXDIRS=_markdown_*/
OUTPUT=context-mkii.pdf context-mkiv.pdf latex-pdftex.pdf \
	latex-luatex.pdf latex-xetex.pdf

# This is the default pseudo-target.
all: $(OUTPUT) clean

# This target typesets the ConTeXt example using the Mark II format.
context-mkii.pdf: context.tex
	cp context.tex context-mkii.tex
	texexec --passon="--shell-escape --8bit" context-mkii.tex

# This target typesets the ConTeXt example using the Mark IV format.
context-mkiv.pdf: context.tex
	cp context.tex context-mkiv.tex
	context context-mkiv.tex

# This target typesets the LaTeX example using the pdfTeX engine.
latex-pdftex.pdf: latex.tex
	cp latex.tex latex-pdftex.tex
	pdflatex --shell-escape latex-pdftex
	pdflatex --shell-escape latex-pdftex

# This target typesets the LaTeX example using the XeTeX engine.
latex-xetex.pdf: latex.tex
	cp latex.tex latex-xetex.tex
	xelatex --shell-escape --8bit latex-xetex
	xelatex --shell-escape --8bit latex-xetex

# This target typesets the LaTeX example using the LuaTeX engine.
latex-luatex.pdf: latex.tex
	cp latex.tex latex-luatex.tex
	lualatex latex-luatex
	lualatex latex-luatex

# This pseudo-target removes any existing auxiliary files and directories.
clean:
	rm -f $(AUXFILES)
	rm -rf $(AUXDIRS)

# This pseudo-target removes any makeable files.
implode: clean
	rm -f $(OUTPUT)
