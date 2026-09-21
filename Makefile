export QUARTO_PYTHON=$(PWD)/.venv/bin/python

.PHONY: all sync render publish clean dist-clean

all: render

sync:
	uv sync

render: sync
	quarto render

publish: sync
	quarto publish gh-pages --no-prompt --no-browser

clean:
	rm -rf _book/

dist-clean: clean
	rm -rf .quarto/