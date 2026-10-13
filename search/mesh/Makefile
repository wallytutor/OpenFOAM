# Ensure you have typst installed:
# sudo apt install pkg-config
# cargo install --locked typst-cli

.PHONY: all render watch clean

all: render

render:
	typst compile main.typ

watch:
	typst watch main.typ

clean:
	rm -f *.pdf
