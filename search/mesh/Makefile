# Ensure you have typst installed:
# sudo apt install pkg-config
# cargo install --locked typst-cli

.PHONY: all render watch clean

all: render

render:
	typst compile --root .. main.typ

watch:
	typst watch --root .. main.typ

clean:
	rm -f *.pdf
