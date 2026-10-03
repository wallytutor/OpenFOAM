export QUARTO_PYTHON="$(PWD)/.venv/bin/python"

# _name_ of Containerfile._name_:
name := foam

# Name the project based on the current user:
project := ${name}-$(shell whoami)

##############################################################################
# Basics
##############################################################################

.PHONY: all
all: render

.PHONY: sync
sync:
	uv sync

##############################################################################
# Documentation
##############################################################################

.PHONY: render
render: sync
	quarto render

.PHONY: publish
publish: sync
	quarto publish gh-pages --no-prompt --no-browser

##############################################################################
# Distributable image
##############################################################################

# Build the container image:
.PHONY: container
container: Containerfile
	podman build -t "${project}" -f Containerfile .

# Dump container to portable .tar file:
$(project).tar: container
	podman save -o "$(project).tar" "localhost/${project}"

# Convert container into apptainer:
$(project).sif: $(project).tar
	apptainer build "$(project).sif" "docker-archive://$(project).tar"

.PHONY: image
image: $(project).sif

##############################################################################
# Cleaning
##############################################################################

.PHONY: clean
clean:
	rm -rf _book/

.PHONY: dist-clean
dist-clean: clean
	rm -rf .quarto/

##############################################################################
# EOF
##############################################################################