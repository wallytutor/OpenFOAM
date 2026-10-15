# Name of container:
name := cantera

# Name the project based on the current user:
project := ${name}-$(shell whoami)

##############################################################################
# Basics
##############################################################################

.PHONY: all
all: image

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
	rm -rf $(project).tar

.PHONY: dist-clean
dist-clean: clean
	rm -rf $(project).sif

##############################################################################
# EOF
##############################################################################