ifndef ENGINE
ENGINE = podman
endif

ifdef MAC
PLATFORM = --platform linux/amd64
endif

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  container_build     : Build the container image"
	@echo "  container_run       : Run the container image"
	@echo "  mac_container_build : Build the container image on MacOS"
	@echo "  mac_container_run   : Run the container image on MacOS"
	@echo "  container_push      : Push the container image to quay.io"
	@echo "  singularity         : Build the singularity image"
	@echo "  singularity_pull    : Pull the singularity image needed to build the singularity file."
	@echo ""
	@echo "Variables:"
	@echo "  ENGINE : Container engine to use (default: podman)"
	@echo "  MAC    : Set to use podman on MacOS (default: unset)"

render:
	quarto render models/2023-08-26-bipartite-ergms_multi.qmd

container_build:
	$(ENGINE) build $(PLATFORM) -t epicenter -f ContainerFile

container_run:
	$(ENGINE) run $(PLATFORM) -it --rm \
		--mount type=bind,source=$(PWD),target=/epicenter \
		--workdir /epicenter epicenter

container_push:
	$(ENGINE) push $(PLATFORM) epicenter quay.io/gvegayon/epicenter:latest

mac_container_build:
	MAC=1 $(MAKE) container_build

mac_container_run:
	MAC=1 $(MAKE) container_run

singularity:
	$(ENGINE) run $(PLATFORM) -it --rm \
		--mount type=bind,source=$(PWD),target=/epicenter \
		--workdir /epicenter \
		quay.io/singularity/singularity:v4.1.0 build epicenter.sif docker://quay.io/gvegayon/epicenter:latest

singularity_run:
	$(ENGINE) run $(PLATFORM) -it --rm \
		--mount type=bind,source=$(PWD),target=/epicenter \
		--workdir /epicenter \
		quay.io/singularity/singularity:v4.1.0 exec epicenter.sif bash

singularity_pull:
	$(ENGINE) pull $(PLATFORM) quay.io/singularity/singularity:v4.1.0

singularity_render_chpc:
	singularity exec --bind=$(PWD):/epicenter --pwd=epicenter \
		epicenter.sif make render