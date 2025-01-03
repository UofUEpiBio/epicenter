ifndef ENGINE
ENGINE = podman
endif


container_build:
	$(ENGINE) build -t epicenter -f ContainerFile

container_run:
	$(ENGINE) run -it --rm --mount type=bind,source=$(PWD),target=/epicenter \
		--workdir /epicenter epicenter

singularity_build:
	singularity build --fakeroot epicenter ContainerFile