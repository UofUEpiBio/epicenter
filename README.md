# bipartite ERGMs to the LFCT dataset

## Running on CHPC

The project can be executed in a CHPC environment. Currently, the `Makefile` and `render.slurm` are set up to be executed in George's CHPC account. The `render.slurm` file is intended to be used with the `sbatch` command:

```bash
sbatch render.slurm
```

Which will load singularity and execute the `singularity_render_chpc` target in the `Makefile`:

```makefile
singularity_render_chpc:
	singularity exec --bind=$(PWD):/epicenter --pwd=epicenter \
		epicenter.sif make render
```

## Singularity and container images

The `epicenter.sif` image was created using a combination of podman/docker and singularity. Because of the weight, the image is not distributed with GitHub. Eventually, it will be shared via a cloud service.

The steps to build the image are:

1. Build the container image using the `container_build` target.

2. Push the container image to a cloud service using the `container_push` target. Currently set up to be used with quay.io/gvegayon/epicenter.

3. Build the sif file using the `singularity_build` target. This can then be copied to the CHPC environment.
