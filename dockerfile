# The build-stage image:
FROM continuumio/miniconda3 AS build

# created using:
# https://pythonspeed.com/articles/conda-docker-image-size/
# https://pythonspeed.com/articles/activate-virtualenv-dockerfile/

# Install the package as normal:
COPY environment.yml .
RUN conda env create -f environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n cycas -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM python:3.9-slim-bullseye AS runtime
# Nextflow needs ps, part of procps
RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# add the app
COPY cycas /app
RUN mkdir -p /app/plots
# When image is run, run the code with the environment
# activated:
SHELL ["/bin/bash", "-c"]
ENV VIRTUAL_ENV=/venv/
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

CMD python /app/cycas.py --help