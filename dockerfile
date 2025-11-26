FROM mambaorg/micromamba:2.4.0

WORKDIR /app
COPY environment.yml environment.yml

USER root
RUN apt-get update && apt-get install -y procps && apt-get clean

RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# Add code form the cycas git submodule 
COPY Cycas/ /app/Cycas
