FROM mambaorg/micromamba:2.4.0

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

USER root
RUN apt-get update && apt-get install -y procps && apt-get clean

RUN micromamba install -y -n base -f /tmp/environment.yml && \
        micromamba clean --all --yes

ENV MAMBA_DOCKERFILE_ACTIVATE=1

RUN ln -s /opt/conda/bin/python /usr/bin/python

# Add code form the cycas git submodule 
COPY Cycas/ /app/Cycas

RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid
