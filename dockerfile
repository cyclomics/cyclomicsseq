FROM continuumio/miniconda3:4.11.0

WORKDIR /app
COPY environment.yml environment.yml
RUN conda env update -n base -f environment.yml

# Add code form the cycas git submodule 
COPY Cycas/cycas/ /cycas
