FROM mambaorg/micromamba
# Root user is required for compatibility with nextflow execution.
USER root
RUN apt update 
RUN apt install -y procps
RUN apt clean autoclean
RUN micromamba install --yes --name base -c bioconda -c conda-forge \
      python=3.12.4 \
      samtools \
      loguru \
      click \
      pysam \
      biopython \
      bokeh \
      tqdm \
      bwa
RUN micromamba clean --all --yes
# add software to container.
RUN mkdir /app
RUN mkdir /app/plots
COPY . /app
WORKDIR /app
# Just show the help message when executing empty container cmd.
CMD python /app/cycas/cycas.py --help
