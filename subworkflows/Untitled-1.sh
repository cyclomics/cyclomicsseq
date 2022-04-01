nextflow run main.nf \
--input_read_dir /home/dami/data/raw_data/MAR7196 \
--output_dir /home/dami/data/nextflow/MAR7196\
--backbone_name BB22\
--backbone_fasta /home/dami/data/backbones/backbones/backbones_db_valid.fasta \
--read_pattern "fastq_pass/*.fastq" \
--reference ~/data/references/Homo_sapiens/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa


export NXF_SINGULARITY_CACHEDIR=/hpc/compgen/users/drebergen/.singularity/
export NXF_HOME=/hpc/compgen/users/drebergen
export NXF_TEMP=/hpc/compgen/users/drebergen/tmp/.nxf/
#export SINGULARITY_TEMPDIR=/hpc/compgen/users/drebergen/tmp/singularity-tmp/
#export SINGULARITY_WORKDIR=/hpc/compgen/users/drebergen/tmp/singularity-tmp/
export SINGULARITY_SCRATCH=/hpc/compgen/users/drebergen/tmp/singularity-scratch
#export SINGULARITY_TMPDIR=/hpc/compgen/users/drebergen/tmp/singularity-tmp/
export SINGULARITY_CACHEDIR=/hpc/compgen/users/drebergen/.singularity/
