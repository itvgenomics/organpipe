FROM continuumio/miniconda3:latest
RUN conda install -c conda-forge -c bioconda \
    mitos=2.1.10 \
    python=3.12 \
    && conda clean -afy
COPY fix_blast.py /opt/conda/lib/python3.12/site-packages/mitos/blast/__init__.py
