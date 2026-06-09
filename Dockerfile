FROM condaforge/miniforge3:latest

WORKDIR /app

RUN conda install -y -c conda-forge -c bioconda \
    padloc python=3.11 && \
    conda clean -afy

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .

RUN mkdir -p data/genomes \
             data/padloc_output \
             data/phage_panel

RUN padloc --db-update

EXPOSE 8000

CMD ["python", "-m", "uvicorn", "app:app", \
     "--host", "0.0.0.0", "--port", "8000"]