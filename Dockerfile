FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /app

RUN apt-get update && apt-get install -y \
    wget curl python3.11 python3-pip \
    r-base r-base-dev \
    hmmer prodigal && \
    apt-get clean

RUN pip3 install padloc

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY . .

RUN mkdir -p data/genomes \
             data/padloc_output \
             data/phage_panel

RUN padloc --db-update

EXPOSE 8000

CMD ["python3", "-m", "uvicorn", "app:app", \
     "--host", "0.0.0.0", "--port", "8000"]
