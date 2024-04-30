FROM ubuntu:latest

# Update and install essential packages:
RUN apt-get update && apt-get install -y wget ca-certificates bash tar
# Install Miniconda:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh \
    && chmod +x /miniconda.sh \
    && /miniconda.sh -b -p /opt/miniconda \
    && rm /miniconda.sh

# Set path:
ENV PATH="/opt/miniconda/bin:$PATH"

COPY alpaca.yml /app/alpaca.yml
COPY run_all_figures.sh /app/run_all_figures.sh
COPY bin /app/bin
COPY input /app/input

# Change working directory:
WORKDIR /app
RUN wget https://zenodo.org/records/11060928/files/alpaca_data.tar.gz && tar -xzvf alpaca_data.tar.gz
# Create the Conda environment:
RUN conda env create -f alpaca.yml

# set commands:
SHELL ["conda", "run", "-n", "alpaca", "/bin/bash", "-c"]
CMD ["conda", "run", "--no-capture-output", "-n", "alpaca", "python3"]

