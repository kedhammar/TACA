FROM python:3.10

# Update pip to latest version
RUN python -m pip install --upgrade pip


# Install Nextflow dependencies
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y git \
    && apt-get install -y curl


RUN git clone https://github.com/SciLifeLab/flowcell_parser.git

RUN cd flowcell_parser \
    && python -m pip install -r requirements.txt \
    && python setup.py install

# Install dependencies
COPY requirements.txt requirements.txt
COPY requirements-dev.txt requirements-dev.txt
RUN python -m pip install -r requirements-dev.txt

RUN mkdir /root/.taca/
COPY tests/data/taca_test_cfg.yaml /root/.taca/taca.yaml