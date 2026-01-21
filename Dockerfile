# NRXN1 Structural Variant Analysis Pipeline
# Docker image for cloud deployment

FROM python:3.10-slim

LABEL maintainer="Nagi Debbah"
LABEL description="NRXN1 Structural Variant Analysis Pipeline"
LABEL version="0.1.0"

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    wget \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install samtools and bcftools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -xjf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.17*

RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 \
    && tar -xjf bcftools-1.17.tar.bz2 \
    && cd bcftools-1.17 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.17*

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /data/input /data/output /data/reference logs

# Set entrypoint
ENTRYPOINT ["python", "main.py"]

# Default command
CMD ["--help"]
