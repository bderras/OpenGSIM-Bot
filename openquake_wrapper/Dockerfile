# =======================================================================
# Dockerfile - Dockerfile for the OpenQuake Scientific Core
# =======================================================================

# 1. Use a standard Python base image.
FROM python:3.9-slim

# 2. Set the working directory inside the container.
WORKDIR /app

# 3. Install system dependencies required by OpenQuake and its libraries.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 4. Install OpenQuake and its scientific dependencies using pip.
# This is the standard installation method.
RUN pip install --no-cache-dir \
    "h5py<3.0.0" \
    openquake.engine \
    matplotlib \
    WeasyPrint

# 5. Copy the local wrapper script(s) into the container.
COPY . .

# 6. Set a default command (optional, as we will call scripts directly).
CMD ["python3", "run_gmm_calculation_v8.py"]