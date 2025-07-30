# Base image with Python
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy requirements first (better Docker cache; rebuilds only if requirements change)
COPY requirements.txt .

# Install system dependencies + Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    python3-pip \
 && pip install --no-cache-dir -r requirements.txt \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Copy the rest of the project
COPY . .

# Default command
CMD ["python", "scripts/mergeMiRPrecursors.py", "--help"]
