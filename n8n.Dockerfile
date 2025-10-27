# =======================================================================
# n8n.Dockerfile - Dockerfile for the n8n Orchestration Service
# =======================================================================

# 1. Use the official n8n base image.
# This is much more efficient and secure than building from a generic Node image.
FROM n8nio/n8n

# 2. Set user to root to install system packages.
USER root

# 3. Install necessary system dependencies.
# We only need 'python3' and 'python3-pip' if we were to run Python scripts
# directly within n8n, but it's good practice to have them.
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# 4. Switch back to the non-root 'node' user for security.
USER node