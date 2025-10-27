#!/bin/bash

# ==============================================================================
# Quick Test Script for the OpenGSIM-Bot Computation Layer
# ==============================================================================
# This script verifies that the OpenQuake Docker container is running and
# can correctly execute a GMM calculation via the Python wrapper.
# ==============================================================================

# --- Configuration ---
# The name of the OpenQuake service as defined in the docker-compose.yml file.
# IMPORTANT: If you change the service name in docker-compose.yml, you must update it here.
OPQUAKE_SERVICE_NAME="openquake"
# Find the running container ID for the specified service
CONTAINER_ID=$(docker-compose ps -q ${OPQUAKE_SERVICE_NAME})

# --- Script Logic ---

echo "--- Running Quick Test for OpenGSIM-Bot's Computation Layer ---"

# 1. Check if the Docker container is running.
# This is a critical prerequisite.
if [ -z "${CONTAINER_ID}" ] || [ ! "$(docker ps -q -f id=${CONTAINER_ID})" ]; then
    echo
    echo "ERROR: The OpenQuake Docker container ('${OPQUAKE_SERVICE_NAME}') is not running."
    echo "Please start all services with 'docker-compose up -d' before running this test."
    echo
    exit 1
fi

echo "Container '${OPQUAKE_SERVICE_NAME}' found. Proceeding with the test..."
echo "Executing a test calculation for AbrahamsonEtAl2014..."
echo

# 2. Execute the test command inside the running container.
# This command simulates the call that the n8n workflow would make.
docker exec ${CONTAINER_ID} python3 /app/openquake_wrapper/run_gmm_calculation_v8.py \
    --module abrahamson_2014 \
    --class_name AbrahamsonEtAl2014 \
    --keys "mag,rrup,vs30,ztor,rake,dip,width,z1pt0,vs30measured,rx,ry0" \
    --values "6.5,25,500,2,90,70,10,50,False,15,0"

# Capture the exit code of the Docker command
EXIT_CODE=$?

echo
echo "--- Test Complete ---"

# 3. Provide a final status message to the user.
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ SUCCESS: The script executed without errors."
    echo "If you see a JSON output above with '\"success\": true', the scientific core is working correctly."
else
    echo "❌ FAILED: The script returned a non-zero exit code (${EXIT_CODE})."
    echo "Please check the error messages above. The container might still be initializing or there could be an issue with the Python script."
fi
echo