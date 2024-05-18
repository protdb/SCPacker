#!/bin/bash
# Container and image base names
CONTAINER_BASE_NAME="scpacker"
IMAGE_BASE_NAME="scpacker_modification"

declare -A DOCKER_PARAMS
# launch parameters for every environment
DOCKER_PARAMS["localhost"]="--shm-size=4G"
DOCKER_PARAMS["standalone"]="--shm-size=4G"
DOCKER_PARAMS["test"]="--shm-size=4G"

declare -A VOLUMES
# volume setup for every environment
VOLUMES["localhost"]="-v /home/gluck/modif_data:/data"
VOLUMES["standalone"]="-v /home/gluck/ptm_test_1:/data"
VOLUMES["test"]="-v /root/workers/test_data:/data"

declare -A COMMANDS
# available commands list
# if passed, overrides command in Dockerfile's CMD