#!/bin/bash
# Container and image base names
CONTAINER_BASE_NAME="scpacker"
IMAGE_BASE_NAME="scpacker_modification"

declare -A DOCKER_PARAMS
# launch parameters for every environment
DOCKER_PARAMS["localhost"]="--shm-size=4G"

declare -A VOLUMES
# volume setup for every environment
VOLUMES["localhost"]="-v /home/gluck/test_md:/data"

declare -A COMMANDS
# available commands list
# if passed, overrides command in Dockerfile's CMD
COMMANDS["bash"]="bash"