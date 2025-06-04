#!/bin/bash

#start R studio server
rstudio-server start

# Start JupyterLab
jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root &

# Keep the container running
tail -f /dev/null
