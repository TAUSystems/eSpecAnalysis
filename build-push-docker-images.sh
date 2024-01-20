#!/bin/bash

# if necessary, reauthenticate:
# aws ecr get-login-password --region us-west-1 | docker login --username AWS --password-stdin 460578213037.dkr.ecr.us-west-1.amazonaws.com

docker build -t image-backend-especanalysis .
docker tag image-backend-especanalysis 460578213037.dkr.ecr.us-west-1.amazonaws.com/tau-image-processing-backend/especanalysis:latest
docker push 460578213037.dkr.ecr.us-west-1.amazonaws.com/tau-image-processing-backend/especanalysis:latest
