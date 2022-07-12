# syntax=docker/dockerfile:1

#Python base image
FROM python:latest 

# Set working direcotry (also creates directory)
WORKDIR ~/app

# Find the package versions running 
# conda list | grep numpy

# Prepare environment
RUN python3 -m pip install numpy==1.21.2
RUN python3 -m pip install pandas==1.3.3
RUN python3 -m pip install pyjanitor==0.22.0
RUN python3 -m pip install statsmodels==0.13.1
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install pyjanitor==0.22.0
RUN python3 -m pip install regex==2021.8.3
RUN python3 -m pip install argparse==1.1

# Copy over test data
RUN mkdir Example
RUN mkdir functions

# copy files over to image
COPY Example Example
COPY functions functions

#load the python script and tell docker to run that script
#when someone tries to execute the container
COPY AdjHE.py .
ENTRYPOINT ["python3", "AdjHE.py"]
