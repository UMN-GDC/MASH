#Python base image
FROM jupyter/scipy-notebook:latest


# Find the package versions running 
# conda list | grep numpy

# Prepare environment
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install pyjanitor==0.22.0
RUN python3 -m pip install regex==2021.8.3
RUN python3 -m pip install argparse==1.1

# Copy over test data
RUN mkdir -p ~/Example
RUN mkdir -p ~/functions


COPY Example/* ~/Example/
COPY functions/* ~/functions/

#load the python script and tell docker to run that script
#when someone tries to execute the container
COPY AdjHE.py .
ENTRYPOINT ["python3", "AdjHE.py"]
