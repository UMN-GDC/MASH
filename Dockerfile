# syntax=docker/dockerfile:1

#Python base image
FROM python:latest 

# Set working direcotry (also creates directory)
WORKDIR ~/app

# Install GCTA and plink
RUN curl -L https://zzz.bwh.harvard.edu/plink/dist/plink-1.08-x86_64.zip ; unzip plink-1.08-x86_64.zip ; rm plink-1.08-x86_64.zip
RUN curl -L https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip -o gcta.zip ; unzip gcta.zip ; rm gcta.zip

# Find the package versions running 
# conda list | grep numpy

# Prepare environment
RUN python3 -m pip install numpy==1.21.2 pandas==1.3.3 pandas-plink==2.2.9
RUN python3 -m pip install scipy==1.7.1 scikit-learn==0.24.2
RUN python3 -m pip install matplotlib==3.4.3 seaborn==0.11.2 
RUN python3 -m pip install pytest==6.2.5
RUN python3 -m pip install pyjanitor==0.22.0
RUN python3 -m pip install statsmodels==0.13.1
RUN python3 -m pip install nibabel==3.2.2 regex==2021.8.3 argparse==1.1 

# Copy over test data
RUN mkdir Estimate
RUN mkdir Simulate

# copy files over to image
COPY Estimate Estimate 
COPY Simulate Simulate 

#load the python script and tell docker to run that script
#when someone tries to execute the container
ENTRYPOINT ["python3", "Estimate.py"]
