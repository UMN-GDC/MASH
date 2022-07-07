#Python base image
FROM jupyter/scipy-notebook:latest

# contains 


# Find the package versions running 
# conda list | grep numpy

# Prepare environment
RUN python3 -m pip install nibabel==3.2.2
RUN python3 -m pip install pyjanitor==0.22.0
RUN python3 -m pip install re==2.2.1
RUN python3 -m pip install json==2.0.9
RUN python3 -m pip install argparse==1.1


# Copy over test data
RUN mkdir /data
COPY Example/* /data


#load the python script and tell docker to run that script
#when someone tries to execute the container
RUN mkdir /functions
COPY functions/* /functions/
COPY AdjHE.py . 
ENTRYPOINT ["python3", "AdjHE.py"]

# Give permissions
RUN chmod 555 -R /code /data AdjHE.py


