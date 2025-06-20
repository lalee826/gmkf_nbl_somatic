# gmkf_nbl_somatic

This is the data and analysis repository for the somatic analyses of the Gabriella Miller Kids Foundation neuroblastoma cohort

### Directories

Data - All necessary data files and tables to perform analyses. Many of these are results files from the analyses that can be reproduced using the scripts provided.  
Docker - Dockerfiles used to create an image that contains all necessary software with appropriate versions to replicate analyses.  
Scripts - All necessary scripts to perform analyses in R and Python.  
Supplemental - Supplementary data tables, figures, and analyses from the manuscript  


### Usage
The easiest way to replicate the analyses presented in the manuscript and additional studies included in this repository is to through Docker.  A dockerfile has been included that can be used to create an image on a local machine. Opening a subsequent container from this image includes all necessary softwares to perform analyses in this repository.  

Installation instructions for standalone Docker engine (recommended for Linux only) can be found here: https://docs.docker.com/engine/  
Installation instructions for Docker Desktop (Windows/Mac/Linux) can be found here: https://docs.docker.com/desktop/  

After installation, navigate to the directory where the Dockerfile and helper shell script are saved and create a Docker image with the following command  
 
 <code> docker build -t (image name) </code>  

 Alternatively, you can pull the docker image from dockerhub 

 <code> docker pull laleeupenn/gmkf:latest </code> 

Once completed, open a Docker container with port 8787 allowing usage of rstudio and port 8888 allowing jupyter notebook  
 
 <code> docker run --rm -ti -e PASSWORD=none -p 8787:8787 -p 8888:8888 (image name) </code>  

You can view the list of open containers and their names with  
 
 <code> docker container list </code>  

When a container is created, navigate to it in a different window using the command:  
 
 <code> docker exec -it (container name) /bin/bash </code>  

Once inside the container, clone this repository  
 
 <code> git clone https://www.github.com/lalee826/gmkf_nbl_somatic </code>  

Rstudio can be opened on your local machine by navigating to localhost:8787 in your web browser. The username/password will be rstudio/rstudio.


Jupyter Notebook can be opened on your local machine by navigating to localhost:8888 in your web browser. You will need to copy/paste the token that is shown in your terminal after opening the Docker container.  

### Raw Data

Raw data is not stored in this repository. The files found in this repository's raw data folder are example files to demonstrate usage of the scripts. The full raw data can be obtained through the Kids First Data Resource Portal (https://portal.kidsfirstdrc.org/).  

Please note that not all raw data files used in analyses from the manuscript can be found in the portal, for example, raw output from Sequenza, Shatterseek, GridSS, and others.
