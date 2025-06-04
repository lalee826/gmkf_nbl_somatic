FROM rocker/tidyverse:3.6.3

#update system packages
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list
RUN apt-get update

#install R devtools
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.0.2.tar.gz', repos=NULL, type='source')"

#install version 1.0.0 of developmental version of dplyr and dependencies
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/testthat/testthat_2.1.1.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('bindrcpp', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.1.1.tar.gz')"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/vctrs/vctrs_0.3.0.tar.gz')"
RUN R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/tidyselect/tidyselect_1.1.0.tar.gz')"
RUN R -e "devtools::install_github('hadley/dplyr@v1.0.0',force=TRUE)"

#install pip for python 3
RUN apt-get update && apt-get install -y python3-pip
RUN python3 -m pip install --upgrade pip setuptools wheel
# Install JupyterLab
RUN pip3 install jupyterlab ipykernel

# Expose port and run
EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]

FROM python:3.7-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    libopenblas-dev \
    liblapack-dev \
    gfortran \
    && pip install --upgrade pip setuptools wheel

#install numpy and scipy from binary
RUN pip install --only-binary :all: numpy==1.16.6
RUN pip install --only-binary :all: scipy==1.1.0

#install pandas,matplotlib
RUN pip3 install pandas==1.1.5
RUN pip3 install matplotlib==3.1.3

#install scikit-learn from binary
RUN pip3 install --only-binary :all: scikit-learn==0.19.2
#install tensorflow
RUN pip3 install tensorflow==2.3.0


