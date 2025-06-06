FROM rocker/tidyverse:3.6.0
WORKDIR /rocker-build/

#Disable GPG checks and allow insecure repositories for archived Debian snapshot
RUN echo 'Acquire::AllowInsecureRepositories "true";' > /etc/apt/apt.conf.d/90insecure \
 && echo 'Acquire::Check-Valid-Until "false";' >> /etc/apt/apt.conf.d/90insecure \
 && echo 'Acquire::AllowDowngradeToInsecureRepositories "true";' >> /etc/apt/apt.conf.d/90insecure \
 && echo 'APT::Get::AllowUnauthenticated "true";' >> /etc/apt/apt.conf.d/90insecure

#Use frozen Debian 9 snapshot for consistent package versions
RUN echo "deb http://snapshot.debian.org/archive/debian/20191101T000000Z stretch main" > /etc/apt/sources.list \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* \
 && apt-get update
 
#Install development packages without pinning specific versions
RUN apt-get install -y --allow-downgrades \
    build-essential \
    dpkg-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    perl \
    libalgorithm-merge-perl \
    libalgorithm-diff-perl \
    libdpkg-perl \
    fakeroot \
    gnupg \
    gnupg-agent \
    libfile-fcntllock-perl \
    curl \
    libx11-dev \
    libxt-dev \
    x11proto-core-dev \
    libpng-dev \
    libjpeg-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
 
#Install R devtools, dependencies, and other packages
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.0.2.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/glue/glue_1.3.2.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/testthat/testthat_2.1.1.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('bindrcpp', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('lifecycle', repos='http://cran.rstudio.com/')" \
 && R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.1.1.tar.gz')" \
 && R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/vctrs/vctrs_0.3.0.tar.gz')" \
 && R -e "devtools::install_url('https://cran.r-project.org/src/contrib/Archive/tidyselect/tidyselect_1.1.0.tar.gz')" \
 && R -e "devtools::install_github('hadley/dplyr@v1.0.0', force=TRUE)" \
 && R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" \
 && R -e "BiocManager::install(version='3.10')" \
 && R -e "install.packages('wordcloud', repos='https://cloud.r-project.org')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lsei/lsei_1.2-0.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_0.55.0.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('futile.logger', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('snow', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('registry', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/bibtex/bibtex_0.4.0.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/pkgmaker/pkgmaker_0.27.tar.gz')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/rngtools/rngtools_1.2.4.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('gridBase', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('foreach', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/NMF/NMF_0.21.0.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Formula/Formula_1.2-3.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('png', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/jpeg/jpeg_0.1-8.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-29.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('acepack', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('gridExtra', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('checkmate', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('htmlwidgets', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/htmlTable/htmlTable_1.12.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('viridis', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.1-1.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('dichromat', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/cpp11/cpp11_0.2.0.tar.gz', repos=NULL, type='source')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_0.8.3.tar.gz', repos=NULL, type='source')" \
 && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" \
 && R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/deconstructSigs/deconstructSigs_1.8.0.tar.gz', repos=NULL, type='source')" \
 && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" \
 && R -e "BiocManager::install('maftools')" \
 && R -e "BiocManager::install('BiocParallel')" \
 && R -e "BiocManager::install('DelayedArray')" \
 && R -e "BiocManager::install('SummarizedExperiment')" \
 && R -e "BiocManager::install('Rsamtools')" \
 && R -e "BiocManager::install('rtracklayer')" \ 
 && R -e "BiocManager::install('BSgenome')" \    
 && R -e "BiocManager::install('GenomicFeatures')" \      
 && R -e "BiocManager::install('VariantAnnotation')" \      
 && R -e "BiocManager::install('NMF')" \      
 && R -e "BiocManager::install('biovizBase')" \
 && R -e "install.packages('shinythemes', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('R.utils', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('extremevalues', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "install.packages('gplots', dependencies=TRUE, repos='http://cran.rstudio.com/')" \
 && R -e "BiocManager::install('GSVA', version='3.10', ask=FALSE)"

#install pip for python 3 and other dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    libffi-dev
RUN python3 -m pip install --upgrade pip setuptools wheel
#Install JupyterLab
RUN pip3 install jupyterlab ipykernel

#Install Python libraries from binary wheels where possible
RUN pip3 install \
    --only-binary :all: numpy==1.16.6 \
    --only-binary :all: scipy==1.1.0 \
    --only-binary :all: scikit-learn==0.19.2 \
    pandas==0.25.3 \
    matplotlib==2.2.5 \
    https://github.com/termcolor/termcolor/archive/refs/tags/1.1.0.zip \
    git+https://github.com/serge-sans-paille/gast.git@0.2.2 \
    tensorflow==1.15.0 \
    jupyterlab ipykernel

#Add custom startup script
COPY start_services.sh /usr/local/bin/start_services.sh
RUN chmod +x /usr/local/bin/start_services.sh

#Expose both RStudio and Jupyter ports
EXPOSE 8787
EXPOSE 8888

#Set password for RStudio user
RUN echo "rstudio:rstudio" | chpasswd

#start services
CMD ["/usr/local/bin/start_services.sh"]
