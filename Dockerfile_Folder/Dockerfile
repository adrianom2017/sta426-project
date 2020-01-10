#Tutorial: http://ropenscilabs.github.io/r-docker-tutorial/
#Build this docker file with the Zheng4eq.filt.wPCA.rds file available in the current directory and run: docker build -t sta426-project -f Dockerfile .
#Run docker image by docker run -e PASSWORD=docker -p 8787:8787 -v sta426_output:/home/rstudio/output sta426-project
#Access Rstudio in browser with port 8787 (mac: localhost/8787, windows: IP/8787)

#Start from rocker/rstudio image version 3.6.1
FROM rocker/rstudio:3.6.1
MAINTAINER Adriano Martinelli adrianom@student.ethz.ch

ENV LD_LIBRARY_PATH /usr/local/lib/R/lib

RUN apt-get update --fix-missing && \
    apt-get install -yq --no-install-recommends \
	libz-dev \
	libxml2-dev \
        bzip2 \
        ca-certificates \
        curl \
        gnupg \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        vim \
        wget && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update -y && apt-get install -y libpng-dev

WORKDIR /tmp

#Install packages, -e makes sure build fails if package cannot be installed
COPY install.R /tmp/install.R
RUN R -f /tmp/install.R

#Add data folder
WORKDIR ../home/rstudio/data
ADD week13_SCE_clustered.rds .

#Add script
WORKDIR ../code
ADD Final_Report.Rmd .

#Prepare output
WORKDIR ../output

#Run script
WORKDIR ..
#Rscript -e 'rmarkdown::render("code/Final_Report.Rmd")'