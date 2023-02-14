# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/geospatial:4.2.2

# required
MAINTAINER Ben Marwick <benmawick@gmail.com>

WORKDIR /CTtps
COPY . /CTtps


RUN  sudo apt-get update -y \
  && sudo apt-get install -y libnlopt-dev libpoppler-cpp-dev libtesseract-dev  tesseract-ocr-eng \
  && R -e "install.packages(c('BiocManager', 'remotes'), repos = c(CRAN = 'https://cloud.r-project.org'))" \
  && R -e "remotes::install_github(c('rstudio/renv', 'quarto-dev/quarto-r'))" \
  # install pkgs we need
  && R -e "renv::restore()" \
  # run all the code
  && R -e "quarto::quarto_render('/CTtps/analysis/paper/paper.qmd')"
