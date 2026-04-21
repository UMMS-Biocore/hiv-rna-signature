FROM rocker/r-ver:4.3.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    gfortran \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('Seurat','harmony','Matrix','dplyr','tibble','tidyr','purrr','stringr','readr','readxl','ggplot2','patchwork','cowplot','viridis','scales','ggvenn'), repos='https://cloud.r-project.org')"

WORKDIR /workspace
COPY . /workspace
RUN chmod +x /workspace/scripts/run_pipeline.sh

ENTRYPOINT ["/workspace/scripts/run_pipeline.sh"]
CMD ["all"]
