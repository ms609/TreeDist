FROM rocker/r-ver:4.4.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra \
    pandoc \
    pandoc-citeproc \
    libglpk40 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('devtools', 'remotes', 'pak', 'covr'), repos='https://cloud.r-project.org/')"

RUN R -e "remotes::install_github('RcppCore/Rcpp')"

RUN R -e "pak::pkg_install(c('fGarch', 'univariateML', 'kdensity'))"

RUN R -e "remotes::install_github('ms609/TreeDistData')"

RUN R -e "pak::pkg_install(c('TreeDist', 'phangorn'))"

WORKDIR /usr/local/src/TreeDist

COPY . /usr/local/src/TreeDist/


RUN R CMD INSTALL .

CMD ["R"]
