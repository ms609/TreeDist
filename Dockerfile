# Use the official R image from the Rocker project
FROM rocker/r-ver:4.4.0

# Install system dependencies
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

# Install R packages
RUN R -e "install.packages(c('devtools', 'remotes', 'pak', 'covr'), repos='https://cloud.r-project.org/')"

# Install latest Rcpp from GitHub
RUN R -e "remotes::install_github('RcppCore/Rcpp')"

# Install your package's dependencies
RUN R -e "pak::pkg_install(c('fGarch', 'univariateML', 'kdensity'))"

# Install TreeDistData from GitHub
RUN R -e "remotes::install_github('ms609/TreeDistData')"

# Install the remaining dependencies for your package
RUN R -e "pak::pkg_install(c('TreeDist', 'phangorn'))"

# Set the working directory
WORKDIR /usr/local/src/TreeDist

# Copy the project files into the Docker image
COPY . /usr/local/src/TreeDist/

# Install the R package (assumes the package is in the working directory)
RUN R CMD INSTALL .

# Command to run when starting the container
CMD ["R"]
