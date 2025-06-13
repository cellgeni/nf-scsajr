FROM r-base:4.4.1

# install system dependencies
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  liblzma-dev libbz2-dev libgsl-dev libcairo2-dev libpng-dev \
  build-essential libcurl4-openssl-dev libfftw3-dev libgfortran5 \
  libgmp3-dev libgraphviz-dev libgtk-3-dev libgtkmm-3.0-dev \
  libxml2-dev libssl-dev libunwind-dev libxt-dev pandoc libsvn1 \
  python3 python-is-python3 python3-dev python3-venv python3-pip \
  git curl default-jdk libv8-dev \
  openjdk-17-jdk openjdk-17-jre \
  cargo libmagick++-dev libhdf5-dev


# install R packages
RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  install.packages(c("devtools","randomcoloR","knitr","doMC","RcppHungarian","hdf5r","reshape"),dependencies=TRUE,upgrade="never");'

RUN R -e 'install.packages("BiocManager"); \
  options(error=function(e)quit(status=2,save="no"),warn=2); \
  BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10","BSgenome.Hsapiens.UCSC.hg38","hexbin","clusterProfiler","org.Hs.eg.db","EBImage"),upgrade="never");'

RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  install.packages(c("fitdistrplus","ggridges","ica","leiden","lmtest","pbapply","RANN","ROCR","scattermore","sctransform","spatstat.explore","spatstat.geom"));'

RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  install.packages(c("Matrix","Seurat"));' 

RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  BiocManager::install(c("rhdf5"));'

RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  devtools::install_github(c("cellgeni/visutils","mamarkevi/plotCoverage","iaaka/sajr"),upgrade = "never")'

RUN R -e 'options(error=function(e)quit(status=2,save="no"),warn=2); \
  devtools::install_github(c("cellgeni/scsajr"),upgrade = "never")'

COPY Dockerfile /docker/
RUN chmod -R 755 /docker

ENTRYPOINT []
