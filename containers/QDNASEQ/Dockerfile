FROM nfcore/base:2.1

LABEL \
  author="Jacob Munro" \
  description="Container for sv-pipe" \
  maintainer="Bahlo Lab"

# set the conda env name
ARG NAME='sv_pipe_qdnaseq'

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# set env variables
ENV TZ=Etc/UTC PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}"
