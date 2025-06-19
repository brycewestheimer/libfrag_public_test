FROM continuumio/anaconda3:5.0.1

# install sys requirements
RUN apt-get update
RUN apt-get install -y build-essential git

# clone in repo
RUN git clone https://github.com/brycewestheimer/libfrag

# update conda
RUN conda update --yes conda

# go into repo
RUN cd libfrag && conda env create -f libfrag-dev-requirements.yml

# activate
RUN echo "source activate libfrag-dev-requirements" > ~/.bash
ENV PATH /opt/conda/envs/libfrag-dev-requirements/bin:$PATH

# activate
#RUN conda activate libfrag-dev-requirements

# make build dir
RUN cd libfrag && mkdir build

# run cmake
RUN  cd libfrag/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=$(conda info --base)/envs/libfrag-dev-requirements

# build and install project
RUN cd libfrag/build && \
    make -j2 install