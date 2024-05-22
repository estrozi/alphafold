# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

ARG CUDA=11.8.0
FROM nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu18.04
# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA
#From https://stackoverflow.com/questions/44683119/dockerfile-replicate-the-host-user-uid-and-gid-to-the-image
ARG GID
ARG UNAME
ARG UID
RUN groupadd -g $GID -o racine
RUN useradd -r -u $UID -g $GID $UNAME

# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        cmake \
        cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
        git \
        hmmer \
        kalign \
        tzdata \
        wget \
        gnuplot \
        tcsh \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && pushd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make -j 4 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/bin \
    && popd \
    && rm -rf /tmp/hh-suite

# Install Miniconda package manager.
ARG MINICONDA_FILE='Miniconda3-py310_24.1.2-0-Linux-x86_64.sh'
RUN wget -q -P /tmp \
#  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    "https://repo.anaconda.com/miniconda/${MINICONDA_FILE}" \
    && bash /tmp/${MINICONDA_FILE} -b -p /opt/conda \
#    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/${MINICONDA_FILE}
#    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

# Install conda packages.
ENV PATH="/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"
#RUN conda install -qy conda==24.1.2 pip python=3.11 \
#    && conda install -y -c nvidia cuda=${CUDA_VERSION} \
#    && conda install -y -c conda-forge openmm=8.0.0 pdbfixer cuda-toolkit=${CUDA_VERSION} \
#    && conda clean --all --force-pkgs-dirs --yes \
#    && conda install cuda -c nvidia/label/cuda-11.8.0
RUN conda install -qy conda==23.5.2 \
    && conda install -y -c conda-forge \
      openmm=7.7.0 \
      cudatoolkit==${CUDA_VERSION} \
      pdbfixer \
      pip \
      python=3.10 \
      && conda install -c nvidia cuda-nvcc \
      && conda clean --all --force-pkgs-dirs --yes

COPY . /app/alphafold
RUN wget -q -P /app/alphafold/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Install pip packages.
RUN pip3 install --upgrade pip --no-cache-dir \
    && pip3 install -r /app/alphafold/requirements.txt --no-cache-dir \
    && pip3 install --upgrade --no-cache-dir \
      jax==0.4.19 \
      jaxlib==0.4.19+cuda11.cudnn86 \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html


# Add SETUID bit to the ldconfig binary so that non-root users can run it.
RUN chmod u+s /sbin/ldconfig.real

# Currently needed to avoid undefined_symbol error.
#RUN ln -sf /usr/lib/x86_64-linux-gnu/libffi.so.7 /opt/conda/lib/libffi.so.7
RUN ln -sf /usr/lib/x86_64-linux-gnu/libffi.so.8 /opt/conda/lib/libffi.so.7

# We need to run `ldconfig` first to ensure GPUs are visible, due to some quirk
# with Debian. See https://github.com/NVIDIA/nvidia-docker/issues/1399 for
# details.
# ENTRYPOINT does not support easily running multiple commands, so instead we
# write a shell script to wrap them up.
WORKDIR /app/alphafold
RUN echo $'#!/bin/bash\n\
ldconfig\n\
python /app/alphafold/run_alphafold.py "$@"' > /app/run_alphafold.sh \
  && chmod 755       /app/run_alphafold.sh
USER $UNAME
ENTRYPOINT ["/app/run_alphafold.sh"]
