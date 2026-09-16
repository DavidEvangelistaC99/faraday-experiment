FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# =========================
# Paquetes del sistema
# =========================
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    autoconf \
    automake \
    libtool \
    pkg-config \
    hdf5-tools \
    libhdf5-dev \
    apache2 \
    apache2-dev \
    wget \
    curl \
    ca-certificates \
    bzip2 \
    git \
    nano \
    && rm -rf /var/lib/apt/lists/*

# =========================
# Instalar Miniforge (Conda): recommended by Madrigal web page
# =========================
ENV CONDA_DIR=/opt/conda

RUN wget -O /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh

ENV PATH=${CONDA_DIR}/bin:$PATH

# =========================
# Crear entorno Conda Python 3.9
# =========================
RUN conda create -y -n schain-madrigal python=3.9 && \
    conda clean -afy

ENV PATH=${CONDA_DIR}/envs/schain-madrigal/bin:${CONDA_DIR}/bin:$PATH

# =========================
# Instalar dependencias Python
# =========================
RUN /bin/bash -c "source ${CONDA_DIR}/etc/profile.d/conda.sh && \
conda activate schain-madrigal && \
conda install -y -c conda-forge \
numpy=1.23 \
scipy=1.8.0 \
matplotlib=3.5.1 \
pandas \
cartopy \
pyzmq \
h5py \
netcdf4 \
&& pip install pip==23.3.1 && \
pip install 'setuptools<60' && \
pip install digital_rf==2.6.7 && \
pip install mysql-connector-python && \
pip install 'click<8.2' && \
pip install mod_wsgi"

# =========================
# Variables de entorno
# =========================
ENV MADROOT=/opt/madrigal

RUN mkdir -p ${MADROOT}

# Copiar código Madrigal
COPY madrigal329/ ${MADROOT}/

# Crear directorio bin para el enlace de python
RUN mkdir -p ${MADROOT}/bin

# Enlazar el Python del entorno Conda
RUN rm -f ${MADROOT}/bin/python && \
    ln -s ${CONDA_DIR}/envs/schain-madrigal/bin/python ${MADROOT}/bin/python

# Directorio de trabajo
WORKDIR ${MADROOT}

# Exponer Apache
EXPOSE 80

CMD ["/bin/bash"]