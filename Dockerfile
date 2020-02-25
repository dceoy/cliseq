FROM ubuntu:latest AS builder

ENV DEBIAN_FRONTEND noninteractive

COPY --from=dceoy/gatk:latest /usr/local /usr/local
COPY --from=dceoy/gatk:latest /opt/conda /opt/conda
COPY --from=dceoy/samtools:latest /usr/local/src/htslib /usr/local/src/htslib
COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bcftools:latest /usr/local/src/bcftools /usr/local/src/bcftools
COPY --from=dceoy/bwa:latest /usr/local/src/bwa /usr/local/src/bwa
COPY --from=dceoy/trim_galore:latest /usr/local/src/FastQC /usr/local/src/FastQC
COPY --from=dceoy/trim_galore:latest /usr/local/src/TrimGalore /usr/local/src/TrimGalore
ADD . /tmp/vcline

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        gcc libbz2-dev libc-dev libcurl4-gnutls-dev libgsl-dev libperl-dev \
        liblzma-dev libssl-dev libz-dev make \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && ln -s /usr/local/src/bwa/bwa /usr/local/bin \
      && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin \
      && ln -s /usr/local/src/TrimGalore/trim_galore /usr/local/bin \
      && cd /usr/local/src/htslib \
      && make install \
      && cd /usr/local/src/samtools \
      && make install \
      && cd /usr/local/src/bcftools \
      && make install

RUN set -e \
      && /opt/conda/envs/gatk/bin/pip install -U --no-cache-dir scipy \
      && /opt/conda/bin/python3 -m pip install -U --no-cache-dir \
        cutadapt /tmp/vcline \
        https://github.com/dceoy/shoper/archive/master.tar.gz

FROM ubuntu:latest

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local /usr/local
COPY --from=builder /opt /opt

RUN set -e \
      && ln -sf /bin/bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates curl openjdk-8-jre \
        libcurl3-gnutls libncurses5 pbzip2 perl pigz python r-base \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
      && echo '. /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc \
      && echo 'conda activate base' >> ~/.bashrc

ENV PATH /opt/conda/envs/gatk/bin:/opt/conda/bin:${PATH}

ENTRYPOINT ["/opt/conda/bin/vcline"]
