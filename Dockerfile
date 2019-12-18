FROM ubuntu:latest AS builder

ENV DEBIAN_FRONTEND noninteractive

COPY --from=dceoy/gatk:latest /usr/local /usr/local
COPY --from=dceoy/samtools:latest /usr/local/src/htslib /usr/local/src/htslib
COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bwa:latest /usr/local/src/bwa /usr/local/src/bwa
COPY --from=dceoy/trim_galore:latest /usr/local/src/TrimGalore /usr/local/src/TrimGalore
COPY --from=dceoy/fastqc:latest /usr/local/src/FastQC /usr/local/src/FastQC
COPY --from=dceoy/picard:latest /usr/local/src/picard /usr/local/src/picard
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        gcc make python3.7-dev python3.7-distutils \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && cd /usr/local/src/htslib \
      && make install \
      && cd /usr/local/src/samtools \
      && make install

RUN set -e \
      && ln -s /usr/local/src/bwa/bwa /usr/local/bin \
      && ln -s /usr/local/src/TrimGalore/trim_galore /usr/local/bin \
      && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin \
      && echo '/usr/bin/java -jar /usr/local/src/picard/picard.jar' \
        > /usr/local/bin/picard \
      && chmod +x /usr/local/bin/picard

RUN set -e \
      && /usr/bin/python3.7 /tmp/get-pip.py \
      && pip install -U --no-cache-dir apache-airflow cutadapt pip \
      && rm -f /tmp/get-pip.py

FROM ubuntu:latest

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local /usr/local

RUN set -e \
      && ln -sf /bin/bash /bin/sh \
      && ln -sf /usr/bin/python3.7 /usr/bin/python3

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates curl default-jre \
        libcurl3-gnutls libncurses5 perl pigz python python3.7 \
        python3.7-distutils r-base \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/bin/bash", "-c"]
