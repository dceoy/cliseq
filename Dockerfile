FROM ubuntu:latest AS builder

ENV DEBIAN_FRONTEND noninteractive

COPY --from=dceoy/gatk:latest /usr/local /usr/local
COPY --from=dceoy/gatk:latest /opt/conda /opt/conda
COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bwa:latest /usr/local/src/bwa /usr/local/src/bwa
COPY --from=dceoy/trim_galore:latest /usr/local/src/FastQC /usr/local/src/FastQC
COPY --from=dceoy/trim_galore:latest /usr/local/src/TrimGalore /usr/local/src/TrimGalore
COPY --from=dceoy/bcftools:latest /usr/local/src/bcftools /usr/local/src/bcftools
COPY --from=dceoy/bedtools:latest /usr/local/src/bedtools2 /usr/local/src/bedtools2
COPY --from=dceoy/manta:latest /opt/manta /opt/manta
COPY --from=dceoy/strelka:latest /opt/strelka /opt/strelka
COPY --from=dceoy/delly:latest /usr/local/bin/delly /usr/local/bin/delly
COPY --from=dceoy/canvas:latest /opt/canvas /opt/canvas
COPY --from=dceoy/msisensor:latest /usr/local/bin/msisensor /usr/local/bin/msisensor
COPY --from=dceoy/lumpy:latest /opt/lumpy-sv /opt/lumpy-sv
COPY --from=dceoy/lumpy:latest /usr/local/src/samblaster /usr/local/src/samblaster
COPY --from=dceoy/lumpy:latest /usr/local/bin/sambamba /usr/local/bin/sambamba
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD . /tmp/vcline

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        g++ gcc libbz2-dev libc-dev libcurl4-gnutls-dev libgsl-dev \
        libperl-dev liblzma-dev libssl-dev libz-dev make python \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && /opt/conda/bin/conda update -n base -c defaults conda \
      && /opt/conda/bin/conda clean -ya \
      && rm -rf /root/.cache/pip

RUN set -e \
      && find \
        /usr/local/src/bwa /usr/local/src/FastQC /usr/local/src/TrimGalore \
        /usr/local/src/samblaster \
        -maxdepth 1 -type f -executable -exec ln -s {} /usr/local/bin \; \
      && cd /usr/local/src/samtools/htslib-* \
      && make install \
      && cd /usr/local/src/samtools \
      && make install \
      && cd /usr/local/src/bcftools \
      && make install \
      && cd /usr/local/src/bedtools2 \
      && make install

RUN set -e \
      && /usr/bin/python /tmp/get-pip.py \
      && /usr/bin/python -m pip install -U --no-cache-dir \
        numpy pip pysam \
      && /opt/conda/bin/python /tmp/get-pip.py \
      && /opt/conda/bin/python -m pip install -U --no-cache-dir \
        cutadapt pip /tmp/vcline \
      && rm -f /tmp/get-pip.py

FROM ubuntu:latest

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local /usr/local
COPY --from=builder /opt /opt
COPY --from=dceoy/gatk:latest /etc/apt /etc/apt
ADD http://mirrors.edge.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb /tmp/libpng12-0.deb

RUN set -e \
      && ln -sf /bin/bash /bin/sh

RUN set -e \
      && mv /etc/apt/sources.list.d/r.list /tmp/ \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https ca-certificates curl gnupg software-properties-common \
      && sed -ne 's/^DISTRIB_RELEASE=\(.*\)$/\1/p' /etc/lsb-release \
        | xargs -i curl -SLO \
          https://packages.microsoft.com/config/ubuntu/{}/packages-microsoft-prod.deb \
      && apt-get -y install ./packages-microsoft-prod.deb /tmp/libpng12-0.deb \
      && rm -f packages-microsoft-prod.deb /tmp/libpng12-0.deb

RUN set -e \
      && mv /tmp/r.list /etc/apt/sources.list.d/ \
      && add-apt-repository universe \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils bsdmainutils ca-certificates curl \
        dotnet-runtime-2.1 gawk openjdk-8-jre libcurl3-gnutls libgsl23 \
        libncurses5 pbzip2 perl pigz python r-base wget \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
      && echo '. /opt/conda/etc/profile.d/conda.sh' >> ~/.bashrc \
      && echo 'source activate gatk' >> /usr/local/src/gatk/gatkenv.rc \
      && echo 'source /usr/local/src/gatk/gatk-completion.sh' \
        >> /usr/local/src/gatk/gatkenv.rc

ENV PYTHONPATH /opt/manta/lib/python:/opt/strelka/lib/python:${PYTHONPATH}
ENV PATH /opt/conda/envs/gatk/bin:/opt/conda/bin:/opt/manta/bin:/opt/strelka/bin:/opt/canvas/bin:/opt/lumpy-sv/bin:${PATH}
ENV BCFTOOLS_PLUGINS /usr/local/src/bcftools/plugins

ENTRYPOINT ["/opt/conda/bin/vcline"]
