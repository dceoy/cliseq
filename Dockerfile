FROM ubuntu:18.04 AS builder

ENV DEBIAN_FRONTEND noninteractive

COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bwa:latest /usr/local/src/bwa /usr/local/src/bwa
COPY --from=dceoy/trim_galore:latest /usr/local/src/FastQC /usr/local/src/FastQC
COPY --from=dceoy/trim_galore:latest /usr/local/src/TrimGalore /usr/local/src/TrimGalore
COPY --from=dceoy/bcftools:latest /usr/local/src/bcftools /usr/local/src/bcftools
COPY --from=dceoy/bedtools:latest /usr/local/src/bedtools2 /usr/local/src/bedtools2
COPY --from=dceoy/gatk:latest /opt/gatk /opt/gatk
COPY --from=dceoy/manta:latest /opt/manta /opt/manta
COPY --from=dceoy/strelka:latest /opt/strelka /opt/strelka
COPY --from=dceoy/delly:latest /usr/local/bin/delly /usr/local/bin/delly
COPY --from=dceoy/canvas:latest /opt/canvas /opt/canvas
COPY --from=dceoy/msisensor:latest /usr/local/bin/msisensor /usr/local/bin/msisensor
COPY --from=dceoy/snpeff:latest /opt/snpEff /opt/snpEff
ADD https://bootstrap.pypa.io/get-pip.py /tmp/get-pip.py
ADD https://raw.githubusercontent.com/dceoy/print-github-tags/master/print-github-tags /usr/local/bin/print-github-tags
ADD . /tmp/vcline

RUN set -e \
      && ln -sf bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates gnupg \
      && echo "deb https://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/" \
        > /etc/apt/sources.list.d/r.list \
      && apt-key adv --keyserver keyserver.ubuntu.com \
        --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        curl g++ gcc libbz2-dev libc-dev libcurl4-gnutls-dev libfreetype6-dev \
        libgsl-dev liblzma-dev libncurses5-dev libperl-dev libpng-dev \
        libssl-dev libz-dev make pkg-config python python3.8-dev \
        python3.8-distutils r-base \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -eo pipefail \
      && /usr/bin/python3.8 /tmp/get-pip.py \
      && grep -ne '- pip:' /opt/gatk/gatkcondaenv.yml \
        | cut -d : -f 1 \
        | xargs -i sed -e '1,{}d' /opt/gatk/gatkcondaenv.yml \
        | cut -d - -f 2- \
        | sed -e 's/\(gatk.*.zip\)/\/opt\/gatk\/\1/' \
        | cut -d = -f 1 \
        | tr -d '\n' \
        | xargs /usr/bin/python3.8 -m pip install -U --no-cache-dir \
          cutadapt pip tensorflow /tmp/vcline \
      && sed -ne 's/^- anaconda::\([^=]\+\).*$/\1/p' \
        /opt/gatk/gatkcondaenv.yml \
        | tr '\n' ' ' \
        | xargs /usr/bin/python3.8 -m pip install -U --no-cache-dir \
      && rm -f /tmp/get-pip.py

RUN set -eo pipefail \
      && chmod +x /usr/local/bin/print-github-tags \
      && print-github-tags --release --tar tensorflow/tensorflow \
        | grep 'archive\/v1\.' \
        | head -1 \
        | xargs -i curl -SL {} -o /tmp/tensorflow.tar.gz \
      && tar xvf /tmp/tensorflow.tar.gz -C /tmp --remove-files \
      && mv /tmp/tensorflow-* /tmp/tensorflow \
      && /usr/bin/python3.8 -m pip install -U --no-cache-dir \
        /tmp/tensorflow/tensorflow/lite/tools/pip_package \
      && rm -rf /tmp/tensorflow

RUN set -e \
      && Rscript /opt/gatk/install_R_packages.R

RUN set -e \
      && cd /usr/local/src/bwa \
      && make clean \
      && make \
      && cd /usr/local/src/samtools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/samtools \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/bcftools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && cd /usr/local/src/bcftools \
      && make clean \
      && ./configure --enable-libgsl --enable-perl-filters \
      && make \
      && make install \
      && cd /usr/local/src/bedtools2 \
      && make clean \
      && make \
      && make install \
      && find \
        /usr/local/src/bwa /usr/local/src/FastQC /usr/local/src/TrimGalore \
        -maxdepth 1 -type f -executable -exec ln -s {} /usr/local/bin \;

FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

COPY --from=builder /usr/local /usr/local
COPY --from=builder /opt /opt
COPY --from=builder /etc/apt /etc/apt
ADD http://mirrors.edge.kernel.org/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb /tmp/libpng12-0.deb

RUN set -e \
      && ln -sf bash /bin/sh

RUN set -ea pipefail \
      && mv /etc/apt/sources.list.d/r.list /tmp/ \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        apt-transport-https apt-utils ca-certificates curl gnupg \
        software-properties-common \
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
        dotnet-runtime-2.1 libcurl3-gnutls libgsl23 libncurses5 openjdk-8-jre \
        pbzip2 perl pigz python python3.8 python3.8-distutils r-base wget \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && ln -sf python3 /usr/bin/python \
      && ln -sf python3.8 /usr/bin/python3

ENV PYTHONPATH /opt/manta/lib/python:/opt/strelka/lib/python:${PYTHONPATH}
ENV PATH /opt/gatk/bin:/opt/manta/bin:/opt/strelka/bin:/opt/canvas/bin:/opt/snpEff/bin:${PATH}
ENV BCFTOOLS_PLUGINS /usr/local/src/bcftools/plugins

ENTRYPOINT ["/usr/local/bin/vcline"]
