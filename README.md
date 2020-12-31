vcline
======

Variant Calling Pipeline for Cancer Genome Sequencing

[![wercker status](https://app.wercker.com/status/7e550a64d29cd12d0c2eef936867f9be/s/master "wercker status")](https://app.wercker.com/project/byKey/7e550a64d29cd12d0c2eef936867f9be)

Installation
------------

```sh
$ pip install -U \
    https://github.com/dceoy/vanqc/archive/main.tar.gz \
    https://github.com/dceoy/vcline/archive/master.tar.gz
```

Dependent commands:

- pigz
- pbzip2
- bgzip
- tabix
- samtools
- bcftools
- java
- gatk
- cutadapt
- fastqc
- trim_galore
- bwa
- python3
- python2
- configManta.py
- configureStrelkaSomaticWorkflow.py
- configureStrelkaGermlineWorkflow.py
- delly
- R
- msisensor-pro
- snpEff (`java -jar /path/to/snpEff.jar`)
- vep

Docker image
------------

Pull the image from [Docker Hub](https://hub.docker.com/r/dceoy/vcline/).

```sh
$ docker image pull dceoy/vcline
```

Usage
-----

1.  Download resource data.

    ```sh
    $ vcline download-resources --dest-dir=/path/to/download/dir
    ```

2.  Create `vcline.yml`.

    ```sh
    $ vcline init
    $ vi vcline.yml     # => edit configurations
    ```

3.  Run the pipeline.

    ```sh
    $ vcline run --workers=2
    ```

Run `vcline --help` for more information.
