vcline
======

Variant Calling Pipeline for Cancer Genome Sequencing

[![Test](https://github.com/dceoy/vcline/actions/workflows/test.yml/badge.svg)](https://github.com/dceoy/vcline/actions/workflows/test.yml)
[![Docker](https://github.com/dceoy/vcline/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/dceoy/vcline/actions/workflows/docker-publish.yml)

Installation
------------

```sh
$ pip install -U https://github.com/dceoy/vcline/archive/master.tar.gz
```

Dependent commands:

- `pigz`
- `pbzip2`
- `bgzip`
- `tabix`
- `samtools`
- `bcftools`
- `java`
- `gatk`
- `cutadapt`
- `fastqc`
- `trim_galore`
- `bwa` or `bwa-mem2`
- `python3`
- `python2`
- `configManta.py`
- `configureStrelkaSomaticWorkflow.py`
- `configureStrelkaGermlineWorkflow.py`
- `delly`
- `R`
- `msisensor-pro`
- `snpEff` (`java -jar /path/to/snpEff.jar`)
- `vep`

Docker image
------------

Pull the image from [Docker Hub](https://hub.docker.com/r/dceoy/vcline/).

```sh
$ docker image pull dceoy/vcline
```

Usage
-----

1.  Download and process resource data.

    ```sh
    $ vcline download --dest-dir=/path/to/download/dir
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
