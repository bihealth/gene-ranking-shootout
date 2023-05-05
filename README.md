# Gene Ranking Shootout

A benchmark for methods that rank genes according to their relevance for a given phenotype (list of HPO terms).

## Methods

The following methods are currently included in the benchmark:

- AMELIE (via web service)
- CADA (via custom Docker/Podman image)
- Phen2Gene (via official Docker/Podman image)
- Phenix algorithm (as implemented in VarFish)

## Installation

Simply install with `pip` (probably inside a conda environment or virtualenv):

```bash
$ git clone https://github.com/bihealth/gene-ranking-shootout.git
$ cd gene-ranking-shootout
$ pip install -e .
$ gene-ranking-shootout --help
```

## Usage

List datasets in the benchmark:

```bash
$ gene-ranking-shootout dataset list
cada_cases_test
```

Show the first entry in a dataset:

```bash
$ gene-ranking-shootout dataset head cada_cases_test --count 2
{"name": "Patient:SCV000281758", "disease_omim_id": "OMIM:617360", "disease_gene_id": "Entrez:8621", "hpo_terms": ["HP:0001508"], "candidate_gene_ids": null}
{"name": "Patient:SCV000864231", "disease_omim_id": "OMIM:132900", "disease_gene_id": "Entrez:4629", "hpo_terms": ["HP:0011499", "HP:0000021"], "candidate_gene_ids": null}
```

Simulate cases based on a dataset.
This will pick a number of cases from the dataset.
Further, it will pick another number of random genes based on the number of rare variants (freq below 0.1% in gnomAD genomes).
The results are written into a JSON file with the cases.

```bash
$ gene-ranking-shootout dataset simulate cada_cases_test \
    /tmp/cases.json \
    --case-count 100 \
    --candidate-genes-count 100
```

You can then run the benchmark on the cases with the different methods:

```
$ gene-ranking-shootout benchmark
Usage: gene-ranking-shootout benchmark [OPTIONS] COMMAND [ARGS]...

  Group for benchmark sub commands.

Options:
  --help  Show this message and exit.

Commands:
  amelie          Benchmark the AMELIE web server.
  phen2gene       Benchmark the Phen2Gene container.
  summarize       Summarize the results.
  varfish-phenix  Benchmark the VarFish implementation of the Phenix...
$ gene-ranking-shootout benchmark amelie /tmp/cases.json /tmp/result-amelie.json
$ gene-ranking-shootout benchmark phen2gene /tmp/cases.json /tmp/result-phen2gene.json
$ gene-ranking-shootout benchmark varfish-phenix http://127.0.0.1:8081/hpo/sim/term-gene /tmp/cases.json /tmp/result-varfish-phenix.json
$ gene-ranking-shootout benchmark cada /tmp/cases.json /tmp/result-cada.json
```

You can also visualize the details of the benchmark results for each result file.

```bash
$ for f in /tmp/result-*.json; do (set -x; gene-ranking-shootout benchmark summarize $f); echo; done
+ gene-ranking-shootout benchmark summarize /tmp/result-amelie.json
    1:   48  ################################
    2:   17  ###########
    3:    7  ####
    4:    3  ##
    5:    6  ####
    6:    3  ##
    7:    2  #
    8:    3  ##
    9:    1  .
   10:    0  

11-..:    9  ######
mssng:    0  

+ gene-ranking-shootout benchmark summarize /tmp/result-cada.json
    1:   17  ###########
    2:    6  ###
    3:    2  #
    4:    4  ##
    5:    3  #
    6:    5  ###
    7:    3  #
    8:    0  
    9:    0  
   10:    2  #

11-..:   58  ######################################
mssng:    0  

+ gene-ranking-shootout benchmark summarize /tmp/result-phen2gene.json
    1:   57  #####################################
    2:    6  ###
    3:    4  ##
    4:    2  #
    5:    2  #
    6:    3  #
    7:    4  ##
    8:    0  
    9:    3  #
   10:    0  

11-..:   19  ############
mssng:    0  
+ gene-ranking-shootout benchmark summarize /tmp/result-varfish-phenix.json
    1:   39  #########################
    2:   10  ######
    3:    7  ####
    4:    4  ##
    5:    2  #
    6:    2  #
    7:    1  .
    8:    1  .
    9:    5  ###
   10:    4  ##

11-..:   25  ################
mssng:    0  
```

## Building CADA Podman Image

There is no public REST API or docker image for CADA (yet).
Here is how to build the needed CADA Podman image:

```bash
# cd docker/cada
# bash build.sh
...
# podman run -it --rm localhost/cada-for-shootout:latest --help
```

## Running Phenix in VarFish

Send an email to the author to get a copy of the necessary data.
Then, run the following in the background.

```bash
$ varfish-server-worker server pheno --path-hpo-dir path/to/varfish-server-worker-db/hpo
```
