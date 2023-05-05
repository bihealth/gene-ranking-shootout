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
cada_clinvar_cases
cada_cases_test
cada_cases_validate
cada_all_cases
cada_collaborator_cases
cada_cases_train
```

Show the first entry in a dataset:

```bash
$ gene-ranking-shootout dataset head cada_cases_test --count 2
{"name": "Patient:SCV000281758", "disease_omim_id": "OMIM:617360", "disease_gene_id": "Entrez:8621", "hpo_terms": ["HP:0001508"], "candidate_gene_ids": null}
{"name": "Patient:SCV000864231", "disease_omim_id": "OMIM:132900", "disease_gene_id": "Entrez:4629", "hpo_terms": ["HP:0011499", "HP:0000021"], "candidate_gene_ids": null}
```

Simulate cases based on one or more datasets.
This will pick a number of cases from the datasets.
Further, it will pick another number of random genes based on the number of rare variants (freq below 0.1% in gnomAD genomes).
The results are written into a JSON file with the cases.
The simulation is randomized with a fixed seed that can be adjusted on the command line if necessary.

```bash
$ gene-ranking-shootout dataset simulate \
    /tmp/cases.json \
    $(gene-ranking-shootout dataset list) \
    --case-count 4714 \
    --candidate-genes-count 100
2023-05-05 11:16:35 | INFO   | Loading data
2023-05-05 11:16:35 | INFO   | ... 4714 cases overall (9428 duplicates)
2023-05-05 11:16:35 | INFO   | Simulating cases
2023-05-05 11:16:37 | INFO   | Wrote 4714 cases
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
$ gene-ranking-shootout benchmark exomiser http://localhost:8081/ phenix /tmp/cases.json /tmp/result-exomiser-phenix.json
$ gene-ranking-shootout benchmark exomiser http://localhost:8081/ phive /tmp/cases.json /tmp/result-exomiser-phive.json
$ gene-ranking-shootout benchmark exomiser http://localhost:8081/ hiphive /tmp/cases.json /tmp/result-exomiser-hiphive.json
$ gene-ranking-shootout benchmark exomiser http://localhost:8081/ hiphive-mouse /tmp/cases.json /tmp/result-exomiser-hiphive-mouse.json
$ gene-ranking-shootout benchmark exomiser http://localhost:8081/ hiphive-human /tmp/cases.json /tmp/result-exomiser-hiphive-human.json
```

You can also visualize the details of the benchmark results for each result file (below for 100 cases).

```bash
$ gene-ranking-shootout benchmark summarize /tmp/result-amelie.json
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

## Running Exomiser

The following are more rough notes than a full manual.
This uses the current Exomiser RES API version 13.2.0 (current at: 2023-05-05).
You will need approximately 75GB of storage for download and extraction and afterwards 49GB.
Probably one could get rid of a lot of the variant-specific data but I did not go into detail here.

```bash
$ wget https://github.com/exomiser/Exomiser/releases/download/13.2.0/exomiser-rest-prioritiser-13.2.0.jar
$ wget https://data.monarchinitiative.org/exomiser/latest/2302_phenotype.zip
$ wget https://data.monarchinitiative.org/exomiser/latest/2302_hg19.zip
$ unzip 2302_phenotype.zip
$ unzip 2302_hg19.zip
$ cat <<EOF > application.properties
exomiser.data-directory=$PWD
exomiser.hg19.data-version=1909
exomiser.phenotype.data-version=2302
exomiser.phenotype.random-walk-file-name=rw_string_10.mv
EOF
$ java -Xmx6G -Xms2G -Dserver.address=0.0.0.0 -Dserver.port=8081 -jar exomiser-rest-prioritiser-13.2.0.jar
```

## Datasets

The following datasets are included at the moment:

- `cada_cases_test.json` - converte from CADA's `cases_test.tsv`
- `cada_cases_train.json` - converte from CADA's `cases_train.tsv`
- `cada_cases_validate.json` - converte from CADA's `cases_validate.tsv`
- `cada_clinvar_cases.json` - converte from CADA's `clinvar_cases.tsv`
- `cada_collaborator_cases.json` - converte from CADA's `collaborator_cases.tsv`

You can conver TSV files with the following structure with `gene-ranking-shootout dataset convert-tsv`.

- Column 1: name for the case; must start with `Patient:` or is ignored.
- Column 2: disease_omim_id; as `OMIM:123456` or `unknown`
- Column 3: disease_gene_id; as `Entrez:123456`
- Column 4: hpo_terms; as semicolon-separated list, e.g., `HP:0001234;HP:0005678`

If a row has less than 4 columns, we assume that column 2 is missing.
All further columns are ignored.
The file should not have a header.
You can find some files in the CADA repository here:

- https://github.com/Chengyao-Peng/CADA/tree/main/src/CADA

The call to `gene-ranking-shootout dataset convert-tsv` should be as follows.

```bash
$ gene-ranking-shootout dataset convert-tsv input.tsv output.json
```

## Some Preliminary Results

The following was generated on 2023/05/05 with all 4714 cases.

```
$ for f in /tmp/result-*.json; do (set -x; gene-ranking-shootout benchmark summarize --bars-top-n 20 $f); echo; done
TODO: MISSING - CADA
TODO: MISSING - AMELIE

+ gene-ranking-shootout benchmark summarize --bars-top-n 20 /tmp/result-phen2gene.json
    1: 2426  ##################################
    2:  470  ######
    3:  209  ##
    4:  125  #
    5:  101  #
    6:   67  .
    7:   51  .
    8:   62  .
    9:   53  .
   10:   41  .
   11:   33  .
   12:   37  .
   13:   42  .
   14:   33  .
   15:   34  .
   16:   28  .
   17:   18  .
   18:   29  .
   19:   17  .
   20:   19  .

21-..:  763  ##########
mssng:    0

+ gene-ranking-shootout benchmark summarize --bars-top-n 20 /tmp/result-varfish-phenix.json
    1: 1709  #######################
    2:  616  ########
    3:  357  ####
    4:  277  ###
    5:  184  ##
    6:  152  ##
    7:  131  #
    8:  118  #
    9:  105  #
   10:   78  #
   11:   71  .
   12:   57  .
   13:   67  .
   14:   64  .
   15:   67  .
   16:   71  .
   17:   56  .
   18:   48  .
   19:   34  .
   20:   48  .

21-..:  403  #####
mssng:    0
```
