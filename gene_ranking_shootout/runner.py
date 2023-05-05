"""Code for running the benchmark."""

from collections import Counter
import csv
import json
import subprocess
import sys
import tempfile
import typing

import cattrs
from loguru import logger
import requests
import tqdm

from gene_ranking_shootout import models


class BarPrinter:
    """Helper for printing bars."""

    def __init__(self, *, bars_top_n=10, total_width=40, outf: typing.TextIO = sys.stdout):
        #: The number of top genes to print bars for.
        self.bars_top_n = bars_top_n
        #: The total display width.
        self.total_width = total_width
        #: The output file.
        self.outf = outf

    def print(self, results: typing.List[models.Result]):
        ranks = [result.rank for result in results]
        above_bars_top_n = len(
            [rank for rank in ranks if rank is not None and rank > self.bars_top_n]
        )
        missing = len([rank for rank in ranks if rank is None])
        counter = Counter(ranks)
        tot_width = self.total_width - 14
        max_value = len(results)

        def gen_bar(value):
            if max_value:
                hash_count = int(tot_width / max_value * value)
            else:
                hash_count = 0
            if value == 0:
                return ""
            elif hash_count == 0:
                return "."
            else:
                return "#" * hash_count

        for i in range(1, self.bars_top_n + 1):
            value = counter.get(i, 0)
            bar = gen_bar(value)
            print(f"  {i:3}: {value:>4}  {bar}", file=self.outf)
        print(file=self.outf)
        bar = gen_bar(above_bars_top_n)
        print(f"{self.bars_top_n + 1}-..: {above_bars_top_n:>4}  {bar}", file=self.outf)
        bar = gen_bar(missing)
        print(f"mssng: {missing:>4}  {bar}", file=self.outf)


class BaseRunner:
    """Base class for the runners."""

    def __init__(self, *, total_width=80, bars_top_n=10):
        #: The total display width.
        self.total_width = total_width
        #: The number of top genes to print bars for.
        self.bars_top_n = bars_top_n

        logger.info("Loading data ...")
        #: The gnomAD counts.
        self.gnomad_data = models.load_gnomad_counts()
        #: Mapping from Entrez gene ID to gene symbol.
        self.entrez_to_symbol = {gene.entrez_id: gene.gene_symbol for gene in self.gnomad_data}
        #: Mapping from gene symbol to Entrez gene ID.
        self.symbol_to_entrez = {gene.gene_symbol: gene.entrez_id for gene in self.gnomad_data}
        logger.info("... done loading data")

    def run(self, path_simulated_json: str, path_results_json: str):
        """Run the benchmark."""
        logger.info("Loading cases ...")
        cases = models.load_cases_json(path_simulated_json)
        logger.info("... done loading {} cases", len(cases))

        logger.info("Running benchmark ...")
        results = []
        for case in tqdm.tqdm(cases):
            result = self.run_ranking(case)
            if result is not None:
                results.append(result)
        logger.info("... done running benchmark")

        logger.info("Writing results ...")
        with open(path_results_json, "wt") as outf:
            json.dump(cattrs.unstructure(results), outf, indent=2)
        logger.info("... done writing results")

        logger.info("Displaying results overview ...")
        self.print_bars(results)
        logger.info("All done. Have a nice day!")

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        """Run the ranking for the given case.

        :param case: The case to run the ranking for.
        :returns: result for the case or ``None`` if no result was found.
        """
        _ = case
        raise NotImplementedError()

    def print_bars(self, results: typing.List[models.Result], outf: typing.TextIO = sys.stdout):
        """Print the bars for the results.

        :param results: The results to print the bars for.
        :param outf: The output file handle to write to.
        """
        BarPrinter(bars_top_n=self.bars_top_n, total_width=40, outf=outf).print(results)


class PhenixVarFishRunner(BaseRunner):
    """Run benchmark for phenix within ``varfish-server-worker``."""

    def __init__(self, base_url: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #: Base URL of Varfish server.
        self.base_url = base_url

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        url_terms = ",".join(case.hpo_terms)
        url_gene_symbols = ",".join(
            [self.entrez_to_symbol[gene_id] for gene_id in case.candidate_gene_ids or []]
        )
        disease_gene_symbol = self.entrez_to_symbol.get(case.disease_gene_id)
        if not disease_gene_symbol:
            logger.warning("Disease gene {} not found in gnomAD data", case.disease_gene_id)
            return None

        url = f"{self.base_url}?terms={url_terms}&gene_symbols={url_gene_symbols},{disease_gene_symbol}"
        # logger.debug("Running query: {}", url)
        result_container = requests.get(url).json()
        # Translate the gene symbols from the result to entrez ids
        result_entrez_ids = []
        for result_entry in result_container["result"]:
            result_entrez_ids.append(self.symbol_to_entrez[result_entry["gene_symbol"]])
        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            logger.error("Disease gene {} not found in results?", case.disease_gene_id)
            return None

        return models.Result(case=case, rank=rank, result_entrez_ids=result_entrez_ids)


class Phen2GeneRunner(BaseRunner):
    """Run benchmark for phen2gene.

    The benchmark will be run using podman, thus podman must be available
    on the system.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #: The name ofthe image to use via podman.
        self.image_name = "docker.io/genomicslab/phen2gene"

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        # Prepare list of all gene symbols.
        gene_symbols = [self.entrez_to_symbol[gene_id] for gene_id in case.candidate_gene_ids or []]
        gene_symbols.append(self.entrez_to_symbol[case.disease_gene_id])

        # Resulting entrez ids are collected here.
        result_entrez_ids = []

        # Run phen2gene with temporary directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            # Write terms and genes to files.
            with open(f"{tmpdir}/terms.txt", "wt") as outf:
                outf.write("\n".join(case.hpo_terms))
            with open(f"{tmpdir}/genes.txt", "wt") as outf:
                outf.write("\n".join(gene_symbols))
            # Run phen2gene using podman.
            cmd = [
                "podman",
                "run",
                "--rm",
                "-v",
                f"{tmpdir}:/code/out",
                "-t",
                self.image_name,
                "-f",
                "/code/out/terms.txt",
                "-l",
                "/code/out/genes.txt",
            ]
            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                logger.error("Error running phen2gene")
                return None

            # Read the output file.
            with open(f"{tmpdir}/output_file.associated_gene_list") as inputf:
                reader = csv.DictReader(inputf, delimiter="\t")
                # Translate the gene symbols from the result to entrez ids.
                for row in reader:
                    result_entrez_ids.append(self.symbol_to_entrez[row["Gene"]])

        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            logger.error("Disease gene {} not found in results?", case.disease_gene_id)
            return None

        return models.Result(case=case, rank=rank, result_entrez_ids=result_entrez_ids)


class AmelieRunner(BaseRunner):
    """Run benchmark with AMELIE web service."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #: URL of the AMELIE API.
        self.api_url = "https://amelie.stanford.edu/api/gene_list_api/"

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        gene_symbols = [self.entrez_to_symbol[gene_id] for gene_id in case.candidate_gene_ids or []]
        disease_gene_symbol = self.entrez_to_symbol.get(case.disease_gene_id)
        if not disease_gene_symbol:
            logger.warning("Disease gene {} not found in gnomAD data", case.disease_gene_id)
            return None
        else:
            gene_symbols.append(disease_gene_symbol)

        payload = {
            "patientName": "query",
            "phenotypes": ",".join(case.hpo_terms),
            "genes": ",".join(gene_symbols),
        }

        response = requests.post(self.api_url, data=payload)

        # Translate the gene symbols from the result to entrez ids.
        result_entrez_ids = []
        try:
            response_json = response.json()
        except json.JSONDecodeError:
            logger.error("Error decoding JSON response: {}", response.text)
            return None
        for row in response_json:
            result_entrez_ids.append(self.symbol_to_entrez[row[0]])

        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            logger.error("Disease gene {} not found in results?", case.disease_gene_id)
            return None

        return models.Result(case=case, rank=rank, result_entrez_ids=result_entrez_ids)


class CadaRunner(BaseRunner):
    """Run benchmark for CADA.

    The benchmark will be run using podman, thus podman must be available
    on the system.  Further, the custom built container image is assumed.
    See the README for details.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #: The name ofthe image to use via podman.
        self.image_name = "localhost/cada-for-shootout:latest"

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        result_entrez_ids = []

        # Run CADA with temporary directory.
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run phen2gene using podman.
            cmd = [
                "podman",
                "run",
                "--rm",
                "-v",
                f"{tmpdir}:/data",
                "-t",
                self.image_name,
                "--hpo_terms",
                ",".join(case.hpo_terms),
                "--out_dir",
                "/data",
            ]
            try:
                subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                logger.error("Error running CADA")
                return None

            # Read the output file.
            with open(f"{tmpdir}/result.txt", "rt") as inputf:
                reader = csv.DictReader(inputf, delimiter="\t")
                # Translate the gene symbols from the result to entrez ids.
                for row in reader:
                    result_entrez_ids.append(row["gene_id"])

        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            logger.error("Disease gene {} not found in results?", case.disease_gene_id)
            return None

        return models.Result(case=case, rank=rank, result_entrez_ids=result_entrez_ids)


class ExomiserRunner(BaseRunner):
    """Run benchmark with Exomiser web service."""

    def __init__(self, base_url, algorithm, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #: Algorithm parameters.
        self.algo_params = {
            "hiphive": ("hiphive", ["human", "mouse", "fish", "ppi"]),
            "hiphive-human": ("hiphive", ["human"]),
            "hiphive-mouse": ("hiphive", ["human", "mouse"]),
            "phenix": ("phenix", []),
            "phive": ("phive", []),
        }
        valid_algos = list(self.algo_params.keys())
        if algorithm not in valid_algos:
            raise ValueError(f"Invalid algorithm {algorithm}, must be one of {valid_algos}")

        #: URL of the exomiser server.
        while base_url.endswith("/"):
            base_url = base_url[:-1]
        self.base_url = base_url
        #: Algorithm to use.
        self.algorithm = algorithm

    def run_ranking(self, case: models.Case) -> typing.Optional[models.Result]:
        gene_ids = [gene_id.replace("Entrez:", "") for gene_id in case.candidate_gene_ids or []]
        gene_ids.append(case.disease_gene_id.replace("Entrez:", ""))

        prio_algorithm, prio_params = self.algo_params[self.algorithm]
        payload = {
            "prioritiser": prio_algorithm,
            "prioritiserParams": ",".join(prio_params),
            "phenotypes": case.hpo_terms,
            "genes": gene_ids,
        }

        url = f"{self.base_url}/exomiser/api/prioritise/"
        response = requests.post(url, json=payload)

        # Translate the gene symbols from the result to entrez ids.
        result_entrez_ids = []
        for entry in response.json()["results"]:
            gene_id = entry["geneId"]
            result_entrez_ids.append(f"Entrez:{gene_id}")

        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            logger.error("Disease gene {} not found in results?", case.disease_gene_id)
            return None

        return models.Result(case=case, rank=rank, result_entrez_ids=result_entrez_ids)
