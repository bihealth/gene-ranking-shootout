"""Code for running the benchmark."""

from collections import Counter
import json
import sys
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
        above_bars_top_n = len([rank for rank in ranks if rank is not None and rank > self.bars_top_n])
        missing = len([rank for rank in ranks if rank is None])
        counter = Counter(ranks)
        tot_width = self.total_width - 14
        max_value = max(list(counter.values()) + [missing, above_bars_top_n])

        def gen_bar(value):
            hash_count = int(tot_width / max_value * value)
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

    def __init__(self, *, bars_top_n=10):
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
            json.dump(
                cattrs.unstructure(results),
                outf,
                indent=2,
            )
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
            logger.error("Disease gene {} not found in resutls?", case.disease_gene_id)
            return None

        return models.Result(
            case=case,
            rank=rank,
            result_entrez_ids=result_entrez_ids,
        )
