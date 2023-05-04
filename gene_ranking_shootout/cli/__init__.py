"""CLI interface"""

from collections import Counter
import csv
import json
import pathlib

import attrs
import cattrs
import click
from loguru import logger
import numpy as np
import requests
import tqdm

from gene_ranking_shootout import models


@click.group()
def main():
    """Main entry point for the CLI interface"""


@main.group()
def benchmark():
    """Group for benchmark sub commands."""


@benchmark.command()
@click.argument("simulated_json")
@click.argument("base_url")
def varfish_phenix(simulated_json, base_url):
    """Benchmark the VarFish implementation of the Phenix algorithm."""
    logger.info("Loading data")
    gnomad_data = load_gnomad_counts()
    entrez_to_symbol = {gene.entrez_id: gene.gene_symbol for gene in gnomad_data}
    symbol_to_entrez = {gene.gene_symbol: gene.entrez_id for gene in gnomad_data}
    cases = models.load_cases_json(simulated_json)
    logger.info(" -> {} cases", len(cases))
    ranks = []
    for case in tqdm.tqdm(cases):
        url_terms = ",".join(case.hpo_terms)
        url_gene_symbols = ",".join(
            [entrez_to_symbol[gene_id] for gene_id in case.candidate_gene_ids]
        )
        disease_gene_symbol = entrez_to_symbol.get(case.disease_gene_id)
        if not disease_gene_symbol:
            logger.warning("Disease gene {} not found in gnomAD data", case.disease_gene_id)
            continue
        url = f"{base_url}?terms={url_terms}&gene_symbols={url_gene_symbols},{disease_gene_symbol}"
        # logger.debug("Running query: {}", url)
        result_container = requests.get(url).json()
        # Translate the gene symbols from the result to entrez ids
        result_entrez_ids = []
        for result_entry in result_container["result"]:
            result_entrez_ids.append(symbol_to_entrez[result_entry["gene_symbol"]])
        # Determine rank for case.
        try:
            rank = result_entrez_ids.index(case.disease_gene_id) + 1
        except ValueError:
            rank = None
        ranks.append(rank)

    above_10 = len([rank for rank in ranks if rank is not None and rank > 10])
    missing = len([rank for rank in ranks if rank is None])
    counter = Counter(ranks)
    tot_width = 40
    max_value = max(list(counter.values()) + [missing, above_10])

    def gen_bar(value):
        hash_count = int(tot_width / max_value * value)
        if value == 0:
            return ""
        elif hash_count == 0:
            return "."
        else:
            return "#" * hash_count

    for i in range(1, 11):
        value = counter.get(i, 0)
        bar = gen_bar(value)
        print(f"  {i:3}: {value:>4}  {bar}")
    print()
    bar = gen_bar(above_10)
    print(f"11-..: {above_10:>4}  {bar}")
    bar = gen_bar(missing)
    print(f"mssng: {missing:>4}  {bar}")
    ranks.sort()


@main.group()
def dataset():
    """Group for dataset sub commands."""


@dataset.command("list")
def list_():
    """Listing of available data sets."""
    # List all files ending in .json in the data directory
    data_dir = pathlib.Path(__file__).parent.parent / "data"
    for path in data_dir.glob("*.json"):
        print(path.stem)


def load_dataset(dataset):
    data_dir = pathlib.Path(__file__).parent.parent / "data"
    path = data_dir / f"{dataset}.json"
    return models.load_cases_json(path)


@dataset.command()
@click.argument("dataset")
@click.option("--count", default=10)
def head(dataset, count):
    """Display first entry in the given dataset."""
    cases = load_dataset(dataset)
    for case in cases[:count]:
        print(json.dumps(cattrs.unstructure(case)))


@dataset.command()
@click.argument("dataset")
@click.argument("out_json")
@click.option("--seed", default=42)
@click.option("--case-count", default=10)
@click.option("--candidate-genes-count", default=19)
def simulate(dataset, out_json, case_count, candidate_genes_count, seed):
    """Simulate cases based on the dataset file."""
    # Load dataset and all genes from ``gnomad_counts.tsv``.
    logger.info("Loading data")
    cases = load_dataset(dataset)
    logger.info("... {} cases overall", len(cases))
    gnomad_counts = load_gnomad_counts()
    gnomad_counts_idx = np.array(range(len(gnomad_counts)))
    gnomad_counts_values = np.array([gene.count for gene in gnomad_counts]).astype(np.float64)
    gnomad_counts_sum = np.sum(gnomad_counts_values)
    gnomad_counts_ps = np.divide(gnomad_counts_values, gnomad_counts_sum)
    logger.info("Simulating cases")
    # Create random number generator.
    rng = np.random.default_rng(seed)
    # Pick cases
    cases_idxs = np.array(range(len(cases)))
    cases_idx = rng.choice(cases_idxs, size=case_count, replace=False)
    selected_cases = [cases[idx] for idx in cases_idx]
    # For each case, pick random candidate genes.
    simulated = []
    for case in selected_cases:
        candidate_gene_idxs = rng.choice(
            gnomad_counts_idx, size=candidate_genes_count + 1, replace=False, p=gnomad_counts_ps
        ).tolist()
        candidate_genes = [gnomad_counts[idx] for idx in candidate_gene_idxs]
        if case.disease_gene_id in candidate_genes:
            candidate_genes.remove(case.disease_gene_id)
        else:
            candidate_genes.pop()
        rng.shuffle(candidate_genes)
        simulated.append(
            attrs.evolve(case, candidate_gene_ids=[gene.entrez_id for gene in candidate_genes])
        )
    # Write result to JSON file.
    with open(out_json, "wt") as f:
        json.dump(cattrs.unstructure(simulated), f, indent=2)
    logger.info("Wrote {} cases", len(simulated))


def load_gnomad_counts():
    return models.load_gnomad_counts_tsv(
        pathlib.Path(__file__).parent.parent / "data" / "gnomad_counts.tsv"
    )


@dataset.command()
@click.argument("tsv_in")
@click.argument("json_out")
def convert_tsv(tsv_in, json_out):
    """Convert from TSV to JSON format."""
    logger.info("Converting from {} to {}", tsv_in, json_out)
    with open(tsv_in, "rt") as inputf:
        with open(json_out, "wt") as outputf:
            cases = []
            header = (
                "name",
                "disease_omim_id",
                "disease_gene_id",
                "hpo_terms",
            )
            reader = csv.reader(inputf, delimiter="\t")
            for row in reader:
                data = dict(zip(header, row))
                cases.append(
                    models.Case(
                        name=str(data["name"]),
                        disease_omim_id=str(data["disease_omim_id"]),
                        disease_gene_id=str(data["disease_gene_id"]),
                        hpo_terms=[str(term) for term in data["hpo_terms"].split(",")],
                    )
                )
            logger.info("Writing {} cases", len(cases))
            json.dump(
                cattrs.unstructure(cases),
                outputf,
                indent=2,
            )


if __name__ == "__main__":
    main()
