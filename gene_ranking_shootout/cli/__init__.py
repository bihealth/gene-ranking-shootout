"""CLI interface"""

import csv
import json
import pathlib
import random

import attrs
import cattrs
import click
from loguru import logger

from gene_ranking_shootout import models


@click.group()
def main():
    """Main entry point for the CLI interface"""


@main.group()
def dataset():
    """Group for dataset sub commands."""


@dataset.command()
def list():
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
@click.option("--candidate-genes-count", default=20)
def simulate(dataset, out_json, case_count, candidate_genes_count, seed):
    """Simulate cases based on the dataset file."""
    # Load dataset and all genes from ``gnomad_counts.tsv``.
    logger.info("Loading data")
    cases = load_dataset(dataset)
    gnomad_counts = models.load_gnomad_counts_tsv(
        pathlib.Path(__file__).parent.parent / "data" / "gnomad_counts.tsv"
    )
    logger.info("Simulating cases")
    # Create random number generator.
    rng = random.Random(seed)
    # Pick cases.
    cases = rng.sample(cases, k=case_count)
    # For each case, pick random candidate genes.
    simulated = []
    for case in cases:
        candidate_genes = rng.sample(gnomad_counts, k=candidate_genes_count + 1)
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
