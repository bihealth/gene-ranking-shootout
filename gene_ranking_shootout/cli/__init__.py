"""CLI interface"""

import csv
import json
import pathlib
import sys
import typing

import attrs
import cattrs
import click
from loguru import logger
import numpy as np

from gene_ranking_shootout import models, runner


@click.group()
def main():
    """Main entry point for the CLI interface"""
    logger.remove()
    fmt = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | "
        "<level>{level: <6}</level> | "
        "<level>{message}</level>"
    )
    logger.add(sys.stderr, format=fmt)


@main.group()
def benchmark():
    """Group for benchmark sub commands."""


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("results_json")
def summarize(results_json, bars_top_n, total_width):
    """Summarize the results."""
    with open(results_json, "rt") as inf:
        results = cattrs.structure(json.load(inf), typing.List[models.Result])
    runner.BarPrinter(bars_top_n=bars_top_n, total_width=total_width).print(results)


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("base_url")
@click.argument("simulated_json")
@click.argument("results_json")
def varfish_phenix(base_url, simulated_json, results_json, bars_top_n, total_width):
    """Benchmark the VarFish implementation of the Phenix algorithm."""
    runner.PhenixVarFishRunner(base_url, bars_top_n=bars_top_n, total_width=total_width).run(
        simulated_json, results_json
    )


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("simulated_json")
@click.argument("results_json")
def phen2gene(simulated_json, results_json, bars_top_n, total_width):
    """Benchmark the Phen2Gene container."""
    runner.Phen2GeneRunner(bars_top_n=bars_top_n, total_width=total_width).run(
        simulated_json, results_json
    )


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("simulated_json")
@click.argument("results_json")
def amelie(simulated_json, results_json, bars_top_n, total_width):
    """Benchmark the AMELIE web server."""
    runner.AmelieRunner(bars_top_n=bars_top_n, total_width=total_width).run(
        simulated_json, results_json
    )


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("simulated_json")
@click.argument("results_json")
def cada(simulated_json, results_json, bars_top_n, total_width):
    """Benchmark the CADA container."""
    runner.CadaRunner(bars_top_n=bars_top_n, total_width=total_width).run(
        simulated_json, results_json
    )


@benchmark.command()
@click.option("--bars-top-n", default=10)
@click.option("--total-width", default=80)
@click.argument("base_url")
@click.argument("algorithm")
@click.argument("simulated_json")
@click.argument("results_json")
def exomiser(base_url, algorithm, simulated_json, results_json, bars_top_n, total_width):
    """Benchmark the Exomiser REST Prioritizer."""
    runner.ExomiserRunner(base_url, algorithm, bars_top_n=bars_top_n, total_width=total_width).run(
        simulated_json, results_json
    )


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
@click.argument("out_json")
@click.argument("datasets", nargs=-1)
@click.option("--seed", default=42)
@click.option("--case-count", default=10)
@click.option("--candidate-genes-count", default=19)
def simulate(out_json, datasets, case_count, candidate_genes_count, seed):
    """Simulate cases based on the dataset file."""
    # Load dataset and all genes from ``gnomad_counts.tsv``.
    logger.info("Loading data")
    cases = []
    skipped = 0
    seen_case_names = set()
    for dataset in datasets:
        for case in load_dataset(dataset):
            if case.name in seen_case_names:
                skipped += 1
                continue
            cases.append(case)
            seen_case_names.add(case.name)
    logger.info("... {} cases overall ({} duplicates)", len(cases), skipped)
    gnomad_counts = models.load_gnomad_counts()
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


@dataset.command()
@click.argument("tsv_in")
@click.argument("json_out")
def convert_tsv(tsv_in, json_out):
    """Convert from TSV to JSON format."""
    logger.info("Converting from {} to {}", tsv_in, json_out)
    with open(tsv_in, "rt") as inputf:
        with open(json_out, "wt") as outputf:
            cases = []
            header = ("name", "disease_omim_id", "disease_gene_id", "hpo_terms")
            reader = csv.reader(inputf, delimiter="\t")
            for row in reader:
                if not row[0].startswith("Patient"):
                    logger.info("Skipping row {}; does not start with 'Patient:'", row)
                    continue
                elif len(row) < 3:
                    logger.info("Skipping row {}; not enough columns", row)
                    continue
                if len(row) == 3:
                    data = dict(zip(header[:1] + header[2:], row))
                else:
                    data = dict(zip(header, row))
                cases.append(
                    models.Case(
                        name=str(data["name"]),
                        disease_omim_id=str(data.get("disease_omim_id", "unknown")),
                        disease_gene_id=str(data["disease_gene_id"]),
                        hpo_terms=[str(term) for term in data["hpo_terms"].split(",")],
                    )
                )
            logger.info("Writing {} cases", len(cases))
            json.dump(cattrs.unstructure(cases), outputf, indent=2)


if __name__ == "__main__":
    main()
