"""Data models"""

import csv
import json
import pathlib
import typing

import attrs


@attrs.frozen()
class Case:
    """A case for benchmarking."""

    #: The case's name.
    name: str
    #: The OMIM code of the disease, e.g., ``"OMIM:251110"``.
    disease_omim_id: str
    #: The disease gene ID, e.g., ```"Entrez:1301"`.
    disease_gene_id: str
    #: The HPO terms for describing the patient, e.g. ``["HP:0001263", "HP:0000256",
    #: "HP:0008414"].
    hpo_terms: typing.List[str]
    #: Candidate gene IDs (e.g., ``["Entrez:1301", "Entrez:1302"]``).
    candidate_gene_ids: typing.Optional[typing.List[str]] = None


def load_cases_json(path):
    """Load ``Case`` objects from JSON file."""
    with open(path, "rt") as f:
        return [Case(**case) for case in json.load(f)]


@attrs.frozen()
class GnomadCounts:
    """Per-gene count data for gnomAD"""

    #: HGNC gene symbol.
    gene_symbol: str
    #: HGNC gene ID.
    hgnc_id: str
    #: Entrez gene ID.
    entrez_id: str
    #: Rare variant count in gnomAD.
    count: int


def load_gnomad_counts_tsv(path):
    """Load ``GnomadCounts`` objects from TSV file (with header)."""
    with open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [GnomadCounts(**row) for row in reader]


def load_gnomad_counts():
    """Helper to load the gnomad counts from the data directory."""
    return load_gnomad_counts_tsv(pathlib.Path(__file__).parent / "data" / "gnomad_counts.tsv")


@attrs.frozen()
class Result:
    """Stores prioritization results."""

    #: The case that was run.
    case: Case
    #: The rank that the disease gene was assigned.
    rank: int
    #: The resulting ranked genes as Entrez IDs.
    result_entrez_ids: typing.List[str]
