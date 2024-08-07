import dataclasses
from collections import namedtuple
from typing import NamedTuple
from typing import List
import csv
import json

CDS_COLOR = "#4b61d1"
RBS_COLOR = "#d14b4b"
PROMOTER_COLOR = "#2ab717"
TERMINATOR_COLOR = "#ff0000"


@dataclasses.dataclass
class CryptResults:
    """Cryptkeeper results class"""

    name: str
    sequence: str
    translation_sites: List[NamedTuple]
    rho_dep_terminators: List[NamedTuple]
    int_terminators: List[NamedTuple]
    promoters: List[NamedTuple]
    annotations: List[NamedTuple]
    burden: float

    def to_csv(self, output_path: str):
        """Export results to csv files"""
        to_csv(self, output_path)

    def to_summary(self, output_path: str):
        """Export results to summary file"""
        to_summary(self, output_path)

    def to_json(self, output_path: str):
        """Export results to json file"""
        to_json(self, output_path)


def to_csv(results: CryptResults, output_path: str):
    """Write results to csv files"""

    def _csv_writer():
        if not dataset:
            return
        with open(writepath, "w", encoding="utf8") as file:
            header = dataset[0]._fields
            writer = csv.writer(file)
            writer.writerow(header)
            for data in dataset:
                writer.writerow(data._asdict().values())

    # Write RBS predictions to CSV
    # Headers: Position, Strand, Start_codon, score1, score2
    # @TODO: Make scores descriptive
    dataset = results.translation_sites
    writepath = output_path + "_rbs.csv"
    _csv_writer()

    # Rho Independent Termination to CSV
    # Headers: Start, end, strand, coef, hairpin score, tail score
    # old system also has hairpin upstream, open, loop, close, and tail
    dataset = results.int_terminators
    writepath = output_path + "_rit.csv"
    _csv_writer()

    # Promoter Sites:
    # Headers: Score, strand, TTSpos, box35pos, box35seq, box10pos, box10seq
    dataset = results.promoters
    writepath = output_path + "_tss.csv"
    _csv_writer()

    # Also need to output Rho Dependent Termination to CSV, which is new
    dataset = results.rho_dep_terminators
    writepath = output_path + "_rdt.csv"
    _csv_writer()


def to_summary(results: CryptResults, output_path: str):
    """Write summary of results to file"""
    with open(output_path, "w") as file:
        file.write("CryptKeeper Results Summary\n")
        file.write("Total Burden: " + str(results.burden) + "\n")
        file.write("Predicted Promoters: " + str(len(results.promoters)) + "\n")
        file.write(
            "Predicted Translation: " + str(len(results.translation_sites)) + "\n"
        )
        file.write(
            "Predicted Rho Independent Terminators: "
            + str(len(results.int_terminators))
            + "\n"
        )
        file.write(
            "Total Rho Dependent Terminators: "
            + str(len(results.rho_dep_terminators))
            + "\n"
        )


def to_json(results: CryptResults, output_path: str):
    """Write results to json file"""
    with open(output_path, "w") as file:
        json.dump(dataclasses.asdict(results), file, indent=4)


def from_json(input_path: str) -> CryptResults:
    """Read results from json file and load them into a CryptResults object"""
    with open(input_path, "r") as file:
        # @TODO: Add sanity checking to json file
        data = json.load(file)
        results = CryptResults(**data)

    # Correct the data type from the json
    expressed_orf = namedtuple(
        "orf_result", "start, end, expression, burden, dG, array, start_codon, strand"
    )
    promocalc_hit = namedtuple(
        "promoter_calculator_result",
        "seq score strand TSSpos box35pos box35seq box10pos box10seq",
    )
    transterm_hit = namedtuple(
        "transterm_result",
        "start end strand conf hairpin_score tail_score seq_upstream seq_hairpin_open seq_hairpin_loop seq_hairpin_close seq_tail",
    )
    rtp_hit = namedtuple(
        "rdtresult",
        """strand,
                                        c_over_g,
                                        start_rut,
                                        end_rut,
                                        term_seq,
                                        downstream_seq,
                                        palindromes,
                                        pause_concensus,
                                        scores""",
    )
    annotations = namedtuple(
        "feature", ["name", "strand", "start", "end", "color", "nest_level"]
    )

    for i, hit in enumerate(results.translation_sites):
        results.translation_sites[i] = expressed_orf(*hit)
    for i, hit in enumerate(results.promoters):
        results.promoters[i] = promocalc_hit(*hit)
    for i, hit in enumerate(results.int_terminators):
        results.int_terminators[i] = transterm_hit(*hit)
    for i, hit in enumerate(results.rho_dep_terminators):
        results.rho_dep_terminators[i] = rtp_hit(*hit)
    for i, hit in enumerate(results.annotations):
        results.annotations[i] = annotations(*hit)

    return results


# Blank EOF line
