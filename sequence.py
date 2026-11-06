"""Loads the DNA sequence shared by both exercises, downloading it from NCBI if needed."""
from pathlib import Path
from urllib.request import urlopen

ACCESSION = "NM_000546.6"  # Homo sapiens tumor protein p53 (TP53) mRNA, 2512 nt
FASTA_PATH = Path(__file__).parent / "data" / "sequence.fasta"
EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text"
)
MIN_LENGTH = 1000
MAX_LENGTH = 3000


def download_fasta():
    FASTA_PATH.parent.mkdir(exist_ok=True)
    with urlopen(EFETCH_URL, timeout=30) as response:
        FASTA_PATH.write_bytes(response.read())


def load_sequence():
    """Returns (description, sequence) and checks the sequence fits the lab's limits."""
    if not FASTA_PATH.exists():
        download_fasta()

    lines = FASTA_PATH.read_text().splitlines()
    description = lines[0].lstrip(">")
    sequence = "".join(lines[1:]).upper()

    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        raise ValueError(f"Expected {MIN_LENGTH}-{MAX_LENGTH} nucleotides, got {len(sequence)}")
    if set(sequence) - set("ACGT"):
        raise ValueError(f"Sequence contains non-ACGT letters: {set(sequence) - set('ACGT')}")
    return description, sequence
