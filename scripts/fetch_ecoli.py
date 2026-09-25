"""Fetch E. coli K-12 MG1655 CDS (NCBI RefSeq FTP) + PaxDb protein abundances."""
import gzip, urllib.request, pathlib

OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "ecoli"
OUT.mkdir(parents=True, exist_ok=True)

cds_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz"
print("downloading CDS...")
with urllib.request.urlopen(cds_url, timeout=120) as r:
    (OUT / "mg1655_cds.fna.gz").write_bytes(r.read())
n = gzip.open(OUT / "mg1655_cds.fna.gz", "rt").read().count(">")
print("CDS records:", n)

pax_url = "https://pax-db.org/downloads/4.2/datasets/paxdb-abundance-files-v4.2/511145/511145-WHOLE_ORGANISM-integrated.txt"
print("downloading PaxDb abundances...")
try:
    with urllib.request.urlopen(pax_url, timeout=120) as r:
        data = r.read()
    (OUT / "paxdb_abundance.tsv").write_bytes(data)
    print("paxdb bytes:", len(data))
except Exception as e:
    print("PAXDB_FETCH_FAILED:", e)
