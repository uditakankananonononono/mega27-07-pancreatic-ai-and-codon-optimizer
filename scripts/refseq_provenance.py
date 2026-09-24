"""RefSeq provenance audit of the ICOR benchmark wild-type sequences (item 07).

The MFE audit and the ICOR head-to-head both score the benchmark's wild-type DNA.
This script verifies those sequences against the CURRENT NCBI RefSeq record named in
each fasta header (efetch rettype=fasta_cds_na -> the joined CDS of the live record),
so version drift, SNPs, or mislabelled transcripts surface instead of being assumed away.

Output: results/refseq_provenance.json. Records cached in data/refseq/.
"""
import json, os, re, time, urllib.request, urllib.parse, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DNA = glob.glob(os.path.join(ROOT, 'data/icor/*/benchmark_sequences/dna'))[0]
CACHE = os.path.join(ROOT, 'data/refseq'); os.makedirs(CACHE, exist_ok=True)
ACC = re.compile(r'^>?([A-Z]{1,3}_?\d+\.\d+)\s')


def read_fasta(p):
    h, s = None, []
    for line in open(p):
        line = line.rstrip()
        if line.startswith('>'):
            h = line
        else:
            s.append(line)
    return h, ''.join(s).upper()


def fetch_cds(acc):
    f = os.path.join(CACHE, acc + '.fasta')
    if not os.path.exists(f):
        url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
               + acc + '&rettype=fasta_cds_na&retmode=text')
        for att in range(4):
            try:
                txt = urllib.request.urlopen(url, timeout=60).read().decode()
                if txt.startswith('>'):
                    break
            except Exception:
                txt = None
            time.sleep(2 * (att + 1))
        if not txt or not txt.startswith('>'):
            return None
        open(f, 'w').write(txt)
        time.sleep(0.4)
    lines = open(f).read().splitlines()
    seqs, cur = [], ''
    for l in lines:
        if l.startswith('>'):
            if cur:
                seqs.append(cur)
            cur = ''
        else:
            cur += l.strip()
    if cur:
        seqs.append(cur)
    return seqs  # one per CDS feature (usually 1)


rows = []
for p in sorted(glob.glob(os.path.join(DNA, '*.fasta'))):
    gene = os.path.basename(p).replace('_dna.fasta', '')
    hdr, wt = read_fasta(p)
    m = ACC.match(hdr or '')
    row = {'gene': gene, 'header': hdr, 'accession': m.group(1) if m else None,
           'benchmark_len': len(wt)}
    if m:
        cds = fetch_cds(m.group(1))
        if cds:
            ref = max(cds, key=len)
            row.update(refseq_cds_len=len(ref), n_cds_features=len(cds),
                       exact_match=bool(ref == wt))
            if len(ref) == len(wt):
                row['identity_frac'] = sum(a == b for a, b in zip(ref, wt)) / len(wt)
            else:
                # check containment either way (UTR/trim differences)
                row['benchmark_in_refseq'] = wt in ref
                row['refseq_in_benchmark'] = ref in wt
        else:
            row['fetch_failed'] = True
    rows.append(row)
    print(gene, row.get('accession'), row.get('exact_match'), flush=True)

ok = [r for r in rows if r.get('exact_match')]
mism = [r for r in rows if r.get('accession') and not r.get('exact_match') and not r.get('fetch_failed')]
out = {
    'n_genes': len(rows),
    'n_with_accession': sum(1 for r in rows if r['accession']),
    'n_no_accession_header': [r['gene'] for r in rows if not r['accession']],
    'n_fetched': sum(1 for r in rows if r.get('accession') and not r.get('fetch_failed')),
    'n_fetch_failed': [r['gene'] for r in rows if r.get('fetch_failed')],
    'n_exact_cds_match': len(ok),
    'n_mismatch': len(mism),
    'mismatches': [{k: v for k, v in r.items() if k != 'header'} for r in mism],
    'note': 'exact_match = benchmark WT DNA identical to the joined CDS of the live NCBI record for the header accession.version; mismatches are reported per gene, not assumed away',
}
json.dump({'summary': out, 'per_gene': rows}, open(os.path.join(ROOT, 'results/refseq_provenance.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
