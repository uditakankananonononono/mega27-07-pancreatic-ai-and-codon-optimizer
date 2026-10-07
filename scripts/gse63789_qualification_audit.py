"""Ingest-only source/field/join qualification; no transfer scoring."""
import argparse,gzip,hashlib,json,re
from pathlib import Path
from collections import Counter
from Bio import SeqIO
ROOT=Path(__file__).resolve().parents[1]
def compute(source):
 source=Path(source)
 if hashlib.sha256(source.read_bytes()).hexdigest()!='1177ee0766cdc664787b1b32044fcbd4a24a55c921ec6844cd74080371c26e93':raise ValueError('source digest')
 raw=(ROOT/'scripts/gse63789_qualification_plan.json').read_bytes();records=[]
 with gzip.open(source,'rt') as f:
  header=f.readline().rstrip().split('\t')
  if header!=['#Name','Length','Sum-mRNA','Sum-FP','FP']:raise ValueError('header')
  for line in f:
   fields=line.rstrip().split('\t')
   if len(fields)!=5:raise ValueError('row fields')
   g,L,rna,fp=fields[:4];L,rna,fp=map(int,[L,rna,fp]);positions=list(map(int,fields[4].split()))
   if L<=0 or min(rna,fp,*positions)<0:raise ValueError('count range')
   records.append({'gene':g,'reported_length':L,'rna_count':rna,'fp_count':fp,'positions':len(positions),'position_sum':sum(positions)})
 names=[r['gene'] for r in records];old=json.loads((ROOT/'results/yeast_translation_audit.json').read_text());cds={}
 for seq in SeqIO.parse(gzip.open(ROOT/'data/yeast/GCF_000146045.2_R64_cds_from_genomic.fna.gz','rt'),'fasta'):
  m=re.search(r'\[locus_tag=([^\]]+)\]',seq.description)
  if m:cds.setdefault(m[1],[]).append(len(seq.seq))
 joined=[r for r in records if r['gene'] in cds and len(cds[r['gene']])==1]
 return {'plan_sha256':hashlib.sha256(raw).hexdigest(),'source_url':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63789/suppl/GSE63789_counts_wt.txt.gz','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'header':header,'rows':len(records),'unique_genes':len(set(names)),'duplicates':len(names)-len(set(names)),'zero_rna':sum(r['rna_count']==0 for r in records),'zero_fp':sum(r['fp_count']==0 for r in records),'position_length_mismatches':sum(r['positions']!=r['reported_length'] for r in records),'position_sum_mismatches':sum(r['position_sum']!=r['fp_count'] for r in records),'unique_cds_joined':len(joined),'unmatched_genes':[r['gene'] for r in records if r['gene'] not in cds],'ambiguous_cds_genes':[r['gene'] for r in records if len(cds.get(r['gene'],[]))>1],'three_times_reported_minus_cds_nt_histogram':dict(sorted(Counter(str(3*r['reported_length']-cds[r['gene']][0]) for r in joined).items())),'length_unit_interpretation':'codon-scale geometry inferred from 5110 exact 3L=CDS_nt matches; do not compare Length directly to nucleotide CDS length; 32 residual mismatches unresolved','length_mismatch_genes':[{**r,'cds_nt':cds[r['gene']][0],'three_times_length_minus_cds_nt':3*r['reported_length']-cds[r['gene']][0]} for r in joined if 3*r['reported_length']!=cds[r['gene']][0]],'old_gse75897_gene_overlap':len(set(names)&set(old['gene_ids'])),'qualification':'conditional: paired gene totals exist; raw count ratios need explicit endpoint/normalization and R63-to-R64 length handling prereg before scoring','metadata_sources':['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63789&targ=self&form=text&view=full','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1557442&targ=self&form=text&view=full','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1557447&targ=self&form=text&view=full'],'metadata_findings':{'samples':['GSM1557442 WT mrna','GSM1557447 WT footprint'],'genome_build':'R63','counts':'total RNA/footprint plus per-position footprint; not documented RPKM columns','filter':'RNA reads>26nt,footprints28..31nt inclusive;uniquely mapped','replicate_scope':'one WT RNA and one WT footprint sample in series; no replicate uncertainty from this pair'},'limits':json.loads(raw)['limits']}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--source',required=True);x=a.parse_args();j=compute(x.source);(ROOT/'results/gse63789_qualification_audit.json').write_text(json.dumps(j,indent=2)+'\n');print({k:v for k,v in j.items() if not k.endswith('histogram') and k not in ['unmatched_genes','metadata_sources']});print('3L-CDS',j['three_times_reported_minus_cds_nt_histogram'])
