"""Inspect landed assay provenance/sequence overlap without measuring transfer performance."""
import collections,hashlib,json,sys
from pathlib import Path
import openpyxl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from codon_optimizer.expr_data import build_expression_dataset
def compute():
 file=ROOT/'data/external/nature16509_data2.xlsx';w=openpyxl.load_workbook(file,read_only=True,data_only=True);sheet=w['E6348.csv'];rr=list(sheet.values);headers=list(rr[0]);assert headers==['name','E','Gz48v','Gave96','Gave48','Gave144','seq'];rows=[dict(zip(headers,r)) for r in rr[1:]];pool=[(lt,seq,y) for lt,gene,seq,y in build_expression_dataset(policy="raw") if len(seq)>=90];train=set(seq.upper() for _,seq,_ in pool);seqs=[r['seq'].upper() for r in rows];counts=collections.Counter(seqs)
 return {'source_article':'https://www.nature.com/articles/nature16509','source_file':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature16509/MediaObjects/41586_2016_BFnature16509_MOESM292_ESM.xlsx','file_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'file_bytes':file.stat().st_size,'sheet':'E6348.csv','headers':headers,'n_rows':len(rows),'expression_categories':dict(collections.Counter(str(r['E']) for r in rows)),'unique_names':len(set(r['name'] for r in rows)),'unique_exact_sequences':len(counts),'duplicate_sequence_excess_rows':sum(n-1 for n in counts.values()),'non_acgt_rows':sum(bool(set(s)-set('ACGT')) for s in seqs),'not_triplet_length_rows':sum(len(s)%3!=0 for s in seqs),'min_max_nt':[min(map(len,seqs)),max(map(len,seqs))],'exact_full_construct_overlap_rows_with_current_paxdb_cds':sum(s in train for s in seqs),'current_paxdb_cds_rows':len(pool),'limits':'E is assay expression category, not protein ppm or log10 continuous abundance. Exact full-construct nonoverlap does not exclude homologous proteins, tags, shared genes or model-selection overlap. No scoring, mRNA foldgain, causal optimization or independent validation claim.'}
if __name__=='__main__':
 j=compute();(ROOT/'results/external_assay_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j,indent=2))
