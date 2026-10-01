"""Retrieve the observed publisher workbook privately; check rights before redistribution."""
import hashlib,urllib.request
from pathlib import Path
URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature16509/MediaObjects/41586_2016_BFnature16509_MOESM292_ESM.xlsx'
SHA='7e32a22549f4565a08284e64f4187c5562181b909bfb382849cd9227fe4ef978'
if __name__=='__main__':
 data=urllib.request.urlopen(URL,timeout=45).read();assert hashlib.sha256(data).hexdigest()==SHA
 p=Path(__file__).resolve().parents[1]/'data/external/nature16509_data2.xlsx';p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data);print('Retrieved hash-verified source. Inspect publisher terms before redistribution:',p)
