#!/usr/bin/env python3
"""Build TAPSy retrieval data from public site pages, curated notes and atlas JSON."""
import json,re,hashlib
from pathlib import Path
from html.parser import HTMLParser
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
DOCS=ROOT/'docs'
BASE='https://henzelerb.github.io/sRNA-TAPS/'
class Text(HTMLParser):
    def __init__(self):
        super().__init__(); self.skip=[]; self.parts=[]
    def handle_starttag(self,tag,attrs):
        if self.skip:
            if tag==self.skip[-1]:self.skip.append(tag)
            return
        if tag in ('script','style','svg','nav','footer','head'):self.skip.append(tag)
        elif tag in ('p','div','li','h1','h2','h3','h4','pre','tr','br'):self.parts.append('\n')
    def handle_endtag(self,tag):
        if self.skip and tag==self.skip[-1]:self.skip.pop()
        elif not self.skip and tag in ('p','div','li','h1','h2','h3','h4','pre','tr'):self.parts.append('\n')
    def handle_data(self,data):
        if not self.skip:self.parts.append(data)
def build():
    chunks=[];manifest={}
    for path in sorted(DOCS.glob('*.html')):
        html=path.read_text()
        if '</body>' not in html:continue
        manifest[path.name]=hashlib.sha256(html.encode()).hexdigest()
        parser=Text();parser.feed(html)
        lines=[' '.join(x.split()) for x in ''.join(parser.parts).splitlines() if x.strip()]
        text='\n'.join(lines)
        title=re.search(r'<title>(.*?)</title>',html,re.S)[1]
        for start in range(0,len(text),1400):
            chunks.append(dict(title=title,url=BASE+path.name,page=path.name,text=text[start:start+1600],kind='page'))
    for row in json.loads((DOCS/'knowledge/literature.json').read_text()):
        chunks.append(dict(row,kind='literature'))
    html=(DOCS/'m5c-atlas.html').read_text()
    rows=json.loads(re.search(r'<script id="site-data"[^>]*>(.*?)</script>',html,re.S)[1])
    fields=('id','biotype','subtype','gene_name','transcript','chrom','position','strand','position_label','mod_rate','coverage','status','literature')
    sites=[{k:r.get(k) for k in fields} for r in rows]
    summary={'total':len(rows),'biotypes':dict(Counter(r['biotype'] for r in rows)),
             'status_by_biotype':{b:dict(Counter(r['status'] for r in rows if r['biotype']==b)) for b in sorted({r['biotype'] for r in rows})}}
    data=dict(schema=1,sources=manifest,chunks=chunks,sites=sites,atlas_summary=summary)
    (DOCS/'knowledge/index.json').write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n')
    print(f'{len(chunks)} chunks, {len(manifest)} pages, {len(sites)} atlas sites')
if __name__=='__main__':build()
