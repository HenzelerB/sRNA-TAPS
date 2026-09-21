// Run from repository root: gjs tests/test_tapsy_knowledge.js
const GLib=imports.gi.GLib;
const BA=imports.byteArray;
function read(path){return BA.toString(GLib.file_get_contents(path)[1]);}
function check(value,message){if(!value)throw new Error(message);}
(0,eval)(read('docs/tapsy-knowledge.js'));
const data=JSON.parse(read('docs/knowledge/index.json'));
let r=TAPSyKnowledge.select(data,'What supports hsa-miR-31-3p position 12?','m5c-atlas.html',null);
check(r.sites.length>0,'miRNA match missing');
check(r.sites[0].transcript==='hsa-miR-31-3p' && r.sites[0].position_label==='nt 12 (mature)','Exact mature position not prioritized');
check(r.sites[0].literature.some(l=>l.pmid==='33980133'),'Expected citation missing');
check(!r.sites.some(s=>s.transcript==='hsa-miR-31-5p'),'Wrong miRNA arm retrieved');
let selected=r.sites[0];
r=TAPSyKnowledge.select(data,'Explain this site','m5c-atlas.html',selected.id);
check(r.sites[0].id===selected.id,'Open detail context lost');
r=TAPSyKnowledge.select(data,'What does NSUN6 modify?','index.html',null);
check(r.chunks.some(c=>c.url.includes('26160102')),'NSUN6 source missing');
r=TAPSyKnowledge.select(data,'How do I install the pipeline?','docs.html',null);
check(r.chunks.some(c=>c.page==='docs.html' && /install/i.test(c.text)),'Documentation retrieval failed');
r=TAPSyKnowledge.select(data,'What are miRNA sites at position 12?','index.html',null);
check(r.sites.length===0,'Generic query falsely matched individual identifiers');
r=TAPSyKnowledge.select(data,'MIR31','index.html',null);
check(r.sites.every(s=>s.gene_name==='MIR31'),'Gene retrieval too broad');
check(r.returned_site_count<=8,'Site context unbounded');
check(data.sites.length===data.atlas_summary.total,'Summary count mismatch');
for (const f of Object.keys(data.sources)) {
 let html=read('docs/'+f);
 check(html.indexOf('tapsy-knowledge.js')<html.indexOf('src="tapsy.js"'),'Script order '+f);
}
print('PASS: exact arm/position, selected site, source citations, page retrieval, negative matches, context limits and page wiring.');
