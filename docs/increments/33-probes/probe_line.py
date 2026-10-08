"""Design probe for increment 33: the Bergen Line from Bane NOR's Banenettverk (NLOD)."""
import xml.etree.ElementTree as ET, collections, pickle, time, sys, glob
import numpy as np, shapely
from shapely.geometry import LineString, MultiLineString, box
ns={'app':'http://skjema.geonorge.no/SOSI/produktspesifikasjon/Banenettverk/1.0','gml':'http://www.opengis.net/gml/3.2'}
t0=time.perf_counter()
lines=collections.defaultdict(list)
for ev,el in ET.iterparse('Samferdsel_0000_Norge_25833_Banenettverk_GML.gml'):
    if el.tag.endswith('}Banelenke'):
        name=el.findtext('.//app:banenavn',namespaces=ns); st=el.findtext('.//app:banestatus',namespaces=ns); med=el.findtext('.//app:medium',namespaces=ns)
        c=np.array(el.find('.//gml:posList',ns).text.split(),float).reshape(-1,2)
        lines[name].append((st,med,LineString(c))); el.clear()
print('parse s',round(time.perf_counter()-t0,1),flush=True)
bb=lines['Bergensbanen']
print('Bergensbanen (status, medium) link counts',collections.Counter((s,m) for s,m,_ in bb),flush=True)
for names in (['Bergensbanen'],['Bergensbanen','Randsfjordbanen','Drammenbanen']):
    sel=[(m,l) for n in names for s,m,l in lines[n] if s=='I']
    km=sum(l.length for _,l in sel)/1e3; tun=sum(l.length for m,l in sel if m=='U')/1e3
    v=sum(len(l.coords) for _,l in sel)
    print(names,'in operation: links',len(sel),'km',round(km,1),'tunnel km',round(tun,1),'vertices',v,'mean seg m',round(km*1e3/(v-len(sel)),1),flush=True)
    ml=MultiLineString([l for _,l in sel])
    for s in (1.0,):
        t=time.perf_counter(); simp=ml.simplify(s,preserve_topology=False); print('  simplify',s,'m: vertices',len(shapely.get_coordinates(simp)),'s',round(time.perf_counter()-t,2),flush=True)
    t=time.perf_counter(); corr=shapely.buffer(ml.simplify(5.0),5000,quad_segs=4); print('  5 km corridor km2',round(corr.area/1e6),'s',round(time.perf_counter()-t,1),flush=True)
    band3=shapely.buffer(ml.simplify(5.0),3000,quad_segs=4); print('  3 km band km2',round(band3.area/1e6),flush=True)
    for r in (100,500,1000,2000):
        print('    band',r,'m km2',round(shapely.buffer(ml.simplify(5.0),r,quad_segs=4).area/1e6),flush=True)
    pickle.dump((ml,corr),open(f"line_{len(names)}.pkl",'wb'))
    tiles=[]
    for tfw in glob.glob('/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/*.tfw'):
        a=[float(x) for x in open(tfw)]; tiles.append(box(a[4]-5,a[5]-50505,a[4]+50505,a[5]+5))
    hit=[t for t in tiles if t.intersects(corr)]
    cov=shapely.union_all(hit)
    print('  DTM10 tiles touching corridor',len(hit),'uncovered corridor km2',round(corr.difference(cov).area/1e6,3),flush=True)
