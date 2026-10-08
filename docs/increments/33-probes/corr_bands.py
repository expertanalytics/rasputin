import numpy as np, pickle, shapely, json
from shapely.geometry import shape
ml,_=pickle.load(open('line_3.pkl','rb'))
dom=shape(json.load(open('corridor_domain.geojson'))['features'][0]['geometry'])
simp=ml.simplify(1.0,preserve_topology=False)
segs=[]
for g in simp.geoms:
    c=np.asarray(g.coords); segs+= [shapely.linestrings(np.stack([c[:-1],c[1:]],1))]
segs=np.concatenate(segs); tree=shapely.STRtree(segs)
print('segments after 1 m simplification',len(segs))
edges=[0,100,250,500,1000,1500,2000,2500,3000,1e9]
rng=np.random.default_rng(1)
def load(T):
    b=open(f'corr_t{T}.ply','rb').read(); h=b.index(b'end_header\n')+len(b'end_header\n'); hdr=b[:h].decode()
    nv=int(hdr.split('element vertex ')[1].split()[0]); nf=int(hdr.split('element face ')[1].split()[0])
    v=np.frombuffer(b,dtype='<f8',count=3*nv,offset=h).reshape(nv,3)
    f=np.frombuffer(b,dtype=np.dtype([('n','u1'),('i','<u4',3)]),count=nf,offset=h+24*nv)['i']
    return v,f
res={}
for T in (20,10,5,2,1):
    v,f=load(T); n=len(f); idx=rng.choice(n,min(n,400000),replace=False); scale=n/len(idx)
    p=v[f[idx].astype(np.int64)][:,:,:2]; c=p.mean(1)
    a=0.5*np.abs((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])-(p[:,2,0]-p[:,0,0])*(p[:,1,1]-p[:,0,1]))
    _,d=tree.query_nearest(shapely.points(c),return_distance=True,all_matches=False)
    k=np.digitize(d,edges)-1
    cnt=np.bincount(k,minlength=9)*scale; ar=np.bincount(k,weights=a,minlength=9)*scale
    res[T]=(cnt,ar); print(T,'triangles',n,'per band',np.round(cnt).astype(int).tolist(),flush=True)
ar1=res[1][1]; print('band areas km2',np.round(ar1/1e6,1).tolist())
Ts=[1,2,5,10,20]; dens={T:res[T][0]/res[T][1] for T in Ts}
mids=np.array([50,175,375,750,1250,1750,2250,2750,4000.])
def est(r0,r1,near=1.0,far=20.0):
    tot=0
    for b in range(9):
        d=mids[b]; t=near if d<=r0 else far if d>=r1 else near+(far-near)*(d-r0)/(r1-r0)
        tot+=np.exp(np.interp(np.log(t),np.log(Ts),[np.log(dens[T][b]) for T in Ts]))*ar1[b]
    return int(tot)
print('estimate linear 0-3000:',est(0,3000),' flat to 100 then linear:',est(100,3000),' flat to 500:',est(500,3000))
