import numpy as np, pickle, shapely
sec,dom=pickle.load(open('section.pkl','rb'))
edges=[0,100,250,500,1000,1500,2000,2500,3000,1e9]
def load(T):
    with open(f'sec_t{T}.ply') as f:
        nv=nf=0
        while True:
            l=f.readline()
            if l.startswith('element vertex'): nv=int(l.split()[2])
            if l.startswith('element face'): nf=int(l.split()[2])
            if l.startswith('end_header'): break
        v=np.loadtxt(f,max_rows=nv); fc=np.loadtxt(f,max_rows=nf,dtype=np.int64)[:,1:4]
    return v,fc
res={}
for T in (20,10,5,2,1):
    v,fc=load(T); p=v[fc][:,:,:2]
    c=p.mean(axis=1); a=0.5*np.abs((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])-(p[:,2,0]-p[:,0,0])*(p[:,1,1]-p[:,0,1]))
    d=shapely.distance(shapely.points(c),sec)
    k=np.digitize(d,edges)-1
    cnt=np.bincount(k,minlength=len(edges)-1); ar=np.bincount(k,weights=a,minlength=len(edges)-1)
    res[T]=(cnt,ar)
    print(T,'triangles',len(fc),'per band',cnt.tolist())
cnt20,ar=res[20]; _,ar1=res[1]
print('band areas km2 (t=1 mesh)',np.round(ar1/1e6,2).tolist())
Ts=np.array([1,2,5,10,20.]); dens={T:res[T][0]/res[T][1] for T in (1,2,5,10,20)}
def tol(d,near=1,far=20,r0=0,r1=3000): return np.where(d<=r0,near,np.where(d>=r1,far,near+(far-near)*(d-r0)/(r1-r0)))
mids=np.array([50,175,375,750,1250,1750,2250,2750,4000.])
for name,(r0,r1) in {'linear 0-3000':(0,3000),'flat 0-100 then linear to 3000':(100,3000)}.items():
    tot=0
    for b in range(len(mids)):
        t=float(tol(mids[b],r0=r0,r1=r1))
        ld=np.interp(np.log(t),np.log(Ts),[np.log(dens[T][b]) for T in (1,2,5,10,20)])
        tot+=np.exp(ld)*ar1[b]
    print(name,'estimated triangles',int(tot))
# step: 1 m within 3000 else 20
st=sum(dens[1][b]*ar1[b] for b in range(8))+dens[20][8]*ar1[8]; print('step 1 m to 3 km, then 20', int(st))
print('density per km2 at 1 m, band 0-100:',round(dens[1][0]*1e6),' band >3000:',round(dens[1][8]*1e6),'; at 20 m band >3000:',round(dens[20][8]*1e6))
