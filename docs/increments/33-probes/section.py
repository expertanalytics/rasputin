import pickle, json, shapely
from shapely.geometry import box, mapping
ml,corr=pickle.load(open('line_1.pkl','rb'))
sec=ml.intersection(box(127000,6720000,150000,6745000))
sec=shapely.line_merge(sec)
print('section km (all links)',round(sec.length/1e3,1))
dom=shapely.buffer(sec.simplify(5.0),5000,quad_segs=4).intersection(box(127000,6700000,150000,6770000)).simplify(20.0)
print('domain km2',round(dom.area/1e6,1),'vertices',len(shapely.get_coordinates(dom)), dom.geom_type)
json.dump({"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":mapping(dom)}]},open('section_domain.geojson','w'))
pickle.dump((sec,dom),open('section.pkl','wb'))
