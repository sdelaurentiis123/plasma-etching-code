import sys, json; sys.path.insert(0,'src'); sys.path.insert(0,'scripts')
import numpy as np
from cascade_funnelling_scan import (make_straight_trench_geometry_3d,
                                     gather_transport_pilot, floor_face_mask,
                                     delivery_split, OPEN_WIDTH_UM)
recs=[]
for ar in (2.0, 12.0):
    geometry, floor_z = make_straight_trench_geometry_3d(
        etched_depth=float(ar)*OPEN_WIDTH_UM, dx=0.01)
    res, _, _, el = gather_transport_pilot(geometry, transport_device="cpu")
    m = floor_face_mask(res, floor_z_um=floor_z)
    d = delivery_split(res, m)
    diag = next((v for v in res.transport.hit_probability.values()
                 if isinstance(v, dict) and "thermalized_rate_per_face" in v), None)
    active = np.asarray(res.active_face_index, dtype=int)
    areas = np.asarray(res.active_face_area, dtype=float)
    pf = np.asarray(diag["thermalized_rate_per_face"], float)
    therm = float(pf[active][m].sum())/max(float(areas[m].sum()),1e-30)
    plasma = sum(v for k,v in d.items() if k.startswith("neutral_"))
    rec=dict(ar=ar, faces=int(m.sum()), e8_source=therm, plasma_neutral=plasma,
             e8_share=therm/max(therm+plasma,1e-30), gather_s=el,
             direct_ion=d["direct_ion"], hot_neutral=d["hot_neutral"])
    recs.append(rec); print(json.dumps(rec), flush=True)
json.dump(recs, open('results/curated/e8_thermalized_return/floor_composition.json','w'), indent=1)
print("DONE")
