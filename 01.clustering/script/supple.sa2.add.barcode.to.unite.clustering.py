"""
Fix unite clustering from L1:
- loss of barcode information under summary picklefile
"""
import os
import pickle
import snapatac2 as sa2


def update_pkl(cid, cll:str = "L1",
               from_dir: str = "L2_dlt2_encoder",
               prefix = "sa2_clustering",
               snap_prefix: str = "nfeat-top_nc50") -> None:
    if not os.path.exists(from_dir):
        raise FileExistsError(f"{from_dir} dos not exist.")
    pkl_fnm = os.path.join(from_dir, f"{prefix}_{cll}_{cid}.pkl")
    if not os.path.exists(pkl_fnm):
        raise FileExistsError(f"{pkl_fnm} does not exist.")
    out_dir = os.path.join(from_dir, "supple.add.barcode")
    os.makedirs(out_dir, exist_ok=True)
    pkl2_fnm = os.path.join(out_dir, f"{prefix}_{cll}_{cid}.pkl")
    
    snap_fnm = os.path.join(from_dir, f"{snap_prefix}_{cid}_unite.h5ad")
    if not os.path.exists(snap_fnm):
        raise FileExistsError(f"{snap_fnm} does not exist.")
    with open(pkl_fnm, 'rb') as f:
        sum = pickle.load(f)
    snap = sa2.read(snap_fnm, 'r')
    if "barcode" in sum.keys():
        print(f"barcode is already in {pkl_fnm}.")
    else:
        sum['barcode'] = snap.obs_names
    snap.close()
    with open(pkl2_fnm, 'wb') as f:
        pickle.dump(sum, f)


cids_unites = list(range(8,37))
list(map(update_pkl, cids_unites))
    
