"""
Descarga los CDF necesarios desde el Solar Orbiter Archive (SOAR).

  - solo_L2_swa-pas-vdf_20220308   (VDF de protones,  ~820 MB)
  - solo_L2_mag-rtn-normal_20220308 (campo magnético, ~11 MB)

SOAR REST:  retrieval_type=LAST_PRODUCT devuelve la última versión del producto.
    https://soar.esac.esa.int/soar-sl-tap/data?retrieval_type=LAST_PRODUCT
        &data_item_id=<id>&product_type=SCIENCE

Uso:  .venv/bin/python src/download_data.py
"""
import os, sys
import requests
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DATA_RAW

SOAR = "https://soar.esac.esa.int/soar-sl-tap/data"
PRODUCTS = {
    "solo_L2_swa-pas-vdf_20220308":     "solo_L2_swa-pas-vdf_20220308_V02.cdf",
    "solo_L2_mag-rtn-normal_20220308":  "solo_L2_mag-rtn-normal_20220308_V02.cdf",
}


def download(data_item_id, filename):
    dest = os.path.join(DATA_RAW, filename)
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        print(f"[skip] {filename} ya existe ({os.path.getsize(dest)/1e6:.1f} MB)")
        return dest
    params = dict(retrieval_type="LAST_PRODUCT", data_item_id=data_item_id,
                  product_type="SCIENCE")
    print(f"[get ] {data_item_id} ...")
    with requests.get(SOAR, params=params, stream=True, timeout=1800) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        done = 0
        with open(dest, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk); done += len(chunk)
                if total:
                    print(f"\r       {done/1e6:7.1f}/{total/1e6:.1f} MB "
                          f"({100*done/total:5.1f}%)", end="", flush=True)
        print()
    print(f"[ok  ] {dest}  ({os.path.getsize(dest)/1e6:.1f} MB)")
    return dest


if __name__ == "__main__":
    os.makedirs(DATA_RAW, exist_ok=True)
    for pid, fn in PRODUCTS.items():
        download(pid, fn)
