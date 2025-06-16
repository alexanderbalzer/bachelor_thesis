import requests
import pandas as pd

URL = "https://gtexportal.org/rest/v2/expression/geneExpression"

GENES = {
    "NACA": "ENSG00000131183",
    "NACAD": "ENSG00000177511",
    "BTF3": "ENSG00000112110"
}

def fetch_expression(gene_symbol, gencode_id):
    params = {
        "gencodeId": gencode_id,
        "datasetId": "gtex_v8",
        "format": "json"
    }
    try:
        r = requests.get(URL, params=params, timeout=10)
        r.raise_for_status()
    except Exception as e:
        print(f"[ERROR] Request failed for {gene_symbol}: {e}")
        return pd.DataFrame()
    try:
        js = r.json()
    except Exception as e:
        print(f"[ERROR] JSON decode failed for {gene_symbol}: {e}")
        print("RAW:", r.text[:200])
        return pd.DataFrame()
    data = js.get("geneExpression", [])
    if not data:
        print(f"[INFO] No expression data found for {gene_symbol}")
        return pd.DataFrame()
    df = pd.DataFrame(data)
    return df[["tissueSiteDetailId", "median"]].set_index("tissueSiteDetailId")

frames = {symbol: fetch_expression(symbol, eid) for symbol, eid in GENES.items()}
df = pd.concat(frames, axis=1)
df.columns = list(GENES.keys())
print(df)
