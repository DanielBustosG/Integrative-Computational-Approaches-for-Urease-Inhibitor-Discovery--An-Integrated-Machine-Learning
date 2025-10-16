# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import os
import re
import time
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations

# =========================================================
# 1) # Configure dataset file paths and human‑readable labels. The script expects
# a SMILES column in each table. For .smi files, the second token on each line
# (if present) is used as the molecule name; otherwise an auto name is created.
# =========================================================
MASTER_PATH = r"C:/Users/Javier/Desktop/Paper_Sofi/master_urease_dataset_unique2.xlsx"
OTHER1_PATH = r"C:/Users/Javier/Desktop/Paper_Sofi/Candidatas2.smi"
OTHER1_LABEL = "Candidate Molecules"
OTHER2_PATH = r"C:/Users/Javier/Desktop/Paper_Sofi/Control.smi"
OTHER2_LABEL = "Control Molecules"
MASTER_LABEL = "Urease inhibitor set"

# List of (path, label) pairs used downstream for loading and labeling.

DATASETS = [
    (MASTER_PATH, MASTER_LABEL),
    (OTHER1_PATH, OTHER1_LABEL),
    (OTHER2_PATH, OTHER2_LABEL),
]

OUT_DIR = os.path.join(os.path.dirname(OTHER1_PATH), "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# =========================================================
# 1.1 Safe Excel saving
# ---------------------------------------------------------
# save_excel_safe(): Robust writer that retries on PermissionError (e.g., file
# open in Excel) and switches writer engine if needed. It also timestamps the
# filename when collisions happen.
# =========================================================
def save_excel_safe(df: pd.DataFrame, path: str, tries: int = 3) -> str:
    path = Path(path)
    stem, suffix = path.stem, path.suffix or ".xlsx"
    for _ in range(tries):
        try:
            df.to_excel(path, index=False, engine="openpyxl")
            return str(path)
        except PermissionError:
            ts = time.strftime("%Y%m%d_%H%M%S")
            path = path.with_name(f"{stem}_{ts}{suffix}")
            time.sleep(0.25)
        except Exception:
            try:
                df.to_excel(path, index=False, engine="xlsxwriter")
                return str(path)
            except PermissionError:
                ts = time.strftime("%Y%m%d_%H%M%S")
                path = path.with_name(f"{stem}_{ts}{suffix}")
                time.sleep(0.25)
    df.to_excel(path, index=False, engine="openpyxl")
    return str(path)

# =========================================================
# 2) Helpers
# Utility functions to locate key columns, read CSVs robustly, normalize keys,
# and pick optional activity metadata when present.
# =========================================================
EXTRA_COLS = ["Activity_Type", "Value_nM", "Units"]

def _find_smiles_col(columns):
    for c in columns:
        if str(c).strip().lower() == "smiles":
            return c
    return None

def _find_name_col(columns):
    candidates = [
        "name", "compound_id", "compound id", "molecule", "molecule_name",
        "molecule name", "title", "preferred_name", "pert_iname", "drug_name", "id"
    ]
    low = {str(c).strip().lower(): c for c in columns}
    for k in candidates:
        if k in low:
            return low[k]
    return None

def _safe_read_csv(path):
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path)

def _normalize_key(s: str) -> str:
    return re.sub(r"[\s_]+", "", str(s).strip().lower())

def _pick_extra_cols(df: pd.DataFrame) -> dict:
    std_keys = {c: _normalize_key(c) for c in EXTRA_COLS}
    inv = {_normalize_key(c): c for c in df.columns}
    out = {}
    for std_name, std_key in std_keys.items():
        candidates = [std_key]
        if std_name == "Activity_Type":
            candidates += ["activitytype"]
        elif std_name == "Value_nM":
            candidates += ["valuenm", "value(nm)", "value_nm"]
        elif std_name == "Units":
            candidates += ["unit", "units"]
        for k in candidates:
            if k in inv:
                out[std_name] = inv[k]
                break
    return out

def load_smiles_table(path):
    ext = os.path.splitext(path)[1].lower()
    base = os.path.splitext(os.path.basename(path))[0]

    def _ensure_name(df, base):
        if "Name" not in df.columns:
            df["Name"] = [f"{base}_{i+1:05d}" for i in range(len(df))]
        else:
            df["Name"] = df["Name"].astype(str).str.strip().replace({"": None})
            df["Name"] = df["Name"].fillna(pd.Series([f"{base}_{i+1:05d}" for i in range(len(df))], index=df.index))
        return df

    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
        col_smi = _find_smiles_col(df.columns)
        if col_smi is None:
            raise ValueError(f"No 'SMILES' column found in: {path}")
        col_name = _find_name_col(df.columns)
        keep = [col_smi] + ([col_name] if col_name else [])
        extra_map = _pick_extra_cols(df)
        for std in EXTRA_COLS:
            if std in extra_map:
                keep.append(extra_map[std])
        df = df[keep].rename(columns={col_smi: "SMILES", (col_name or "Name"): "Name"})
        if extra_map:
            df = df.rename(columns={v: k for k, v in extra_map.items()})
        df = _ensure_name(df, base)

    elif ext in [".csv", ".tsv"]:
        df = _safe_read_csv(path)
        col_smi = _find_smiles_col(df.columns)
        if col_smi is None:
            raise ValueError(f"No 'SMILES' column found in: {path}")
        col_name = _find_name_col(df.columns)
        keep = [col_smi] + ([col_name] if col_name else [])
        extra_map = _pick_extra_cols(df)
        for std in EXTRA_COLS:
            if std in extra_map:
                keep.append(extra_map[std])
        df = df[keep].rename(columns={col_smi: "SMILES", (col_name or "Name"): "Name"})
        if extra_map:
            df = df.rename(columns={v: k for k, v in extra_map.items()})
        df = _ensure_name(df, base)

    elif ext == ".smi":
        smiles, names = [], []
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for idx, line in enumerate(f, 1):
                parts = line.strip().split()
                if not parts:
                    continue
                smi = parts[0]
                name = parts[1] if len(parts) > 1 else f"{base}_{idx:05d}"
                smiles.append(smi)
                names.append(name)
        df = pd.DataFrame({"SMILES": smiles, "Name": names})
    else:
        raise ValueError(f"Unsupported file extension: {path}")

    df = df[df["SMILES"].notna()].copy()
    df["SMILES"] = df["SMILES"].astype(str).str.strip()
    df = df[df["SMILES"] != ""]
    for c in EXTRA_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    return df[["SMILES", "Name"] + EXTRA_COLS]

# =========================================================
# 3) Load datasets
# ---------------------------------------------------------
# Read, de‑duplicate by SMILES, and merge all datasets. Also prepare per‑label
# SMILES sets for quick counts and overlap reasoning. Save the combined table.
# =========================================================
frames, sets_by_label = [], {}
for path, label in DATASETS:
    df = load_smiles_table(path)
    df = df.drop_duplicates(subset=["SMILES"]).reset_index(drop=True)
    df["Label"] = label
    frames.append(df)
    sets_by_label[label] = set(df["SMILES"])

df_all = pd.concat(frames, ignore_index=True)
combined_xlsx = os.path.join(OUT_DIR, "chemspace_combined_SMILES_labels.xlsx")
save_excel_safe(df_all, combined_xlsx)

print("Summary per dataset:")
for label in sets_by_label:
    print(f"  - {label:<22}: {len(sets_by_label[label])} molecules")

# =========================================================
# 4) ChemPlot visualization 
# ---------------------------------------------------------
# Build a ChemPlot Plotter from SMILES + class labels. Choose the reducer
# (PCA/UMAP/t‑SNE). Then export a static PNG and an interactive HTML plot.
# Post‑process the static figure to tweak marker aesthetics.
# =========================================================
from chemplot import Plotter
import matplotlib.pyplot as plt
import os


REDUCER = "pca" # one of: 'pca', 'umap', 'tsne'
REMOVE_OUTLIERS = False # optional outlier removal inside ChemPlot


# Create the Plotter with categorical targets given by dataset Label.
cp = Plotter.from_smiles(
df_all["SMILES"].tolist(),
target=df_all["Label"].tolist(),
target_type="C", # 'C' = categorical targets
sim_type="structural" # structural similarity (vs property‑based)
)


# Run the chosen dimensionality reducer.
if REDUCER.lower() == "umap":
    cp.umap()
elif REDUCER.lower() == "tsne":
    cp.tsne()
else:
    cp.pca()

labels = [label for _, label in DATASETS]
suffix = REDUCER.upper()

png_name = f"chemspace_{'_vs_'.join(labels)}_{suffix}.png"
static_png = os.path.join(OUT_DIR, png_name)
cp.visualize_plot(
    kind="scatter",
    remove_outliers=REMOVE_OUTLIERS,
    is_colored=True,
    filename=static_png,
    title="Chemical Space"
)
print(f" Static PNG saved → {static_png}")

html_name = f"chemspace_{'_vs_'.join(labels)}_{suffix}.html"
interactive_html = os.path.join(OUT_DIR, html_name)
cp.interactive_plot(
    kind="scatter",
    remove_outliers=REMOVE_OUTLIERS,
    is_colored=True,
    filename=interactive_html,
    show_plot=False,
    title="Chemical Space"
)
print(f"Interactive HTML saved → {interactive_html}")

# --- Post‑render: tweak marker edges/size/alpha on the static figure ---
fig = plt.gcf()
ax = plt.gca()
for coll in ax.collections:
    coll.set_edgecolor("black")
    coll.set_linewidth(0.5)
    coll.set_alpha(0.9)
    coll.set_sizes([220])  # 🔹 tamaño de los círculos aumentado
fig.savefig(static_png, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print(f" PNG recolored and finalized → {static_png}")

# =========================================================
# 5) Nearest neighbors in ChemPlot coordinates (TOP‑K) ChemPlot (TOP-K)
# For each candidate molecule (OTHER1_LABEL), compute Euclidean distance in the
# 2D embedding to all molecules from other datasets, select top‑K neighbors, and
# export a long‑format table including activity fields when the neighbor belongs
# to the MASTER set.
# =========================================================
TOP_K = 5

def _pick_xy_columns(df_plot: pd.DataFrame):
    cols = {c.lower(): c for c in df_plot.columns}
    for x_key, y_key in [
        ("x", "y"),
        ("x1", "y1"),
        ("pc1", "pc2"),
        ("dim 1", "dim 2"),
        ("dim1", "dim2"),
    ]:
        if x_key in cols and y_key in cols:
            return cols[x_key], cols[y_key]
    num_cols = [c for c in df_plot.columns if np.issubdtype(df_plot[c].dtype, np.number)]
    if len(num_cols) >= 2:
        return num_cols[0], num_cols[1]
    raise RuntimeError("Could not find X/Y columns in ChemPlot df_plot.")

def get_chemplot_coords_or_fallback(cp, df_all):
    df_plot = getattr(cp, "df_plot", None)
    if isinstance(df_plot, pd.DataFrame) and len(df_plot) == len(df_all):
        xcol, ycol = _pick_xy_columns(df_plot)
        coords_df = df_all[["SMILES", "Name", "Label"] + EXTRA_COLS].copy()
        coords_df["X"] = df_plot[xcol].to_numpy()
        coords_df["Y"] = df_plot[ycol].to_numpy()
        return coords_df

    print("Could not read 2D coords from ChemPlot; using PCA over ECFP4 as fallback.")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from sklearn.decomposition import PCA

    def mol_from_smiles(s):
        try:
            return Chem.MolFromSmiles(s)
        except Exception:
            return None

    N_BITS = 2048
    df_tmp = df_all.copy()
    df_tmp["Mol"] = df_tmp["SMILES"].apply(mol_from_smiles)
    df_tmp = df_tmp[df_tmp["Mol"].notna()].copy()
    df_tmp["FP"] = df_tmp["Mol"].apply(lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=N_BITS))
    X = np.vstack([np.array(fp.ToList(), dtype=np.int8) for fp in df_tmp["FP"]])
    coords = PCA(n_components=2, random_state=42).fit_transform(X)

    coords_df = df_tmp[["SMILES", "Name", "Label"] + EXTRA_COLS].copy()
    coords_df["X"] = coords[:, 0]
    coords_df["Y"] = coords[:, 1]
    return coords_df

coords_df = get_chemplot_coords_or_fallback(cp, df_all).reset_index(drop=True)
queries = coords_df[coords_df["Label"] == OTHER1_LABEL].copy()
pool = coords_df[coords_df["Label"] != OTHER1_LABEL].copy()

rows_long = []
for _, q in queries.iterrows():
    dx = pool["X"].to_numpy() - q["X"]
    dy = pool["Y"].to_numpy() - q["Y"]
    dist = np.sqrt(dx * dx + dy * dy)
    k = min(TOP_K, np.isfinite(dist).sum())
    idx = np.argpartition(dist, k)[:k]
    idx = idx[np.argsort(dist[idx])]
    for rank, i in enumerate(idx, start=1):
        nn = pool.iloc[i]
        row = {
            "Query_Name": q["Name"],
            "Query_SMILES": q["SMILES"],
            "NN_Rank": rank,
            "Neighbor_Name": nn["Name"],
            "Neighbor_SMILES": nn["SMILES"],
            "Neighbor_Origin": nn["Label"],
            "Dist2D": float(dist[i]),
        }
        if nn["Label"] == MASTER_LABEL:
            row["NN_Activity_Type"] = nn.get("Activity_Type", pd.NA)
            row["NN_Value_nM"] = nn.get("Value_nM", pd.NA)
            row["NN_Units"] = nn.get("Units", pd.NA)
        else:
            row["NN_Activity_Type"] = pd.NA
            row["NN_Value_nM"] = pd.NA
            row["NN_Units"] = pd.NA
        rows_long.append(row)

nn_long_df = pd.DataFrame(rows_long)
safe_label = OTHER1_LABEL.replace(" ", "_")
nn_long_xlsx = os.path.join(OUT_DIR, f"nearest_neighbors_top{TOP_K}_long_{safe_label}_{suffix}_2D.xlsx")
nn_long_xlsx = save_excel_safe(nn_long_df, nn_long_xlsx)

print(f"\nSaved NN (long, ChemPlot 2D): {nn_long_xlsx}")
print(f"   Candidates: {len(queries)} | Rows NN: {len(nn_long_df)} "
      f"(expected ≈ {len(queries)*TOP_K} if there is enough pool)")
