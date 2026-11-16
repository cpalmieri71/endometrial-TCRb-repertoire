"""Shared helpers: configuration, patient codes, logistic regression.

Every script reads paths from a JSON configuration file (default: config/config.json,
override with --config or the ENDO_CONFIG environment variable). Relative paths in the
configuration are resolved against "project_root". See config/config.example.json.

Patient identifiers never appear in the code. Published patient codes (P1-P13) are read
from a private key file (private/patient_key.csv) that must not be committed.
"""
import argparse, csv, json, os
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Config:
    def __init__(self, path):
        with open(path) as f:
            self.raw = json.load(f)
        self.root = self.raw["project_root"]

    def path(self, key):
        v = self.raw[key]
        return v if os.path.isabs(v) else os.path.join(self.root, v)

    @property
    def repseq(self): return self.path("repseq_dir")
    @property
    def results(self): return self.path("results_dir")
    @property
    def private(self): return self.path("private_dir")
    @property
    def out(self): return self.path("output_dir")

    def result(self, *parts): return os.path.join(self.results, *parts)

    def out_dir(self, sub):
        d = os.path.join(self.out, sub)
        os.makedirs(d, exist_ok=True)
        return d


def load_config(argv=None, description=None):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--config", default=os.environ.get("ENDO_CONFIG", os.path.join(REPO, "config", "config.json")))
    args, _ = ap.parse_known_args(argv)
    return Config(args.config)


def patient_codes(cfg):
    """{internal_id (int): 'Pn'} from the private key file."""
    with open(os.path.join(cfg.private, "patient_key.csv"), newline="") as f:
        return {int(r["internal_id"]): r["patient_code"] for r in csv.DictReader(f)}


def curated_annotations(cfg):
    """Private per-patient curated annotations (nodal status, MSI status, POLE/TP53 variants,
    molecular class), keyed by internal identifier."""
    with open(os.path.join(cfg.private, "curated_annotations.csv"), newline="", encoding="utf-8-sig") as f:
        return {int(r["internal_id"]): r for r in csv.DictReader(f)}


def repseq_files(cfg):
    """{(internal_id, 'BIOPSY'|'PBMC'): path} for the 26 MiXCR clonotype tables."""
    import glob, re
    out = {}
    for f in sorted(glob.glob(os.path.join(cfg.repseq, "*.tsv"))):
        m = re.match(r"^(\d+)_(PBMC|biopsia|biopsy)\.tsv$", os.path.basename(f))
        if m:
            out[(int(m.group(1)), "PBMC" if m.group(2) == "PBMC" else "BIOPSY")] = f
    return out


def irls(X, y, tol=1e-10, maxit=100):
    """Logistic regression by iteratively reweighted least squares.
    Returns coefficients, fitted probabilities and the Fisher information matrix."""
    beta = np.zeros(X.shape[1])
    for _ in range(maxit):
        p = 1 / (1 + np.exp(-(X @ beta)))
        W = np.clip(p * (1 - p), 1e-9, None)
        H = (X.T * W) @ X
        step = np.linalg.solve(H + 1e-10 * np.eye(len(beta)), X.T @ (y - p))
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    p = 1 / (1 + np.exp(-(X @ beta)))
    return beta, p, H
