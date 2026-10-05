
from pathlib import Path

import pandas as pd
from scipy.interpolate import CubicSpline


def load_metric_splines(csv_path: Path):
    df = pd.read_csv(csv_path, index_col=0)
    u = df.index.values
    f = df["f"].values
    h = df["h"].values
    return CubicSpline(u, f), CubicSpline(u, h)
