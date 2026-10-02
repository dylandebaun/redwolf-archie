#!/usr/bin/env python3
import argparse
import glob
import os
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def natural_key(value: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", value)]


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def infer_feature_names(ncols: int):
    # Expected layout for window stats files here:
    # [feature columns] sample chrom start end
    # feature columns are based on haplotypes h:
    # 1..h IFS, h+1..2h+1 distance distribution,
    # then 7 fixed features: mean,var,skew,kurt,min_dist,S*,private_snps
    # So total feature columns = 2h+8
    if ncols < 8:
        return [f"col{i+1}" for i in range(ncols)]

    feature_cols = ncols - 4
    names = [f"col{i+1}" for i in range(ncols)]

    if feature_cols >= 10 and (feature_cols - 8) % 2 == 0:
        h = (feature_cols - 8) // 2
        idx = 0
        for i in range(1, h + 1):
            names[idx] = f"ifs_{i}"
            idx += 1
        for i in range(1, h + 2):
            names[idx] = f"dist_bin_{i}"
            idx += 1
        fixed = [
            "dist_mean",
            "dist_var",
            "dist_skew",
            "dist_kurt",
            "min_dist_ref",
            "s_star",
            "private_snps",
        ]
        for label in fixed:
            if idx < feature_cols:
                names[idx] = label
                idx += 1

    names[-4] = "sample"
    names[-3] = "chrom"
    names[-2] = "start"
    names[-1] = "end"
    return names


def load_files(paths):
    frames = []
    skipped = []

    for path in paths:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            skipped.append((path, "missing_or_empty"))
            continue
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                header=None,
                compression="gzip",
                on_bad_lines="skip",
            )
        except Exception as exc:  # noqa: BLE001
            skipped.append((path, f"read_error: {exc}"))
            continue

        if df.empty:
            skipped.append((path, "no_rows"))
            continue

        names = infer_feature_names(df.shape[1])
        df.columns = names
        df["__source_file"] = os.path.basename(path)
        frames.append(df)

    if not frames:
        raise SystemExit("No readable stats files found.")

    merged = pd.concat(frames, ignore_index=True)
    return merged, skipped


def add_genome_coordinates(df: pd.DataFrame):
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["mid"] = (df["start"] + df["end"]) / 2.0

    chroms = sorted(df["chrom"].dropna().unique().tolist(), key=natural_key)
    offsets = {}
    running = 0.0
    for chrom in chroms:
        offsets[chrom] = running
        chrom_max = df.loc[df["chrom"] == chrom, "end"].max()
        if pd.notna(chrom_max):
            running += float(chrom_max)

    df["genome_x"] = df["mid"] + df["chrom"].map(offsets).astype(float)
    return df


def get_feature_columns(df: pd.DataFrame):
    exclude = {"start", "end", "mid", "genome_x", "sample", "chrom"}
    feats = []
    for col in df.columns:
        if col in exclude or col.startswith("__"):
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            s = pd.to_numeric(df[col], errors="coerce")
            if s.notna().mean() < 0.5:
                continue
            if s.nunique(dropna=True) <= 1:
                continue
            feats.append(col)
    return feats


def downsample(df: pd.DataFrame, max_points: int):
    if len(df) <= max_points:
        return df
    return df.sample(n=max_points, random_state=42).sort_values("genome_x")


def add_binned_trend(ax, x, y, bins=300):
    if len(x) < bins:
        return
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < bins:
        return

    edges = np.linspace(x.min(), x.max(), bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    idx = np.digitize(x, edges) - 1
    med = np.full(bins, np.nan)
    for i in range(bins):
        vals = y[idx == i]
        if len(vals) > 0:
            med[i] = np.median(vals)
    ax.plot(centers, med, color="black", lw=1.5, alpha=0.9)


def plot_feature(df: pd.DataFrame, feature: str, out_dir: str, max_points: int, dpi: int):
    d = df[["genome_x", "chrom", feature]].copy()
    d[feature] = pd.to_numeric(d[feature], errors="coerce")
    d = d.dropna()
    if d.empty:
        return

    d = downsample(d, max_points)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.scatter(d["genome_x"], d[feature], s=4, alpha=0.4, linewidths=0)
    add_binned_trend(ax, d["genome_x"], d[feature])

    chroms = sorted(d["chrom"].unique().tolist(), key=natural_key)
    xticks = []
    labels = []
    for chrom in chroms:
        sub = d[d["chrom"] == chrom]
        if len(sub) == 0:
            continue
        xticks.append((sub["genome_x"].min() + sub["genome_x"].max()) / 2.0)
        labels.append(chrom)

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("Genome position (concatenated chromosomes)")
    ax.set_ylabel(feature)
    ax.set_title(f"ArchIE window signal: {feature}")
    ax.grid(alpha=0.2)
    fig.tight_layout()

    out_path = os.path.join(out_dir, f"{sanitize_name(feature)}.png")
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_overview(df: pd.DataFrame, features, out_dir: str, max_points: int, dpi: int):
    if not features:
        return
    chosen = [
        f
        for f in [
            "s_star",
            "min_dist_ref",
            "private_snps",
            "dist_mean",
            "dist_var",
            "dist_skew",
            "dist_kurt",
            "ifs_1",
        ]
        if f in features
    ]
    if not chosen:
        chosen = features[:8]

    n = len(chosen)
    fig, axes = plt.subplots(n, 1, figsize=(14, 2.2 * n), sharex=True)
    if n == 1:
        axes = [axes]

    d0 = downsample(df, max_points)
    for ax, feat in zip(axes, chosen):
        d = d0[["genome_x", "chrom", feat]].copy()
        d[feat] = pd.to_numeric(d[feat], errors="coerce")
        d = d.dropna()
        if d.empty:
            continue
        ax.scatter(d["genome_x"], d[feat], s=3, alpha=0.35, linewidths=0)
        add_binned_trend(ax, d["genome_x"], d[feat], bins=250)
        ax.set_ylabel(feat)
        ax.grid(alpha=0.2)

    chroms = sorted(d0["chrom"].dropna().unique().tolist(), key=natural_key)
    xticks = []
    labels = []
    for chrom in chroms:
        sub = d0[d0["chrom"] == chrom]
        if len(sub) == 0:
            continue
        xticks.append((sub["genome_x"].min() + sub["genome_x"].max()) / 2.0)
        labels.append(chrom)

    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right")
    axes[-1].set_xlabel("Genome position (concatenated chromosomes)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "overview.png"), dpi=dpi)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot ArchIE window_stats chr*.stats.tsv.gz files.")
    parser.add_argument(
        "--input-glob",
        default="ArchIE/data/window_stats/chr*.stats.tsv.gz",
        help="Glob for chromosome stats files (default: ArchIE/data/window_stats/chr*.stats.tsv.gz)",
    )
    parser.add_argument(
        "--out-dir",
        default="ArchIE/data/window_stats/plots",
        help="Directory for output PNGs.",
    )
    parser.add_argument(
        "--features",
        default="",
        help="Comma list of columns to plot (e.g. s_star,min_dist_ref or col1,col2). Default plots all detected numeric features.",
    )
    parser.add_argument("--max-points", type=int, default=200000, help="Max points per plot after downsampling.")
    parser.add_argument("--dpi", type=int, default=150, help="PNG resolution.")
    args = parser.parse_args()

    paths = sorted(glob.glob(args.input_glob), key=natural_key)
    if not paths:
        raise SystemExit(f"No files matched: {args.input_glob}")

    os.makedirs(args.out_dir, exist_ok=True)

    df, skipped = load_files(paths)
    df = add_genome_coordinates(df)

    features = get_feature_columns(df)
    if args.features.strip():
        requested = [x.strip() for x in args.features.split(",") if x.strip()]
        features = [f for f in requested if f in df.columns and f in features]
        if not features:
            raise SystemExit("Requested features not found in detected numeric feature columns.")

    for feat in features:
        plot_feature(df, feat, args.out_dir, max_points=args.max_points, dpi=args.dpi)

    plot_overview(df, features, args.out_dir, max_points=args.max_points, dpi=args.dpi)

    print(f"Read files: {len(paths)}")
    print(f"Used rows: {len(df)}")
    print(f"Plotted features: {len(features)}")
    print(f"Output directory: {args.out_dir}")
    if skipped:
        warnings.warn(f"Skipped {len(skipped)} files (in-progress/malformed).")
        for path, reason in skipped[:20]:
            print(f"SKIPPED\t{reason}\t{path}")


if __name__ == "__main__":
    main()
