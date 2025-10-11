#!/usr/bin/env python3
"""
plot_ina226.py  Plot INA226 CSV as 3 separate figures (Voltage, Current, Power)

Options:
  --csv PATH                 (required)
  --rails R1 R2 ...          (optional subset of rails; default all)
  --voltage-unit {mV,V}      (default: mV)
  --current-unit {mA,A}      (default: mA)
  --power-unit   {uW,mW,W}   (default: uW)
  --omit "s1:e1,s2:e2,..."   (optional; relative seconds to remove)
                              e.g., --omit 200:500

Example:
./plot_ina226.py --csv /nfs/versal/rootfs/home/root/app/rails.csv
./plot_ina226.py --csv /nfs/versal/rootfs/home/root/app/rails.csv --power-unit W
./plot_ina226.py --csv /nfs/versal/rootfs/home/root/app/rails.csv --omit 40:1950 --power-unit W
"""

#!/usr/bin/env python3
"""
plot_ina226.py  Plot INA226 CSV as 3 separate figures (Voltage, Current, Power)
with BROKEN X-AXIS for any omitted spans including trimming left/right.

Options:
  --csv PATH                 (required)
  --rails R1 R2 ...          (optional subset of rails; default all)
  --voltage-unit {mV,V}      (default: mV)
  --current-unit {mA,A}      (default: mA)
  --power-unit   {uW,mW,W}   (default: uW)
  --omit "s1:e1,s2:e2,..."   (optional; relative seconds to remove)
                              e.g., --omit 0:100,600:1900,4000:999999  (can trim ends too)
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SUFFIXES = {"mV": "_mV", "mA": "_mA", "uW": "_uW"}

# ---------- helpers ----------


def detect_rails(columns):
    rails = set()
    for col in columns:
        for suf in SUFFIXES.values():
            if col.endswith(suf):
                rails.add(col[: -len(suf)])
    return sorted(rails)


def load_csv(path):
    df = pd.read_csv(path)
    if "ts_ms" not in df.columns:
        raise ValueError("CSV must contain a 'ts_ms' column (epoch milliseconds).")
    df["time"] = pd.to_datetime(df["ts_ms"], unit="ms", errors="coerce")
    df = df.set_index("time").drop(columns=["ts_ms"])
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def subset_metric(df, rails, metric_suffix):
    cols, labels = [], []
    for r in rails:
        col = r + metric_suffix
        if col in df.columns:
            cols.append(col)
            labels.append(r)
    sub = df[cols].copy() if cols else pd.DataFrame(index=df.index)
    sub.columns = labels
    return sub


def apply_units(df_metric, metric, unit):
    # Base units are mV, mA, uW
    if metric == "voltage":
        if unit == "mV":
            return df_metric, "Voltage (mV)"
        if unit == "V":
            return df_metric / 1000.0, "Voltage (V)"
    elif metric == "current":
        if unit == "mA":
            return df_metric, "Current (mA)"
        if unit == "A":
            return df_metric / 1000.0, "Current (A)"
    elif metric == "power":
        if unit == "uW":
            return df_metric, "Power (µW)"
        if unit == "mW":
            return df_metric / 1000.0, "Power (mW)"
        if unit == "W":
            return df_metric / 1_000_000.0, "Power (W)"
    raise ValueError("Invalid unit selection.")


# ---------- omit ? kept segments (broken axis) ----------


def parse_omit(spec):
    if not spec:
        return []
    spans = []
    for part in (p.strip() for p in spec.split(",") if p.strip()):
        if ":" not in part:
            continue
        a, b = part.split(":", 1)
        try:
            s = float(a)
            e = float(b)
            if math.isfinite(s) and math.isfinite(e) and e > s:
                spans.append((s, e))
        except ValueError:
            pass
    if not spans:
        return []
    spans.sort()
    merged = []
    cs, ce = spans[0]
    for s, e in spans[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged


def build_kept_segments(tmin, tmax, omit_spans):
    """
    From total [tmin, tmax] and merged omit_spans, produce a list of kept [a,b] segments.
    Zero-width segments are dropped. Spans are clipped to [tmin,tmax].
    """
    if not omit_spans:
        return [(tmin, tmax)]
    segs = []
    cur = tmin
    for s, e in omit_spans:
        s = max(s, tmin)
        e = min(e, tmax)
        if e <= s:
            continue
        if cur < s:
            segs.append((cur, s))
        cur = max(cur, e)
    if cur < tmax:
        segs.append((cur, tmax))
    # drop degenerate
    return [(a, b) for (a, b) in segs if b > a]


def draw_broken_spines_between(ax_left, ax_right):
    # Hide touching spines
    ax_left.spines["right"].set_visible(False)
    ax_right.spines["left"].set_visible(False)
    ax_left.yaxis.tick_left()
    ax_right.yaxis.tick_right()
    # Diagonal marks
    d = 0.015
    kwargs = dict(transform=ax_left.transAxes, color="k", clip_on=False)
    ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    kwargs.update(transform=ax_right.transAxes)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)
    ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)


def plot_broken_multi(rel_t, df_metric, kept_segments, ylabel, title):
    """
    Create a multi-panel broken axis (N panels) for kept segments,
    with a single legend docked outside on the right.
    """
    if df_metric.empty or not kept_segments:
        return None

    # width ratios proportional to segment durations (avoid 0)
    widths = [max(b - a, 1e-6) for (a, b) in kept_segments]
    n = len(kept_segments)

    # Make a bit wider to leave space for the outside legend
    fig_width = 14 if n <= 3 else 16
    fig, axes = plt.subplots(
        1,
        n,
        sharey=True,
        figsize=(fig_width, 4.5),
        gridspec_kw={"width_ratios": widths, "wspace": 0.05},
    )
    if n == 1:
        axes = [axes]

    # Plot each kept segment; collect handles/labels once
    all_handles = []
    all_labels = []
    for i, (ax, (a, b)) in enumerate(zip(axes, kept_segments)):
        mask = (rel_t >= a) & (rel_t <= b)
        # plot each rail explicitly
        for col in df_metric.columns:
            y = df_metric[col].values
            (h,) = ax.plot(rel_t[mask], y[mask], label=col)
            # collect handle/label from first panel only
            if i == 0:
                all_handles.append(h)
                all_labels.append(col)
        ax.set_xlim(a, b)
        if i == 0:
            ax.set_ylabel(ylabel)
            ax.set_title(title)
        if i == n - 1:
            ax.set_xlabel("Time (s)")  # actual time on each panel

    # Put broken marks between neighbor panels
    for i in range(n - 1):
        draw_broken_spines_between(axes[i], axes[i + 1])

    # Make room on the right for the legend
    fig.subplots_adjust(right=0.80)  # leave ~20% for legend column

    # One combined legend outside the plots
    fig.legend(
        all_handles,
        all_labels,
        loc="center left",
        bbox_to_anchor=(0.82, 0.5),
        frameon=False,
        ncol=1,
    )

    fig.tight_layout(rect=[0.0, 0.0, 0.80, 1.0])  # keep plots inside left 80%
    return fig


# ---------- main ----------


def main():
    ap = argparse.ArgumentParser(
        description="Plot INA226 CSV with multi-panel broken x-axis for omitted spans."
    )
    ap.add_argument("--csv", required=True, help="Path to INA226 CSV file.")
    ap.add_argument(
        "--rails",
        nargs="*",
        default=None,
        help="Optional subset of rails to plot (default: all).",
    )
    ap.add_argument("--voltage-unit", choices=["mV", "V"], default="mV")
    ap.add_argument("--current-unit", choices=["mA", "A"], default="mA")
    ap.add_argument("--power-unit", choices=["uW", "mW", "W"], default="uW")
    ap.add_argument(
        "--omit",
        default=None,
        help='Relative-time spans to remove, e.g. "0:100,600:1900,5000:1e9" (can trim ends too)',
    )
    args = ap.parse_args()

    df = load_csv(args.csv)
    rails_all = detect_rails(df.columns)
    if not rails_all:
        raise SystemExit("No rail columns found with suffixes _mV/_mA/_uW.")
    rails = rails_all if args.rails is None else args.rails

    v_df = subset_metric(df, rails, SUFFIXES["mV"])
    i_df = subset_metric(df, rails, SUFFIXES["mA"])
    p_df = subset_metric(df, rails, SUFFIXES["uW"])
    if v_df.empty and i_df.empty and p_df.empty:
        raise SystemExit("No matching columns for chosen rails.")

    # relative time (seconds)
    t0 = df.index[0]
    rel_sec = (df.index - t0).total_seconds().values
    tmin, tmax = float(rel_sec.min()), float(rel_sec.max())

    # omit ? kept segments
    omit_spans = parse_omit(args.omit)
    kept = build_kept_segments(tmin, tmax, omit_spans)
    # if no omit given or omit cuts everything down to one span, it'll still work:
    if not kept:
        # nothing left  avoid crash
        kept = [(tmin, tmax)]

    # apply units
    v_df, v_ylabel = apply_units(v_df, "voltage", args.voltage_unit)
    i_df, i_ylabel = apply_units(i_df, "current", args.current_unit)
    p_df, p_ylabel = apply_units(p_df, "power", args.power_unit)

    # 3 separate figures, each with N-panel broken x-axis
    plot_broken_multi(rel_sec, v_df, kept, v_ylabel, "INA226 Voltage vs Time")
    plot_broken_multi(rel_sec, i_df, kept, i_ylabel, "INA226 Current vs Time")
    plot_broken_multi(rel_sec, p_df, kept, p_ylabel, "INA226 Power vs Time")

    plt.show()


if __name__ == "__main__":
    main()
