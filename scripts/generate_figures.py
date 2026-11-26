"""Generate 3 publication-quality figures for Paper 87 (EQUIANGULAR-D18-20-CERT-1.0)."""

from __future__ import annotations

import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    out = Path("paper/figures")
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.labelsize": 10,
        "axes.titlesize": 10.5,
        "legend.fontsize": 8.5,
        "figure.dpi": 220,
    })

    # Fig 1: Settled real equiangular line maxima N(d) across d=15..24 & nested Leech sections
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.3))
    dims = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    exact_N = [36, 40, 48, 57, 72, 90, 126, 176, 276, 276]
    prev_ub = [36, 40, 48, 59, 74, 94, 126, 176, 276, 276]
    ax1.plot(dims, prev_ub, ls="--", marker="s", ms=5, color="#999999", label="Prior upper bound [GS24]")
    ax1.plot(dims, exact_N, lw=2.0, marker="o", ms=5.5, color="#1f4e79", label="Exact $N(d)$ (settled)")
    ax1.scatter([18, 19, 20], [57, 72, 90], color="#c0392b", s=65, zorder=6, label="Newly settled $d=18,19,20$")
    ax1.set_xlabel(r"Dimension $d$")
    ax1.set_ylabel(r"Max equiangular lines $N(d)$")
    ax1.set_title(r"(a) Exact Maxima $N(d)$ Through $d \leq 24$")
    ax1.grid(alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    sec_dims = [23, 22, 21, 20]
    axes_cnt = [276, 176, 126, 90]
    ghw = [0, 100, 150, 186]
    x = np.arange(len(sec_dims))
    w = 0.36
    b1 = ax2.bar(x - w / 2, axes_cnt, width=w, color="#1f4e79", label=r"Retained axes $276 - d_j(C)$")
    b2 = ax2.bar(x + w / 2, ghw, width=w, color="#d95f02", label=r"Gen. Hamming wt $d_j(C)$")
    for bars in (b1, b2):
        for b in bars:
            h = b.get_height()
            if h > 0:
                ax2.text(b.get_x() + b.get_width() / 2, h + 4, f"{int(h)}", ha="center", va="bottom", fontsize=8)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"d={d}\n(j={23-d})" for d in sec_dims])
    ax2.set_ylim(0, 315)
    ax2.set_ylabel("Count")
    ax2.set_title("(b) Nested Sharp Leech Sections (Sec. 4.6)")
    ax2.grid(axis="y", alpha=0.25)
    ax2.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    fig.savefig(out / "fig1_maxima_and_leech_sections.png", dpi=220)
    plt.close(fig)

    # Fig 2: Theorem 5.1 sharp fixed-angle intervals d_+ - d_- + 1 vs odd q
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.3))
    L7 = 72.0 / (7.0 * math.sqrt(10.0) + 2.0)
    L23 = 600.0 / (23.0 * math.sqrt(26.0) + 2.0)

    def interval_len(q_val: int, k: int) -> int:
        M, b, Lk = (28, 10, L7) if k == 7 else (276, 56, L23)
        Q = q_val + 1
        t = 2 * ((Q + 2 * b - 1) // (2 * b))
        r = (b * t - Q) // 2
        dm = M * t - k
        dp = int(math.floor(Lk * Q)) - k - 1
        return (dp - dm + 1) if (r < t and dm <= dp) else 0

    qs7 = list(range(319, 1201, 2))
    lens7 = [interval_len(q, 7) for q in qs7]
    ax1.plot(qs7, lens7, lw=1.2, color="#1f4e79", label=r"$k=7$: slope $L_7 - 28/10 \approx 0.1831$")
    ax1.axvline(319, ls=":", color="#c0392b", lw=1.4, label=r"Threshold $q \geq 319$")
    ax1.set_xlabel(r"Odd denominator $q$")
    ax1.set_ylabel(r"Interval size $d_+ - d_- + 1$")
    ax1.set_title(r"(a) Theorem 5.1 ($k=7$, $N_{1/q}(d)=d+7$)")
    ax1.grid(alpha=0.25)
    ax1.legend(loc="upper left", frameon=True)

    qs23 = list(range(5599, 9001, 2))
    lens23 = [interval_len(q, 23) for q in qs23]
    ax2.plot(qs23, lens23, lw=0.9, color="#2e6b3e", label=r"$k=23$: slope $L_{23} - 276/56 \approx 0.1017$")
    ax2.axvline(5599, ls=":", color="#c0392b", lw=1.4, label=r"Threshold $q \geq 5599$")
    ax2.set_xlabel(r"Odd denominator $q$")
    ax2.set_ylabel(r"Interval size $d_+ - d_- + 1$")
    ax2.set_title(r"(b) Theorem 5.1 ($k=23$, $N_{1/q}(d)=d+23$)")
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper left", frameon=True)
    fig.tight_layout()
    fig.savefig(out / "fig2_fixed_angle_intervals.png", dpi=220)
    plt.close(fig)

    # Fig 3: Lemma 4.8 evaluation-code weight spectrum & Theorem 5.2 Paley conference family
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.3))
    wts = [100, 112, 128, 132, 144, 148, 164, 176, 276]
    y_pos = np.ones(len(wts))
    ax1.stem(wts, y_pos, linefmt="#1f4e79", markerfmt="o", basefmt=" ")
    ax1.axvline(100, ls="--", color="#c0392b", lw=1.4, label=r"$d_1(C)=100$ ($y^2=5/3$)")
    ax1.axvline(186, ls=":", color="#2e6b3e", lw=1.6, label=r"$d_3(C)=186$ (Prop. 4.9)")
    ax1.set_yticks([])
    ax1.set_ylim(0, 1.35)
    ax1.set_xlim(80, 290)
    ax1.set_xlabel("Allowed nonzero codeword weight in $C$")
    ax1.set_title("(a) Leech Evaluation Code Weights (Lemma 4.8)")
    ax1.grid(axis="x", alpha=0.25)
    ax1.legend(loc="upper right", frameon=True)

    ps = [5, 13, 17, 29, 37]
    ds = [3 * (p + 1) // 2 for p in ps]
    ns = [2 * (p + 1) for p in ps]
    ax2.plot(ds, ns, lw=2.0, marker="o", color="#8e44ad", label=r"Paley $t=2$: $N = \frac{4}{3} d$")
    ax2.scatter(ds[:3], ns[:3], color="#d95f02", s=58, zorder=6, label="Certified C13–C15 ($p=5,13,17$)")
    for p, d, n in zip(ps[:3], ds[:3], ns[:3]):
        ax2.annotate(f"p={p} ({d},{n})", (d, n), textcoords="offset points", xytext=(7, -4), fontsize=8.0)
    ax2.set_ylim(5, 82)
    ax2.set_xlabel(r"Dimension $d = 3(p+1)/2$")
    ax2.set_ylabel(r"Exact line count $2(p+1)$")
    ax2.set_title(r"(b) Quadratic Irrational Cosines (Thm. 5.2)")
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper left", frameon=True)
    fig.tight_layout()
    fig.savefig(out / "fig3_golay_weights_and_paley.png", dpi=220)
    plt.close(fig)
    print("Generated 3 figures for Paper 87.")


if __name__ == "__main__":
    main()
