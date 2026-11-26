"""Generate all 16 canonical certificates for Paper 87 (EQUIANGULAR-D18-20-CERT-1.0)."""

from __future__ import annotations

import json
from pathlib import Path
from src.equiangular_dim18_20_verifier import SCHEMA, verify_certificate


def main() -> None:
    out = Path("certificates")
    out.mkdir(parents=True, exist_ok=True)
    certs = [
        # Mode 1: Dim 18 construction (57 lines) & Lemma 1.3 paired lattice reductions
        {
            "schema": SCHEMA,
            "id": "C01_dim18_prop23_57_lines",
            "mode": "dim18_construction_and_reduction",
            "submode": "prop23_57_lines",
            "num_lines": 57,
            "ambient_dim": 18,
            "inv_cosine_q": 5,
        },
        {
            "schema": SCHEMA,
            "id": "C02_lemma13_dim18_A2",
            "mode": "dim18_construction_and_reduction",
            "submode": "lemma13_paired_lattice",
            "span_dim_r": 18,
            "inv_cosine_q": 5,
            "marking_rank_T": 2,
            "claimed_w2": 6,
            "claimed_x2": 4,
            "claimed_xw": 3,
            "embeds_in_rank24": True,
        },
        {
            "schema": SCHEMA,
            "id": "C03_lemma13_dim19_A1",
            "mode": "dim18_construction_and_reduction",
            "submode": "lemma13_paired_lattice",
            "span_dim_r": 19,
            "inv_cosine_q": 5,
            "marking_rank_T": 1,
            "claimed_w2": 6,
            "claimed_x2": 4,
            "claimed_xw": 3,
            "embeds_in_rank24": True,
        },
        {
            "schema": SCHEMA,
            "id": "C04_lemma13_dim20_leech",
            "mode": "dim18_construction_and_reduction",
            "submode": "lemma13_paired_lattice",
            "span_dim_r": 20,
            "inv_cosine_q": 5,
            "marking_rank_T": 0,
            "claimed_w2": 6,
            "claimed_x2": 4,
            "claimed_xw": 3,
            "embeds_in_rank24": True,
        },
        # Mode 2: Extended binary Golay code G_24, Lemma 4.8 weights, nested Leech sharp sections
        {
            "schema": SCHEMA,
            "id": "C05_golay24_lemma48_weights",
            "mode": "leech_evaluation_code_dim20",
            "submode": "golay_and_lemma48",
            "claimed_nonzero_weights": [100, 112, 128, 132, 144, 148, 164, 176, 276],
        },
        {
            "schema": SCHEMA,
            "id": "C06_leech_section_dim22_176",
            "mode": "leech_evaluation_code_dim20",
            "submode": "nested_sharp_sections",
            "section_dim": 22,
            "claimed_num_axes": 176,
            "claimed_generalized_hamming_weight": 100,
        },
        {
            "schema": SCHEMA,
            "id": "C07_leech_section_dim21_126",
            "mode": "leech_evaluation_code_dim20",
            "submode": "nested_sharp_sections",
            "section_dim": 21,
            "claimed_num_axes": 126,
            "claimed_generalized_hamming_weight": 150,
        },
        {
            "schema": SCHEMA,
            "id": "C08_leech_section_dim20_90",
            "mode": "leech_evaluation_code_dim20",
            "submode": "nested_sharp_sections",
            "section_dim": 20,
            "claimed_num_axes": 90,
            "claimed_generalized_hamming_weight": 186,
        },
        # Mode 3: Theorem 5.1 sharp fixed-angle intervals (all 6 table rows + explicit Gram for q=19)
        {
            "schema": SCHEMA,
            "id": "C09_thm51_q19_k7_gram",
            "mode": "fixed_angle_intervals_thm51",
            "inv_cosine_q": 19,
            "excess_k": 7,
            "claimed_t": 2,
            "claimed_r": 0,
            "claimed_d_minus": 49,
            "claimed_d_plus": 51,
            "verify_explicit_seidel_gram": True,
        },
        {
            "schema": SCHEMA,
            "id": "C10_thm51_q111_k23",
            "mode": "fixed_angle_intervals_thm51",
            "inv_cosine_q": 111,
            "excess_k": 23,
            "claimed_t": 2,
            "claimed_r": 0,
            "claimed_d_minus": 529,
            "claimed_d_plus": 539,
            "verify_explicit_seidel_gram": False,
        },
        {
            "schema": SCHEMA,
            "id": "C11_thm51_q321_k7",
            "mode": "fixed_angle_intervals_thm51",
            "inv_cosine_q": 321,
            "excess_k": 7,
            "claimed_t": 34,
            "claimed_r": 9,
            "claimed_d_minus": 945,
            "claimed_d_plus": 952,
            "verify_explicit_seidel_gram": False,
        },
        {
            "schema": SCHEMA,
            "id": "C12_thm51_q5601_k23",
            "mode": "fixed_angle_intervals_thm51",
            "inv_cosine_q": 5601,
            "excess_k": 23,
            "claimed_t": 102,
            "claimed_r": 55,
            "claimed_d_minus": 28129,
            "claimed_d_plus": 28155,
            "verify_explicit_seidel_gram": False,
        },
        # Mode 4: Theorem 5.2 Paley conference matrices + Proposition 5.3 N_{1/7}(10)=12, N_{1/7}(11)=13
        {
            "schema": SCHEMA,
            "id": "C13_thm52_paley_p5",
            "mode": "conference_and_cosine_1_7",
            "submode": "thm52_paley_conference",
            "prime_p": 5,
            "repeat_t": 2,
            "claimed_dim": 9,
            "claimed_lines": 12,
        },
        {
            "schema": SCHEMA,
            "id": "C14_thm52_paley_p13",
            "mode": "conference_and_cosine_1_7",
            "submode": "thm52_paley_conference",
            "prime_p": 13,
            "repeat_t": 2,
            "claimed_dim": 21,
            "claimed_lines": 28,
        },
        {
            "schema": SCHEMA,
            "id": "C15_thm52_paley_p17",
            "mode": "conference_and_cosine_1_7",
            "submode": "thm52_paley_conference",
            "prime_p": 17,
            "repeat_t": 2,
            "claimed_dim": 27,
            "claimed_lines": 36,
        },
        {
            "schema": SCHEMA,
            "id": "C16_prop53_cosine_1_7_d10_d11",
            "mode": "conference_and_cosine_1_7",
            "submode": "prop53_cosine_1_7",
            "dim_d": 11,
            "claimed_lines": 13,
        },
    ]
    for c in certs:
        verify_certificate(c)
        (out / f"{c['id']}.json").write_text(json.dumps(c, indent=2) + "\n", encoding="utf-8")
    print(f"Generated and verified {len(certs)} certificates in {out}")


if __name__ == "__main__":
    main()
