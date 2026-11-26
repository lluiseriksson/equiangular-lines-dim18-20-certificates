"""Comprehensive unit, canonical, forgery, and AST zero-float tests for Paper 87."""

from __future__ import annotations

import ast
import json
from fractions import Fraction
from pathlib import Path
import pytest

from src.equiangular_dim18_20_verifier import (
    SCHEMA,
    exact_psd_and_rank,
    verify_certificate,
    verify_file,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
CERT_DIR = REPO_ROOT / "certificates"
CERT_FILES = sorted(CERT_DIR.glob("*.json"))


def test_sixteen_canonical_certificates_exist() -> None:
    assert len(CERT_FILES) == 16


@pytest.mark.parametrize("cert_path", CERT_FILES, ids=[p.stem for p in CERT_FILES])
def test_canonical_certificate_passes(cert_path: Path) -> None:
    res = verify_file(cert_path)
    assert res["file"] == cert_path.name


def test_zero_float_ast_audit() -> None:
    src_path = REPO_ROOT / "src" / "equiangular_dim18_20_verifier.py"
    tree = ast.parse(src_path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, float):
            pytest.fail(f"float literal found: {node.value}")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "float":
            pytest.fail("float() call found")
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert alias.name not in ("math", "cmath", "numpy", "scipy")
        if isinstance(node, ast.ImportFrom):
            assert node.module not in ("math", "cmath", "numpy", "scipy")


def test_prop53_dim10_also_passes() -> None:
    res = verify_certificate({
        "schema": SCHEMA,
        "id": "T_prop53_d10",
        "mode": "conference_and_cosine_1_7",
        "submode": "prop53_cosine_1_7",
        "dim_d": 10,
        "claimed_lines": 12,
    })
    assert res["gram_rank"] == 10


def test_thm51_additional_table_rows() -> None:
    for q, k, t, r, dm, dp in [
        (5601, 7, 562, 9, 15729, 16703),
        (6001, 23, 108, 23, 29785, 30167),
    ]:
        res = verify_certificate({
            "schema": SCHEMA,
            "mode": "fixed_angle_intervals_thm51",
            "inv_cosine_q": q,
            "excess_k": k,
            "claimed_t": t,
            "claimed_r": r,
            "claimed_d_minus": dm,
            "claimed_d_plus": dp,
            "verify_explicit_seidel_gram": False,
        })
        assert res["interval"] == [dm, dp]


def test_psd_rejects_indefinite_matrix() -> None:
    psd, _ = exact_psd_and_rank([
        [Fraction(1), Fraction(2)],
        [Fraction(2), Fraction(1)],
    ])
    assert not psd


def test_reject_wrong_schema() -> None:
    c = json.loads((CERT_DIR / "C01_dim18_prop23_57_lines.json").read_text(encoding="utf-8"))
    c["schema"] = "BAD"
    with pytest.raises(ValueError, match="schema"):
        verify_certificate(c)


def test_reject_forged_dim18_count() -> None:
    c = json.loads((CERT_DIR / "C01_dim18_prop23_57_lines.json").read_text(encoding="utf-8"))
    c["num_lines"] = 58
    with pytest.raises(ValueError):
        verify_certificate(c)


def test_reject_forged_lemma13_even_q() -> None:
    c = json.loads((CERT_DIR / "C02_lemma13_dim18_A2.json").read_text(encoding="utf-8"))
    c["inv_cosine_q"] = 6
    with pytest.raises(ValueError, match="odd"):
        verify_certificate(c)


def test_reject_forged_lemma48_weights() -> None:
    c = json.loads((CERT_DIR / "C05_golay24_lemma48_weights.json").read_text(encoding="utf-8"))
    c["claimed_nonzero_weights"] = [100, 112, 128]
    with pytest.raises(ValueError, match="claimed_nonzero_weights"):
        verify_certificate(c)


def test_reject_forged_leech_section_axes() -> None:
    c = json.loads((CERT_DIR / "C08_leech_section_dim20_90.json").read_text(encoding="utf-8"))
    c["claimed_num_axes"] = 91
    with pytest.raises(ValueError, match="expected 91"):
        verify_certificate(c)


def test_reject_forged_thm51_dplus() -> None:
    c = json.loads((CERT_DIR / "C09_thm51_q19_k7_gram.json").read_text(encoding="utf-8"))
    c["claimed_d_plus"] = 52
    with pytest.raises(ValueError, match="mismatch"):
        verify_certificate(c)


def test_reject_nonprime_paley() -> None:
    c = json.loads((CERT_DIR / "C13_thm52_paley_p5.json").read_text(encoding="utf-8"))
    c["prime_p"] = 9
    with pytest.raises(ValueError, match="not prime"):
        verify_certificate(c)


def test_reject_forged_paley_lines() -> None:
    c = json.loads((CERT_DIR / "C14_thm52_paley_p13.json").read_text(encoding="utf-8"))
    c["claimed_lines"] = 29
    with pytest.raises(ValueError, match="mismatch"):
        verify_certificate(c)


def test_reject_forged_prop53_lines() -> None:
    c = json.loads((CERT_DIR / "C16_prop53_cosine_1_7_d10_d11.json").read_text(encoding="utf-8"))
    c["claimed_lines"] = 14
    with pytest.raises(ValueError, match="relative bound"):
        verify_certificate(c)


def test_reject_bool_as_int_injection() -> None:
    c = json.loads((CERT_DIR / "C02_lemma13_dim18_A2.json").read_text(encoding="utf-8"))
    c["marking_rank_T"] = True
    with pytest.raises(ValueError, match="exact int"):
        verify_certificate(c)


def test_reject_int_as_bool_injection() -> None:
    c = json.loads((CERT_DIR / "C09_thm51_q19_k7_gram.json").read_text(encoding="utf-8"))
    c["verify_explicit_seidel_gram"] = 1
    with pytest.raises(ValueError, match="bool"):
        verify_certificate(c)
