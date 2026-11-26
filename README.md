# Formal Certification of Real Equiangular Lines in Dimensions 18–20 and Sharp Fixed-Angle Bounds (`EQUIANGULAR-D18-20-CERT-1.0`)

**Author:** Lluis Eriksson  
**License:** `MIT` (code, Lean 4 proofs, certificates, tests); `CC-BY-4.0` (manuscript, figures)  
**Target Preprint:** Valjakas Djaljapayan, Zhen-Qi Liao, and Wei-Hsuan Yu (in collaboration with OpenAI's GPT-6), *"Real equiangular lines in dimensions 18–20 and sharp fixed-angle bounds"*, `arXiv:2611.34855v1` (25 November 2026).

## Verification Summary
- **Constructive Lean 4.35.0 kernel** (`lean/EquiangularDim1820Cert.lean`): 18/18 theorems (`0 sorry`, `0` custom axioms, `0 Classical.choice`, `0` warnings).
- **Zero-float exact rational Python verifier** (`src/equiangular_dim18_20_verifier.py`): AST-audited zero-float implementation across 4 modes (Proposition 2.3 explicit 57-line $\mathbb{R}^{18}$ Gram PSD/rank-18 witness & Lemma 1.3 paired even-lattice reductions, extended binary Golay code $G_{24}$ + Lemma 4.8 evaluation-code weights + nested sharp Leech sections $(22,176)/(21,126)/(20,90)$, Theorem 5.1 exact quadratic-surd fixed-angle intervals & Seidel-block Gram matrices, Theorem 5.2 symmetric Paley conference matrices & Proposition 5.3 $N_{1/7}(10)=12, N_{1/7}(11)=13$).
- **Certificates & Tests:** 16/16 canonical JSON certificates (`certificates/C01`–`C16`), 32/32 pytest unit + forgery + AST checks.

## Replay Commands
```bash
(cd lean && lake build)
python -m src.equiangular_dim18_20_verifier certificates
pytest -q
```
