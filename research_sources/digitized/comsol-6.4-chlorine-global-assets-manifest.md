# COMSOL 6.4 chlorine-global-model asset manifest

## Retrieval route and retention boundary

- Official model page:
  `https://www.comsol.com/model/chlorine-discharge-global-model-47611`
- The page's public `models/get-the-files` route supplied the model bundle and
  text assets on 2026-08-07.
- The raw COMSOL files are **not redistributed in this repository**. This
  manifest pins identity and inspected support only; any future use must
  retrieve the assets from the official page and verify the hashes.

## Identity ledger

| official asset | SHA-256 | inspected content |
|---|---|---|
| `chlorine_global_model.mph` | `c61b8f244437fcf2d9c48798a12f7ebc4d4d3c5bca2d45a20ef264a78575edce` | ZIP-contained `dmodel.xml`; 38 forward volume rows including two elastic rows, plus eight explicit reverse rows |
| `chlorine_global_model_variables_1.txt` | `997840c6f744e510ab03d0bfc91b39716f60a76aa489710ea9c098fc3e3c73f4` | level/event variable file |
| `chlorine_global_model_variables_2.txt` | `06caf3c8054b9333a989d2551bb581c86189fe38bdc3276e43b2963546e30734` | auxiliary variable file |
| `Cl2_mom_xsec.txt` | `6bf31f960206ce09ab7603df2e6fbbd4688340a737f911e4822684986118782d` | 49 numeric records, `0.01983--29.17 eV` |
| `Cl_mom_xsec.txt` | `6c1391eda01d90ac504e73e3f5306b510e8a5d721860b03b50cdba39edd29661` | 28 numeric records, `0--25 eV` |

The 6.2, 6.3, and 6.4 downloads of both momentum-transfer tables are
byte-identical. The Cl2 file has 48 newline characters but 49 numeric records
because its final record has no terminating newline; row counts must therefore
parse records rather than use `wc -l`.

## Raw-model findings frozen here

- Reverse rows 24--31 multiply forward fits by `exp(deltaE/Te)` and expose no
  statistical-weight factor.
- `ev1/ev2/ev3 = 0.07/0.14/0.21 eV`.
- `eCl12/eCl52 = 1.35/10.17 eV` are used directly as atomic excitation gaps,
  despite Figure 10 placing ground Cl at `1.25 eV`.
- `eionCl = 14.25 eV` and `ediss = 4 eV` are used as event energies.
- Forward reaction 20 uses `exp(-13.29/Te)`, diverging from the primary
  paper's printed `exp(-13.19/Te)`.

These observations define an exact implementation-reproduction target. They
do not make the bundled elastic tables or energy variables evaluated physical
data.
