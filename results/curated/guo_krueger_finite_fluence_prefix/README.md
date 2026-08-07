# Guo/Kwon–Krüger finite-fluence feature prefix

Status: **time and space gates passed at 0.5 s; a 60 s no-fit forecast is
authorized numerically, but the Krüger boundary remains nonpredictive**.

This gate advances the source-fixed Guo/Kwon translating-layer state from a
bare-oxide initial condition instead of replacing every feature face with its
steady composition at the first step. No yield, flux, translating-layer
capacity, or endpoint-depth parameter was adjusted.

## Time convergence at 10 nm

| nominal step (s) | steps | depth (nm) | rate (nm/s) | mask opening (nm) |
|---:|---:|---:|---:|---:|
| 0.25 | 2 | 5.62791033 | 11.2558207 | 85.0892158 |
| 0.125 | 4 | 5.91041642 | 11.8208328 | 85.4576138 |
| 0.0625 | 8 | 5.96472993 | 11.9294599 | 85.6848977 |

The observed order is 2.37890. The Richardson limit is 5.97765742 nm and the
fine result is 0.2163% from that limit. The medium-to-fine change is 0.9106%;
the coarse-to-medium change is still 4.7798%. A long trajectory should
therefore extend the 0.125 s checkpoint rather than the cheaper 0.25 s one.

## Space convergence at 0.0625 s

| spacing | steps | depth (nm) | rate (nm/s) | mask opening (nm) |
|---:|---:|---:|---:|---:|
| 10 nm | 8 | 5.96472993 | 11.9294599 | 85.6848977 |
| 5 nm | 8 | 6.13490241 | 12.2698048 | 85.4152643 |

The fine-minus-coarse depth difference is +2.8530%, inside the preregistered
5% spatial gate. Both final profile images were inspected at native
resolution: the floor and mouth are connected and symmetric, and neither
contains a disconnected gas pocket or a topology/plotting artifact.

Material ledgers close exactly in all four trajectories. The maximum
neutral-radiosity relative balance residual is 1.804e-12.

The 10 nm fine-prefix rate is 13.24% below the 13.75 nm/s mean required to
reach 825 nm in 60 s; the 5 nm rate is 10.77% below. These are early-time
observations, **not linear endpoint forecasts**. Aspect-ratio transport,
surface-state evolution, and mouth evolution remain nonlinear.

## Evidence boundary

Passing the numerical gate does not identify the physical boundary. Krüger
publishes only an aggregate energetic-ion row and combined IEAD, not the
ion-species mixture or stable C4F6 wafer flux. C2F3/C3F4 require declared
topology transfers outside Guo's printed neutral list, most of Krüger's IEAD
is above the Guo/Yin regression board, and the 2.5 nm translating-layer
capacity is a source-bounded cross-chemistry transfer. The mask model also
contains unmeasured density and reduced-film parameters.

The authorized next result is therefore a **no-fit published-boundary
sensitivity forecast**, not a Tier-A prediction. It answers whether this
fully declared transfer happens to match 825 nm; it cannot prove that the
unpublished Krüger boundary has been reconstructed.

## Receipts

`audit.json` records every configuration hash, raw audit/profile SHA-256,
implementation revision, convergence calculation, and physical blocker. The
full local artifacts are checksum-bound at:

- `/private/tmp/krueger_guo_transient_dx10`
- `/private/tmp/krueger_guo_transient_dt125_dx10`
- `/private/tmp/krueger_guo_transient_dt0625_dx10`
- `/private/tmp/krueger_guo_transient_dt0625_dx5`

Rebuild the receipt with:

```bash
python scripts/audit_guo_krueger_finite_fluence_prefix.py
```
