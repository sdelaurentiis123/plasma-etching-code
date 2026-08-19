# christophorou-olthoff-rao-1996-cfx

**NIST-evaluated electron ionization of free CFx radicals**

- **Citation:** L. G. Christophorou, J. K. Olthoff, and M. V. V. S. Rao,
  “Electron Interactions with CF4,” *Journal of Physical and Chemical
  Reference Data* **25**, 1341--1388 (1996).
- **DOI:** `10.1063/1.555986`
- **Official NIST PDF:**
  `https://www.nist.gov/system/files/documents/srd/jpcrd512.pdf`
- **PDF SHA256:**
  `381e368840e28c84bb03eb4e684691e5c330c8c4f592bc666ae13bc292a39244`
- **Executable tables:**
  `data/experimental/christophorou_olthoff_rao_1996_cfx/`
- **Status:** PRIMARY NIST FULL TEXT + TABLES 31--33 VISUALLY AUDITED AT
  240 DPI

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| C1 | Table 31 tabulates parent `CF3+`, `CF2+`, and `CF+` production from their free radicals through 200 eV, with 18%, 16%, and 15% overall uncertainties. | Supplies measured Maxwellian rates instead of copied CFx daughter-ionization fits. |
| C2 | Table 32 tabulates dissociative `CF3 -> CF2+ + F` and `CF3 -> CF+ + 2F` curves with 20% overall uncertainty. | Makes light-ion production depend on the evolved radical inventory and EEDF rather than direct C4F6 branching alone. |
| C3 | Table 33 tabulates net `CF2 -> CF+` and states that the curve has two onsets: single-ion `CF+ + F` and ion-pair `CF+ + F+`. | The net `CF+` curve is usable for CF+ production; exact electron/F+ bookkeeping requires a later source-resolved split. |
| C4 | Tables 32 and 33 give `F+` only at 70 eV, with 30% uncertainty. | The anchors are retained as evidence and forbidden from becoming invented energy-dependent curves. |
| C5 | The NIST review emphasizes that CF3, CF2, CF and their ions are abundant/reactive CF4-discharge fragments. | Supports a secondary-fragment network; it does not supply C4F6 parent fragmentation or a reactor state. |

## Executable decision

`build_nist_1996_cfx_ionization_network()` exposes exactly six measured,
atom/charge-conserving electron-impact channels. No Benck current and no
Krueger depth selects a rate. The provider closes one major topology gap, but
it does not authorize an absolute C4F6/Ar reactor flux or feature depth.
