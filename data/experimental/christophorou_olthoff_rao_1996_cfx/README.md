# NIST-evaluated CFx secondary-ionization curves

Tables 31--33 of Christophorou, Olthoff, and Rao, *J. Phys. Chem. Ref.
Data* **25**, 1341--1388 (1996), DOI `10.1063/1.555986`, reproduce the
Tarnovsky/Becker measurements of electron-impact ionization of free `CF`,
`CF2`, and `CF3` radicals. The six energy-resolved channels are transcribed
here in SI units. The official NIST scan has SHA256
`381e368840e28c84bb03eb4e684691e5c330c8c4f592bc666ae13bc292a39244`.

The journal pages containing Tables 31--33 were rendered at 240 dpi and
inspected at original resolution. Run the checksum and committed-data replay:

```bash
python scripts/extract_nist_cfx_ionization.py \
  --source-pdf /path/to/jpcrd512.pdf \
  --render-directory /path/to/rendered/pages \
  --check
```

Table 31 reports 15%, 16%, and 18% overall uncertainty for `CF`, `CF2`, and
`CF3`, respectively. Tables 32 and 33 report 20% for the `CF3` dissociative
branches and 16% for the net `CF2 -> CF+` curve. The `F+` values at 70 eV are
single-energy anchors and are not promoted into rate curves. In particular,
the Table-33 `CF+` signal contains single-ion and ion-pair onsets; retaining it
as net measured `CF+` production is valid for the measured light-ion current,
but a later charge-resolved `F+` closure must split those branches.

These tables close secondary light-ion production. They do not supply C4F6
neutral fragmentation, heavy-ion branching, residence/wall losses, absolute
reactor flux, a Krueger boundary, or feature depth.
