# NIST C4F6 electron-ionization fragmentation board

This directory records a vision-audited numerical transcription of the NIST
Chemistry WebBook SRD 69 electron-ionization mass spectrum for
1,3-hexafluorobutadiene (`C4F6`, CAS `685-63-2`, NIST MS `5987`). The source
pixels are not committed because the WebBook explicitly prohibits downloading
the spectrum. The local 800 x 600 printable image is checksum-locked in the
manifest, and every red stick is replayed by
`scripts/digitize_nist_c4f6_mass_spectrum.py`.

The base peak is `C3F3+`; intact `C4F6+` is about 44 relative-intensity units.
Together those two heavy channels carry 70% of the monoisotopic stick
intensity. Direct `CF+`, `CF2+`, and `CF3+` are only 15.23, 1.59, and 4.77.
Therefore Lan--Jeon's aggregate C4F6 ionization cross section cannot honestly
be mapped only onto the three light `CFx+` ions.

This is a molecular fragmentation prior, not a plasma composition. Benck's
absolute C4F6/Ar reactor board shows a much larger `CF2+/CF+` ratio and a large
unresolved ion current; Kim's wafer-facing CCP spectrum also directly observes
`C3F3+` and heavier ions. A reactor model must retain direct heavy products,
secondary neutral-fragment ionization, ion-neutral chemistry, wall return, and
mass-resolved transport before it can claim a C4F6 wafer boundary.

Replay with the locally held image:

```bash
python scripts/digitize_nist_c4f6_mass_spectrum.py \
  --source-image /path/to/nist_c4f6_mass_spectrum.png \
  --overlay /tmp/nist_c4f6_mass_spectrum_overlay.png
```

The numbers must not be used as an absolute ion flux, transplanted into
Krueger's reactor, or fitted to the 825 nm endpoint.
