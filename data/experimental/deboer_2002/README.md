# de Boer 2002 Figure 9 direct experimental evidence

Primary source: M. J. de Boer et al., “Guidelines for etching silicon MEMS structures using
fluorine high-density plasmas at cryogenic temperatures,” *Journal of Microelectromechanical
Systems* **11**, 385–401 (2002), DOI `10.1109/JMEMS.2002.800928`, Figure 9.

Author-university copy used for extraction:
`https://ris.utwente.nl/ws/files/6683265/guidelines-etching-boer.pdf`.

- source PDF SHA-256: `45c245a9b19671f532945155dc16c3e00d35464eb8e49480a09f90a90498ff6c`
- PDF page: 13 (printed page 397)
- extraction command: `pdfimages -f 13 -l 13 -png guidelines-etching-boer.pdf deboer-fig`
- extracted Figure 9 image: image index 6, 1500 x 911 px
- extracted image SHA-256: `0f78ae30e5cc2e128f4fdb84217551fe350bd7696966c6ea40233f70a9a765c4`

The copyrighted PDF and image are not vendored. The CSV retains their checksums, every marker center,
and both linear axis transforms so the digitization is replayable by a lawful source holder. A
conservative 0.30 µm opening and 0.50 µm depth digitization bound exceeds two source-image pixels.
The paper does not report a statistical measurement uncertainty for these markers, so it remains
`not_reported` and cannot be silently replaced by the digitization bound.

The original preregistered split was intentionally small and fixed before the first direct Figure 9
engine score:

- each curve’s widest-opening marker is a boundary/open-rate anchor and is not scored;
- the narrowest 12.5 minute marker is the sole fluorine-sticking calibration condition;
- the other 12 markers were held-out width/time transfers for that first score.

That score has now occurred.  Consequently **none of the Figure 9 markers are held out anymore**.
The old split remains in the checksummed CSV and frozen protocol solely to replay the historical
test.  All 16 markers are development evidence for subsequent mechanism repair; a different,
preregistered experiment is required for any new validation claim.

## Boundary-input completeness

The Figure 9 caption reports 600 W ICP power, 10 mtorr, 90 sccm SF6, 3.0 W CCP with a -30 V
self-bias, a SiO2 mask, and 5% exposed area.  It does **not** report the oxygen flow, electrode
temperature, mask thickness/profile, ion flux, IEDF, or IADF for this figure.  The parent paper also
does not report statistical uncertainty for the plotted depths.  Those omissions prevent an
absolute predictive-chemistry claim from Figure 9 alone and must not be filled by hidden fit
parameters.

The mechanism diagnosis is nevertheless unusually strong.  The directly related primary study by
Jansen et al., *Microelectronic Engineering* 35, 45-50 (1997), DOI
`10.1016/S0167-9317(96)00142-6`, used horizontal-trench controls to show that radical depletion was
small in this SF6/O2 regime, that inhibitor depletion would produce the wrong (inverse-lag) sign,
and that angular/electrostatic ion depletion controls the observed RIE lag.  Development after the
first score therefore corrects the ion boundary input and kinetic transport before adding any new
reaction channel.  No photon flux or photon-assisted yield is reported, so no photon channel is
introduced.

This direct dataset supersedes any description of `1.00 / 0.43 / 0.29 / 0.20` at AR
`0 / 10 / 20 / 40` as digitized de Boer experiment. That older curve is the Blauw/Clausing model
evaluated with a fitted sticking coefficient. It remains useful as a model cross-check, but it is not
raw Figure 9 evidence.
