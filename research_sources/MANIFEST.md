# Research sources archive (durable — tmp/ is gitignored and wipeable)

Text extractions of fetched primary sources + digitized figure data. Raw PDFs
(211MB, tmp/pdfs/) are NOT committed except chang_thesis.pdf (load-bearing for
the angular curve + thresholds; MIT DSpace 1721.1/50356, open bitstream).
Refetch routes are in the RESEARCH_* docs that cite each source (curl -k for
TLS-blocked; the HAL Anubis proof-of-work solver was at scratchpad/getpdf.py —
recreate from RESEARCH_LER_EXPERIMENTAL_GATES notes if needed).

## thesis_extracts/ (grep-ready, line refs in RESEARCH docs point here)
- HG_jap97.txt
- HG_jvstb97.txt
- hamilton_2018_cl2_dissociation.txt
- kemaneci_2014_chlorine_global.txt
- chantry_1987_wall_diffusion.txt
- economou_chlorine_neutral_transport.txt
- economou_chlorine_charged_transport.txt
- guha_2008_cl_wall_interactions.txt
- malyshev_1998_lam_chlorine_dissociation.txt
- stafford_2009_cl_wall_recombination.txt
- stafford_2010_cl_wall_recombination.txt
- wang_olthoff_1999_ion_energy.txt
- lan_jeon_2014_c4f6_verified_excerpt.txt
- casey_2021_pt_verified_excerpt.txt
- an_2026_nnp_etch_verified_excerpts.txt
- arts_2021_apr_angular_verified_excerpt.txt
- Konina_Kseniia_PhD_Thesis_2024.txt
- Lanham_Steven_PhD_Thesis_2022.txt
- Qu_Chenhui_PhD_Thesis_2020.txt
- coatings2023_bowing_narrowing.txt
- deboer-2002.txt
- huang_thesis.txt
- huard_chad_phd_thesis.txt
- jeong_2023.txt
- krueger-2024.txt
- krueger_thesis.txt
- lam_shen_lill_jjap2023.txt
- logue_michael_phd_thesis.txt
- mask_geometry_micromachines_2023.txt
- nanomaterials2024_necking.txt
- nanomaterials2024_necking_layout.txt
- skku_acl_2013.txt
- song_sangheon_phd_thesis.txt
- tian_peng_phd_thesis.txt
- tuwien_rodrigues_2023_fc_silica.txt
- wang_mingmei_phd_thesis.txt
- zhang_yiting_phd_thesis.txt

## digitized/
- krueger_fig7a_simulated_aperture.csv / fig7b_experimental (aperture vs depth,
  600dpi digitization; source of the 38.8@271 / 39.0@200 neck targets)
- hamilton_2018_cl2_* (exact liborigin extraction of the official CC-BY OPJ:
  eight state cross sections, independent total, and Figure-5 reference rates)
- Wang–Olthoff 1999 Figure 9 is under
  `data/experimental/wang_olthoff_1999/`: 24 absolute/species-resolved
  chlorine-ion markers, with a 600-dpi PIL audit and exact NIST PDF/render
  hashes.
- Stafford et al. 2010 Figure 8 is under
  `data/experimental/stafford_2010/`: 39 conditioned-wall chlorine
  recombination markers, with a 600-dpi human-vision/PIL audit and exact
  IUPAC PDF/render hashes.
- Malyshev et al. 1998 Figures 3, 7--8, and 11 are under
  `data/experimental/malyshev_1998_lam/`: 62 measured electron-temperature
  markers, 38 measured Cl2-dissociation markers, and 27 volume-averaged
  electron-density markers from a commercial Lam Alliance reactor. The
  replays extract and checksum the publisher PDF's native one-bit figure
  rasters before the PIL marker audits; fitted curves, inseparable marker
  overlaps, and the documented dissociation anomaly are excluded.
- extract_mouth_profiles.py (regenerator)

## library/ — the internal literature library (START HERE)

`LIBRARY.md` indexes every source this project has cited, grouped by topic, with a
reverse index mapping each petch constant/law/decision to the source that fixed it.
One `library/<bibkey>.md` per source carries citation, DOI, retrieval route, fetch
status, and the claims table (verbatim rows relocated from the RESEARCH_/RESULTS_
docs, each pointing back at `doc:line`).

Two binding conventions: (1) fetch => extract + library entry in the SAME commit;
(2) provenance strings name bibkeys so a constant is one grep from its evidence.
