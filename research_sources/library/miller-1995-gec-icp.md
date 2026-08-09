# miller-1995-gec-icp

**Primary GEC inductively coupled source geometry and diagnostic reference**

- **Citation:** P. A. Miller, G. A. Hebner, K. E. Greenberg, P. D. Pochan,
  and B. P. Aragon, “An Inductively Coupled Plasma Source for the Gaseous
  Electronics Conference RF Reference Cell,” *Journal of Research of the
  National Institute of Standards and Technology* **100**, 427--439 (1995).
- **DOI:** `10.6028/jres.100.032`
- **Official NIST PDF:**
  `https://nvlpubs.nist.gov/nistpubs/jres/100/4/j14mil.pdf`
- **PDF SHA-256:**
  `e67774edeefae2fb3c50c8343479519a88d6020f3022961b0ded6708cad5ceb9`
- **Status:** PRIMARY NIST FULL TEXT; FIGURE 1 VISUALLY AUDITED AT 300 DPI

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| M1 | Figure 1 gives a 165.1 mm lower-electrode extension, 40.5 mm window-to-electrode gap, 15 mm probe height, 9.5 mm silica-window thickness, 114.3 mm upper clear span, and a 248 mm-ID standard cell body. | Exact geometry receipt for the Wise Figure-3 board. The spatial plasma domain is not silently replaced by the 248 mm body because the electrode/window assembly and side ports make the boundary non-cylindrical outside the central gap. |
| M2 | The source uses a five-turn planar coil of 3 mm copper tubing; the hand-wound coil outside diameter is approximately 10 cm. | Coil-side field support and an uncertainty boundary: the winding geometry is approximate and generator power does not determine the plasma-side tangential field. |
| M3 | The radial Langmuir-probe and microwave-interferometer paths are midway between the upper and lower electrodes in the Miller hardware paper, while the dimensioned Figure-1 scan line is 15 mm above the lower electrode. | Diagnostic-plane provenance. Wise Figure 3 inherits the reference-cell diagnostic data through its reference 15; exact experiment-to-drawing alignment is retained as a source boundary. |
| M4 | Reported plasma power subtracts vacuum-measured coil/hardware resistive loss from matched-network input power; the authors state that plasma-induced changes in the hardware current distribution leave uncertainty in that correction. | The 180 W Wise value is a diagnostic-conditioned plasma-power quantity, not a generator-setpoint conversion law. |
| M5 | Across the source measurements, typical electron-density full width at half maximum was 7--9 cm and the density peaked on axis. | Independent qualitative/interval check on the Wise Figure-3 digitization and on future spatial reactor solutions; not a wafer-flux or depth validation. |

## Use decision

The dimensioned geometry is executable as an independent GEC spatial board.
No PDF bytes are packaged.  The exact official-PDF and 300-dpi render hashes,
the visually audited dimensions, and their measurement boundary are stored in
`data/experimental/wise_1996_gec_icp/gec_icp_geometry.json`.

The drawing does not measure the coil-side complex electric field, absorbed-
power distribution, wall recombination state, or species-resolved wafer flux.
Those quantities cannot be chosen from the Wise density curve or from a
feature depth.
