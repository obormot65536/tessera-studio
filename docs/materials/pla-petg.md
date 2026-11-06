# Built-in PLA and PETG: sources and evidence

This file records where every number in the two printed library materials comes from
(`packages/materials/src/tessera/materials/library.py`, issue #61, ADR-010 action item 3). The
maintainer approves each number in the pull request that adds or changes it (hard rule 6).

Axes: 1 along the raster, 2 across the road in the layer, 3 the build direction. Units are MPa
and t/mm³. Poisson ratios are the major ones, nu_ij = -eps_j / eps_i under stress along i.

Each entry is named after its source, because the same polymer measures very differently on
different printers and settings. Some examples from the sources below:
- PLA through-layer tensile strength: 5 MPa in [2], about 35 MPa in [1].
- PLA through-layer stiffness: 2000 MPa in [2], about 3800 MPa in [1].
- PETG in-layer tensile strength: 30 to 55 MPa across desktop studies.

## `pla-fff-gonabadi2020`

The base source is [2]. It was printed on an Ender-3 at 0.15 mm layers, 220 °C (its Sec 3.1
says 200 °C), 60 mm/s, with 100% infill, a 1.2 mm shell and a raft. Tension was tested to ASTM D638
type IV and shear to ASTM D5379 (Iosipescu), at 2 mm/min, with strain measured by DIC and five or
more specimens per condition. [2] did not test across the road in the layer (axis 2) or in
compression. Those values come from [1]: 0.2 mm layers, 220 °C, unidirectional, machined from
dense blocks, strain gauges, 2.5e-4 1/s, single tests.

| Value | | Source and location | Notes |
|---|---|---|---|
| e1 | 3500 | [2] p701, Fig 6 (0° on-edge) | |
| e2 | 4040 | [1] Table 6 | From [1], a different printer |
| e3 | 2000 | [2] p701, Fig 6 (upright) | |
| nu12 | 0.32 | [2] p701 (0° flat) | [1] gives 0.34 |
| nu13 | 0.35 | [2] p701 (0° on-edge) | Cross-check: the upright 0.2 is nu31, and 0.2 × 3500 / 2000 = 0.35 |
| nu23 | 0.37 | [1] Table 6 | From [1] |
| g12 | 1270 | [2] Table 7 (Iosipescu, 0° flat) | [1] gives 1500, from off-axis tension |
| g13 | 1210 | [2] Table 7 (Iosipescu, 0° on-edge) | |
| g23 | 950 | [2] Table 7 (Iosipescu, upright) | **Assumed plane.** The test shears the layer interface but does not say whether that is 13 or 23 |
| density | 1.22e-9 | [2] p698 | Filament, not a printed part |
| t1 | 55 | [2] p701, abstract | Peak; [2] gives no yield at 0°. [1]'s 0.2% yield is 54.84 |
| t2 | 45.16 | [1] Table 4 row 1 | 0.2% yield; from [1] |
| t3 | 5 | [2] p701, abstract | Peak, equal to failure (brittle, 0.45% strain). Process-sensitive (see above) |
| s12 | 30 | [2] Table 7 (0° flat) | "Yield"; [2] does not define the offset |
| s13 | 17.5 | [2] Table 7 (upright) | **Assumed:** the interlayer value, used because it is lower than the 23 MPa across-road value (0° on-edge) |
| s23 | 17.5 | [2] Table 7 (upright) | **Assumed plane**, as for g23 |
| c1 | 95.79 | [1] Table 4 row 9 | 0.2% yield; from [1] |
| c2 | 94.87 | [1] Table 4 row 11 | 0.2% yield; from [1] |
| c3 | — | not published | Taken equal to t3, and Tessera says so |

Strengths are first-failure values (`tessera.model.strength`): yield where the source gives one,
otherwise peak. The values from [1] were measured at a strain rate 20 times higher than its slow
set, which was about 38% lower. They are not conservative for slow or sustained loads.

## `petg-pellet-romeijn2022`

The only source is [3]. It was printed on a pellet-fed screw extruder on a robot arm: 5 mm roads,
2 mm layers, 220 °C, unheated bed, 38 °C enclosure. Specimens were heat-treated 10 min at 100 °C
and tested to ASTM D3039 at 2 mm/min with DIC. **This is not desktop filament printing.** Desktop
PETG is likely closer to isotropic in the layer: a desktop study quoted in [5] gives E 1620 and
1560 MPa at 0° and 90°. Treat e2 in particular as specific to this process.

| Value | | Source and location | Notes |
|---|---|---|---|
| e1 | 1898.7 | [3] Table 5 | CoV 3.0% |
| e2 | 1272.4 | [3] Table 5 | CoV 28% |
| e3 | 1894.8 | [3] Table 5 | CoV 3.1% |
| nu12 | 0.401 | [3] Table 6 | |
| nu13 | 0.462 | [3] Table 6, Eq 6 | **Derived:** measured nu31 = 0.461; nu13 = nu31 × E1 / E3 |
| nu23 | 0.373 | [3] Table 6 | CoV 24% |
| g12 | 748.8 | [3] Table 8 | From 45° off-axis tension (Eqs 2-4), not a shear test |
| g13 | 611.7 | [3] Table 8 | As g12 |
| g23 | 685.5 | [3] Table 8 | As g12 |
| density | 1.2742e-9 | [3] Sec 4.1.2 | Printed material, 2 samples |
| strength | none | | **No source gives PETG through-layer tension or interlayer shear.** The maintainer decided on #61 to ship stiffness only, so Tessera reports "not checked for strength" for this material |

## References

1. Song Y, Li Y, Song W, Yee K, Lee K-Y, Tagarielli VL. Measurements of the mechanical response
   of unidirectional 3D-printed PLA. Materials & Design 123 (2017) 154-164.
   doi:10.1016/j.matdes.2017.03.051. Accepted manuscript, Imperial College Spiral repository.
2. Gonabadi H, Yadav A, Bull SJ. The effect of processing parameters on the mechanical
   characteristics of PLA produced by a 3D FFF printer. Int J Adv Manuf Technol 111 (2020)
   695-709. doi:10.1007/s00170-020-06138-4.
3. Romeijn T, Behrens M, Paul G, Wei D. Instantaneous and long-term mechanical properties of
   Polyethylene Terephthalate Glycol (PETG) additively manufactured by pellet-based material
   extrusion. Additive Manufacturing (2022) 103145. doi:10.1016/j.addma.2022.103145.
4. Popa CF, Galatanu SV, Marsavina L. Prediction of the tensile strength of FDM specimens based
   on Tsai Hill criteria. Frattura ed Integrità Strutturale 73 (2025) 153-165.
   doi:10.3221/IGF-ESIS.73.11. Desktop PETG t1 55.21 and t2 48.89 MPa (Table 2). Not used: it
   has no through-layer or reliable shear data, and it comes from a different process than [3].
5. Algarni M, Ghazali S. Comparative study of the sensitivity of PLA, ABS, PEEK, and PETG's
   mechanical properties to FDM printing process parameters. Crystals 11 (2021) 995.
   doi:10.3390/cryst11080995. A review; used only for cross-checks.

## Reviewed and not used

| Paper | Why not |
|---|---|
| Mehir et al., Biopolymer Applications Journal 5(1) (2026) 40-52 (PLA) | Strain from crosshead travel, so moduli are less than half those above. No Poisson ratios or strengths |
| Hernández-Alvarez et al., J. Manuf. Mater. Process. 10 (2026) 64 (PETG) | Flat bars only, 2 specimens each. Moduli from crosshead travel. Tables 6 and 9 disagree |
| Radkovský, Bečvář, MM Science Journal 2026/03, doi:10.17973/MMSJ.2026_03_2025148 (PLA, PETG) | 15% infill, so the stresses describe a sparse structure |
| Almajdoub, Ahmed, Al-Mutawassit J. Basic Appl. Sci. 2(1) (2026) 1-18 (PLA, PETG, PA) | Cross-ply layups only. Printer and temperature not given. Moduli from crosshead travel |

## Still wanted

- PETG from desktop filament: through-layer tension (upright specimens), interlayer shear
  (Iosipescu, V-notched rail or short-beam), and all nine elastic constants from one study.
  With these, PETG would get strengths and a desktop entry.
- PLA: one study that covers the in-layer and through-layer directions, with the interlayer
  shear plane named (13 or 23), and through-layer compression (c3).
