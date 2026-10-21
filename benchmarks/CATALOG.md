# Validation catalog

This catalog plans the validation suite for Tessera Studio. Each case is a problem with a
closed-form or published reference answer. Together the cases show that a Tessera result can be
trusted before a maker relies on it. The target is a quality baseline comparable to the Ansys
Verification Manual and the NAFEMS standard benchmarks, built only from sources a GPL-3.0
repository may reuse. Most cases come from the classical literature (Timoshenko's books first,
then Roark, Lekhnitskii, Hertz and Johnson), plus published NAFEMS, Code_Aster and CalculiX test
problems used as cross-checks. The catalog is a plan, not the source of truth. A case becomes
binding only when its module exists in `benchmarks/cases/` and its reference file exists in
`benchmarks/references/`.

**Reference values are maintainer-approved.** Every number a test compares against lives in
`benchmarks/references/<name>.toml`, which only the maintainer edits (CLAUDE.md rule 3; a hook
blocks Claude). The values and tolerances below are *proposals*. Rule 6 applies: a reference
value, material number or tolerance becomes binding only after the maintainer approves it
explicitly in the PR.

## Legend

| Prefix | Domain | Roadmap phase |
|---|---|---|
| `TS` | Linear statics, isotropic | v0.1 |
| `TO` | Orthotropic and region-wise materials (FDM) | v0.5 |
| `TC` | Contact and assemblies: bonded/tied, frictional contact, bolts with preload, pins | v0.5 |
| `TM` | Modal (natural frequencies) | v1.x |
| `TN` | Geometric nonlinearity / large deflection | v2 |
| `TT` | Thermal conduction and thermal stress | v3 |
| `TF` | CFD | v4 |

(v1.0 adds the GUI and needs no new physics benchmarks. Its "3 reference parts" gate is a
separate acceptance test, not part of this catalog.)

**Status.** `passing`: the case exists in `benchmarks/cases/` and its reference file is
approved. `awaiting reference`: the case exists but its reference file is not approved yet, so
it skips. `planned`: not written yet.

**Verification of the reference value** (what this draft checked; see "How values were checked"):

- `verified`: value or formula seen in the cited public page (search-engine excerpt of that page).
- `computed`: number computed for this catalog from a formula whose source is cited. Recompute it
  in the case module (`test_reference_matches_formula`) before approval.
- `unverified`: written from memory or not found in an accessible source. The maintainer must
  check it against the book before it goes into a `.toml` file.
- `TBD`: value not stated on purpose. It comes from the cited source when the case is written.

**Tolerance.** These are proposals for the gate mesh. The v0.1 milestone gate is 2 %. Patch
tests should be exact up to round-off, because C3D10 reproduces linear strain fields exactly, and
use `1e-6` relative. Every case also runs a mesh-convergence sweep with at least three
refinements, as `fea-verifier` requires.

## How values were checked (October 2026)

Direct page fetches were blocked in the drafting environment. The proxy refused every
documentation host. Each `verified` mark below therefore rests on a web-search excerpt of the
cited URL, not on a full read of the page or book. No book was opened. Book chapter and section
numbers are marked `unverified` unless a public page quoted them. Before a reference goes into a
`.toml` file, the maintainer (or `fea-verifier` with page access) should confirm it in the
primary source.

## Summary

| ID | Phase | Problem | Status |
|---|---|---|---|
| TS-01 | v0.1 | End-loaded cantilever, tip deflection | passing |
| TS-02 | v0.1 | Strip with central hole in tension, Kt (Howland) | passing |
| TS-03 | v0.1 | Uniform tension patch test (isotropic) | planned |
| TS-04 | v0.1 | Bar hanging under its own weight | passing |
| TS-05 | v0.1 | Simply supported beam, uniform load | planned |
| TS-06 | v0.1 | Thick-walled cylinder, internal pressure (Lamé) | passing |
| TS-07 | v0.1 | Thick-walled sphere, internal pressure (Lamé) | awaiting reference |
| TS-08 | v0.1 | Circular plate, uniform pressure, clamped and simply supported | planned |
| TS-09 | v0.1 | Torsion of circular and square bars (Saint-Venant) | passing |
| TS-10 | v0.1 | Thin ring under diametral compression (curved beam) | planned |
| TS-11 | v0.1 | NAFEMS LE1 elliptic membrane | planned |
| TS-12 | v0.1 | NAFEMS LE10 thick plate under pressure (3D) | planned |
| TO-01 | v0.5 | Orthotropic patch tests: 3 axial + 3 shear | planned |
| TO-02 | v0.5 | Off-axis (rotated orientation) tension | planned |
| TO-03 | v0.5 | Orthotropic plate with hole (Lekhnitskii) | planned |
| TO-04 | v0.5 | Shear-flexible orthotropic cantilever | planned |
| TO-05 | v0.5 | Region-wise two-material beam (perimeter + infill) | planned |
| TC-01 | v0.5 | Tied interface patch test, non-matching meshes | planned |
| TC-02 | v0.5 | Tied cantilever split in two parts | planned |
| TC-03 | v0.5 | Frictionless contact patch test | planned |
| TC-04 | v0.5 | Hertz sphere on flat | planned |
| TC-05 | v0.5 | Shrink (interference) fit of two cylinders | planned |
| TC-06 | v0.5 | Frictional block: stick/slip threshold | planned |
| TC-07 | v0.5 | Single bolt with preload, then external axial load | planned |
| TC-08 | v0.5 | Bolted lap joint in shear, slip load | planned |
| TC-09 | v0.5 | Beam on two pins | planned |
| TM-01 | v1.x | Cantilever bending frequencies | planned |
| TM-02 | v1.x | Clamped circular plate frequency | planned |
| TM-03 | v1.x | NAFEMS FV52 solid square plate | planned |
| TN-01 | v2 | Elastica: cantilever with a large tip load | planned |
| TN-02 | v2 | Cantilever rolled up by an end moment | planned |
| TN-03 | v2 | Post-buckled column (elastica of a buckled bar) | planned |
| TT-01 | v3 | Free and fully constrained thermal expansion | planned |
| TT-02 | v3 | Bimetal strip (Timoshenko 1925) | planned |
| TT-03 | v3 | NAFEMS LE11 cylinder/taper/sphere, temperature load | planned |
| TT-04 | v3 | Steady radial conduction and thermal stress in a thick cylinder | planned |
| TF-01 | v4 | Laminar pipe flow (Hagen-Poiseuille) | planned |

Count: v0.1 12 (5 passing), v0.5 14, v1.x 3, v2 3, v3 4, v4 1. Total 37.

---

## v0.1: linear statics, isotropic (TS)

### TS-01 Cantilever, tip deflection (passing)

- **Problem:** steel prism 200 x 10 x 10 mm, root face clamped, 100 N uniform tip shear.
- **Quantity:** mean u_z of the tip-face nodes.
- **Reference:** Timoshenko beam, delta = PL^3/(3EI) + PL/(kappa G A), with Cowper's kappa.
  Value 1.526724 mm. Sources: Timoshenko & Gere, *Mechanics of Materials*; Timoshenko & Goodier,
  *Theory of Elasticity* 3rd ed. sec. 21; Cowper, J. Appl. Mech. 33 (1966). Full derivation in
  `benchmarks/cases/test_cantilever.py`.
- **Licence:** formula only, no reuse issue.
- **Equivalents:** none recorded.
- **Features:** face supports, face force (both exist).
- **Tolerance:** 2 % (approved, `cantilever_tip_deflection.toml`).

### TS-02 Strip with a central hole in tension (passing)

- **Problem:** W = 50, d = 10, t = 1 mm, 100 MPa gross tension, eighth model.
- **Quantity:** net-section Kt at the mid-thickness hole-edge node.
- **Reference:** Howland's strip solution as tabulated in Timoshenko & Goodier sec. 35:
  Kt_g = 3.14 at d/W = 0.2, so Kt_n = 2.512. Details in `benchmarks/cases/test_plate_hole.py`.
  The approved reference file marks it **PROVISIONAL**: the 3.14 is still to be confirmed in the
  book.
- **Licence:** formula and tabulated fact only.
- **Equivalents:** none recorded.
- **Features:** directional (roller) supports (`Support.dofs`, exists).
- **Tolerance:** 2 % (approved, `plate_hole_kt.toml`).

### TS-03 Uniform tension patch test

- **Problem:** box 10 x 10 x 10 mm with rollers on three orthogonal faces and uniform traction
  sigma on the opposite face. Run with a distorted, unstructured mesh.
- **Quantity:** every node and integration point gives sigma_xx = sigma, other stresses 0,
  eps_xx = sigma/E, eps_yy = eps_zz = -nu sigma/E.
- **Reference:** exact (linear field). Patch-test idea: MacNeal & Harder, "A proposed standard
  set of problems to test finite element accuracy", Finite Elements in Analysis and Design 1
  (1985) 3-20 (`unverified`: not looked up).
- **Licence:** no reuse issue.
- **Equivalents:** none recorded.
- **Features:** directional supports, face force (exist).
- **Tolerance:** 1e-6 relative.

### TS-04 Bar hanging under its own weight

- **Problem:** prism of length L hanging from its top face under gravity g. The case uses a
  roller (u_z = 0) top face: a clamp adds a -0.6 % Poisson bias unrelated to gravity (see the
  case docstring).
- **Quantity:** tip elongation delta = rho g L^2 / (2E); axial stress rho g (L - x).
- **Reference:** elementary (Timoshenko, *Strength of Materials* Part I, chapter on tension;
  section `unverified`). The case mainly exercises ADR-002 density units (t/mm^3, g = 9810
  mm/s^2): a unit slip gives an error of 10^3 or more.
- **Licence:** no reuse issue.
- **Equivalents:** none recorded.
- **Features:** gravity/body load, material density.
- **Tolerance:** 2 % (the clamped top end adds a local Poisson effect).

### TS-05 Simply supported beam, uniform load

- **Problem:** prismatic beam on two line supports (one pinned line, one roller line), with
  uniform pressure on the top face.
- **Quantity:** midspan deflection 5qL^4/(384EI), plus a Timoshenko shear term qL^2/(8 kappa G A).
- **Reference:** Timoshenko & Gere, *Mechanics of Materials* (deflection tables appendix;
  `unverified`); Roark's *Formulas for Stress and Strain*, beam table (`unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** pressure load, supports on edges/lines (or narrow faces), roller dofs.
- **Tolerance:** 2 %. Choose L/h >= 20 so that support indentation stays negligible.

### TS-06 Thick-walled cylinder under internal pressure (Lamé)

- **Problem:** long cylinder a = 10, b = 20 mm, internal pressure p, open ends. Model a slice
  with symmetry rollers.
- **Quantity:** sigma_theta and sigma_r at r = a and r = b, and the radial displacement u_r(a).
  sigma_theta(a) = p (a^2 + b^2)/(b^2 - a^2), i.e. 5p/3 for b/a = 2 (`computed`).
- **Reference:** Lamé. Timoshenko, *Strength of Materials* Part II, 3rd ed. (1956), art. 42,
  as cited by Ansys VM25 (`verified`: art. number seen in the VM25 public page excerpt). Also
  Timoshenko & Goodier (section `unverified`). Formula `verified` on the Wikipedia "Cylinder
  stress" page excerpt.
- **Licence:** formula only. Do not copy VM25 text or inputs.
- **Equivalents:** Ansys VM25 (thick-walled cylinder; VM38 is its elastic-plastic sibling,
  art. 70); Code_Aster SSLV04 (V3.04.004, hollow cylinder under internal pressure, Lamé
  reference, 3D, plane-stress and axisymmetric models).
- **Features:** pressure on curved faces, rollers on cut planes, cylindrical result
  components (post-processing).
- **Tolerance:** 2 % on sigma_theta(a), 1 % on u_r(a).

### TS-07 Thick-walled sphere under internal pressure (Lamé) (awaiting reference)

- **Problem:** hollow steel sphere a = 10, b = 20 mm, internal pressure p = 10 MPa, outer
  surface free. One octant with symmetry rollers (u_x = 0 on x = 0, u_y = 0 on y = 0,
  u_z = 0 on z = 0), which are exact for the point-symmetric field. Case:
  `benchmarks/cases/test_lame_sphere.py`.
- **Quantity:** sigma_t(a) = p (b^3 + 2a^3) / (2 (b^3 - a^3)) = 5p/7 = 7.142857 MPa
  (`computed`), measured as the mean of sigma_kk over the inner-surface nodes on the plane
  normal to k (k = x, y, z); u_r(a) = p a^4 / (E (b^3 - a^3)) [(1 - 2 nu) + (1 + nu) b^3 /
  (2 a^3)] = 3.809524e-4 mm (`computed`), mean over all inner-surface nodes. sigma_r(a) = -p is
  not gated (nodal extrapolation gives -1.5 % at the gate mesh).
- **Reference:** Timoshenko & Goodier, *Theory of Elasticity*, spherical vessel (section
  `unverified`); Ugural & Fenster, *Advanced Strength and Applied Elasticity*, thick-walled
  spheres (section `unverified`); Air Force Stress Analysis Manual "Thick pressure vessels" (seen
  as a search result only, formula not read). Formula `unverified`: it is the standard result,
  but no source was opened. The case checks it internally (sigma_r(a) = -p, sigma_r(b) = 0,
  sigma_r + 2 sigma_t constant, u_r consistent with Hooke's law).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** pressure on a doubly curved face, rollers on three cut planes (both exist).
- **Tolerance:** 2 % on sigma_t(a), 1 % on u_r(a), as TS-06 (proposed, `lame_sphere.toml`).
- **Measured** (CalculiX 2.21; the gate mesh is 1.25 mm, 8 elements through the
  wall, about 14 s):

  | h [mm] | nodes | sigma_t(a) | u_r(a) |
  |---|---|---|---|
  | 2.5 | 2 561 | +1.39 % | -0.084 % |
  | 1.667 | 7 066 | +0.66 % | -0.023 % |
  | 1.25 | 15 483 | +0.48 % | -0.009 % |
  | 1.0 | 27 785 | +0.23 % | -0.005 % |
  | 0.833 | 45 103 | +0.15 % | -0.003 % |

  The sweep in the case runs 2.5 / 1.667 / 1.25. The hoop mean converges from above. Single
  nodes scatter from -0.5 % to +1.6 % at 1.25 mm (up to +3.6 % at 2.5 mm), so the measure
  averages them.

### TS-08 Circular plate under uniform pressure

- **Problem:** disc of radius a, thickness h (a/h = 20), pressure q. Variant A: outer face
  clamped. Variant B: simply supported, i.e. a ring support (u_z = 0 on the bottom outer edge).
- **Quantity:** centre deflection. A: w0 = q a^4 / (64 D). B: w0 = (5 + nu)/(1 + nu) q a^4 /
  (64 D), with D = E h^3 / (12 (1 - nu^2)). Thin-plate formulas `verified` (structville.com
  derivation, tge.readthedocs.io page excerpts).
- **Reference:** Timoshenko & Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed. (1959),
  chapter 3, symmetrical bending of circular plates (section numbers `unverified`). The same
  chapter gives a correction for transverse shear in thicker plates. Its formula is `unverified`
  and must be taken from the book. At a/h = 20 the shear correction is still of the order of
  1 %, so the case should either use it or use a/h >= 40.
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** pressure, edge/ring support (variant B), thin-body meshing (at least 3
  quadratic elements through the thickness).
- **Tolerance:** 2 % (A), 3 % (B, the ring-support idealisation is ambiguous in 3D).

### TS-09 Torsion of a circular and a square bar

- **Problem:** (a) circular shaft of radius R. (b) square bar of side a. One end clamped, torque
  T on the other end (applied as a tangential traction or a moment on a face).
- **Quantity:** twist per unit length theta and maximum shear stress.
  (a) theta = T/(G J), J = pi R^4 / 2, tau_max = T R / J.
  (b) J = k1 a^4 with k1 = 0.1406, tau_max = T / (k2 a^3) with k2 = 0.2082 (`computed` from the
  Saint-Venant series; tabulated 0.141 / 0.208 `verified` in the Atlantis Press and CliffsNotes
  excerpts).
- **Reference:** Timoshenko & Goodier, *Theory of Elasticity*, ch. 10, torsion of a
  rectangular bar (sec. 109 per Karami's course notes; `unverified`); Timoshenko, *Strength of
  Materials* Part I (circular shaft).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** moment/torque load on a face (or a distributed tangential traction), kinematic
  coupling to a reference point.
- **Tolerance:** 2 % on theta (measured away from both ends); 3 % on tau_max (midside of the
  square, stress recovery at the surface).

### TS-10 Thin ring under diametral compression

- **Problem:** circular ring of mean radius R and rectangular section with R/h >= 10, pinched by
  two opposite forces P. Quarter model with symmetry rollers.
- **Quantity:** shortening of the loaded diameter, (pi/4 - 2/pi) P R^3 / (E I) =
  0.14878 P R^3 / (E I) (`computed`; 0.149 `verified` in a Scribd lab-report excerpt, which is a
  weak source).
- **Reference:** Timoshenko, *Strength of Materials* Part I, curved bars / thin ring, by
  Castigliano (section `unverified`); Roark's *Formulas for Stress and Strain*, ring table
  (`unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** rollers, point/small-face force, curved-beam bending.
- **Tolerance:** 3 % (curvature and shear corrections at R/h = 10 are about 1 %).

### TS-11 NAFEMS LE1, elliptic membrane

- **Problem:** quarter plate between an inner ellipse (semi-axes 2 x 1 m) and an outer ellipse
  (3.25 x 2.75 m), thickness 0.1 m. Outward pressure 10 MPa on the outer edge; E = 210 GPa,
  nu = 0.3. Symmetry on both axes. In Tessera units: mm, with the sizes scaled.
- **Quantity:** sigma_yy at D = (2, 0), the inner-ellipse end on the x axis. Target
  **92.7 MPa** (`verified`: COMSOL, Cast3M and Abaqus benchmark pages, which say the value is
  analytical). The outer-ellipse semi-axes are `unverified`: they match the wall widths 1.25 /
  1.75 shown in the COMSOL figure.
- **Reference:** NAFEMS, *The Standard NAFEMS Benchmarks*, TNSB Rev. 3, October 1990 (NAFEMS
  R0015), test LE1, as cited by Abaqus and Altair.
- **Licence:** the NAFEMS publication is copyrighted. Geometry and target values are facts and
  are widely reproduced; describe the problem in our own words and do not copy the figure.
- **Equivalents:** Abaqus Benchmarks LE1; Altair OS-V 0010; COMSOL "elliptic membrane".
- **Features:** pressure on a curved edge face, rollers, plane-stress-like thin solid (mesh at
  least 2 elements through the thickness), local refinement at D.
- **Tolerance:** 2 % (published fine-mesh quadratic results: Abaqus CPS8 92.81 MPa, +0.12 %).

### TS-12 NAFEMS LE10, thick plate under pressure (true 3D)

- **Problem:** the LE1 planform, thickness 0.6 m, 1 MPa pressure on the upper face.
  u_y = 0 on face DCD'C', u_x = 0 on face ABA'B', u_x = u_y = 0 on the outer face BCB'C', and
  u_z = 0 **only on the mid-plane line EE' of the outer face**.
- **Quantity:** sigma_yy at D, on the top surface at the inner-ellipse x-axis point. Target
  **-5.38 MPa** (compressive; `verified` in the Abaqus page, Altair OS-V 0060 and FeenoX
  excerpts. Abaqus and Altair quote it as a magnitude, "5.38").
- **Revision issue:** the support differs between transcriptions. The original TNSB Rev. 3
  constrains u_z only on the mid-plane line of the outer face, and -5.38 MPa belongs to that
  support. Some later guides constrain u_z over the whole outer face, and with that support the
  converged value is about -5.24 to -5.25 MPa, 2.4-2.6 % lower. This was reported via an ESRD
  NAFEMS benchmarks guide and GitHub issue andeplane/fem-lab#183. Both were seen only as search
  excerpts, so treat the -5.25 value as `unverified`. Tessera must model the **line**
  constraint (needs edge/line supports) and compare with -5.38. With a whole-face constraint
  the case would pass or fail for the wrong reason.
- **Reference:** NAFEMS TNSB Rev. 3 (1990), test LE10.
- **Licence:** as TS-11.
- **Equivalents:** Abaqus Benchmarks LE10; Altair OS-V 0060; COMSOL "thick plate"; FeenoX
  parametric LE10 (open source, about -5.380 MPa reported); Code_Aster runs reported by SimScale.
- **Features:** pressure, rollers, **supports on an edge (line)**, local refinement at D.
- **Tolerance:** 3 %. Published C3D10 results scatter (Abaqus fine C3D10 -5.77 MPa, +7 %;
  Altair Tet10 fine 5.03 MPa), so the gate needs a refined mesh near D and a convergence sweep.

---

## v0.5: orthotropic and region-wise materials (TO)

FDM parts are modelled as orthotropic with axes 1 = raster/extrusion direction, 2 = in-layer
transverse direction, 3 = build (Z) direction. The tests use made-up but valid engineering
constants. They do not use material data a user would trust, so rule 6 does not apply to them.
Example set: E1 = 3000, E2 = 2400, E3 = 1800, G12 = 1000, G13 = 700, G23 = 650 MPa,
nu12 = 0.35, nu13 = nu23 = 0.30 (check positive definiteness in the case).

### TO-01 Orthotropic patch tests (3 axial + 3 shear)

- **Problem:** unit cube with a distorted mesh. (a) Uniaxial traction along 1, 2 and 3 in turn,
  with rollers on the three opposite faces. (b) Pure shear in planes 12, 13 and 23, imposed as a
  linear displacement field on all boundary nodes (displacement patch test).
- **Quantity:** strains eps = S sigma with the orthotropic compliance S (axial: eps_ii =
  sigma/E_i, eps_jj = -nu_ij sigma/E_i; shear: gamma_ij = tau/G_ij). Also check that the
  reciprocity nu_ij/E_i = nu_ji/E_j holds in the material conversion.
- **Reference:** exact. Compliance form: Lekhnitskii, *Theory of Elasticity of an Anisotropic
  Body* (Holden-Day 1963 / Mir 1981), ch. 1; Jones, *Mechanics of Composite Materials*, 2nd ed.
  (1999), ch. 2 (both locations `unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** orthotropic material, material orientation, prescribed non-zero displacements.
- **Tolerance:** 1e-6 relative (v0.5 gate: "orthotropic patch tests exact").

### TO-02 Off-axis tension (rotated material orientation)

- **Problem:** as TO-01(a), but with axis 1 rotated by theta = 30 deg and 45 deg about Z from
  the load direction x.
- **Quantity:** eps_xx = sigma/E_x(theta) with 1/E_x = cos^4/E1 + (1/G12 - 2 nu12/E1) sin^2
  cos^2 + sin^4/E2, and the shear-extension coupling gamma_xy != 0 (sign checks the rotation
  sense).
- **Reference:** Jones, *Mechanics of Composite Materials*, ch. 2, transformed compliance;
  Lekhnitskii, *Anisotropic Plates* (Gordon & Breach 1968) (locations `unverified`; formula
  written from memory, so `unverified` until checked).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** per-region orientation as an angle or a local frame; rollers that leave the
  shear coupling free (one corner fixed, minimal rollers).
- **Tolerance:** 1e-6 relative.

### TO-03 Orthotropic plate with a circular hole (Lekhnitskii)

- **Problem:** thin wide plate (W/d >= 10, so finite-width effects are small) under uniaxial
  tension along axis 1, with a central hole. Eighth model as in TS-02.
- **Quantity:** Kt at the hole edge on the net section:
  Kt_inf = 1 + sqrt(2 (sqrt(E1/E2) - nu12) + E1/G12).
  It reduces to 3 for isotropy (`computed`). With E1 = 3000, E2 = 2000, G12 = 1000, nu12 = 0.35
  it gives 3.179 (`computed`).
- **Reference:** Lekhnitskii, *Anisotropic Plates* (1968), circular hole in an orthotropic plate
  (section `unverified`). The formula form was confirmed only as garbled search excerpts (NASA
  NTRS 19860005880; cdmHUB "Orthotropic Stress Concentration Factor" tool), so it is
  `unverified` until checked in the book. A finite-width correction (Tan, *Stress
  Concentrations in Laminated Composites*, 1994) is needed if W/d < 10 (`unverified`).
- **Licence:** formula only.
- **Equivalents:** TS-02 (isotropic limit). No NAFEMS/Code_Aster/Ansys case recorded.
- **Features:** orthotropic material + orientation, rollers, local refinement.
- **Tolerance:** 3 % (infinite-plate idealisation).

### TO-04 Shear-flexible orthotropic cantilever

- **Problem:** the TS-01 geometry with axis 1 along the span and a deliberately low G13. Use
  L/h = 5 so the shear term is about 20-30 % of the deflection.
- **Quantity:** tip deflection P L^3/(3 E1 I) + P L/(kappa G13 A).
- **Reference:** Timoshenko beam (as TS-01) with orthotropic moduli. The shear coefficient for an
  orthotropic section is `unverified`: choose kappa = 5/6 and keep its effect inside the
  tolerance, or take it from Cowper (1966).
- **Licence:** formula only.
- **Equivalents:** TS-01.
- **Features:** orthotropic material, orientation.
- **Tolerance:** 3 % (3D end effects grow at L/h = 5; set it from an independent run).

### TO-05 Region-wise two-material beam (perimeter shell + infill)

- **Problem:** rectangular cantilever split into an outer skin region (stiff, "perimeters") and
  a core region (soft, "infill"), each with its own orthotropic constants and orientation, and a
  conforming shared interface.
- **Quantity:** tip deflection with the transformed-section bending stiffness
  EI_eff = sum(E1_k I_k) (plus the shear term); also the axial stiffness of the same bar in
  tension, E A_eff = sum(E1_k A_k).
- **Reference:** composite (two-material) beams, Timoshenko, *Strength of Materials* Part I
  (section `unverified`); Timoshenko & Gere, *Mechanics of Materials*, composite beams
  (`unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** region-wise materials in one part (print-setup regions), per-region orientation.
- **Tolerance:** 2 % axial, 3 % bending.

---

## v0.5: contact and assemblies (TC)

### TC-01 Tied (bonded) interface patch test, non-matching meshes

- **Problem:** two blocks stacked along z and tied on their common face, meshed independently
  with different sizes so that the interface nodes do not match. Uniform tension along z, then
  along x (in-plane).
- **Quantity:** uniform stress sigma everywhere (max deviation), and a displacement field with
  no jump across the interface.
- **Reference:** exact (constant-stress patch). The tie must pass it to be consistent.
- **Licence:** no reuse issue.
- **Equivalents:** NAFEMS R0081 CGS-1 is a *contact* patch test (see TC-03). No tied case
  recorded.
- **Features:** bonded/tied connection (CalculiX `*TIE`), assemblies of several parts.
- **Tolerance:** 1e-4 relative for the z load (CalculiX ties of quadratic faces are not exactly
  patch-test-consistent, so ask `calculix-specialist` to confirm the achievable level); 1e-6
  for matching meshes.

### TC-02 Cantilever split into two tied parts

- **Problem:** the TS-01 beam cut at mid-span into two parts tied together, with non-matching
  meshes.
- **Quantity:** tip deflection.
- **Reference:** TS-01 value (1.526724 mm).
- **Licence:** none.
- **Equivalents:** TS-01.
- **Features:** tie, assemblies.
- **Tolerance:** 2 %, and within 0.5 % of the single-part TS-01 result at the same mesh size.

### TC-03 Frictionless contact patch test

- **Problem:** two blocks in frictionless contact with non-matching meshes, pressed together by
  a uniform pressure.
- **Quantity:** uniform contact pressure equal to the applied pressure; no stress oscillation.
- **Reference:** exact. NAFEMS R0081, Feng & Prinja, *Benchmark tests for finite element
  modelling of contact, gapping and sliding* (2001), CGS-1, is a contact patch test (`verified`
  that CGS-1 exists and is a patch test, via the NAFEMS R0081 page and a NAFEMS FENET 2004
  presentation excerpt; its values are `TBD`). Taylor & Papadopoulos (1991) on the contact patch
  test (`unverified`).
- **Licence:** NAFEMS R0081 is a paid, copyrighted publication. Use only the problem idea and
  the exact solution.
- **Equivalents:** NAFEMS CGS-1.
- **Features:** surface-to-surface contact (CalculiX `*CONTACT PAIR`, mortar or penalty).
- **Tolerance:** 1 % with mortar contact. Penalty node-to-surface may not pass exactly;
  `calculix-specialist` decides which formulation is gated.

### TC-04 Hertz sphere on a flat

- **Problem:** elastic sphere of radius R pressed on an elastic half-space (large block) by a
  force F. Axisymmetric quarter model with symmetry rollers.
- **Quantity:** contact radius a = (3 F R / (4 E*))^(1/3); peak pressure p0 = 3F/(2 pi a^2);
  approach delta = a^2 / R, with 1/E* = (1 - nu1^2)/E1 + (1 - nu2^2)/E2 (formulas `verified`
  in MIT OCW lecture notes and the Arizona optomechanics tutorial excerpts, consistent with
  Johnson ch. 4). Example, steel on steel, R = 10 mm, F = 100 N: a = 0.1866 mm,
  p0 = 1371 MPa, delta = 3.48 um (`computed`; uses `steel-generic` as in TS-01, so no new
  material data needs approval).
- **Reference:** Johnson, *Contact Mechanics*, CUP 1985, ch. 4 (Hertz theory); Timoshenko &
  Goodier, *Theory of Elasticity* 3rd ed., art. 140, pp. 409-413 (`verified` as the citation
  used in Ansys VM63, seen in its public page excerpt).
- **Licence:** formula only.
- **Equivalents:** Ansys VM63 (sphere on a rigid plane); Code_Aster SSNV104 (V6.04.104, two
  hemispheres, R = 50 mm, E = 20 000 MPa, nu = 0.3, 2 mm approach each side; reference centre
  sigma_zz = -2798.3 MPa as quoted by SimScale's validation page, `verified` via the search
  excerpt only); NAFEMS R0081 CGS-3 (Hertz cylinder on a rigid foundation, values `TBD`).
- **Features:** contact, strong local refinement under the contact (a/h_mesh >= 5), a body
  large enough to approximate a half-space (>= 20 a).
- **Tolerance:** 3 % on a, 5 % on p0, 3 % on F(delta).

### TC-05 Shrink (interference) fit of two cylinders

- **Problem:** hub on a shaft with radial interference delta, both long (plane strain-like
  slice).
- **Quantity:** contact pressure p from the Lamé compatibility (for equal materials:
  p = E delta (b^2 - c^2)(c^2 - a^2) / (2 c^3 (b^2 - a^2)), with radii a < c < b). Hoop stress
  at the hub bore.
- **Reference:** Timoshenko, *Strength of Materials* Part II, compound cylinders / shrink fits
  (`unverified`; the formula is written from memory, so `unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** contact with initial overlap (interference) or with an adjustment step,
  rollers.
- **Tolerance:** 2 %.

### TC-06 Frictional block: stick/slip threshold

- **Problem:** block pressed onto a fixed base by a normal force N, then a tangential force Q
  ramped up.
- **Quantity:** no relative slip for Q < mu N; slip begins at Q = mu N. The contact traction
  ratio tau/p equals mu on the slipping area.
- **Reference:** Coulomb friction (exact). Johnson, *Contact Mechanics* ch. 7 for partial slip
  (`unverified`).
- **Licence:** no reuse issue.
- **Equivalents:** Code_Aster SSNV128 (plate pressed on a rigid plane with friction,
  Lagrangian / penalty / continuous methods; `verified` via the V6.04.128 search excerpt).
  NAFEMS R0081 CGS-4 (sliding wedge on springs with friction, values `TBD`).
- **Features:** frictional contact, load stepping (nonlinear contact solve).
- **Tolerance:** 2 % on the slip load.

### TC-07 Single bolt with preload, then external axial load

- **Problem:** two plates clamped by one bolt (solid or beam-connector bolt) with preload Fp,
  then an external separating force FA on the plates, below the separation load.
- **Quantity:** (1) bolt force = Fp with no external load (equilibrium, exact); (2) bolt force
  increase = Phi FA with Phi = k_b/(k_b + k_c).
- **Reference:** (1) exact. (2) k_c from a pressure-cone model: VDI 2230 Part 1 (2015) or
  Shigley's *Mechanical Engineering Design*, ch. 8 (both `unverified`; copyrighted, cite the
  method only). The cone model is an engineering approximation, so (2) is an informative check,
  not a sharp benchmark.
- **Licence:** VDI 2230 and Shigley are copyrighted. Use only formulas and cite them.
- **Equivalents:** Code_Aster SSND103 (bolt as a discrete element with preload, DIS_BILI_ELAS;
  `verified` that it exists via the search excerpt).
- **Features:** bolt connection with preload, contact between the plates (and under the
  head/nut, or tied), load steps (preload, then service load).
- **Tolerance:** (1) 0.5 %; (2) 15 %.

### TC-08 Bolted lap joint in shear (friction grip)

- **Problem:** two overlapping plates joined by n bolts, each preloaded to Fp, friction
  coefficient mu at the faying surface, pulled apart in shear.
- **Quantity:** slip load F_slip = mu n_s n Fp (n_s = number of friction interfaces) and the load
  transferred before slip.
- **Reference:** Coulomb friction with the clamp force (exact for a uniform pressure
  distribution); slip-resistance form as in EN 1993-1-8 (`unverified`, copyrighted, cite only).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** bolts with preload, frictional contact, nonlinear load stepping. This is the
  "bolted-joint" case in the v0.5 gate.
- **Tolerance:** 3 % on F_slip.

### TC-09 Beam on two pins

- **Problem:** the TS-05 beam with pin holes at both ends on fixed pins (pin connectors free to
  rotate about the pin axis, one with axial sliding), with a midspan load.
- **Quantity:** midspan deflection equal to the simply supported value P L^3/(48 E I) (plus
  shear); zero end moments (pins transmit no moment).
- **Reference:** Timoshenko & Gere, *Mechanics of Materials* (`unverified` location).
- **Licence:** formula only.
- **Equivalents:** TS-05.
- **Features:** pin (revolute) connection, assemblies.
- **Tolerance:** 3 % (local hole flexibility adds compliance; keep L/h >= 20 and the holes
  small).

---

## v1.x: modal (TM)

### TM-01 Cantilever bending frequencies

- **Problem:** the TS-01 prism (L/h = 20), clamped root, density rho.
- **Quantity:** first three bending frequencies f_n = (beta_n L)^2 / (2 pi L^2)
  sqrt(E I / (rho A)), with beta_n L = 1.875104, 4.694091, 7.854757 (`computed` as roots of
  cos x cosh x = -1; `verified` in the arXiv 2303.18126 and SolidWorks tutorial excerpts). Also
  the first torsional and axial modes.
- **Reference:** Timoshenko, Young & Weaver, *Vibration Problems in Engineering*, 4th ed.
  (1974), transverse vibration of prismatic bars (section `unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** modal analysis (`*FREQUENCY`), density, mode identification.
- **Tolerance:** 1 % on mode 1, 2 % on modes 2-3 (shear and rotary inertia lower the higher
  modes; use Timoshenko-beam values if they move by more than 0.5 %).

### TM-02 Clamped circular plate, fundamental frequency

- **Problem:** the TS-08A disc, a/h = 40.
- **Quantity:** omega_1 = (lambda^2 / a^2) sqrt(D / (rho h)) with lambda^2 = 10.2158
  (`computed` from J0(l) I1(l) + I0(l) J1(l) = 0, lambda = 3.19622).
- **Reference:** Leissa, *Vibration of Plates*, NASA SP-160 (1969), clamped circular plate
  (public domain US government work; `verified` that the report exists via the DTIC copy; the
  table value is `unverified`); Timoshenko, Young & Weaver (`unverified`).
- **Licence:** NASA SP-160 is a US government publication, free to cite and reuse.
- **Equivalents:** none recorded.
- **Features:** modal analysis, thin-plate meshing.
- **Tolerance:** 2 %.

### TM-03 NAFEMS FV52, simply supported solid square plate

- **Problem:** 10 x 10 x 1 m solid plate, u_z = 0 on the bottom edges, E = 200 GPa,
  nu = 0.3, rho = 8000 kg/m^3.
- **Quantity:** frequencies of modes 4, 5/6, 7: **45.897, 109.44, 167.89 Hz** (modes 1-3 are
  in-plane rigid-body modes). `verified` in the Altair OS-V 0455 (2017/2023) and Ansys Benchmark
  C5 public page excerpts. Some Altair 2021/2022 pages list 265.81 Hz for mode 7, so that mode
  is disputed: gate only modes 4-6 until the original is checked. The Ansys case is labelled
  "FV52 (modified)" and differs in its support line (bottom edges, not the mid-plane).
- **Reference:** NAFEMS TNSB Rev. 3 (1990), FV52 (also NAFEMS R0015).
- **Licence:** as TS-11.
- **Equivalents:** Ansys Benchmark C5 (FV52 modified); Altair OS-V 0455.
- **Features:** modal, edge supports, density.
- **Tolerance:** 1 % (modes 4-6).

---

## v2: geometric nonlinearity (TN)

### TN-01 Elastica: cantilever under a large tip load

- **Problem:** slender cantilever (L/h >= 50), tip force P that stays vertical (dead load),
  alpha = P L^2 / (E I) = 1, 2, 5, 10.
- **Quantity:** tip deflections over L (vertical / horizontal):
  alpha = 1: 0.30172 / 0.05643; alpha = 2: 0.49346 / 0.16064; alpha = 5: 0.71379 / 0.38763;
  alpha = 10: 0.81061 / 0.55500 (`computed` by shooting on the elastica ODE
  theta'' = -alpha cos theta, RK4 + bisection).
- **Reference:** Bisshopp & Drucker, "Large deflection of cantilever beams", Quart. Appl. Math.
  3 (1945) 272-276 (citation `verified` via the AMS QAM 1973 follow-up page excerpt; the
  tabulated values themselves are `unverified`, so compare the computed values with the paper's
  table before approval).
- **Licence:** formula only. The 1945 paper is copyrighted (AMS); use only its numbers.
- **Equivalents:** none recorded for a straight tip-loaded elastica (Ansys VM136 is the buckled
  bar, see TN-03).
- **Features:** geometrically nonlinear static analysis (`NLGEOM`), load stepping.
- **Tolerance:** 1 %.

### TN-02 Cantilever rolled up by an end moment

- **Problem:** straight cantilever with an end moment M (as a force couple or a coupled
  moment) increased until M L / (E I) = 2 pi.
- **Quantity:** tip position x = (E I / M) sin(theta), y = (E I / M)(1 - cos(theta)) with
  theta = M L / (E I). At 2 pi the beam closes into a full circle.
- **Reference:** exact (constant curvature M/(E I)); elementary beam theory, Timoshenko & Gere
  (`unverified` location).
- **Licence:** none.
- **Equivalents:** Ansys VM26 is the shell-plate version under an end moment (`verified` that
  it exists via the public page excerpt).
- **Features:** NLGEOM, follower moment (rotating couple), kinematic coupling.
- **Tolerance:** 1 % at theta = pi/2 and pi; 2 % at 2 pi.

### TN-03 Post-buckled column (elastica of a buckled bar)

- **Problem:** column fixed at the base and free at the top, axial load P/P_cr from 1.015 to
  about 1.9, triggered by a small lateral perturbation.
- **Quantity:** tip lateral and axial displacement and the tip rotation versus P/P_cr.
- **Reference:** Timoshenko & Gere, *Theory of Elastic Stability*, 2nd ed., p. 78 (`verified` as
  the citation used by Ansys VM136, public page excerpt). Values `TBD` from the book or an
  elliptic-integral computation.
- **Licence:** formula only. Do not copy VM136 tables.
- **Equivalents:** Ansys VM136 "Large Deflection of a Buckled Bar (the Elastica)"; Code_Aster
  SSNL502 (cantilever buckling in large displacements, `verified` via the search excerpt).
- **Features:** NLGEOM, arc-length or displacement control near P_cr.
- **Tolerance:** 2 %.

---

## v3: thermal and thermal stress (TT)

### TT-01 Free and fully constrained thermal expansion

- **Problem:** block heated uniformly by Delta T. (a) Free (minimal statically determinate
  supports). (b) Clamped on all faces normal to x, y and z.
- **Quantity:** (a) eps = alpha Delta T in all directions and zero stress. (b) sigma =
  -E alpha Delta T / (1 - 2 nu) hydrostatic (all directions constrained); a bar constrained only
  axially gives sigma = -E alpha Delta T. With orthotropic alpha_i (FDM), (a) gives
  eps_i = alpha_i Delta T.
- **Reference:** exact. Timoshenko & Goodier, *Theory of Elasticity*, thermal stress chapter
  (`unverified` location).
- **Licence:** none.
- **Equivalents:** none recorded.
- **Features:** temperature field input, thermal expansion (isotropic and orthotropic alpha).
- **Tolerance:** 1e-6 relative.

### TT-02 Bimetal strip (Timoshenko 1925)

- **Problem:** two bonded layers (thicknesses t1, t2, moduli E1, E2, expansion alpha1, alpha2),
  heated by Delta T, free.
- **Quantity:** curvature. Timoshenko's formula (as usually quoted):
  1/rho = 6 (alpha2 - alpha1) Delta T (1 + m)^2 / (h [3 (1 + m)^2 + (1 + m n)(m^2 + 1/(m n))]),
  with m = t1/t2, n = E1/E2, h = t1 + t2. For m = n = 1 it gives 3 Delta alpha Delta T / (2h).
  Formula `unverified`: written from memory, and a patent excerpt reproduces a formula of this
  kind but it was not checked against the paper.
- **Reference:** S. Timoshenko, "Analysis of bi-metal thermostats", J. Opt. Soc. Am. 11 (1925)
  233-255 (citation `verified` via patent excerpts that cite it).
- **Licence:** formula only. The 1925 paper is copyrighted (Optica).
- **Equivalents:** none recorded. This is the thermal analogue of a two-region print (TO-05).
- **Features:** thermal expansion, region-wise materials or tied parts.
- **Tolerance:** 2 % (L/h >= 20; end effects decay within about one thickness).

### TT-03 NAFEMS LE11, solid cylinder/taper/sphere with a temperature load

- **Problem:** axisymmetric solid (cylinder, taper, sphere) with temperature
  Delta T = sqrt(x^2 + y^2) + z, E = 210 GPa, nu = 0.3, alpha = 2.3e-4 /deg C.
- **Quantity:** sigma_zz at point A: target **-105 MPa** (`verified` in the Abaqus and Altair
  OS-V 0070 excerpts; Cast3M quotes the magnitude 105 MPa).
- **Reference:** NAFEMS TNSB Rev. 3 (1990), LE11.
- **Licence:** as TS-11.
- **Equivalents:** Abaqus Benchmarks LE11; Altair OS-V 0070; Cast3M test case.
- **Features:** spatially varying temperature field (analytical function on the nodes),
  thermal expansion.
- **Tolerance:** 5 % (fine quadratic meshes published at -1.7 % (C3D20) to -6 % (Tet10); tighten
  after Tessera's own convergence study).

### TT-04 Steady radial conduction and thermal stress in a thick cylinder

- **Problem:** long hollow cylinder with T(a) = T_i and T(b) = T_o, steady conduction, then
  thermoelastic stress.
- **Quantity:** T(r) = T_i + (T_o - T_i) ln(r/a)/ln(b/a) (exact); hoop stress at r = a and r = b
  from the closed-form thermal-stress solution (`TBD` from the source).
- **Reference:** Timoshenko & Goodier, *Theory of Elasticity*, thermal stresses in a long
  circular cylinder (section `unverified`).
- **Licence:** formula only.
- **Equivalents:** none recorded.
- **Features:** heat conduction solve, temperature boundary conditions, coupled
  (sequential) thermal-stress step.
- **Tolerance:** 1 % on T, 2 % on stresses.

---

## v4: CFD (TF)

### TF-01 Laminar pipe flow (Hagen-Poiseuille)

- **Problem:** fully developed laminar flow in a circular pipe (Re < 100) with a pressure drop
  Delta p over length L.
- **Quantity:** u_max = Delta p R^2 / (4 mu L) and the parabolic profile; flow rate
  Q = pi R^4 Delta p / (8 mu L).
- **Reference:** exact. Any fluid-mechanics text, e.g. White, *Fluid Mechanics* (`unverified`
  location).
- **Licence:** none.
- **Equivalents:** none recorded.
- **Features:** incompressible laminar flow solver, inlet/outlet pressure, no-slip walls.
- **Tolerance:** 1 %.

---

## Licence notes

- **Classical books** (Timoshenko, Timoshenko & Goodier, Timoshenko & Woinowsky-Krieger,
  Timoshenko/Young/Weaver, Roark, Lekhnitskii, Johnson) are copyrighted. Formulas, tabulated
  coefficients and reference values are facts and may be used with a citation. Do not copy text,
  figures or tables verbatim into the repository.
- **NAFEMS** benchmark publications (TNSB Rev. 3 / R0015; R0081 for contact) are copyrighted
  and sold. Geometry and target values are facts and widely reproduced. Write our own problem
  description, cite the publication, and do not copy figures or text.
- **Code_Aster** documentation (including the V-series validation manual) carries "Document
  diffusé sous licence GNU FDL" on the versions found (v11-v15 era pages and mirrors, some with
  "Copyright 20xx EDF R&D"). The notice on the **current** (v16/v17) documentation was not
  verified, because code-aster.org could not be reached. The GFDL is not compatible with
  GPL-3.0 for copied text, so cite Code_Aster cases and use their numbers, but do not paste
  their prose or figures into GPL files. Code_Aster `.comm` inputs are GPL (the code is GPL
  since 2001); check the exact version before copying any.
- **CalculiX** test suite (`ccx_<version>.test.tar.bz2` from dhondt.de). Packagers disagree on
  the licence: Debian and openSUSE say **GPL-2.0-only**, MSYS2 says GPL-2.0-or-later. A
  GPL-2.0-only file cannot be included in a GPL-3.0 project. Until the source headers or COPYING
  file are checked, **do not copy CalculiX `.inp` test decks into the repo**: reproduce the
  problems independently and cite them.
- **Ansys Verification Manual:** no text or input files used. VM numbers are named only where a
  public ansyshelp page appeared in search results (VM25, VM26, VM38, VM63, VM136, Benchmark
  C5).
- **NASA SP-160** (Leissa) is a US government work and free to reuse.

## Sources

Books and papers (cited, not opened in this draft; locations marked `unverified` above need a
check):

- S. P. Timoshenko, *Strength of Materials*, Parts I and II, 3rd ed., Van Nostrand 1955-56
  (Part II art. 42 and art. 70 as cited by Ansys VM25 / VM38).
- S. P. Timoshenko & J. N. Goodier, *Theory of Elasticity*, 3rd ed., McGraw-Hill 1970 (sec. 21,
  35; art. 140 pp. 409-413 as cited by Ansys VM63; torsion ch. 10).
- S. P. Timoshenko & S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed., McGraw-Hill
  1959, ch. 3.
- S. P. Timoshenko, D. H. Young & W. Weaver, *Vibration Problems in Engineering*, 4th ed.,
  Wiley 1974.
- S. P. Timoshenko & J. M. Gere, *Theory of Elastic Stability*, 2nd ed., McGraw-Hill 1961, p. 78;
  *Mechanics of Materials*.
- S. P. Timoshenko, "Analysis of bi-metal thermostats", J. Opt. Soc. Am. 11 (1925) 233-255.
- W. C. Young, R. G. Budynas, A. M. Sadegh, *Roark's Formulas for Stress and Strain*.
- S. G. Lekhnitskii, *Anisotropic Plates*, Gordon & Breach 1968; *Theory of Elasticity of an
  Anisotropic Body*.
- R. M. Jones, *Mechanics of Composite Materials*, 2nd ed., Taylor & Francis 1999.
- K. L. Johnson, *Contact Mechanics*, Cambridge University Press 1985.
- K. E. Bisshopp & D. C. Drucker, Quart. Appl. Math. 3 (1945) 272-276.
- G. R. Cowper, J. Appl. Mech. 33 (1966) 335-340.
- A. W. Leissa, *Vibration of Plates*, NASA SP-160 (1969).
- NAFEMS, *The Standard NAFEMS Benchmarks*, TNSB Rev. 3, 1990 (R0015).
- Q. Feng & N. K. Prinja, *Benchmark Tests for Finite Element Modelling of Contact, Gapping and
  Sliding*, NAFEMS R0081, 2001.

Web pages seen (as search-engine excerpts, October 2026):

- NAFEMS LE1: https://doc.comsol.com/6.4/doc/com.comsol.help.models.sme.elliptic_membrane/elliptic_membrane.html ;
  https://www-cast3m.cea.fr/html/CasTestsCastem/node8.html ;
  https://abaqus.uclouvain.be/English/SIMACAEBMKRefMap/simabmk-c-le1.htm ;
  https://help.altair.com/hwsolvers/os/topics/solvers/os/nafems_test_problem_le1_r.htm
- NAFEMS LE10: https://abaqus-docs.mit.edu/2017/English/SIMACAEBMKRefMap/simabmk-c-le10.htm ;
  https://help.altair.com/hwsolvers/os/topics/solvers/os/nafems_test_problem_le10_r.htm ;
  https://doc.comsol.com/6.4/doc/com.comsol.help.models.sme.thick_plate/thick_plate.html ;
  https://www.seamplex.com/feenox/tests/nafems/le10/ ;
  https://www.seamplex.com/fino/cases/012-nafems-le10 ;
  https://www.esrd.com/wp-content/uploads/dlm_uploads/Benchmarks-Guide-Standard-NAFEMS-Benchmarks-Linear-Elastic-Tests.pdf ;
  https://github.com/andeplane/fem-lab/issues/183 ;
  https://www.simscale.com/docs/validation-cases/thick-plate-under-pressure/
- NAFEMS LE11: https://docs.software.vt.edu/abaqusv2024/English/SIMACAEBMKRefMap/simabmk-c-le11.htm ;
  https://help.altair.com/hwsolvers/os/topics/solvers/os/nafems_test_problem_le11_r.htm ;
  https://www-cast3m.cea.fr/html/CasTestsCastem/node9.html
- NAFEMS FV52: https://2017.help.altair.com/2017/hwsolvers/simply-supported_solid_square_plate.htm ;
  https://help.altair.com/hwsolvers/os/topics/solvers/os/simply_supported_soid_square_plate_fv52_r.htm ;
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_vm/Hlp_V_VMC5.html
- NAFEMS R0081 (contact): https://www.nafems.org/publications/resource_center/r0081/ ;
  https://nafems.org/downloads/FENet_Meetings/Glasgow_UK_Oct_2004/FENET_Glasgow_Oct2004_DLE_Becker.pdf
- Ansys VM (named only): https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_vm/Hlp_V_VM25.html ;
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_vm/Hlp_V_VM38.html ;
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_vm/Hlp_V_VM63.html ;
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_vm/Hlp_V_VM136.html ;
  https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ans_vm/Hlp_V_VM26.html
- Code_Aster: SSLV04 https://ericca.uqtr.ca/fr14.4/man_v/v3/v3.04.004.pdf ;
  SSNV104 https://ericca.uqtr.ca/fr13.6/man_v/v6/v6.04.104.pdf and
  https://simscale.com/docs/validation-cases/hertzian-contact-between-two-spheres ;
  SSNV128 https://ericca.uqtr.ca/fr12.4/man_v/v6/v6.04.128.pdf ;
  SSND103 https://ericca.uqtr.ca/fr13.6/man_v/v6/v6.08.103.pdf ;
  SSNL502 https://ericca.uqtr.ca/fr11.7/man_v/v6/v6.02.502.pdf ;
  licence notices https://ericca.uqtr.ca/fr14.4/man_u/u0/u0.00.01.pdf ,
  https://code-aster.org/V2/doc/v13/fr/man_u/u0/u0.00.01.pdf
- CalculiX licence: https://metadata.ftp-master.debian.org/changelogs/main/c/calculix-ccx/stable_copyright ;
  https://src.opensuse.org/pool/ccx/src/branch/devel/ccx.spec ;
  https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-calculix-ccx
- Formulas: https://en.wikipedia.org/wiki/Cylinder_stress (Lamé);
  https://structville.com/2020/11/elastic-analysis-of-circular-plates.html and
  https://tge.readthedocs.io/en/latest/content/mechanical/bending_forces.html (circular plates);
  https://www.atlantis-press.com/article/126000183.pdf and
  https://www.cliffsnotes.com/study-notes/598479 (torsion coefficients);
  https://ocw.iti.hr/courses/mechanical-engineering/2-76-multi-scale-system-design-fall-2004/readings/lec12cont_read.pdf and
  https://wp.optics.arizona.edu/optomech/wp-content/uploads/sites/53/2016/12/Tutorial_Taylor_Brian.pdf (Hertz);
  https://arxiv.org/pdf/2303.18126 and
  https://help.solidworks.com/2025/English/simconntutorial/Dynamic_Analysis_of_a_Cantilever_Beam.htm (cantilever roots);
  https://www.ams.org/qam/1973-30-04/S0033-569X-1973-99716-X (Bisshopp & Drucker citation);
  https://apps.dtic.mil/sti/tr/pdf/ADA307623.pdf (NASA SP-160);
  https://ntrs.nasa.gov/api/citations/19860005880/downloads/19860005880.pdf and
  https://cdmhub.org/resources/orthoscf (Lekhnitskii hole);
  https://www.scribd.com/document/475261208/Lab-Report-9-docx (ring coefficient, weak source)
