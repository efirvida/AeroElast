# Removal / migration manifest: NuMAD + windio2 meshers and windIO 1.x (gmsh-official-mesher S9)

Derived read-only at HEAD b5d5c49 (explorer mwnocwpc-b-eyl5, grep + reads, no tests run; "?" = judged from grep
hits only). Executed by tasks T10-T17 of odd/tasks/gmsh-official-mesher.md. Every DELETE is recorded in that file's
Log with the commit that removes it. Physics guards are MIGRATED, never dropped; anchors are never re-banded
(AGENTS rules 1-2): a reference number that moves on v2/gmsh is a finding to log.

## A. src

### A1. NuMAD meshing cluster (src/aeroelast/models/blade/numad/mesh_gen/)
After b5d5c49 the only importer outside the cluster is core/mesh/generators.py, inside dead methods (1191, 1248, 1284).

| file | remaining src importers | disposition |
|---|---|---|
| mesh_gen.py (`__init__.py:1` imports get_shell_mesh) | generators.py:1191/1248/1284 (dead) | DELETE; fix `__init__.py` |
| mesh2d, shell_region, local_refine, segment2d, surface, boundary2d, spatial_grid_list2d/3d, element_utils | each other, mesh_gen | DELETE |
| mesh_tools.py | windio2_basis.py:12 (get_direction_cosines) | MOVE get_direction_cosines into windio2_basis, DELETE |
| windio2_mesh.py | core/mesh/windio2.py:13 (generate_mesh), gmsh_skin.py:1974 (Pchip, _centroid), :2040 (RawMesh) | MOVE Pchip, _centroid, RawMesh (and map_mesh_materials only if still used) to a gmsh-side module; DELETE the rest |
| windio2_collapse.py, collapse_opt.py, windio2_smooth.py | windio2_mesh path only | DELETE |
| windio2_basis.py | gmsh_skin:1973/2039 | KEEP (shell_bases, map_stack, physical_fiber_angle) |
| gmsh_skin.py | gmsh_blade.py:213 | KEEP; fix docstrings naming windio2_mesh (1949, 1955, 1964, 2021, 2036) |

span_domain already lives in objects/windio2_laminate.py:58. tests/support/v2_blade.py:96 hashes windio2_*.py and
collapse_opt.py: switch the digest to gmsh_skin + windio2_basis (invalidates cached npz).

### A2. Dead code
- core/mesh/generators.py BladeMesh: _insert_web_snap_stations (1177), _insert_web_end_stations (1235),
  _insert_row_stations (1274), _limit_bay_length (1302), _refine_high_gradient_sections (1321),
  _deduplicate_and_create_mesh (1366) -> DELETE (~1177-1486). Callers deleted with them: tests/test_blade_bay_length.py,
  test_blade_span_decoupling.py:245, software/mesh/test_windio2_public_mesh.py (21, 155-157), tools/run_s1_sectional.py:122
  (migrate that call). Remove excel_file / airfoil_dir from the BladeMesh signature.
- models/blade/model.py: generate_blocks_sections (164, unreachable body to ~600, splineZi uses) DELETE;
  BladeModel.generage_mesh (891) and BladeModel n_samples: check callers, delete.
- excel_file plumbing: runner.py:716-868 and :1326, cli/run_bem_fsi.py:299-338 and :546, cli/run_aero_fsi.py:132,
  cli/aeroelast.py:59-96 (--excel-file), standalone.py:125 -> DELETE.

### A3. windIO 1.x / excel reading
| item | disposition |
|---|---|
| io/excel_to_blade.py, objects/blade.py read_excel (5, 65-66, 103-119) | DELETE; drop openpyxl if unused |
| io/xml_to_airfoil, airfoil_to_xml | KEEP (yaml_to_blade imports airfoil_to_xml) |
| models/blade/windio.py | KEEP normalize_windio (2.x legacy view read by yaml_to_blade:81, aerodynamics.py:660) |
| windio.py v1 branch (startswith("1.") early return; is_windio2 default "1.0") | ValueError naming the v2 file, after every v1 input is gone (check tests/exp_aero_damp_off/IEA-15-240-RWT.yaml) |
| yaml_to_blade.py v1 else-branches (50, 137, 190; hub_outer_shape_bem fallback 87-88; _add_stations v1 axis 213-222) | DELETE after v1 is rejected |

### A4. core/mesh/windio2.py geometry_mesh
Imports windio2_mesh.generate_mesh. No live src consumer (BladeMesh geometry_only raises, T7). Tests:
test_windio2_public_mesh, test_windio2_precision.py:106-131. DELETE geometry_mesh with the windio2 mesher; keep
WindIO2Metadata / ElementMetadata only if gmsh_blade uses them (grep). Check cli/mesh.py for geometry_only.

## B. tests

### B1. DELETE
- software/mesh: test_blade_local_refine, test_blade_mesh, test_blade_mesh_quality, test_collapse_opt,
  test_windio2_collapse, test_windio2_smooth, test_windio2_mesh, test_windio2_web_mesh (move its
  physical_fiber_angle frame check to the gmsh path if it guards the frame), test_windio2_precision,
  test_windio2_public_mesh.
- tests/: test_blade_mesher_minor, test_blade_span_decoupling, test_blade_bay_length, test_blade_web_position,
  test_blade_span_interpolation, test_blade_chord_grid, test_blade_le_band_width, test_blade_te_band_width,
  test_blade_side_assignment, test_blade_ply_order, test_blade_skin_reference (not ..._default), test_blade_mass_mapping,
  test_sharpy_bem_numad_rotor_contract.
- v1 twins (v2 twin in tests/validation/blade/v2/): test_iea15mw_s1_sectional, test_iea15mw_s2_static_prescribed
  (incl. UTD rows), test_iea15mw_s4_rotating_modal, test_iea15mw_v03_static_gravity,
  test_iea15mw_v05_structural_properties, validation/blade/test_iea15mw_s7_torsion,
  validation/blade/test_iea15mw_v02_natural_frequencies. Check first what they import from tools (run_s1/s2/s4).
- software/test_audit_blade_input (with tools/audit_blade_input.py).
- software/mesh/test_yaml_materials_optional.py:57 (IEA34 v1): MIGRATE if it guards optional materials, else DELETE (?).

### B2. MIGRATE (target: upstream v2 yaml, gmsh chord_divisions 20 span_aspect 1, skin_reference mid)
| test | guard | note |
|---|---|---|
| validation/blade/test_blade_rated_twist | rated twist sign/magnitude | v1 references: re-measure |
| test_blade_twist_mechanism, test_blade_twist_anchor_beam | twist mechanism, beam anchor | re-measure |
| test_blade_section_frame_tear | section frame tear | re-measure |
| test_blade_section_torsion_stiffness | section GJ | re-measure |
| test_sectional_tangent_frame | tangent frame | calls BladeMesh directly |
| test_blade_deloading_vs_reference | deloading | re-measure |
| validation/rotor/test_bem_fsi_deformed_geometry | deformed geometry | |
| test_multicell_shear_flow (619, 1256) | multicell shear flow | |
| test_blade_ccx_parity, test_ccx_fsi_parity_parked_v50 | CalculiX parity | CalculiX rerun |
| validation/parity/test_blade_iea15mw_mesh_convergence | mesh convergence | via blade_val.YAML |
| validation/blade/test_blade_iea15mw_validation | hosts _to_rust_mesh, SPAN_DIRECTION (imported by test_multicell_shear_flow:1248) | split helpers out first |
| software/mesh/test_turbine_mesh (_YAML line 12) | tower/nacelle/hub | re-point to v2 |
| software/mesh/test_windio2_regions:77, test_windio2_frames:225 | v1-vs-v2 comparisons | DELETE the v1 halves (?) |
| test_aeroelast_cli_mesh_export:94/108/130, software/cli/test_cli_mesh:12, test_rotor_mesh_export:229, test_rotor_aero_surface:159 | CLI / rotor mesh | re-point |
| test_vlm_bem_validation, test_aero_backend_contract_coherence, test_iea15mw_v01_rotor_performance (V-01), validation/bem/test_bem_polars, test_bem_openfast_parity, test_bem_prebend_polar_interp, test_force_projection_load_frame | BEM/VLM via load_blade_aero | iea_blade_yaml fixture -> v2; V-01 / BEM parity may move: highest re-measure risk |

### B3. KEEP
test_blade_mesh_gmsh_adapter, test_gmsh_*, test_blade_no_webs, test_blade_skin_reference_default,
validation/blade/v2/*, software/test_migrate_blademesh_configs.py, test_windio2_iea22, test_windio2_iea34,
test_windio2_material_mapping, v2 halves of test_windio2_regions/frames, support/test_v2_blade,
support/test_v2_sectional, test_ch6_dlc_prepared_cases.

### B4. Support helpers
- tests/conftest.py:169-251: iea_blade_yaml -> v2 upstream; iea_numad_blade re-point or delete; iea_blade_xlsx DELETE.
- tests/support/v2_blade.py: MESHERS ("gmsh",), default "gmsh" (75, 91, 106, 154, 256); drop the windio2 digest
  glob and _generate_raw windio2 branch (157-163); drop V2_BLADE_MESHER and mesher= in callers.
- tests/support/v2_fidelity.py:22 docstring; tests/support/v2_sectional.py check for windio2 calls (?).

## C. tools
Imported by tests (MIGRATE, never delete): run_s1_sectional (laminate_a11, native_bay_centres), run_s2_static_cases,
run_s4_rotating_modal, run_s5_oneway, run_s6_oneway_dynamic, run_s7_torsion, l3_mode_shapes.
MIGRATE: run_blade_structural_matrix, l3_torque_sensitivity (drop numad/numad-mid/windio2 choices), mesh_quality
(drop windio2 import 225), run_rigid_bem_check, s2_ring_metric_blade, l3_box_energy, l3_s1_blade, l3_s4_blade,
l3_bredt_mesh, l3_twist_integral (IEA-3.4 -> IEA34 upstream), t6_ccx_centrifugal_buckling, t7_panel_stiffness_parity,
t10_ccx_arbiter, s2_bd_equivalence, s4_modal_convergence, v02_modal_convergence, blade_ccx_convergence,
yaw_sweep_plots, yaml_layup_independent, l5_t9_psi_osc, l5_twist_sense_chain (if still used).
DELETE (one-off diagnostics; results already in odd/ and logs): the diagnose_* naming v1/xlsx/BladeMesh
(newmark_only, clt_compliance_steiner, bending_compliance_2x2, te_layup_dump, twist_sense_frame, edge_dip_locator,
edge_shear_split, edge_chord_ei, sign_chain, mode_beam_projection, shear_centre, fsi_deformed, zhou_loads_reverse,
damping_matrix, rotating_assembly, gc_twist_chord, wall_flow_static_ab, rotating_prestress, s1_gj_tube, s1_gj_station,
section_shear, section_moment_realization, vlm_sectional, yaw_aero); s1_midsurface_counterfactual, s1_ei_partition,
s1_clt_steiner_ei, s1_reference_sweep; probe_s2_edge_convergence, run_s2e_diagnostic, run_s8a_gj_circumferential,
yaml_layup_vs_shell, blade_mass_diagnostic, blade_gravity_kinematic_check, blade_gravity_nl_increment,
blade_gravity_beam_shortening, debug_kg_magnitudes, audit_blade_input, validate_vlm_vs_bem_sweep, l3_t12_reconcile,
l3_w3_lc3_integral. Final `grep -l` before deleting (~67 files name v1).

## D. Data and scripts
DELETE: tests/IEA-15-240-RWT.yaml (+ its .ech), tests/IEA-15-240-RWT_rigid.yaml, tests/IEA-3.4-130-RWT.yaml,
tests/NuMAD_utd_iea15mw.xlsx, tests/archive/numad_xlsx_configs/, tests/exp_aero_damp_off/IEA-15-240-RWT.yaml (check
readers), tests/IEA15MW/frontiersin_2025_yaw_anim/solid_corotational_rigid_yaw_0.yaml, tests/IEA15MW/run_yaw0_anim_rigid.srm.
EDIT: tests/run_step1c_mp_off.srm:60, comments in tests/IEA15MW/submit_ch6_corotational.sh:10, submit_yaw_batch.sh:83,
submit_ch6_matrix.py:35; README.md 194-200 (CLI example -> upstream v2; keep pyNuMAD attribution while numad/objects
and io stay vendored); docs/validation/contract.yaml (~20 rows), levels.yaml:110, 144. AGENTS.md:75 (rule 7) and
.pi/skills/fem-shell/SKILL.md:100 name the v1 yaml: edit at the end, outside a working session's prompt prefix.
Do not touch tests/IEA15MW/validation papers/ or tests/IEA34/reference.

## E. Anchors whose test runs on v1/NuMAD only
| anchor | v2/gmsh replacement | status |
|---|---|---|
| S-0, S-3, box EI, D-Tube, UL elastica | none needed | unaffected |
| S-1 | v2/test_v2_iea15mw_s1_sectional | at risk: "S-1 cannot be read on gmsh" (campaign L8) unless the v2 twin reads gmsh |
| S-2 | v2/..._s2_static_prescribed (windio2 default) | switch to gmsh; UTD rows retire |
| S-4, S-7, V-02, V-03, V-05 | v2 twins | switch to gmsh |
| S-5 | tools/run_s5_oneway (gmsh branch at 164, v1 at 174) | migrate the tool |
| S-6 | tools/run_s6_oneway_dynamic (v1 at 99) | MIGRATE (NaN open) |
| V-01 | none (v1 polars via load_blade_aero) | MIGRATE fixture; numbers may move |

## Execution order (commits of ~400 changed lines)
1. Move shared helpers (get_direction_cosines; Pchip, _centroid, RawMesh), fix gmsh_skin imports/docstrings; v2_blade
   gmsh-only + digest (invalidates cache).
2. src dead code: generators.py dead methods, generate_blocks_sections, excel_file plumbing.
3. Delete NuMAD-only tests and v1 twins (B1); then the NuMAD cluster, windio2 remainder, geometry_mesh, excel reader.
4. Reject v1 in windio.py; delete v1 branches of yaml_to_blade (after every v1 input is gone).
5. Test migrations (B2) in 3 commits: aero fixture (V-01), blade structural guards, parity/CLI/rotor; moved numbers logged.
6. Tools: migrate (2 commits), delete groups (2 commits).
7. Data and scripts (D).
8. Full suite on SLURM, known-red re-record.
Riskiest: step 5 (v1-measured references move). Second: step 1 (cache rebuild, large meshes > 8 GB).
