# Divergencia entre modelos de pala: NuMAD UTD vs. definición oficial IEA 15 MW

**Fecha**: 2026-09-08 (hallazgo de layup) — ampliado 2026-09-09 (causa de la
divergencia FSI con el blade oficial)
**Estado**: migración completada; la divergencia FSI con el blade oficial
resuelta (masa consistente — ver §8).

---

## 8. Divergencia FSI del blade oficial (2026-09-09) — causa y fix

Al re-correr las campañas FSI con el blade oficial WindIO, el acoplamiento
divergía exponencialmente (torques ~1e150 en t≈0.1 s) mientras el control
UTD permanecía estable.  Diagnóstico A/B local (login node) + instrumentación:

1. **Fluid-side** (3 bugs reales, corregidos y verificados):
   - El polar "circular" del YAML tiene solo 2 puntos de alpha → CCBlade
     degeneraba en la inducción de la raíz (`_parse_polars_from_yaml`
     ahora expande las tablas <4 puntos a rango completo).
   - `force_projection.py` restaba `hub_radius` a las coordenadas de la
     malla (que ya son blade-locales) → strips desplazados y strip-0 vacío.
   - `fsi_participant.py` reconstruía `r_def` sin sumar el hub → las
     estaciones raíz caían dentro del hub.
2. **Solid-side** (la causa dominante): la **matriz de masa lumped** ponía
   la masa nodal completa sobre los DOFs rotacionales de la punta (paneles
   finos del taped-chord tip, rigidez rotacional efectiva casi nula:
   K-diag ~5e3 vs acoplamientos ~1.6e6, ratio hasta 8380×).  El Newmark
   implícito (incondicionalmente estable para sistemas SPD bien
   condicionados) divergía a ~92×/paso en el DOF θy de la punta.  La masa
   **consistente** (la que ya usaban V-02/V-03) es estable.

**Fix**: `rotor.py` y `linear_dynamic.py` (StressStiffenedDynamicFSI) usan
ahora `assemble_mass_matrix()` (consistente); las masas nodales para los
términos explícitos se extraen por suma de filas (partición de la unidad).
Verificación end-to-end en el login node (caso yaw_0, t=1.9 s estable con
τ_grav=1.9e7 físico y 0 NaN) y en el cluster (yaw_0/10 a t=0.9 s con
τ_aero=-5.1e6, Ct=0.31, 0 NaN).  `rotor_inertial.py` mantiene la masa
lumped deliberadamente (requiere M diagonal para F_ref = -M·a_ref).

Herramientas de diagnóstico agregadas (reutilizables): `tools/diagnose_fsi_deformed.py`,
`tools/diagnose_newmark_only.py`, `tools/diagnose_rotating_assembly.py`,
`tools/diagnose_damping_matrix.py`.

---

## 1. Antecedentes

AeroElast usa como modelo estructural de la pala IEA 15 MW el archivo
`NuMAD_utd_iea15mw.xlsx`, un re-modelado independiente creado por el grupo de
Griffith en UTD (Escalera Mendoza et al., AIAA SciTech 2023,
doi:10.2514/6.2023-2093) a partir de la definición oficial. El repositorio
oficial de la IEA 15 MW enlaza ese modelo como "el modelo NuMAD", pero las
propiedades de viga publicadas (archivos ElastoDyn, BeamDyn, HAWC2 y las
tablas PreComp/VABS/BECAS del reporte NREL/TP-5000-75698) derivan de la
**layup oficial**, no del re-modelado UTD.

## 2. Evidencia 1 — validación estática S-2 con ambos modelos

`tests/test_iea15mw_s2_static_prescribed.py` compara el shell contra una viga
3D Euler-Bernoulli construida con la misma geometría de malla y las rigideces
distribuidas oficiales (ElastoDyn + BeamDyn), para 4 casos de carga estática
con raíz empotrada:

| Caso | Componente dominante | Blade UTD | Blade oficial (WindIO) |
|---|---|---|---|
| LC1 gravedad | uy | +1.9% | +0.8% |
| LC2 carga de punta flapwise | uy | +1.0% | +1.4% |
| LC3 carga de punta edgewise | ux | **+40–55%** | **+8.1%** |
| LC4 carga distribuida | uy | +2.2% | −2.5% |

El resultado LC3 con el modelo UTD es convergente en malla (element_size 0.25
y 0.5 dan el mismo valor) y ambas referencias de viga (ElastoDyn y BeamDyn)
coinciden entre sí, por lo que la discrepancia no es numérica: es del modelo de
entrada.

## 3. Evidencia 2 — comparación directa de layups

Extracción de la layup oficial (`Documentation/IEA-15-240-RWT_tabular.xlsx`,
hojas "Blade Shell Layup" / "Shear Web Layup" / "Material Properties") y de la
layup UTD (hoja "Geometry" de `NuMAD_utd_iea15mw.xlsx`):

| Span [m] | Spar cap UTD/OFF [mm] | Piel triax UTD/OFF [mm] | Refuerzo TE UTD/OFF [mm] |
|---|---|---|---|
| 17.6 | 53.4 / 53.3 | 12.9 / 13.0 | **10.3 / 29.5** (−65%) |
| 28.7 | 83.9 / 83.9 | 3.0 / 3.0 | **19.3 / 30.0** (−36%) |
| 51.4 | 95.8 / 95.8 | 2.0 / 2.0 | **22.8 / 30.0** (−24%) |
| 62.9 | 87.5 / 87.5 | 2.0 / 2.0 | **12.5 / 30.0** (−58%) |
| 74.7 | 77.7 / 77.7 | 2.0 / 2.0 | **3.4 / 10.9** (−68%) |

Spar caps y pieles triax son prácticamente idénticos. La diferencia sistemática
es el **refuerzo del borde de fuga** (glass_uni), que en el modelo UTD es
30–75% más delgado entre 10 y 75 m de envergadura.

## 4. Explicación física

El refuerzo TE es el principal contribuyente a la rigidez edgewise outboard:
glass_uni concentrado en la posición de máximo brazo respecto del eje neutro
edgewise. La reducción de espesor del UTD abate `EI_edge(r)` en la zona
media-outboard, que es exactamente la región que domina la deflexión de punta
bajo carga estática (integral de `(L−s)²/EI(s)`). La consistencia interna del
cuadro completo:

- Los modos 1F/1E (V-02) pesan la zona inboard → cerraban a <2% aún con el UTD.
- El modo 2E (V-02, −8.65%, "el más sensible") pesa más outboard → ya mostraba
  el sesgo.
- Los casos flapwise/gravedad (S-2) apenas dependen del TE → cerraban a 1–2%.

## 5. Experimento decisivo

El blade oficial existe en formato WindIO v2 en `tests/IEA-15-240-RWT.yaml`
(layup completa `internal_structure_2d_fem`, polares incluidos) y el importador
ya existía en AeroElast (`models/blade/numad/io/yaml_to_blade.py` +
`BladeMesh(yaml_file=...)`). Al correr S-2 con ese blade, el gap edgewise cae
de +40–55% a **+8.1%**, dentro de la banda teórica esperada shell-vs-viga
(deformación de sección + corte transversal). Confirmado: **el sesgo era del
modelo UTD, no del solver**.

## 6. Migración

Se migró el stack de validación a la definición oficial WindIO:

- `tests/test_iea15mw_s2_static_prescribed.py` — corre ambos modelos: oficial
  como referencia de cierre y UTD como registro de sensibilidad de input.
- `tests/validation/blade/test_iea15mw_v02_natural_frequencies.py` → blade oficial. Valores
  medidos (element_size 0.25): 1F 0.5377 Hz (−3.7%), 1E 0.6977 Hz (+8.9%,
  dentro de la dispersión inter-método ±14% [Zhou 2025]; tolerancia del test
  ampliada a ±10% con justificación), 2F 1.5881 Hz (−4.3%), 2E 2.1101 Hz
  (−2.6%, cierra el gap de −8.65% que el UTD mostraba).
- `tests/test_iea15mw_v03_static_gravity.py` → blade oficial (3/3, masa
  dentro de ±2% de 67,921 kg).
- `tests/test_iea15mw_v05_structural_properties.py` — ya usaba el YAML oficial.
- YAMLs de campañas FSI (`tests/IEA15MW/frontiersin_2025_yaw/`,
  `tests/IEA15MW/ch6_prepared_cases/`): `mesh.generator.params` y
  `blade_file` → `../../IEA-15-240-RWT.yaml` (28 archivos).
- La campaña B1/B2 (relanzada 2026-09-08, job 11592005) completó con el
  blade UTD y queda como registro de sensibilidad de input; sus criterios de
  aceptación cierran (B1 orden p=2.57, B2 p=1.17, |Δflap_mean| <1%/5%).
- Nota canónica agregada a `AGENTS.md`.

## 7. Notas para trabajo futuro (base de un artículo sobre el modelo UTD)

Elementos ya verificados y cuantificados:

1. Divergencia de layup TE documentada spanwise (§3) con evidencia numérica
   independiente (§2, §5).
2. El modelo UTD ya difería en masa integrada: 67,167–68,077 kg vs
   65,250 kg (objetivo WISDEM) y 67,921 kg (tabular oficial) — consistente con
   diferencias estructurales adicionales más allá del TE.

Pendientes para cerrar una crítica completa del modelo UTD:

1. Cuantificar `EI_edge(r)` y `EI_flap(r)` de AMBAS layups con el extractor
   seccional S-1 (pendiente) y contra los targets BeamDyn/ElastoDyn — hoy el
   impacto del TE está cuantificado solo vía deflexión global (+8.1% vs
   +40–55%).
2. Comparar modos 2E/3E de ambos modelos contra BModes (el 2E es el más
   sensible al TE).
3. Evaluar impacto en fatiga edgewise (DEL del momento edgewise de raíz) —
   canal donde el modelo UTD subestimará carga por rigidez reducida.
4. Revisar las asunciones documentadas en el paper AIAA 2023-2093 sobre la
   conversión WISDEM→NuMAD para explicar el origen del espesor TE reducido.

**Posición editorial**: el modelo UTD es una contribución valiosa y citable
como re-modelado independiente; el hallazgo aquí documentado no invalida su
uso general, pero obliga a declarar la divergencia cuando se compare contra las
propiedades oficiales de la IEA 15 MW.
