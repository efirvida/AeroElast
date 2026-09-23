# Validación seccional del shell contra las tablas de viga: cadena de evidencia

Este documento registra la investigación (2026-09-10) sobre la diferencia
observada entre las propiedades seccionales extraídas del modelo shell y las
tablas de viga publicadas en los archivos OpenFAST. La conclusión: **la
diferencia es la discrepancia entre las dos representaciones estructurales
que publica el propio yaml WindIO, no un error del shell**.

## 1. La observación

Comparación mid-span (0.2–0.85R) del shell oficial (WindIO, element_size 0.5)
contra las tablas de viga:

| Propiedad | Shell vs tabla | Fuente de la tabla |
|-----------|----------------|--------------------|
| Masa m(r) | +4–8 % | ElastoDyn BMassDen / yaml M[0,0] |
| EI_flap   | 10–25 % (alto sistemático, mezcla flap/edge por el twist) | BeamDyn K 6×6 |
| EI_edge   | ~10 % | BeamDyn K 6×6 |
| **GJ**    | **+17 %** (estable en 4 estimadores) | BeamDyn K 6×6 |
| EA        | +20 % (CLT directa del layup) | BeamDyn K 6×6 |

El GJ se verificó con cuatro estimadores independientes del twist rate
(splines s=1e-3/1e-2/1e-1 y diferencias locales): 17.5–24.5 %, por lo que el
sesgo no es un artefacto del método de extracción.

## 2. Verificación de la cadena de entrada

El yaml WindIO (`tests/IEA-15-240-RWT.yaml`) publica **dos representaciones
estructurales lado a lado**:

1. `elastic_properties_mb.six_x_six` — las matrices 6×6 de viga (K y M), y
2. `internal_structure_2d_fem` — el layup 2D (webs + capas) del que se
   construye el shell.

Verificado por comparación directa:

- `six_x_six.stiff_matrix` ≡ el archivo `IEA-15-240-RWT_BeamDyn_blade.dat`
  (matrices K): diferencia máxima **3.9e-12** — idénticas a precisión de
  máquina.
- Twist: idéntico (diferencia 7.8e-15; el campo del yaml está en radianes).
- Chord: idéntico salvo la última estación (r/R=0.995 — convención de punta,
  1.47 m vs 0.82 m).
- Masa: +4.1 % (el M[0,0] del yaml vs el BMassDen de ElastoDyn — la misma
  diferencia del V-03).

Conclusión: **la referencia de viga de OpenFAST es exactamente la
representación 6×6 del mismo yaml** — no hay modelos distintos entre
repositorios; la diferencia vive dentro del yaml, entre sus dos
representaciones.

## 3. Jueces independientes

### 3.1 Estimado Bredt-Batho (sección cerrada) y CLT directa

Para una sección hueca de pared delgada, la rigidez torsional tiene solución
analítica (Bredt-Batho): GJ = 4A²/∮(ds/t), con A el área encerrada. Evaluada
con el A66 del layup del propio yaml (CLT por capa), el estimado trackea al
shell (20 % mutuo, limitado por el mapeo por-vértice del estimado) y ambos
quedan por encima de las 6×6 — es decir, el shell es **CLT-consistente con
su layup**.

El EA se resolvió de la misma manera (2026-09-10): la extracción por
respuesta global (carga axial en punta → derivada del desplazamiento
tangente-proyectado) queda contaminada por el pre-bend — el bending
inducido por la excentricidad axial del pre-bend proyecta sobre la tangente
con magnitud comparable a la señal axial. El método inmune es la **CLT
directa**: EA(z) = Σ sobre los elementos de la estación de A11(layup)·(ancho
spanwise), con A11 = Σ Q̄11_k·h_k. Resultado: **+20 %** vs el BeamDyn — la
misma familia que el GJ (+17 %). El método de respuesta queda verificado en
el cajón (test_box_ea_matches_analytic, ratio 1.000) — la contaminación es
específica de la geometría pre-doblada.

### 3.2 Benchmark analítico del cajón de sección cerrada

`tests/test_box_torsion_bending_benchmark.py` — cajón rectangular de pared
delgada (w=h=1 m, t=0.01 m, L=10 m, isotrópico E=70 GPa) con soluciones
cerradas exactas:

- **Torsión**: GJ numérico = 2.6320e8 N·m² vs el exacto G·2w²h²t/(w+h) =
  2.632e8 — **coincidencia a la precisión del test** (tol 2 %). El camino de
  corte de membrana del MITC4 (incluido el esquema SRI) es correcto para la
  torsión de secciones cerradas.
- **Flexión**: EI por curvatura en el mid-span = ratio **1.000** vs el
  I_x = wt·h²/2 + t·h³/6 analítico. (Nota metodológica: la deflexión en
  punta no sirve de medida — las cargas puntuales distorsionan la sección
  localmente; las esquinas dan +2.7 % y el mid-width se abulta. La
  curvatura es inmune.)

## 4. Conclusión

Con la torsión y la flexión del elemento validadas analíticamente, la banda
observada contra las tablas de viga queda **probada como la diferencia entre
las dos representaciones del yaml** (las 6×6 generadas por una herramienta de
sección — PreComp/BECAS — y el layup 2D FEM que modela el shell). El shell
no es la fuente del error.

## 4b. La física de ElastoDyn verificada en el fuente (2026-09-13)

Revisión del código fuente de ElastoDyn (github.com/OpenFAST/openfast,
modules/elastodyn, ElastoDyn.f90) para complementar la comparación con
evidencia a nivel de formulación, no solo de resultados:

1. **La pala es un modelo modal de 3 DOFs** (2 flap + 1 edge) con formas
   modales **polinomiales** (el orden 6; la función `SHP` evalúa
   φ(z) = Σ cᵢ·z^(i+1) y sus derivadas). La rigidez generalizada es
   K = Σ `StiffBF × ΔR × φ²` sobre los elementos de la viga 1D.

2. **El stiffening centrífugo es 1D, basado en la tensión**:
   ```
   ElmntStff = FMomAbvNd(K,J) * p%DRNodes(J) * p%RotSpeed**2
   KBFCent   = KBFCent + ElmntStff * Shape'^2
   ```
   donde `FMomAbvNd` es el primer momento de masa por encima del nodo —
   exactamente la tensión N(z) = ω²∫ρ(s)·s ds de la viga 1D. El K_Cent de
   ElastoDyn = ∫N(z)(φ')² dz, la rigidez geométrica 1D clásica.

3. **El twist** entra solo como las segundas derivadas de las shapes
   torsadas (un acople geométrico en las formas modales), y **no existe
   DOF de torsión**.

4. La torre incluye el **destiffening gravitatorio** (el P-Δ de la
   gravedad sobre la torre — `KTFAGrav`).

**Diferencias con Aeroelast**: ElastoDyn usa la rigidez geométrica 1D de
tensión (la N(z)(φ')²); Aeroelast usa la matriz de stress inicial 3D
completa (el tensor σ × los tres bloques de gradientes, con los frames de
sección rotados por el twist). La consecuencia observada — el ~0-2% de
stiffening estático en ElastoDyn (la tensión ponderada por las formas
polinomiales pesa poco en la deflexión estática) frente al ~13% del shell
(la formulación 3D responde más) — queda así explicada a nivel de fuente.
El modal (+8.5% en ElastoDyn) es la misma tensión ponderada por las
pendientes de los polinomios. La ausencia de torsión y la truncación modal
explican también por qué `IPDefl1` no es comparable con el u_x del shell.

## 5. Consecuencias para la estrategia de validación

1. **Claims absolutos del shell**: benchmarks analíticos (E-01/E-02 y el
   cajón torsional de este documento). Cero dependencia de tablas externas.
2. **Comparación contra las tablas de viga**: se reporta como chequeo de
   consistencia entre las dos representaciones del yaml, con las bandas
   documentadas (GJ +17 %, EI 10–25 %, masa +4–8 %) y esta explicación.
3. **Tests en el repo**: `test_iea15mw_s1_sectional.py` (incluido el GJ con
   tolerancia 30 %) y `test_box_torsion_bending_benchmark.py` son guards de
   regresión alrededor de las bandas documentadas — no claims de acuerdo
   exacto con las tablas.
4. **Dirimir cuál representación es la "verdadera"** requeriría un tercer
   juez externo (BECAS/VABS sobre el layup 2D) o datos físicos de la pala —
   no bloquea el paper.

## 6. Re-medición post-fix del twist (2026-09-18)

El fix de signo del twist (`rotorspin = -1`, 2026-09-15) cambió la geometría
de la pala (twist nose-up como el yaml) y con eso la extracción de GJ por
respuesta estática pasó de 18.9 % a 31.1 % medio |rel| — fuera de la banda
física documentada (+17 %).

**Causa**: el par de punta se aplicaba como fuerzas globales-y (par sobre z
global) y el twist se medía como rotación sobre z global. En la pala
pre-doblada (sweep) el par global-z tiene componente de flexión, y la
rotación global-z mezcla flexión con torsión — el mismo artefacto ya
documentado para EI_flap y el canal EA de respuesta global.

**Fix**: el par de punta ahora se aplica sobre la **tangente local**
pre-doblada de la sección de punta y la rotación se mide proyectada sobre
esa tangente (`tools/run_s1_sectional.py`, `ea_gj_from_static_response`).

**Resultado post-fix**: GJ mid-span +23.9 % medio |rel| (rango −12 a +55 %),
alineado con el juez CLT/Bredt del layup (~+20 %) y con la familia física
+17–25 % de las dos representaciones del yaml. El test S-1 queda con
tolerancia 30 % como guarda de regresión.
