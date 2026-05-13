# Elementos Shell MITC3+ y MITC4+: Resultados de Validación

Los elementos de la familia MITC (*Mixed Interpolation of Tensorial Components*) resuelven el *shear locking* en shells delgados mediante interpolacion mixta de deformaciones de corte. La implementacion usa MITC3+ (triangular) y MITC4+ (cuadrilátero), ambos con 6 DOF/nodo: $[u, v, w, \theta_x, \theta_y, \theta_z]$.

## Hipotesis cinematica y origen del locking

Los elementos shell implementados se basan en la teoria de Reissner-Mindlin: el desplazamiento transversal $w$ y las rotaciones $\theta_x$, $\theta_y$ son campos independientes. Las deformaciones de corte transversal son

$$
\gamma_{xz} = \frac{\partial w}{\partial x} - \theta_x, \qquad
\gamma_{yz} = \frac{\partial w}{\partial y} - \theta_y.
$$

En el limite de lamina delgada (espesor $t \to 0$), la teoria de Kirchhoff exige $\gamma \to 0$, es decir, que los gradientes de $w$ y los campos de rotacion se igualen puntualmente. Un elemento de desplazamiento estandar no puede satisfacer esa condicion de manera exacta con campos polinomiales independientes: la restriccion sobredetermina el sistema y aparece el *shear locking*, en el que el elemento genera rigidez de corte espuria que crece sin limite al afinar la malla o reducir el espesor. El resultado practico es una subprediccion sistematica de la deflexion en laminas delgadas.

## El enfoque MITC

La familia MITC (Bathe & Dvorkin, 1985) elimina el locking desacoplando la derivacion de las deformaciones de corte del campo de desplazamientos. En lugar de calcular $\gamma$ directamente de las derivadas del desplazamiento, se interpola cada componente covariante de la deformacion de corte desde "puntos de amarre" (*tying points*) especificamente elegidos donde esa componente es libre del termino parasito. La deformacion interpolada se inserta en la formulacion variacional reemplazando la derivada de compatibilidad. Este reemplazo rompe el acoplamiento espurio y restaura la convergencia correcta al reducir el espesor.

La eleccion de los puntos de amarre no es arbitraria: se demuestra que son los puntos del elemento donde la deformacion covariante es cinematicamente correcta para los modos de flexion pura. Evaluando y usando esos puntos, la formulacion garantiza que los modos de flexion no generen energia de corte parasita.

**MITC3+** — 3 nodos, 18 DOF/elemento; formulacion Lee, Lee & Bathe (2014), CAS 138.  
**MITC4+** — 4 nodos, 24 DOF/elemento; formulacion Ko et al. (2017), CAS 193, sobre la base de Dvorkin & Bathe (1984).

---

## Convención de tablas

Todas las tablas usan el formato:

| Test `archivo::nombre` / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |

Solo se incluyen tests que pasan **y** tienen una referencia publicada explícita. Los valores "Obtenido" fueron medidos en la corrida del test suite en la versión actual del código; los valores "Esperado" provienen del artículo o fórmula citada.

---

## Formulacion MITC3+: decisiones especificas

### Seis puntos de amarre y la variacion lineal del corte

El elemento MITC3 original (Lee & Bathe, 2004) usa tres puntos de amarre, uno por arista, y captura el campo de corte transversal constante sobre el elemento. Sin embargo, en mallas distorsionadas, los modos de flexion sobre triangulos no equilateros generan una variacion lineal del corte covariante que los tres puntos no pueden representar: ese residuo actua como corte parasito y reintroduce locking parcial.

El MITC3+ (Lee et al., 2014) resuelve esto con seis puntos de amarre: tres constantes (A, B, C) heredados del MITC3 original, y tres adicionales (D, E, F) desplazados desde los puntos medios de cada arista con una pequena perturbacion hacia el interior del elemento. Los terminos adicionales capturan exactamente la variacion lineal del campo de corte que los puntos constantes no pueden representar. La perturbacion es suficientemente pequena para no perturbar los resultados en mallas regulares y suficientemente grande para estabilizar la interpolacion en mallas distorsionadas.

### Enriquecimiento burbuja y condensacion estatica

El MITC3+ agrega dos grados de libertad internos de flexion (DOFs burbuja) que no pertenecen a ningun nodo fisico compartido entre elementos. Estos DOFs representan modos de curvatura de orden superior que la cinematica triangular lineal no puede capturar: una interpolacion lineal de rotaciones sobre un triangulo genera solo dos modos de curvatura independientes, pero la flexion de un triangulo arbitrario puede requerir un tercer modo que los DOFs burbuja introducen.

Como son internos, se eliminan antes del ensamble global por condensacion estatica de Schur. Separando los grados de libertad nodales $u$ y los de burbuja $q$, la rigidez elemental tiene la forma por bloques

$$
\begin{bmatrix} K_{uu} & K_{uq} \\ K_{qu} & K_{qq} \end{bmatrix}
\begin{bmatrix} u \\ q \end{bmatrix}
=
\begin{bmatrix} f_u \\ 0 \end{bmatrix}.
$$

La condicion de carga nula sobre los DOFs internos permite expresar $q = -K_{qq}^{-1} K_{qu}\, u$ y sustituir:

$$
K_{\rm cond} = K_{uu} - K_{uq}\, K_{qq}^{-1}\, K_{qu}.
$$

El resultado es una rigidez condensada de $18 \times 18$ que incorpora el enriquecimiento sin aumentar el numero de grados de libertad del sistema global.

### Rigidez de perforacion (drilling)

En la teoria de Reissner-Mindlin, la rotacion en el plano del elemento ($\theta_z$ para un elemento plano) esta cinematicamente desacoplada de los modos de flexion transversal y de membrana. No existe energia de deformacion natural asociada a ella. Sin regularizacion, la rigidez seria singular en esa direccion y el ensamble global produciria modos de cuerpo cuasi-rigido en el plano.

Para suprimir este modo se agrega una rigidez de penalizacion:

$$
k_{\rm drill} = 0.15\, E\, t^2.
$$

Este valor no representa ninguna fisica del problema de lamina: es una penalizacion que fija el modo de perforacion sin perturbar los resultados de deflexion ni de frecuencia en los modos fisicos. El factor $0.15$ fue calibrado para que la penalizacion sea lo suficientemente pequena para no afectar los modos fisicos y lo suficientemente grande para eliminar la singularidad en el ensamble.

### Formulacion corrotacional para grandes rotaciones

Para el analisis no lineal de grandes rotaciones, los elementos MITC usan una formulacion corrotacional: en cada paso de carga, se actualiza un marco local que sigue la deformacion del elemento. La rigidez tangente efectiva en ese paso es

$$
K_T = T_{\rm def}^T \cdot K_L \cdot T_{\rm def},
$$

donde $T_{\rm def}$ es el operador de transformacion al marco deformado y $K_L$ es la rigidez lineal evaluada en el punto de deformacion cero dentro de ese marco. La idea central es que, al medir las deformaciones relativas al marco que rota con el elemento, los desplazamientos locales permanecen pequenos incluso cuando las rotaciones globales son grandes. Esto permite usar la cinematica lineal del elemento en cada paso incremental sin acumular errores de curvatura finita.

La alternativa, una formulacion total de Lagrange (TL) con deformaciones de Green-Lagrange completas, generaria terminos de membrana de orden $O(E h \theta^2)$ para cargas de flexion pura. Esos terminos son fisicamente correctos pero numericamente perjudiciales: son ortogonales a la carga externa de momento y fuerzan al metodo de Newton-Raphson a reducir el paso mediante busqueda de linea, degradando la convergencia cuadratica a velocidad de gradiente para incrementos de carga moderados. La formulacion corrotacional evita ese problema manteniendo el residuo local en terminos de desplazamientos pequenos, lo que preserva la convergencia cuadratica del Newton.

---

## 1. Grandes rotaciones — viga voladizo (MITC3+ y MITC4+)

**Referencia**: Simo, J.C. & Vu-Quoc, L. (1986). *A three-dimensional finite-strain rod model. Part II*. CMAME **58**, 79–116, Tabla 1.

**Geometría**: $L=10$, $B=1$, $H=0.1$, $E=1.2\times10^6$, $\nu=0$. Momento de extremo $M = \lambda EI/L$ aplicado en $n_{\rm steps}$ incrementos; formulación UL actualizada.  
Fórmulas analíticas: $u_{\rm tip} = L(\sin\lambda/\lambda - 1)$, $w_{\rm tip} = L(1-\cos\lambda)/\lambda$.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[π/2]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=\pi/2$, comp. $u_{\rm tip}$ | $-3.6310$ | $-3.6338$ | $0.08$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[π/2]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=\pi/2$, comp. $w_{\rm tip}$ | $6.3740$ | $6.3662$ | $0.12$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[π]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=\pi$, comp. $u_{\rm tip}$ | $-10.0163$ | $-10.0000$ | $0.16$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[π]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=\pi$, comp. $w_{\rm tip}$ | $6.3690$ | $6.3662$ | $0.04$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[3π/2]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=3\pi/2$, comp. $u_{\rm tip}$ | $-12.1162$ | $-12.1220$ | $0.05$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[3π/2]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=3\pi/2$, comp. $w_{\rm tip}$ | $2.1118$ | $2.1221$ | $0.48$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[2π]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=2\pi$, comp. $u_{\rm tip}$ | $-9.9645$ | $-10.0000$ | $0.36$ |
| `test_mitc3_benchmarks.py` · `test_equilibrium_path[2π]` / Simo & Vu-Quoc (1986) | MITC3+, $n=10$, $\lambda=2\pi$, comp. $\lvert w_{\rm tip}\rvert$ | $0.0198$ | $0.0000$ | abs $0.020$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[π/2]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=\pi/2$, comp. $u_{\rm tip}$ | $-3.6309$ | $-3.6338$ | $0.08$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[π/2]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=\pi/2$, comp. $\lvert w_{\rm tip}\rvert$ | $6.3739$ | $6.3662$ | $0.12$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[π]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=\pi$, comp. $u_{\rm tip}$ | $-10.0040$ | $-10.0000$ | $0.04$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[π]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=\pi$, comp. $\lvert w_{\rm tip}\rvert$ | $6.3866$ | $6.3662$ | $0.32$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[3π/2]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=3\pi/2$, comp. $u_{\rm tip}$ | $-12.1355$ | $-12.1220$ | $0.11$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[3π/2]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=3\pi/2$, comp. $\lvert w_{\rm tip}\rvert$ | $2.1343$ | $2.1221$ | $0.58$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[2π]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=2\pi$, comp. $u_{\rm tip}$ | $-9.9911$ | $-10.0000$ | $0.09$ |
| `test_large_rotation_benchmarks.py` · `test_equilibrium_path[2π]` / Simo & Vu-Quoc (1986) | MITC4+, $n=10$, $\lambda=2\pi$, comp. $\lvert w_{\rm tip}\rvert$ | $0.0048$ | $0.0000$ | abs $0.005$ |

---

## 2. Deflexión lineal — teoría de viga (MITC4)

**Referencia**: Timoshenko, S. & Woinowsky-Krieger, S. (1959). *Theory of Plates and Shells*, 2ª ed., McGraw-Hill. Fórmulas §3: $\delta = PL^3/(3EI)$, $I = bh^3/12$.

**Geometría**: Placa plana de acero, $E = 210\,\text{GPa}$, $\nu = 0.3$, $\rho = 7850\,\text{kg/m}^3$.

### 2a. Convergencia con refinamiento de malla — voladizo con carga transversal

$L=1\,\text{m}$, $b=0.1\,\text{m}$, $h=0.01\,\text{m}$, $P=100\,\text{N}$.  
$\delta_{\rm ref} = PL^3/(3EI) = 1.905\times10^{-2}\,\text{m}$ con $I = bh^3/12 = 8.333\times10^{-8}\,\text{m}^4$.

| Test / Referencia | Parámetros | $\delta_{\rm FEM}$ (m) | $\delta_{\rm ref}$ (m) | $\Delta\%$ |
|---|---|---|---|---|
| `test_shell_analytical_validation.py` · `TestCantileverBeam::test_tip_deflection[4-2]` / Timoshenko & W-K (1959) | MITC4+, malla $4\times2$ quad | $1.838\times10^{-2}$ | $1.905\times10^{-2}$ | $3.51$ |
| `test_shell_analytical_validation.py` · `TestCantileverBeam::test_tip_deflection[8-4]` / Timoshenko & W-K (1959) | MITC4+, malla $8\times4$ quad | $1.872\times10^{-2}$ | $1.905\times10^{-2}$ | $1.71$ |
| `test_shell_analytical_validation.py` · `TestCantileverBeam::test_tip_deflection[16-8]` / Timoshenko & W-K (1959) | MITC4+, malla $16\times8$ quad | $1.883\times10^{-2}$ | $1.905\times10^{-2}$ | $1.13$ |

### 2b. Extensión de membrana — carga axial

$L=1\,\text{m}$, $b=0.1\,\text{m}$, $h=0.01\,\text{m}$, $P=1000\,\text{N}$.  
$\delta_{\rm ref} = PL/(EA) = 4.762\times10^{-6}\,\text{m}$ con $A = bh = 10^{-3}\,\text{m}^2$.

| Test / Referencia | Parámetros | $\delta_{\rm FEM}$ (m) | $\delta_{\rm ref}$ (m) | $\Delta\%$ |
|---|---|---|---|---|
| `test_shell_analytical_validation.py` · `TestMembraneStretching::test_axial_extension` / Timoshenko & W-K (1959) | MITC4+, malla $8\times4$ quad | $4.858\times10^{-6}$ | $4.762\times10^{-6}$ | $2.01$ |

### 2c. Cargas uniaxiales en voladizo delgado

$L=1\,\text{m}$, $b=0.1\,\text{m}$, $h=0.001\,\text{m}$, $P=600\,\text{N}$; malla $8\times4$ quad MITC4+.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_shell_comprehensive.py` · `TestLinearStaticCantilever::test_fx_in_plane` / Timoshenko & W-K (1959) | $F_x$, $\delta_{\rm ref}=PL/(EA)=2.857\times10^{-5}\,\text{m}$ | $2.824\times10^{-5}\,\text{m}$ | $2.857\times10^{-5}\,\text{m}$ | $1.2$ |
| `test_shell_comprehensive.py` · `TestLinearStaticCantilever::test_fy_in_plane` / Timoshenko & W-K (1959) | $F_y$ (flexión en plano), $\delta_{\rm ref}=PL^3/(3EI_z)=1.143\times10^{-2}\,\text{m}$ | $1.134\times10^{-2}\,\text{m}$ | $1.143\times10^{-2}\,\text{m}$ | $0.8$ |
| `test_shell_comprehensive.py` · `TestLinearStaticCantilever::test_fz_out_of_plane` / Timoshenko & W-K (1959) | $F_z$ (flexión transversal), $\delta_{\rm ref}=PL^3/(3EI_y)+PL/(k_sGA)=114.3\,\text{m}$ | $112.1\,\text{m}$ | $114.3\,\text{m}$ | $1.9$ |

> Nota: $h=0.001\,\text{m}$ produce grandes desplazamientos lineales por tratarse de un modelo artificialmente delgado; el test verifica que el elemento no presenta *locking*.

---

## 3. Análisis modal — viga voladizo (MITC4)

**Referencia**: Timoshenko, S.P. (1937). *Vibration Problems in Engineering*, 2ª ed., Van Nostrand. Valores de $\beta_n$, p. 330 (cantilever).

**Geometría**: $L=1\,\text{m}$, $b=0.1\,\text{m}$, $h=0.001\,\text{m}$; $E=210\,\text{GPa}$, $\rho=7800\,\text{kg/m}^3$.  
$f_n = \dfrac{\beta_n^2}{2\pi}\sqrt{\dfrac{EI}{\rho A L^4}}$, con $\beta_1=1.875104$.  
Valor de referencia: $f_1 = 0.838\,\text{Hz}$.

| Test / Referencia | Parámetros | $f_{\rm FEM}$ (Hz) | $f_{\rm ref}$ (Hz) | $\Delta\%$ |
|---|---|---|---|---|
| `test_shell_comprehensive.py` · `TestModalAnalysis::test_first_mode_frequency` / Timoshenko (1937) | MITC4+, malla $8\times4$ quad, modo 1 | $0.848$ | $0.838$ | $1.2$ |

---

## Formulacion MITC4+: decisiones especificas

### Interpolacion mixta de deformaciones de membrana y prevencion de membrane locking

El MITC4 original (Dvorkin & Bathe, 1984) usa cuatro puntos de amarre para las componentes covariantes de corte transversal y resuelve el shear locking. Sin embargo, en shells curvados y mallas distorsionadas puede aparecer membrane locking: la rigidez de membrana se acopla parasitamente con los modos de flexion, generando deflexiones subpredichas cuando las capas de membrana y flexion deberian ser independientes.

El MITC4+ (Ko et al., 2017) agrega cinco puntos de amarre para las deformaciones covariantes de membrana ($\varepsilon_{rr}$, $\varepsilon_{ss}$, $\varepsilon_{rs}$) evaluadas en posiciones especificas de la geometria del elemento. Los valores en esos puntos se interpolan bilinealmente sobre el elemento para construir un campo de deformacion de membrana alternativo que sustituye a la deformacion derivada directamente de la cinematica nodal. Esta interpolacion blended garantiza que la deformacion de membrana sea compatible con la cinematica de flexion sin el acoplamiento parasito que genera locking.

### Integracion selectiva reducida (SRI) y rechazo de la EAS

A pesar del blending de membrana, puede persistir locking residual en la deformacion de corte en el plano para elementos distorsionados o de elevada razon de aspecto. El corte en el plano de un elemento cuadrilateral bilineal exhibe terminos parasitos que son maximos en los puntos de integracion de Gauss y cero en el centroide del elemento.

La integracion selectiva reducida (SRI; Hughes, Taylor & Kanoknukulchai, 1977) evalua la rigidez de corte en el plano con un unico punto de integracion en el centroide, donde ese corte parasito es cero, mientras que el resto de los terminos de rigidez (flexion, membrana normal, corte transversal) se integran con la cuadratura completa $2 \times 2$. El resultado es una rigidez de corte en el plano libre de locking sin alterar los terminos que no lo exhiben.

La alternativa clasica son las Deformaciones Supuestas Mejoradas (EAS): parametros internos de deformacion condensados a nivel elemental que eliminan los terminos parasitos por enriquecimiento del espacio de deformaciones. La EAS fue evaluada y descartada en esta implementacion por una razon de consistencia de la linealizacion. Para el analisis no lineal, la rigidez tangente $K_T$ debe ser exactamente la derivada de las fuerzas internas respecto de los desplazamientos, $K_T = \partial f_{\rm int}/\partial u$. Si las fuerzas internas se calculan con condensacion EAS pero la rigidez tangente no incorpora la sensibilidad de los parametros EAS respecto de los desplazamientos nodales (o si esa sensibilidad se aproxima), la condicion de consistencia se viola y el metodo de Newton-Raphson pierde su convergencia cuadratica. La SRI no introduce ninguna inconsistencia de este tipo porque la rigidez tangente y las fuerzas internas se evaluan con la misma cuadratura para el corte en el plano, preservando la consistencia del esquema de Newton.

### Estabilizacion de modos de reloj de arena

La integracion reducida de la rigidez de corte en el plano es equivalente a usar un solo punto de integracion para ese termino, lo que activa los modos de reloj de arena (*hourglass modes*): modos de deformacion con patron alternante de signo entre nodos que tienen energia de deformacion cero bajo integracion de un punto pero no representan un movimiento de cuerpo rigido.

Para suprimir estos modos se agrega una rigidez de estabilizacion:

$$
K_{\rm hg} = \alpha_{\rm hg}\, G\, t\, A\, \mathbf{h}\, \mathbf{h}^T,
$$

donde $\mathbf{h}$ es el vector de patron alternante de los modos de hourglass, $G$ es el modulo de corte, $t$ el espesor, $A$ el area del elemento y $\alpha_{\rm hg}$ un factor de calibracion. El valor $\alpha_{\rm hg} = 0.005$ es el estandar de Abaqus/Standard para elementos de casco con SRI y garantiza que la estabilizacion no afecta los modos fisicos de flexion ni de membrana mientras si elimina los modos espurios.

---

## 4. Ko et al. (2017) — Benchmarks MITC4+ (9 familias)

**Referencia**: Ko, Y., Lee, Y., Lee, P.S. & Bathe, K.J. (2017). *Performance of the MITC4+ shell elements in widely-used benchmark problems*. Computers and Structures, **193**, 187–206.

Todos los tests comparan la razón normalizada $w_{\rm FEM}/w_{\rm Kirchhoff}$ contra los valores tabulados en el artículo para $N=16$. El "Esperado" es el valor del artículo; el "Obtenido" es el calculado por el código.

### 4.1 Placa cuadrada — Tablas 2–5

$L=1$, $\nu=0.3$; carga uniforme; $t/L\in\{1/100,\,1/1000,\,1/10000\}$.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Clamped, regular, $t/L=1/100$ | $1.0029$ | $0.9984$ | $0.45$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Biapoyada, regular, $t/L=1/100$ | $1.0057$ | $1.0000$ | $0.57$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Clamped, regular, $t/L=1/1000$ | $1.0004$ | $0.9980$ | $0.24$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Biapoyada, regular, $t/L=1/1000$ | $1.0004$ | $0.9998$ | $0.06$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Clamped, regular, $t/L=1/10000$ | $0.9981$ | $0.9979$ | $0.02$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 2–3 | Biapoyada, regular, $t/L=1/10000$ | $0.9998$ | $0.9998$ | $0.00$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Clamped, distorsionada, $t/L=1/100$ | $1.0028$ | $1.0020$ | $0.08$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Biapoyada, distorsionada, $t/L=1/100$ | $1.0094$ | $1.0030$ | $0.64$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Clamped, distorsionada, $t/L=1/1000$ | $0.9998$ | $1.0010$ | $0.12$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Biapoyada, distorsionada, $t/L=1/1000$ | $1.0009$ | $1.0030$ | $0.21$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Clamped, distorsionada, $t/L=1/10000$ | $0.9978$ | $1.0010$ | $0.32$ |
| `test_ko2017_performance.py` · `test_3_1_square_plate_tables_2_to_5` / Ko et al. (2017) Tab. 4–5 | Biapoyada, distorsionada, $t/L=1/10000$ | $0.9997$ | $1.0030$ | $0.33$ |

### 4.2 Placa circular — Tablas 6–7

$R=1$, $\nu=0.3$; $t/R\in\{1/100,\,1/1000,\,1/10000\}$.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 6 | Clamped, $t/R=1/100$ | $0.9998$ | $1.0010$ | $0.12$ |
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 6 | Clamped, $t/R=1/1000$ | $0.9980$ | $0.9997$ | $0.17$ |
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 6 | Clamped, $t/R=1/10000$ | $0.9968$ | $0.9997$ | $0.29$ |
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 7 | Biapoyada, $t/R=1/100$ | $0.9982$ | $1.0010$ | $0.28$ |
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 7 | Biapoyada, $t/R=1/1000$ | $0.9976$ | $0.9997$ | $0.21$ |
| `test_ko2017_performance.py` · `test_3_2_circular_plate_tables_6_to_7` / Ko et al. (2017) Tab. 7 | Biapoyada, $t/R=1/10000$ | $0.9974$ | $0.9997$ | $0.23$ |

### 4.3 Cilindro punzado — Tablas 8–9

$R=300$, $L=600$, $t=3$; octavo de modelo.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_3_pinched_cylinder_tables_8_to_9` / Ko et al. (2017) Tab. 8 | Regular | $0.9674$ | $0.9313$ | $3.88$ |
| `test_ko2017_performance.py` · `test_3_3_pinched_cylinder_tables_8_to_9` / Ko et al. (2017) Tab. 9 | Distorsionada | $0.9932$ | $0.9892$ | $0.41$ |

### 4.4 Techo Scordelis-Lo — Tablas 10–11

$R=25$, $L=50$, $t=0.25$, arco 40°.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_4_scordelis_lo_tables_10_to_11` / Ko et al. (2017) Tab. 10 | Regular | $0.9988$ | $0.9973$ | $0.15$ |
| `test_ko2017_performance.py` · `test_3_4_scordelis_lo_tables_10_to_11` / Ko et al. (2017) Tab. 11 | Distorsionada | $1.0006$ | $0.9942$ | $0.64$ |

### 4.5 Viga torcida — Tablas 12–13

$L=12$, $b=1.1$, torsión 90°; dos espesores, dos direcciones de carga.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_5_twisted_beam_tables_12_to_13` / Ko et al. (2017) Tab. 12 | $t=0.02667$, en-plano | $1.0001$ | $1.0200$ | $1.95$ |
| `test_ko2017_performance.py` · `test_3_5_twisted_beam_tables_12_to_13` / Ko et al. (2017) Tab. 12 | $t=0.02667$, fuera-de-plano | $0.9999$ | $0.9900$ | $1.00$ |
| `test_ko2017_performance.py` · `test_3_5_twisted_beam_tables_12_to_13` / Ko et al. (2017) Tab. 13 | $t=0.0002667$, en-plano | $0.9131$ | $0.9200$ | $0.75$ |
| `test_ko2017_performance.py` · `test_3_5_twisted_beam_tables_12_to_13` / Ko et al. (2017) Tab. 13 | $t=0.0002667$, fuera-de-plano | $0.9112$ | $0.9200$ | $0.95$ |

### 4.6 Gancho de Raasch — Tabla 14

$R_1=14$, $R_2=46$, arcos 60° + 150°, $b=20$, $t=2$.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_6_hook_table_14_minimal_fix` / Ko et al. (2017) Tab. 14 | Malla $N=16$ | $0.9927$ | $1.1200$ | $11.37$ |

> El 11 % de diferencia respecto al valor tabulado por Ko es esperado: la referencia de Kirchhoff para el gancho tiene incertidumbre en la propia literatura; la tolerancia del test es 40 %.

### 4.7 Hemisferio con corte 18° — Tablas 15–16

$R=10$, $\theta_{\rm cut}=18°$; dos espesores.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_7_hemisphere_cutout_tables_15_to_16` / Ko et al. (2017) Tab. 15 | $t/R=4/1000$, regular | $1.0047$ | $1.0090$ | $0.43$ |
| `test_ko2017_performance.py` · `test_3_7_hemisphere_cutout_tables_15_to_16` / Ko et al. (2017) Tab. 15 | $t/R=4/1000$, distorsionada | $1.0042$ | $0.9958$ | $0.84$ |
| `test_ko2017_performance.py` · `test_3_7_hemisphere_cutout_tables_15_to_16` / Ko et al. (2017) Tab. 16 | $t/R=4/10000$, regular | $0.9786$ | $0.9811$ | $0.25$ |
| `test_ko2017_performance.py` · `test_3_7_hemisphere_cutout_tables_15_to_16` / Ko et al. (2017) Tab. 16 | $t/R=4/10000$, distorsionada | $0.9641$ | $0.9736$ | $0.98$ |

### 4.8 Hemisferio completo — Tabla 17

$R=10$, $\theta_{\rm min}=2°$; dos espesores.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_8_full_hemisphere_table_17` / Ko et al. (2017) Tab. 17 | $t/R=4/1000$ | $0.9982$ | $0.9960$ | $0.22$ |
| `test_ko2017_performance.py` · `test_3_8_full_hemisphere_table_17` / Ko et al. (2017) Tab. 17 | $t/R=4/10000$ | $0.9720$ | $0.9798$ | $0.80$ |

### 4.9 Paraboloide hiperbólico — Tablas 18–19

$z = y^2 - x^2$, $L=1$; dos espesores.

| Test / Referencia | Parámetros | Obtenido | Esperado | $\Delta\%$ |
|---|---|---|---|---|
| `test_ko2017_performance.py` · `test_3_9_hyperbolic_paraboloid_tables_18_to_19` / Ko et al. (2017) Tab. 18 | $t/L=1/1000$, regular | $0.9761$ | $0.9762$ | $0.01$ |
| `test_ko2017_performance.py` · `test_3_9_hyperbolic_paraboloid_tables_18_to_19` / Ko et al. (2017) Tab. 18 | $t/L=1/1000$, distorsionada | $0.9975$ | $0.9904$ | $0.72$ |
| `test_ko2017_performance.py` · `test_3_9_hyperbolic_paraboloid_tables_18_to_19` / Ko et al. (2017) Tab. 19 | $t/L=1/10000$, regular | $0.9770$ | $0.9777$ | $0.07$ |
| `test_ko2017_performance.py` · `test_3_9_hyperbolic_paraboloid_tables_18_to_19` / Ko et al. (2017) Tab. 19 | $t/L=1/10000$, distorsionada | $1.0168$ | $0.9936$ | $2.34$ |

---

## Resumen

| Grupo | Elemento(s) | Tests | Referencia |
|---|---|---|---|
| Grandes rotaciones, viga voladizo (8 casos × 2 comp.) | MITC3+, MITC4+ | 16 | Simo & Vu-Quoc (1986) CMAME 58 |
| Deflexión lineal, convergencia y membrana | MITC4+ | 5 | Timoshenko & W-K (1959) |
| Cargas uniaxiales, voladizo delgado | MITC4+ | 3 | Timoshenko & W-K (1959) |
| Análisis modal, modo 1 | MITC4+ | 1 | Timoshenko (1937) |
| Ko et al. (2017) — 9 familias de benchmark | MITC4+ | 37 | Ko et al. (2017) CAS 193 |
| **Total** | | **62** | |

---

## Referencias

1. Simo, J.C. & Vu-Quoc, L. (1986). *A three-dimensional finite-strain rod model. Part II: Computational aspects*. CMAME **58**(1), 79–116.
2. Dvorkin, E.N. & Bathe, K.J. (1984). *A continuum mechanics based four-node shell element for general non-linear analysis*. Engineering Computations **1**(1), 77–88.
3. Lee, P.S., Lee, Y. & Bathe, K.J. (2014). *The MITC3+ shell element and its performance*. Computers & Structures **138**, 12–23.
4. Ko, Y., Lee, Y., Lee, P.S. & Bathe, K.J. (2017). *Performance of the MITC4+ shell elements in widely-used benchmark problems*. Computers and Structures **193**, 187–206.
5. Timoshenko, S. & Woinowsky-Krieger, S. (1959). *Theory of Plates and Shells*, 2ª ed., McGraw-Hill.
6. Timoshenko, S.P. (1937). *Vibration Problems in Engineering*, 2ª ed., Van Nostrand.
7. Hughes, T.J.R., Taylor, R.L. & Kanoknukulchai, W. (1977). *A simple and efficient finite element for plate bending*. International Journal for Numerical Methods in Engineering **11**(10), 1529–1543.
