# Nomenclatura y convencion comun

Este documento centraliza la notacion usada en los cinco informes de validacion teorica de AeroElast. Cada informe mantiene su propia tabla local orientada al lector aislado; esta referencia resuelve ambiguedades entre informes y hace explicitas las convenciones que se comparten o que se especializan en cada contexto.

## Mecanica estructural (informes 01–04)

| Simbolo | Definicion | Informes | Observacion |
| --- | --- | --- | --- |
| $u$ | Vector de desplazamientos estructurales nodales | 01, 02 | En 03 denota desplazamientos en el marco rotante; en 04 se usa $u_e$ para el desplazamiento elastico global. |
| $\dot{u}$ | Vector de velocidades estructurales nodales | 01–04 | |
| $\ddot{u}$ | Vector de aceleraciones estructurales nodales | 01–04 | |
| $u_e$ | Desplazamiento elastico en coordenadas globales | 04 | Reservado para la formulacion inertial. |
| $M$ | Matriz de masa estructural | 01–04 | En 01 es masa consistente; en la ruta FSI de alta performance (02, 03) es masa lumped por suma de filas. |
| $K$ | Matriz de rigidez elastica lineal | 01–04 | En 04 se escribe $K(\theta)$ para indicar dependencia de la geometria rotada. |
| $C$ | Matriz de amortiguamiento de Rayleigh | 02–04 | Ausente en el problema modal (01). En 04 se escribe $C(\theta)$. |
| $K_G$ | Matriz de rigidez geometrica por prestress | 02, 03 | Se ensambla solo con la parte tensil del estado de membrana ($\sigma^+$). |
| $K_{SP}$ | Matriz de spin softening | 03 | Aproximacion lineal de la dependencia centrifuga con el desplazamiento transversal. |
| $G_{cor}$ | Matriz giroscopica de Coriolis | 03 | $G_{cor} = 2 M \widetilde{\Omega}$. |
| $T$ | Operador booleano de reduccion al subespacio libre | 01–04 | Elimina grados de libertad con condicion de Dirichlet homogenea. |
| $u_r$ | Desplazamiento reducido (subespacio libre) | 01–04 | $u = T u_r$. |
| $K_r,\, M_r,\, C_r$ | Matrices reducidas al subespacio libre | 01–04 | $(\cdot)_r = T^T (\cdot) T$. |
| $f_{fsi}$ | Fuerza nodal recibida desde el participante fluido | 02–04 | |
| $f_g$ | Fuerza gravitacional nodal | 02–04 | Puede escribirse en marco rotante ($f_g^{loc}$) o global ($f_g^{glob}$). |
| $f_{cf}$ | Fuerza centrifuga base | 03 | Evaluada en posicion no deformada cuando $K_{SP}$ esta activo. |
| $f_{euler}$ | Fuerza de Euler por aceleracion angular no nula | 03 | |
| $\eta_m,\, \eta_k$ | Coeficientes de amortiguamiento de Rayleigh | 02–04 | $C = \eta_m M + \eta_k K$. |
| $\Delta t$ | Ancho de ventana temporal | 02–04 | |
| $\beta,\, \gamma$ | Parametros del esquema de Newmark | 02–04 | $\beta = 0.25$, $\gamma = 0.5$ (aceleracion promedio constante). |
| $a_0 \dots a_5$ | Coeficientes derivados del esquema de Newmark | 02–04 | |
| $K_{\mathrm{eff}}$ | Rigidez efectiva del sistema de Newmark | 02–04 | Incluye contribuciones inerciales, de amortiguamiento y geometricas. |

## Analisis modal (informe 01)

| Simbolo | Definicion |
| --- | --- |
| $\phi_i$ | Forma modal del modo $i$ |
| $\lambda_i$ | Autovalor modal ($\lambda_i = \omega_i^2$) |
| $\omega_i$ | Frecuencia angular natural del modo $i$ |
| $f_i$ | Frecuencia natural del modo $i$ en Hz |
| $m_i$ | Masa modal generalizada del modo $i$ |
| $\Gamma_i$ | Factor de participacion modal del modo $i$ |
| $m_{\mathrm{eff},i}$ | Masa efectiva modal del modo $i$ |

## Cinematica y dinamica angular del rotor (informes 03–04)

| Simbolo | Definicion | Informes |
| --- | --- | --- |
| $\theta$ | Angulo acumulado del rotor | 03, 04 |
| $\Delta \theta^n$ | Incremento angular de la ventana $n$ | 03, 04 |
| $\omega$ | Velocidad angular del rotor | 03, 04 |
| $\omega_{step}$ | Velocidad angular efectiva fija dentro de una ventana | 03, 04 |
| $\alpha$ | Aceleracion angular del rotor | 03, 04 |
| $\hat{n}$ | Eje unitario de rotacion | 03, 04 |
| $I$ | Momento de inercia total respecto del eje de rotacion | 03, 04 |
| $R(\theta)$ | Matriz de rotacion rigida del rotor al angulo $\theta$ | 03, 04 |
| $\widetilde{\Omega}$ | Operador antisimetrico de la velocidad angular | 03 |
| $\tau_{aero}$ | Torque aerodinamico respecto del eje de rotacion | 03, 04 |
| $\tau_g$ | Torque gravitatorio respecto del eje de rotacion | 03, 04 |
| $\tau_{shaft}$ | Torque de eje (control o generador) | 03, 04 |
| $\tau_{driving}$ | Torque motriz: $\tau_{aero} + \tau_g$ | 03, 04 |
| $a_{ref}$ | Aceleracion de la referencia rigidamente rotada | 04 |
| $F_{ref}$ | Carga equivalente de referencia: $-M a_{ref}$ | 04 |
| $\hat{x}$ | Posicion de la referencia rigidamente rotada | 04 |
| $x$ | Posicion total del punto material | 04 |
| $X_0$ | Configuracion de referencia no rotada | 04 |

## Mecanica de la rigidez geometrica (informes 02–03)

| Simbolo | Definicion |
| --- | --- |
| $B_G$ | Matriz cinematica geometrica (gradientes de desplazamiento fuera del plano) |
| $\sigma^+$ | Parte tensil del tensor de tensiones de membrana (direcciones principales positivas) |
| $N_{\alpha\beta}^+$ | Resultante de membrana tensil por unidad de longitud |
| $\mathbf{P}$ | Matriz de vectores principales del tensor de tensiones de membrana |
| $N_1,\, N_2$ | Valores principales de las resultantes de membrana |

## Aerodinamica BEM (informe 05)

| Simbolo | Definicion |
| --- | --- |
| $r_k$ | Posicion radial de la estacion $k$ |
| $r_k^{def}$ | Radio deformado de la estacion $k$ |
| $c_k$ | Cuerda de la estacion $k$ |
| $\vartheta_k$ | Twist aerodinamico total de la estacion $k$ |
| $\Delta \vartheta_k$ | Incremento de twist elastico |
| $a_k$ | Factor de induccion axial |
| $a'_k$ | Factor de induccion tangencial |
| $W_k$ | Velocidad relativa en la estacion $k$ |
| $\phi_k$ | Angulo de inflow |
| $\alpha_k$ | Angulo de ataque de la estacion $k$ |
| $\theta_{twist,k}$ | Twist geometrico de la estacion $k$ |
| $\theta_{pitch}$ | Angulo de paso (pitch) global de la pala |
| $N_{p,k}$ | Carga normal por unidad de longitud en la estacion $k$ |
| $T_{p,k}$ | Carga tangencial por unidad de longitud en la estacion $k$ |
| $M_{p,k}$ | Momento de pitching por unidad de longitud en la estacion $k$ |
| $F_i$ | Fuerza nodal proyectada sobre el nodo $i$ de la malla estructural |
| $\psi$ | Azimut del rotor |
| $\Delta \psi$ | Incremento de azimut por ventana convergida |
| $V_\infty$ | Velocidad libre del flujo incidente |
| $\rho$ | Densidad del fluido |
| $C_{l,k},\, C_{d,k},\, C_{m,k}$ | Coeficientes de sustentacion, resistencia y momento de la estacion $k$ |
| $\mathcal{I}_k$ | Conjunto de nodos de la malla estructural asociados a la franja $k$ |
| $\hat{e}_s$ | Direccion unitaria del vano de la pala |
| $\hat{c}_k$ | Direccion unitaria de la cuerda de la estacion $k$ |
| $x_{AC,k}$ | Posicion del centro aerodinamico de la franja $k$ |

## Simbolos sobreargados: advertencias explicitas

Los tres simbolos que se usan con significados distintos segun el contexto son:

### $u$ — desplazamiento estructural

| Informe | Significado |
| --- | --- |
| 01, 02 | Vector global de desplazamientos nodales. |
| 03 | Desplazamiento elastico expresado en el marco rotante. |
| 04 | Se evita deliberadamente; se usa $u_e$ para el desplazamiento elastico global. |

### $M$ — masa

| Informe | Tipo de masa |
| --- | --- |
| 01 | Masa consistente (acoplamiento inercial completo). |
| 02, 03, 04 | Masa lumped por suma de filas en la ruta FSI de alta performance. |

Las frecuencias obtenidas con masa consistente (modal) no son directamente comparables con un analisis de Newmark que use masa lumped, aunque las diferencias en los modos bajos son generalmente pequenas. Cuando se comparan resultados modales con transitorios FSI, conviene verificar que ambos usen el mismo tipo de masa.

### $\alpha$ — aceleracion y angulo de ataque

| Informe | Significado |
| --- | --- |
| 03, 04 | Aceleracion angular del rotor ($\mathrm{rad/s}^2$). |
| 05 | Angulo de ataque de la estacion $k$ ($\alpha_k$, en radianes). |

En los informes de rotor no aparece un angulo de ataque como variable primaria. En el informe BEM no aparece la aceleracion angular como simbolo autonomo. El contexto elimina la ambiguedad, pero conviene no mezclar ecuaciones de ambos informes sin aclaracion previa.

## Operadores y supraindices comunes

| Notacion | Significado |
| --- | --- |
| $(\cdot)^{loc}$ | Cantidad expresada en el marco rotante local |
| $(\cdot)^{glob}$ | Cantidad expresada en el marco global inercial |
| $(\cdot)^+$ | Parte positiva (filtro tensil): valores propios negativos llevados a cero |
| $(\cdot)_r$ | Cantidad reducida al subespacio de grados de libertad libres |
| $(\cdot)_n$ | Evaluacion en el instante de tiempo $n$ |
| $\widetilde{v}$ | Operador antisimetrico (producto vectorial) asociado al vector $v$ |
| $\hat{v}$ | Vector unitario en la direccion de $v$ |
