# A3 — Statistical stability (yaw=0°, T_w=10s, corotational K=6)

Window: [40.0, 100.0] s, sub-window length = 10.0 s.
Reports mean ± IC95 across non-overlapping sub-windows. The corotational
rows are the evidence used by the standalone corotational report; any
additional solver rows are kept only for cross-solver diagnostics.

**Caveat para `std` y `pp`**: estas columnas reportan IC95 de la dispersión
*dentro de una ventana de 10 s*. Eso es lo correcto si interesa saber si la
amplitud típica observada en un fragmento de 10 s es robusta. La `std` global
sobre la ventana completa es un objeto distinto: no se mide bien por
promedio de stds locales, porque el promedio de stds locales subestima la std
global cuando la señal contiene escalas de período mayor a 10 s. La columna
`mean` no sufre este sesgo.

| Signal | Stat | Unit | Corot mean ± IC95 | Inert mean ± IC95 | Δ (inert-corot) | Overlap? |
|---|---|---|---|---|---:|:---:|
| edge | mean | m | -0.5071 ± 0.1914 | -0.6381 ± 0.3471 | -0.1310 | sí |
| edge | std | m | 1.0662 ± 0.0635 | 0.8570 ± 0.5661 | -0.2092 | sí |
| flap | mean | m | 12.8072 ± 0.0476 | 12.8408 ± 0.1889 | +0.0335 | sí |
| flap | pp | m | 2.0624 ± 0.4247 | 1.9950 ± 1.6212 | -0.0674 | sí |
| flap | std | m | 0.4894 ± 0.1000 | 0.4761 ± 0.3820 | -0.0133 | sí |
| power | mean | MW | 14.7053 ± 0.0247 | 14.6840 ± 0.0299 | -0.0213 | sí |
| power | pp | MW | 0.6747 ± 0.0095 | 0.5965 ± 0.3905 | -0.0783 | sí |
| thrust | mean | MN | 2.0098 ± 0.0058 | 2.0041 ± 0.0068 | -0.0057 | sí |
| torque | mean | MN·m | 18.5994 ± 0.0312 | 18.5725 ± 0.0379 | -0.0269 | sí |
