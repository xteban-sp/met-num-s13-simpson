🔗 **App publicada:** https://xteban-sp.github.io/met-num-s13-simpson/

# Integración numérica · Sesión 13

Regla de Simpson 1/3 compuesta para estimar el consumo energético de un centro de datos, en web y en Python.

Universidad Peruana Unión · Ingeniería de Sistemas · Fernando Istaña Flores
Curso: Análisis Numérico / Métodos Numéricos · Docente: Jorge Luis Manrique Plasencia

---

## Qué hay aquí

| Archivo | Contenido |
|---|---|
| `index.html` | La app web. Un solo archivo, sin dependencias: se abre con doble clic o desde GitHub Pages |
| `python/simpson.py` | Pide a, b y n, valida con una excepción que n sea entero par, `simpson_13(f, a, b, n)` vectorizada con NumPy, cota del error y tabla de convergencia |
| `python/salida_simpson.txt` | Salida con a = 0, b = 6, n = 6 |
| `python/salida_n_impar.txt` | Salida con n = 5: la validación detiene la ejecución |
| `python/simpson_s13.png` | P(t) con las tres parábolas de n = 6 |
| `python/convergencia_s13.png` | Error real y cotas contra n (log-log, pendiente −4) |
| `informes/` | El informe de la GAA de la Sesión 13 en Word y en PDF |

## La app web

- **Cualquier integral**: f(t) en t o en x (radianes), intervalo [a, b] y n par; si n es impar lo rechaza y lo explica.
- **Tabla** con los coeficientes 1, 4, 2, …, 4, 1, sumas de impares y pares, y la fórmula con los números reemplazados.
- **Gráfica** con las parábolas de cada par de subintervalos sombreadas.
- **Convergencia**: valor, cota con la M del enunciado, cota con la M calculada, error real y orden observado.
- **Cota del error** con la M dada y con el máximo real de |f⁽⁴⁾|, y el n mínimo para una tolerancia.

## Ejecutar la versión en Python

```bash
cd python
python simpson.py
```

Necesita `numpy`, y `matplotlib` para las gráficas.

## Caso resuelto — energía del datacenter

P(t) = 10 + 5t·sen(πt/6) + e^(0.2t) kW durante 6 h; E = ∫₀⁶ P(t) dt.

| i | tᵢ | P(tᵢ) | coef. |
|---|---|---|---|
| 0 | 0 | 11.000000 | 1 |
| 1 | 1 | 13.721403 | 4 |
| 2 | 2 | 20.152079 | 2 |
| 3 | 3 | 26.822119 | 4 |
| 4 | 4 | 29.546049 | 2 |
| 5 | 5 | 25.218282 | 4 |
| 6 | 6 | 13.320117 | 1 |

```
h = 1     I ≈ (1/3)·[11 + 4(65.761804) + 2(49.698128) + 13.320117] = 128.9212 kWh
exacto    E = 60 + 180/π + (e^1.2 − 1)/0.2 = 128.8964 kWh      error real 0.0248 kWh
```

**Sobre la cota de error:** con la M ≈ 0.45 de la guía sale |Eₜ| ≤ 0.015 kWh, que es *menor* que el error real
(0.0248), así que esa M no puede ser el máximo de |P⁽⁴⁾|. El máximo real en [0, 6] es M ≈ 3.4332 (en t ≈ 5.07),
y con él |Eₜ| ≤ 0.1144 kWh, que sí cubre el error.

| n | h | I_Simpson (kWh) | cota (M = 0.45) | cota (M = 3.4332) | error real |
|---|---|---|---|---|---|
| 2 | 3 | 131.608592 | 1.215 | 9.2697 | 2.7122 |
| 4 | 1.5 | 129.027508 | 0.0759 | 0.5794 | 0.1311 |
| 6 | 1 | 128.921195 | 0.015 | 0.1144 | 0.0248 |
| 12 | 0.5 | 128.897878 | 9.4×10⁻⁴ | 7.2×10⁻³ | 1.5×10⁻³ |
| 24 | 0.25 | 128.896458 | 5.9×10⁻⁵ | 4.5×10⁻⁴ | 9.4×10⁻⁵ |
