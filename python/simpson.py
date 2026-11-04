"""
Análisis Numérico · Unidad 3 · Sesión 13
Integración numérica: regla de Simpson 1/3 compuesta. Energía consumida por un centro de datos.

Potencia de los servidores durante 6 horas de carga crítica (kW):
    P(t) = 10 + 5t·sen(pi·t/6) + e^(0.2·t)
Energía total (kWh):
    E = integral de 0 a 6 de P(t) dt

Fernando Istaña Flores · Universidad Peruana Unión
"""

import math
import sys

import numpy as np


# --------------------------------------------------------------------------
# Modelo
# --------------------------------------------------------------------------
def P(t):
    """Potencia en kW. Acepta un número o un arreglo de NumPy (el seno trabaja en radianes)."""
    return 10 + 5 * t * np.sin(np.pi * t / 6) + np.exp(0.2 * t)


def P4(t):
    """Cuarta derivada de P(t)."""
    a = np.pi / 6
    return 5 * (a**4 * t * np.sin(a * t) - 4 * a**3 * np.cos(a * t)) + 0.2**4 * np.exp(0.2 * t)


def energia_exacta(a=0.0, b=6.0):
    """Primitiva: 10t + 5[sen(at)/a² − t·cos(at)/a] + e^(0.2t)/0.2, con a = pi/6."""
    w = math.pi / 6
    F = lambda t: 10 * t + 5 * (math.sin(w * t) / w**2 - t * math.cos(w * t) / w) + math.exp(0.2 * t) / 0.2
    return F(b) - F(a)


M_GUIA = 0.45   # cota de |P''''| que da el enunciado


# --------------------------------------------------------------------------
# Módulo de cálculo
# --------------------------------------------------------------------------
class NImparError(ValueError):
    """n impar: Simpson 1/3 necesita que las parábolas tomen los subintervalos de dos en dos."""


def validar_n(n):
    if not isinstance(n, (int, np.integer)) or isinstance(n, bool):
        raise NImparError("n debe ser un número entero (se recibió %r)." % (n,))
    if n < 2 or n % 2 != 0:
        raise NImparError("n = %d no es válido: Simpson 1/3 necesita un número entero PAR de subintervalos (n ≥ 2)." % n)


def coeficientes(n):
    """Vector de coeficientes 1, 4, 2, 4, …, 2, 4, 1 de longitud n + 1."""
    c = np.full(n + 1, 2.0)
    c[1::2] = 4.0
    c[0] = c[-1] = 1.0
    return c


def simpson_13(f, a, b, n):
    """
    Regla de Simpson 1/3 compuesta, vectorizada:

        h = (b − a) / n,   x = a, a + h, …, b
        I ≈ (h/3) · ( c · f(x) ),   c = [1, 4, 2, 4, …, 2, 4, 1]
    """
    validar_n(n)
    x = np.linspace(a, b, n + 1)
    h = (b - a) / n
    return h / 3 * np.dot(coeficientes(n), f(x))


def cota_error(a, b, n, M):
    """|Et| ≤ (b − a)·h⁴/180 · M."""
    h = (b - a) / n
    return (b - a) * h**4 / 180 * M


# --------------------------------------------------------------------------
# Entrada y validación
# --------------------------------------------------------------------------
def pedir(mensaje):
    """input() que, si los datos llegan redirigidos (no desde el teclado), los repite en pantalla."""
    valor = input(mensaje)
    if not sys.stdin.isatty():
        print(valor)
    return valor


def leer_datos():
    """Pide a, b y n. Si n no es un entero par, detiene la ejecución con un mensaje de error."""
    try:
        a = float(pedir("  a (inicio del intervalo, h) = "))
        b = float(pedir("  b (fin del intervalo, h)    = "))
        texto_n = pedir("  n (subintervalos, par)      = ").strip()
    except EOFError:
        print("  (sin entrada: se usan los datos de la guía a = 0, b = 6, n = 6)")
        return 0.0, 6.0, 6
    print()
    try:
        n = int(texto_n)
    except ValueError:
        raise NImparError("n debe ser un número entero (se recibió %r)." % texto_n)
    if b <= a:
        raise ValueError("b debe ser mayor que a.")
    validar_n(n)
    return a, b, n


if __name__ == "__main__":
    print("=" * 74)
    print("ENTRADA DE DATOS")
    print("=" * 74)
    try:
        a, b, n = leer_datos()
    except ValueError as error:
        print()
        print("  ERROR: %s" % error)
        print("  Ejecución detenida.")
        sys.exit(1)
    print("  a = %g, b = %g, n = %d  (n par: válido)" % (a, b, n))

    exacta = energia_exacta(a, b)
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    c = coeficientes(n)
    y = P(x)

    print()
    print("=" * 74)
    print("PARTE 1 · SIMPSON 1/3 COMPUESTO CON n = %d" % n)
    print("=" * 74)
    print("  h = (b − a)/n = (%g − %g)/%d = %g h" % (b, a, n, h))
    print()
    print("   i     t_i    5t·sen(pi·t/6)    e^(0.2t)         P(t_i)   coef   coef·P")
    for i in range(n + 1):
        s = 5 * x[i] * math.sin(math.pi * x[i] / 6)
        print("  %2d  %6.2f   %14.6f   %10.6f   %12.6f   %4d   %10.6f" % (i, x[i], s + 0.0, math.exp(0.2 * x[i]), y[i], c[i], c[i] * y[i]))
    impares = y[1:-1:2].sum()
    pares = y[2:-1:2].sum()
    I = simpson_13(P, a, b, n)
    print()
    print("  Vector de coeficientes c = %s" % np.array2string(c.astype(int), separator=", "))
    print("  Suma de impares = %.6f     Suma de pares = %.6f" % (impares, pares))
    print("  I = (h/3)·[P(t0) + 4·%.6f + 2·%.6f + P(tn)]" % (impares, pares))
    print("    = (%g/3)·%.6f = %.6f kWh" % (h, np.dot(c, y), I))

    print()
    print("=" * 74)
    print("COTA DEL ERROR TEÓRICO")
    print("=" * 74)
    tt = np.linspace(a, b, 600001)
    M_real = float(np.max(np.abs(P4(tt))))
    t_max = float(tt[np.argmax(np.abs(P4(tt)))])
    e_real = abs(exacta - I)
    print("  |Et| ≤ (b − a)·h⁴/180 · M")
    print("  Con M = %.2f (guía):     |Et| ≤ %g·%g⁴/180·%.2f = %.6f kWh" % (M_GUIA, b - a, h, M_GUIA, cota_error(a, b, n, M_GUIA)))
    print("  Valor exacto: E = 60 + 180/pi + (e^1.2 − 1)/0.2 = %.6f kWh" % exacta)
    print("  Error real = |%.6f − %.6f| = %.6f kWh" % (exacta, I, e_real))
    if e_real > cota_error(a, b, n, M_GUIA):
        print("  OJO: el error real supera esa cota, así que M = %.2f no es el máximo de |P''''|." % M_GUIA)
    print("  Máximo real de |P''''| en [%g, %g]: %.6f (en t = %.4f)" % (a, b, M_real, t_max))
    print("  Con M = %.4f:          |Et| ≤ %.6f kWh  (≥ error real: correcto)" % (M_real, cota_error(a, b, n, M_real)))

    print()
    print("=" * 74)
    print("ANÁLISIS DE CONVERGENCIA")
    print("=" * 74)
    print("  %4s %8s %14s %16s %18s %14s" % ("n", "h", "I_Simpson", "|Et| (M=0.45)", "|Et| (M=%.4f)" % M_real, "error real"))
    filas = []
    for k in (2, 4, 6, 12, 24):
        v = simpson_13(P, a, b, k)
        filas.append((k, (b - a) / k, v, abs(exacta - v)))
        print("  %4d %8.4f %14.6f %16.6e %18.6e %14.6e" % (k, (b - a) / k, v, cota_error(a, b, k, M_GUIA), cota_error(a, b, k, M_real), abs(exacta - v)))
    print()
    print("  Al duplicar n (12 → 24) el error real se divide entre %.1f ≈ 2⁴ = 16: orden h⁴." % (filas[3][3] / filas[4][3]))

    print()
    print("=" * 74)
    print("PRUEBA DE LA VALIDACIÓN (n = 5)")
    print("=" * 74)
    try:
        simpson_13(P, a, b, 5)
    except NImparError as error:
        print("  ERROR: %s" % error)

    # ----------------------------------------------------------------------
    # Gráficas: P(t) con las parábolas de n = 6, y la convergencia
    # ----------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker

        t = np.linspace(a, b, 600)
        fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=200)
        for j in range(0, n, 2):
            xs, ys = x[j:j + 3], y[j:j + 3]
            q = np.polyfit(xs, ys, 2)
            tq = np.linspace(xs[0], xs[-1], 80)
            ax.fill_between(tq, np.polyval(q, tq), color="#c8890a", alpha=0.18, linewidth=0,
                            label="área de las parábolas" if j == 0 else None)
            ax.plot(tq, np.polyval(q, tq), "--", color="#c8890a", linewidth=1.3, label="parábolas (n = %d)" % n if j == 0 else None)
            ax.vlines(xs[0], 0, ys[0], color="#c8890a", linewidth=0.8)
        ax.vlines(x[-1], 0, y[-1], color="#c8890a", linewidth=0.8)
        ax.plot(t, P(t), color="#14366b", linewidth=2, label="P(t) = 10 + 5t·sen(πt/6) + e^(0.2t)")
        ax.scatter(x[0::2], y[0::2], color="#101c33", zorder=3, s=34, label="nodos pares")
        ax.scatter(x[1::2], y[1::2], color="#c8890a", zorder=3, s=40, label="nodos impares (coef. 4)")
        ax.set_xlabel("Tiempo t (h)")
        ax.set_ylabel("Potencia P(t) (kW)")
        ax.set_title("Simpson 1/3 compuesto, n = %d: %.4f kWh  ·  exacto: %.4f kWh" % (n, I, exacta))
        ax.set_ylim(0, max(y) * 1.15)
        ax.grid(color="#dfe5ef", linewidth=0.8)
        ax.set_axisbelow(True)
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        ax.legend(frameon=False, fontsize=7.5, loc="lower center", ncol=2)
        fig.tight_layout()
        fig.savefig("simpson_s13.png")

        fig, ax = plt.subplots(figsize=(7.4, 3.8), dpi=200)
        ns = [2, 4, 6, 8, 12, 16, 24, 48]
        errs = [abs(exacta - simpson_13(P, a, b, k)) for k in ns]
        ax.loglog(ns, errs, "o-", color="#14366b", linewidth=2, label="error real")
        ax.loglog(ns, [cota_error(a, b, k, M_real) for k in ns], "--", color="#c8890a", label="cota con M = %.3f (pendiente −4)" % M_real)
        ax.loglog(ns, [cota_error(a, b, k, M_GUIA) for k in ns], ":", color="#a32330", label="cota con M = 0.45 (guía)")
        ax.set_xticks(ns)
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("Número de subintervalos n")
        ax.set_ylabel("Error (kWh)")
        ax.set_title("Convergencia de Simpson 1/3 compuesto")
        ax.grid(color="#dfe5ef", linewidth=0.8, which="both")
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5)
        fig.tight_layout()
        fig.savefig("convergencia_s13.png")
        print()
        print("  Gráficas guardadas en simpson_s13.png y convergencia_s13.png")
    except ImportError:
        print()
        print("  (matplotlib no está instalado: se omitieron las gráficas)")
