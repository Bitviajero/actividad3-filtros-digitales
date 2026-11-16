"""
Actividad Formativa 3 - Implementación y evaluación de filtros digitales
=====================================================================
Curso: Análisis de Fourier (Unidad 3)
Perfil: Ingeniería en Desarrollo de Software

Objetivo:
  Diseñar e implementar filtros digitales pasa bajos, pasa altos y pasa
  banda (IIR Butterworth, IIR Chebyshev y FIR con ventana) y evaluar su
  efectividad sobre una señal con ruido.

Estructura (como pide la actividad):
  1. Definición de la señal de entrada
  2. Diseño del filtro
  3. Aplicación del filtro a la señal
  4. Visualización de los resultados

Requisitos:
  pip install numpy scipy matplotlib

Ejecución:
  python filtros_digitales.py

Nota:
  Se usa filtfilt (fase cero) para que la comparación visual sea limpia.
  Para aplicaciones en tiempo real se usaría lfilter (causal, con retardo).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, cheby1, firwin, freqz, filtfilt

# ---------------------------------------------------------------------------
# Configuración general
# ---------------------------------------------------------------------------
FS = 1000.0          # Frecuencia de muestreo (Hz)
T = 2.0              # Duración de la señal (s)
N = int(FS * T)      # Número de muestras
t = np.linspace(0, T, N, endpoint=False)

# Frecuencias de las componentes útiles (bien separadas para ver el efecto)
F_LOW = 50.0         # Componente "baja"
F_MID = 150.0        # Componente "media"
F_HIGH = 300.0       # Componente "alta"

RUIDO_AMP = 0.8      # Amplitud del ruido blanco

# Carpeta para guardar figuras (útil para el repositorio)
CARPETA_FIG = "figuras"
os.makedirs(CARPETA_FIG, exist_ok=True)

plt.rcParams.update({"figure.figsize": (10, 6), "font.size": 10})


# ---------------------------------------------------------------------------
# Funciones auxiliares
# ---------------------------------------------------------------------------
def espectro(senal, fs=FS):
    """Espectro de amplitud unilateral de una señal.

    Devuelve (frecuencias, magnitud_normalizada).
    """
    n = len(senal)
    X = np.fft.fft(senal)
    freqs = np.fft.fftfreq(n, 1 / fs)
    mascara = freqs >= 0
    return freqs[mascara], np.abs(X[mascara]) / n


def plot_respuesta(b, a, titulo, ax=None):
    """Grafica la respuesta en frecuencia (magnitud en dB)."""
    if ax is None:
        _, ax = plt.subplots()
    w, h = freqz(b, a, worN=2048, fs=FS)
    ax.plot(w, 20 * np.log10(np.abs(h) + 1e-12))
    ax.set_xlabel("Frecuencia (Hz)")
    ax.set_ylabel("Magnitud (dB)")
    ax.set_title(titulo, fontsize=9)
    ax.set_ylim(-80, 5)
    ax.grid(True, alpha=0.3)
    ax.axvline(F_LOW, color="green", ls=":", alpha=0.6)
    ax.axvline(F_MID, color="orange", ls=":", alpha=0.6)
    ax.axvline(F_HIGH, color="red", ls=":", alpha=0.6)
    return ax


def guardar_figura(nombre):
    """Guarda la figura actual en la carpeta de figuras."""
    ruta = os.path.join(CARPETA_FIG, nombre)
    plt.savefig(ruta, dpi=120, bbox_inches="tight")
    print(f"  Figura guardada: {ruta}")


def atenuacion_db(b, a, freqs):
    """Atenuación (dB) del filtro en las frecuencias dadas."""
    _, h = freqz(b, a, worN=freqs, fs=FS)
    return 20 * np.log10(np.abs(h) + 1e-12)


# ===========================================================================
# 1. DEFINICIÓN DE LA SEÑAL DE ENTRADA
# ===========================================================================
def generar_senal():
    """Crea la señal compuesta (útil) y la señal con ruido."""
    senal_util = (
        np.sin(2 * np.pi * F_LOW * t)
        + 0.7 * np.sin(2 * np.pi * F_MID * t)
        + 0.5 * np.sin(2 * np.pi * F_HIGH * t)
    )
    ruido = RUIDO_AMP * np.random.randn(N)
    senal_ruidosa = senal_util + ruido
    return senal_util, ruido, senal_ruidosa


# ===========================================================================
# 2. DISEÑO DE FILTROS
# ===========================================================================
def disenar_filtros():
    """Diseña filtros FIR e IIR para PB, PA y BP.

    Devuelve dict: filtros[tipo][metodo] = (b, a, descripcion)
    """
    orden_iir = 6
    ripple_cheby = 1.0      # dB de rizado
    numtaps_fir = 201       # coeficientes FIR

    fc_lp = 80.0            # pasa bajos
    fc_hp = 220.0           # pasa altos
    fc_bp = [100.0, 200.0]  # pasa banda

    filtros = {"pasa_bajos": {}, "pasa_altos": {}, "pasa_banda": {}}

    # --- Pasa bajos ---
    filtros["pasa_bajos"]["butterworth"] = (
        *butter(orden_iir, fc_lp, btype="low", fs=FS),
        f"Butterworth LP (orden {orden_iir}, fc={fc_lp} Hz)",
    )
    filtros["pasa_bajos"]["chebyshev"] = (
        *cheby1(orden_iir, ripple_cheby, fc_lp, btype="low", fs=FS),
        f"Chebyshev I LP (orden {orden_iir}, rp={ripple_cheby} dB)",
    )
    filtros["pasa_bajos"]["fir"] = (
        firwin(numtaps_fir, fc_lp, fs=FS), [1.0],
        f"FIR LP (ventana, {numtaps_fir} taps, fc={fc_lp} Hz)",
    )

    # --- Pasa altos ---
    filtros["pasa_altos"]["butterworth"] = (
        *butter(orden_iir, fc_hp, btype="high", fs=FS),
        f"Butterworth HP (orden {orden_iir}, fc={fc_hp} Hz)",
    )
    filtros["pasa_altos"]["chebyshev"] = (
        *cheby1(orden_iir, ripple_cheby, fc_hp, btype="high", fs=FS),
        f"Chebyshev I HP (orden {orden_iir}, rp={ripple_cheby} dB)",
    )
    filtros["pasa_altos"]["fir"] = (
        firwin(numtaps_fir, fc_hp, pass_zero=False, fs=FS), [1.0],
        f"FIR HP (ventana, {numtaps_fir} taps, fc={fc_hp} Hz)",
    )

    # --- Pasa banda ---
    filtros["pasa_banda"]["butterworth"] = (
        *butter(orden_iir, fc_bp, btype="band", fs=FS),
        f"Butterworth BP (orden {orden_iir}, fc={fc_bp} Hz)",
    )
    filtros["pasa_banda"]["chebyshev"] = (
        *cheby1(orden_iir, ripple_cheby, fc_bp, btype="band", fs=FS),
        f"Chebyshev I BP (orden {orden_iir}, rp={ripple_cheby} dB)",
    )
    filtros["pasa_banda"]["fir"] = (
        firwin(numtaps_fir, fc_bp, pass_zero=False, fs=FS), [1.0],
        f"FIR BP (ventana, {numtaps_fir} taps, fc={fc_bp} Hz)",
    )

    return filtros


# ===========================================================================
# 3. APLICACIÓN DE FILTROS
# ===========================================================================
def aplicar_filtro(b, a, senal):
    """Aplica un filtro con fase cero (filtfilt).

    filtfilt aplica el filtro hacia adelante y hacia atrás, eliminando el
    desfase. Es ideal para análisis offline. Para tiempo real usar lfilter.
    """
    return filtfilt(b, a, senal)


# ===========================================================================
# 4. VISUALIZACIÓN DE RESULTADOS
# ===========================================================================
def figura_senal_original(senal_util, ruido, senal_ruidosa):
    """Señal útil, ruido, señal compuesta y su espectro."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    axes[0, 0].plot(t, senal_util, color="green")
    axes[0, 0].set_title("Señal útil (sin ruido)")
    axes[0, 0].set_xlabel("Tiempo (s)"); axes[0, 0].set_ylabel("Amplitud")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(t, ruido, color="gray")
    axes[0, 1].set_title("Ruido blanco")
    axes[0, 1].set_xlabel("Tiempo (s)"); axes[0, 1].set_ylabel("Amplitud")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(t, senal_ruidosa, color="steelblue")
    axes[1, 0].set_title("Señal con ruido (entrada al filtro)")
    axes[1, 0].set_xlabel("Tiempo (s)"); axes[1, 0].set_ylabel("Amplitud")
    axes[1, 0].grid(True, alpha=0.3)

    f, X = espectro(senal_ruidosa)
    axes[1, 1].plot(f, X, color="steelblue")
    axes[1, 1].set_title("Espectro de la señal con ruido")
    axes[1, 1].set_xlabel("Frecuencia (Hz)"); axes[1, 1].set_ylabel("|X(f)|")
    axes[1, 1].set_xlim(0, FS / 2)
    axes[1, 1].grid(True, alpha=0.3)
    for fref, color in [(F_LOW, "green"), (F_MID, "orange"), (F_HIGH, "red")]:
        axes[1, 1].axvline(fref, color=color, ls=":", alpha=0.6)
        axes[1, 1].text(fref + 5, max(X) * 0.85, f"{fref:.0f} Hz",
                        color=color, fontsize=8, rotation=90)

    fig.suptitle("Señal de entrada: útil + ruido", fontsize=13)
    fig.tight_layout()
    guardar_figura("01_senal_entrada.png")
    plt.show()


def figura_tipo_filtro(tipo, filtros_tipo, senal_ruidosa):
    """Para un tipo de filtro (PB, PA o BP) muestra:
       - respuesta en frecuencia de los 3 métodos
       - señal filtrada en tiempo (3 métodos)
       - espectro antes/después (3 métodos)
       - tabla de atenuación en las componentes
    """
    nombres = {"pasa_bajos": "PASA BAJOS",
               "pasa_altos": "PASA ALTOS",
               "pasa_banda": "PASA BANDA"}
    titulo = nombres[tipo]

    # --- Figura A: respuestas en frecuencia ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, (_, (b, a, desc)) in zip(axes, filtros_tipo.items()):
        plot_respuesta(b, a, desc, ax=ax)
    fig.suptitle(f"Respuesta en frecuencia — Filtro {titulo}", fontsize=13)
    fig.tight_layout()
    guardar_figura(f"02_respuesta_{tipo}.png")
    plt.show()

    # --- Figura B: comparación en tiempo ---
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for ax, (_, (b, a, desc)) in zip(axes, filtros_tipo.items()):
        y = aplicar_filtro(b, a, senal_ruidosa)
        ax.plot(t, senal_ruidosa, color="gray", alpha=0.4, label="Con ruido")
        ax.plot(t, y, color="crimson", linewidth=1.5, label=desc)
        ax.set_ylabel("Amplitud")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Tiempo (s)")
    fig.suptitle(f"Señal filtrada en tiempo — {titulo}", fontsize=13)
    fig.tight_layout()
    guardar_figura(f"03_tiempo_{tipo}.png")
    plt.show()

    # --- Figura C: espectros antes y después ---
    f_in, X_in = espectro(senal_ruidosa)
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    for ax, (_, (b, a, desc)) in zip(axes, filtros_tipo.items()):
        y = aplicar_filtro(b, a, senal_ruidosa)
        f_out, X_out = espectro(y)
        ax.plot(f_in, X_in, color="gray", alpha=0.5, label="Antes")
        ax.plot(f_out, X_out, color="crimson", label="Después")
        ax.set_title(desc, fontsize=9)
        ax.set_ylabel("|X(f)|")
        ax.set_xlim(0, FS / 2)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Frecuencia (Hz)")
    fig.suptitle(f"Espectro antes/después — {titulo}", fontsize=13)
    fig.tight_layout()
    guardar_figura(f"04_espectro_{tipo}.png")
    plt.show()

    # --- Reporte numérico de atenuación ---
    print(f"\n--- Atenuación (dB) del filtro {titulo} en las componentes ---")
    print(f"{'Método':<15}{F_LOW:>10.0f} Hz{F_MID:>12.0f} Hz{F_HIGH:>12.0f} Hz")
    for metodo, (b, a, _) in filtros_tipo.items():
        at = atenuacion_db(b, a, [F_LOW, F_MID, F_HIGH])
        print(f"{metodo:<15}{at[0]:>10.2f}{at[1]:>12.2f}{at[2]:>12.2f}")


# ===========================================================================
# PROGRAMA PRINCIPAL
# ===========================================================================
def main():
    np.random.seed(42)  # reproducibilidad

    print("=" * 60)
    print("ACTIVIDAD 3 - FILTROS DIGITALES")
    print("=" * 60)

    # 1. Señal de entrada
    print("\n[1] Generando señal de entrada...")
    senal_util, ruido, senal_ruidosa = generar_senal()
    print(f"    Fs = {FS} Hz, N = {N} muestras, duración = {T} s")
    print(f"    Componentes: {F_LOW} Hz, {F_MID} Hz, {F_HIGH} Hz")
    print(f"    Amplitud del ruido: {RUIDO_AMP}")
    figura_senal_original(senal_util, ruido, senal_ruidosa)

    # 2. Diseño de filtros
    print("\n[2] Diseñando filtros...")
    filtros = disenar_filtros()
    for tipo, metodos in filtros.items():
        print(f"    {tipo}: {list(metodos.keys())}")

    # 3 y 4. Aplicación y visualización
    for tipo, filtros_tipo in filtros.items():
        print(f"\n[3-4] Procesando filtro {tipo}...")
        figura_tipo_filtro(tipo, filtros_tipo, senal_ruidosa)

    print("\n" + "=" * 60)
    print("Proceso completado. Figuras en la carpeta 'figuras/'.")
    print("=" * 60)


if __name__ == "__main__":
    main()