# Actividad 3 — Implementación y evaluación de filtros digitales

**Curso:** Análisis de Fourier — Unidad 3
**Carrera:** Ingeniería en Desarrollo de Software
**Autor:** Alfonso Rios Belmares
**Fecha:** 15 de Noviembre de 2026

---

## 1. Descripción

Este repositorio contiene el código y el análisis de la **Actividad Formativa 3**,
cuyo objetivo es diseñar, implementar y evaluar **filtros digitales pasa bajos,
pasa altos y pasa banda** aplicados a una señal con ruido, usando **Python**
(`numpy`, `scipy.signal`, `matplotlib`).

Se comparan tres familias de filtros para cada tipo:

- **IIR Butterworth** (respuesta máximamente plana en banda de paso).
- **IIR Chebyshev I** (transición más abrupta a costa de rizado).
- **FIR con ventana** (fase lineal, estabilidad garantizada).

La efectividad se evalúa comparando el **espectro de amplitud (FFT)** y la
**señal en el tiempo** antes y después del filtrado, además de cuantificar la
**atenuación en dB** en cada componente de la señal.

---

## 2. Estructura del repositorio
actividad3-filtros-digitales/
├── README.md
├── requirements.txt
├── filtros_digitales.py
└── figuras/ (se genera al ejecutar)


---

## 3. Requisitos

- Python 3.9 o superior
- Librerías:

```bash
pip install -r requirements.txt

## 4. Ejecución
python filtros_digitales.py
