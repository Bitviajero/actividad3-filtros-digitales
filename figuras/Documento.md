# UNIVERSIDAD CIUDADANA DE NUEVO LEÓN
## Ingeniería en Desarrollo de Software
### Análisis de Fourier — Unidad 3

# ACTIVIDAD FORMATIVA 3: IMPLEMENTACIÓN Y EVALUACIÓN DE FILTROS DIGITALES

**Autor:** [Tu nombre completo]
**Matrícula:** [Tu matrícula]
**Profesor:** [Nombre del profesor]
**Fecha:** [Mes Año]

---

## Índice
1. Introducción
2. Desarrollo de la actividad
   2.1 Marco teórico
   2.2 Metodología y herramientas
   2.3 Generación de la señal de prueba
   2.4 Filtro pasa bajos
   2.5 Filtro pasa altos
   2.6 Filtro pasa banda
   2.7 Comparación entre filtros IIR y FIR
3. Conclusión
4. Referencias

---

## 1. Introducción

El procesamiento digital de señales es una de las áreas de mayor aplicación en la ingeniería moderna. Desde la eliminación de ruido en señales biomédicas hasta el filtrado de audio en aplicaciones móviles, el diseño de filtros digitales permite separar la información útil de las componentes no deseadas de una señal. En el contexto de la asignatura de Análisis de Fourier, esta actividad tiene como propósito aplicar los conceptos de espectro de frecuencia, transformada discreta de Fourier y respuesta en frecuencia al diseño y evaluación de filtros digitales de tipo pasa bajos, pasa altos y pasa banda.

El objetivo específico de la actividad es diseñar e implementar filtros digitales utilizando Python, aplicarlos a una señal compuesta por varias sinusoides y ruido blanco, y evaluar su efectividad comparando la señal y su espectro antes y después del filtrado. Para ello se emplean tres familias de filtros: dos de tipo IIR (Butterworth y Chebyshev I) y una de tipo FIR (diseñada con ventana). Esta comparación permite observar el compromiso entre selectividad, fase y costo computacional que caracteriza a cada familia.

El documento se organiza de la siguiente manera: en la sección 2 se presenta el marco teórico mínimo, la metodología, la señal de prueba y el análisis de cada tipo de filtro con sus resultados; en la sección 3 se exponen las conclusiones; y en la sección 4 se listan las referencias bibliográficas y técnicas utilizadas.

---

## 2. Desarrollo de la actividad

### 2.1 Marco teórico

**Contenido de frecuencia y espectro.** Toda señal de tiempo discreto $x[n]$ puede describirse en términos de su contenido de frecuencia mediante la transformada de Fourier de tiempo discreto (DTFT). Para el análisis numérico se emplea la transformada discreta de Fourier (DFT), que corresponde a muestrear la DTFT en $N$ frecuencias equiespaciadas:

$$ X_k = \sum_{n=0}^{N-1} x[n] e^{-j2\pi kn/N}, \quad k = 0, 1, \dots, N-1 $$

La DFT se calcula de forma eficiente mediante el algoritmo de la transformada rápida de Fourier (FFT), cuyo costo es del orden de $N\log_2 N$ en lugar de $N^2$. En este trabajo se utiliza la FFT para obtener el espectro de amplitud de la señal antes y después del filtrado. La frecuencia asociada a cada índice $k$ es $f_k = k F_s/N$, donde $F_s$ es la frecuencia de muestreo.

**Filtros digitales y respuesta en frecuencia.** Un filtro digital es un sistema lineal e invariante en el tiempo (LTI) que modifica el contenido espectral de una señal. Su comportamiento se describe mediante la respuesta en frecuencia $H(f)$, que indica cuánto se atenúa o amplifica cada componente de frecuencia. Un filtro se caracteriza por tres regiones: banda de paso (frecuencias que se conservan), banda de transición (frecuencias intermedias) y banda de rechazo (frecuencias que se atenúan). La atenuación se expresa en decibelios como $20\log_{10}|H(f)|$.

**Filtros IIR y FIR.** Los filtros de respuesta al impulso infinita (IIR) se describen mediante una ecuación de diferencias con realimentación y ofrecen alta selectividad con un orden bajo, aunque su fase no es lineal y pueden presentar inestabilidad numérica. Los filtros de respuesta al impulso finita (FIR) no tienen realimentación, garantizan estabilidad y permiten fase lineal, pero requieren un mayor número de coeficientes para lograr la misma selectividad. En este trabajo se emplean tres familias: Butterworth (IIR, respuesta máximamente plana en la banda de paso), Chebyshev I (IIR, transición más abrupta a costa de rizado en la banda de paso) y FIR con ventana (fase lineal).

**Tipos de filtro.** El filtro pasa bajos conserva las frecuencias por debajo de una frecuencia de corte $f_c$ y atenúa las superiores. El pasa altos hace lo contrario: conserva las frecuencias por encima de $f_c$. El pasa banda conserva un intervalo $[f_1, f_2]$ y atenúa las frecuencias fuera de él.

### 2.2 Metodología y herramientas

La actividad se desarrolló en Python 3.9 o superior, utilizando las siguientes librerías:

- **NumPy**: generación de la señal, cálculo de la FFT y operaciones numéricas.
- **SciPy (`scipy.signal`)**: diseño de filtros (`butter`, `cheby1`, `firwin`), cálculo de la respuesta en frecuencia (`freqz`) y aplicación de los filtros (`filtfilt`).
- **Matplotlib**: generación de las gráficas.

El código se organiza en cuatro secciones que corresponden a las etapas de la actividad:

1. **Definición de la señal de entrada**: generación de la señal útil, el ruido y la señal compuesta, y cálculo de su espectro.
2. **Diseño del filtro**: construcción de los coeficientes $(b, a)$ para cada familia y tipo de filtro.
3. **Aplicación del filtro a la señal**: uso de `filtfilt` para filtrar sin introducir desfase.
4. **Visualización de los resultados**: respuesta en frecuencia, señal en tiempo y espectro antes/después, además de una tabla de atenuación en decibelios.

El repositorio del proyecto se encuentra disponible en GitHub: [URL de tu repositorio]. Contiene el script `filtros_digitales.py`, el archivo `requirements.txt`, el `README.md` y la carpeta `figuras/` con las gráficas generadas.

Se utilizó `filtfilt` en lugar de `lfilter` porque aplica el filtro hacia adelante y hacia atrás, eliminando el desfase y permitiendo una comparación visual limpia entre la señal original y la filtrada. Para aplicaciones en tiempo real sería necesario usar `lfilter` y considerar el retardo introducido.

### 2.3 Generación de la señal de prueba

Se generó una señal de prueba compuesta por tres sinusoides de 50 Hz, 150 Hz y 300 Hz, sumadas a ruido blanco gaussiano de amplitud 0.8. La señal se muestreó a $F_s = 1000$ Hz durante 2 segundos, lo que produce $N = 2000$ muestras. La elección de $F_s$ cumple el teorema de Nyquist para la componente de mayor frecuencia (300 Hz), ya que $F_s > 2 \cdot 300 = 600$ Hz.

La Figura 1 muestra, de arriba a abajo y de izquierda a derecha: la señal útil sin ruido, el ruido blanco, la señal compuesta que se usará como entrada a los filtros, y el espectro de amplitud de esta última. En el espectro se distinguen claramente los picos correspondientes a las tres componentes (50, 150 y 300 Hz) sobre un fondo de ruido distribuido en todo el rango de frecuencias.

![Figura 1. Señal de entrada al sistema de filtrado.](figuras/01_senal_entrada.png)
*Figura 1. Señal de entrada al sistema de filtrado: (a) señal útil, (b) ruido blanco, (c) señal compuesta y (d) espectro de amplitud.*

Esta señal permite evaluar los tres tipos de filtro de forma clara:
- Un **pasa bajos** con frecuencia de corte de 80 Hz debe conservar la componente de 50 Hz y atenuar las de 150 Hz y 300 Hz.
- Un **pasa altos** con frecuencia de corte de 220 Hz debe conservar la componente de 300 Hz y atenuar las de 50 Hz y 150 Hz.
- Un **pasa banda** con frecuencias de corte de 100 Hz y 200 Hz debe conservar la componente de 150 Hz y atenuar las otras dos.

### 2.4 Filtro pasa bajos

Se diseñaron tres filtros pasa bajos con frecuencia de corte de 80 Hz: un Butterworth de orden 6, un Chebyshev I de orden 6 con rizado de 1 dB, y un FIR con ventana de 201 coeficientes. La Figura 2 muestra la respuesta en frecuencia de los tres filtros.

En la figura se observa que los tres filtros conservan la componente de 50 Hz (atenuación cercana a 0 dB) y atenúan progresivamente las componentes de 150 Hz y 300 Hz. El filtro de Chebyshev presenta una transición más abrupta que el de Butterworth, mientras que el FIR, debido a su mayor número de coeficientes, presenta una transición más suave pero con fase lineal.

![Figura 2. Respuesta en frecuencia de los filtros pasa bajos.](figuras/02_respuesta_pasa_bajos.png)
*Figura 2. Respuesta en frecuencia de los filtros pasa bajos (Butterworth, Chebyshev I y FIR).*

La Figura 3 compara la señal original con la salida de cada filtro en el dominio del tiempo. Se observa que la señal filtrada por los tres métodos conserva la componente de baja frecuencia y elimina las variaciones rápidas asociadas al ruido y a las componentes de 150 Hz y 300 Hz.

![Figura 3. Señal filtrada en el dominio del tiempo por los tres filtros pasa bajos.](figuras/03_tiempo_pasa_bajos.png)
*Figura 3. Señal filtrada en el dominio del tiempo por los tres filtros pasa bajos.*

La Figura 4 muestra el espectro de amplitud antes y después del filtrado. En los tres casos se aprecia que el pico de 50 Hz se conserva, mientras que los picos de 150 Hz y 300 Hz se atenúan significativamente. La Tabla 1 cuantifica esta atenuación en decibelios para cada método.

![Figura 4. Espectro de amplitud antes y después del filtrado pasa bajos.](figuras/04_espectro_pasa_bajos.png)
*Figura 4. Espectro de amplitud antes y después del filtrado pasa bajos.*

**Tabla 1. Atenuación en dB del filtro pasa bajos en cada componente.**

| Método | 50 Hz (Paso) | 150 Hz (Rechazo) | 300 Hz (Rechazo) |
| :--- | :--- | :--- | :--- |
| Butterworth | ~ 0 dB | ~ -35 dB | < -80 dB |
| Chebyshev I | ~ -1 dB | ~ -50 dB | < -80 dB |
| FIR (ventana) | ~ 0 dB | ~ -60 dB | < -80 dB |

La tabla confirma que los tres filtros cumplen su función: la componente de 50 Hz prácticamente no se atenúa, mientras que las de 150 Hz y 300 Hz presentan una atenuación considerable. El filtro de Chebyshev logra mayor atenuación en 150 Hz que el de Butterworth con el mismo orden, a costa de un rizado de 1 dB en la banda de paso. El FIR, con 201 coeficientes, ofrece una atenuación comparable pero con fase lineal.

### 2.5 Filtro pasa altos

Se diseñaron tres filtros pasa altos con frecuencia de corte de 220 Hz, con los mismos parámetros que en el caso anterior. La Figura 5 muestra su respuesta en frecuencia.

![Figura 5. Respuesta en frecuencia de los filtros pasa altos.](figuras/02_respuesta_pasa_altos.png)
*Figura 5. Respuesta en frecuencia de los filtros pasa altos.*

A diferencia del pasa bajos, en este caso la banda de paso se ubica por encima de 220 Hz, por lo que el filtro debe conservar la componente de 300 Hz y atenuar las de 50 Hz y 150 Hz. La Figura 6 muestra la señal filtrada en el tiempo, donde se observa que la componente de alta frecuencia (300 Hz) domina la salida, mientras que las variaciones lentas asociadas a 50 Hz y 150 Hz desaparecen.

![Figura 6. Señal filtrada en el dominio del tiempo por los tres filtros pasa altos.](figuras/03_tiempo_pasa_altos.png)
*Figura 6. Señal filtrada en el dominio del tiempo por los tres filtros pasa altos.*

La Figura 7 presenta el espectro antes y después del filtrado. Se confirma que el pico de 300 Hz se conserva y que los picos de 50 Hz y 150 Hz se atenúan. La Tabla 2 cuantifica la atenuación.

![Figura 7. Espectro de amplitud antes y después del filtrado pasa altos.](figuras/04_espectro_pasa_altos.png)
*Figura 7. Espectro de amplitud antes y después del filtrado pasa altos.*

**Tabla 2. Atenuación en dB del filtro pasa altos en cada componente.**

| Método | 50 Hz (Rechazo) | 150 Hz (Rechazo) | 300 Hz (Paso) |
| :--- | :--- | :--- | :--- |
| Butterworth | < -80 dB | ~ -25 dB | ~ 0 dB |
| Chebyshev I | < -80 dB | ~ -30 dB | ~ 0 dB |
| FIR (ventana) | ~ -70 dB | ~ -40 dB | ~ 0 dB |

Los tres filtros atenúan correctamente las componentes de baja frecuencia y conservan la de 300 Hz. Al igual que en el caso anterior, el comportamiento es consistente con la respuesta en frecuencia diseñada: la banda de rechazo se ubica ahora en las frecuencias bajas.

### 2.6 Filtro pasa banda

Se diseñaron tres filtros pasa banda con frecuencias de corte de 100 Hz y 200 Hz. Este caso es el más selectivo de los tres, ya que requiere dos frecuencias de corte y, por tanto, mayor orden para lograr la misma pendiente que los casos anteriores. La Figura 8 muestra la respuesta en frecuencia.

![Figura 8. Respuesta en frecuencia de los filtros pasa banda.](figuras/02_respuesta_pasa_banda.png)
*Figura 8. Respuesta en frecuencia de los filtros pasa banda.*

La banda de paso se ubica entre 100 Hz y 200 Hz, por lo que el filtro debe conservar la componente de 150 Hz y atenuar las de 50 Hz y 300 Hz. La Figura 9 muestra la señal filtrada en el tiempo.

![Figura 9. Señal filtrada en el dominio del tiempo por los tres filtros pasa banda.](figuras/03_tiempo_pasa_banda.png)
*Figura 9. Señal filtrada en el dominio del tiempo por los tres filtros pasa banda.*

La Figura 10 presenta el espectro antes y después del filtrado. Se observa que solo el pico de 150 Hz se conserva, mientras que los de 50 Hz y 300 Hz se atenúan. La Tabla 3 cuantifica la atenuación.

![Figura 10. Espectro de amplitud antes y después del filtrado pasa banda.](figuras/04_espectro_pasa_banda.png)
*Figura 10. Espectro de amplitud antes y después del filtrado pasa banda.*

**Tabla 3. Atenuación en dB del filtro pasa banda en cada componente.**

| Método | 50 Hz (Rechazo) | 150 Hz (Paso) | 300 Hz (Rechazo) |
| :--- | :--- | :--- | :--- |
| Butterworth | < -80 dB | ~ 0 dB | ~ -60 dB |
| Chebyshev I | < -80 dB | ~ 0 dB | ~ -70 dB |
| FIR (ventana) | ~ -70 dB | ~ 0 dB | ~ -70 dB |

El filtro pasa banda es el que exige mayor selectividad, ya que debe atenuar componentes en ambos extremos del espectro. Los resultados confirman que los tres métodos logran conservar la componente central y atenuar las laterales. El filtro FIR, al tener fase lineal, preserva mejor la forma de onda de la componente de 150 Hz, mientras que los IIR introducen un desfase que, en este caso, se compensa con `filtfilt`.

### 2.7 Comparación entre filtros IIR y FIR

A partir de los resultados de las tres secciones anteriores, es posible establecer una comparación entre las familias de filtros empleadas. La Tabla 4 resume las diferencias observadas.

**Tabla 4. Comparación entre filtros IIR y FIR.**

| Criterio | IIR (Butterworth / Chebyshev) | FIR (ventana) |
| :--- | :--- | :--- |
| Selectividad | Alta con orden bajo | Requiere muchos coeficientes |
| Fase | No lineal | Lineal |
| Estabilidad | Puede ser inestable en órdenes altos | Siempre estable |
| Costo computacional | Bajo | Alto (más coeficientes) |
| Rizado en banda de paso | Butterworth: nulo; Chebyshev: controlado | Nulo |
| Uso típico | Filtrado en tiempo real | Análisis offline, audio |

Los resultados confirman lo esperado teóricamente. El filtro de Butterworth ofrece la respuesta más plana en la banda de paso y es el más sencillo de interpretar. El de Chebyshev I logra una transición más abrupta con el mismo orden, a costa de un rizado de 1 dB. El FIR, con 201 coeficientes, ofrece una transición más suave pero garantiza fase lineal, lo cual es importante cuando se desea preservar la forma de onda de la señal filtrada. En términos de costo computacional, el FIR requiere aproximadamente el doble de operaciones por muestra que un IIR de orden 6.

Cabe destacar que el truncamiento implícito en el diseño del FIR y en el análisis mediante FFT produce fuga espectral, fenómeno estudiado en la Unidad 3. Esto explica por qué las sinusoides no aparecen como líneas perfectas en el espectro, sino como picos con lóbulos laterales. Dicha fuga no invalida los resultados, pero debe considerarse al interpretar las gráficas. Adicionalmente, el filtro FIR presenta lóbulos laterales en la banda de rechazo debido al truncamiento de su respuesta al impulso, un fenómeno análogo al de Gibbs estudiado en la teoría de series de Fourier.

---

## 3. Conclusión

La actividad permitió aplicar de forma integrada los conceptos de la Unidad 3: espectro de frecuencia, transformada discreta de Fourier, muestreo y fuga espectral, en el contexto del diseño y evaluación de filtros digitales.

Se diseñaron e implementaron filtros pasa bajos, pasa altos y pasa banda en tres familias distintas (Butterworth, Chebyshev I y FIR con ventana), y se evaluó su efectividad comparando la señal en el tiempo y su espectro de amplitud antes y después del filtrado. Los resultados confirman que los tres tipos de filtro cumplen su objetivo: atenúan las componentes no deseadas y conservan las deseadas, con diferencias en la forma de la transición, el rizado y la fase.

La comparación entre familias permite concluir que no existe un filtro "mejor" en abstracto, sino que la elección depende de los requisitos de la aplicación: si se requiere alta selectividad con bajo costo, un IIR (Butterworth o Chebyshev) es adecuado; si se requiere fase lineal y estabilidad garantizada, un FIR es preferible, a costa de un mayor número de coeficientes.

La FFT resultó ser una herramienta fundamental para la evaluación, ya que permitió cuantificar en decibelios la atenuación lograda en cada componente y verificar que el comportamiento observado corresponde al diseñado. El uso de `filtfilt` facilitó la comparación visual al eliminar el desfase introducido por los filtros IIR.

Finalmente, el trabajo permitió reforzar la relación entre la teoría de Fourier y la práctica del procesamiento digital de señales, y evidenció la importancia de elegir correctamente la frecuencia de muestreo para evitar el aliasing y de interpretar el espectro considerando fenómenos como la fuga espectral por truncamiento.

---

## 4. Referencias

Kamen, E. W., & Heck, B. S. (2008). *Fundamentos de señales y sistemas usando la Web y MATLAB* (3.ª ed.). Pearson Educación.

SciPy Developers. (2024). *scipy.signal* [Documentación]. https://docs.scipy.org/doc/scipy/reference/signal.html

NumPy Developers. (2024). *numpy.fft* [Documentación]. https://numpy.org/doc/stable/reference/routines.fft.html

Matplotlib Development Team. (2024). *Matplotlib documentation*. https://matplotlib.org/stable/

[Tu usuario]. (2025). *actividad3-filtros-digitales* [Repositorio de GitHub]. https://github.com/[tu-usuario]/actividad3-filtros-digitales