# Teoría del análisis espectral de la VDF de protones de Solar Orbiter / SWA-PAS

*Documento de referencia a nivel de posgrado para la reconstrucción de la Función de Distribución de Velocidad (VDF, *Velocity Distribution Function*) de protones del viento solar medida por el sensor SWA-PAS de Solar Orbiter, su expresión en coordenadas alineadas al campo magnético, y su descomposición sobre una base ortonormal de funciones de Hermite mediante cuadratura de Gauss, hasta el espectro bidimensional $\log_{10}(c_{mn}^2/c_{00}^2)$.*

Este documento acompaña a la implementación contenida en `src/` (módulos `config.py`, `cdf_io.py`, `physics.py`, `hermite.py`, `plots.py` y el orquestador `run_analysis.py`). A lo largo del texto se enlaza cada bloque teórico con la función concreta que lo realiza, de modo que la teoría y el código formen una sola pieza. El caso de estudio reproducido es el instante **2022-03-08 14:45:22** (Sección 10).

---

> ## Tres correcciones de atribución (declaradas una sola vez, al inicio)
>
> Estas correcciones provienen de la investigación bibliográfica primaria y deben tenerse presentes en todo el documento:
>
> **(i) No existe un artículo independiente "Livi et al. 2023 de calibración en vuelo de PAS".** La referencia *Livi, S., et al. 2023, A&A, 676, A36* (DOI `10.1051/0004-6361/202346304`) corresponde a **"First results from the Solar Orbiter Heavy Ion Sensor"**, es decir el sensor **SWA-HIS**, un instrumento distinto de la misma suite SWA. La documentación instrumental y de calibración de **PAS** está en **Owen et al. (2020)** (construcción y calibración, su Tabla 8), la **Guía de Usuario de Datos L2 de Fedorov (2020)** y **Louarn et al. (2024)** (momentos del viento solar y calibración cruzada en vuelo). Cítense estas y no un inexistente "Livi-PAS-2023".
>
> **(ii) El artículo de descomposición Hermite-Hermite de Larosa et al. (2025) es de Parker Solar Probe (PSP), no de Solar Orbiter.** Larosa et al. (2025) (arXiv:2512.01492) aplica la descomposición HH con girotropía "cartesiana" (eje perpendicular reflejado, malla $50\times50$) a VDFs de protones del instrumento **SPAN-i de PSP**. El trabajo espectral análogo en **Solar Orbiter** es **Coburn et al. (2024)** (electrones de SWA-EAS, base Hermite-Laguerre, malla $60\times60$). Este proyecto adopta la metodología HH de Larosa pero aplicada a datos de Solar Orbiter / SWA-PAS, lo cual es legítimo siempre que se atribuya correctamente el origen de cada elemento.
>
> **(iii) Las pendientes específicas "$m^{-1.3}$ (paralela), $n^{-3.1}$ (perpendicular)" y el estimador de pendiente local $\alpha_m=(c_{m+2,0}^2-c_{m,0}^2)/(2c_{00}^2)$ son convenciones de trabajo de este curso/proyecto**, no citas textuales de la literatura. Son consistentes *en espíritu* con la cascada anisótropa de Pezzi et al. (2018) (paralela somera, perpendicular abrupta) y con los cortes direccionales de Larosa et al. (2025), y se tratan como tales en la Sección 9. La desigualdad robusta y citable es $|\text{pendiente}_\parallel| < |\text{pendiente}_\perp|$ (cascada paralela extendida, perpendicular suprimida).

---

## Índice

1. [Introducción y contexto físico](#1-introducción-y-contexto-físico)
2. [El instrumento SWA-PAS y el producto L2 `swa-pas-vdf`](#2-el-instrumento-swa-pas-y-el-producto-l2-swa-pas-vdf)
3. [Reconstrucción del espacio de velocidades](#3-reconstrucción-del-espacio-de-velocidades)
4. [Marco alineado al campo y momentos del plasma](#4-marco-alineado-al-campo-y-momentos-del-plasma)
5. [Modelos paramétricos y sesgo humano](#5-modelos-paramétricos-y-sesgo-humano)
6. [Descomposición polinomial: la analogía de Fourier](#6-descomposición-polinomial-la-analogía-de-fourier)
7. [Cuadratura gaussiana](#7-cuadratura-gaussiana)
8. [Interpolación sobre la malla de nodos](#8-interpolación-sobre-la-malla-de-nodos)
9. [El espectro 2D de Hermite](#9-el-espectro-2d-de-hermite)
10. [Resultados de esta reproducción](#10-resultados-de-esta-reproducción)
11. [Referencias](#11-referencias)

**Notación.** $f(\mathbf v)$ es la densidad en el espacio de fases (unidades $\mathrm{s^3\,m^{-6}}$ en SI, o $\mathrm{s^3\,cm^{-6}}$ en CGS); $m$ y $q$ son la masa y carga de la especie (protón: $m=m_p$, $q=e$); $\mathbf B$ el campo magnético y $\hat{\mathbf b}=\mathbf B/|\mathbf B|$ su versor; $k_B$ la constante de Boltzmann. Los índices $i,j,k$ recorren los *bins* de energía, elevación y azimut respectivamente. Todos los sistemas de referencia son dextrógiros.

---

## 1. Introducción y contexto físico

### 1.1 El viento solar como plasma débilmente colisional y magnetizado

El viento solar es un plasma **débilmente colisional** y **magnetizado** que fluye de forma supersónica desde la corona. La escala física esencial es la siguiente comparación de tiempos: el tiempo de colisión protón-protón, $\tau_{pp}$, es típicamente comparable o **mayor** que el tiempo de expansión del flujo (el tiempo que tarda el plasma en recorrer una escala apreciable de distancia heliocéntrica). En consecuencia, las colisiones binarias **no alcanzan a relajar la distribución de velocidades hacia el equilibrio termodinámico local** antes de que las condiciones macroscópicas cambien apreciablemente (Marsch 2006; Verscharen, Klein & Maruca 2019).

Esta es la diferencia fundamental con un gas ordinario. En un gas colisional, el Teorema H garantiza que la VDF relaja rápidamente a una **Maxwelliana** (la solución de máxima entropía); toda la física queda capturada por unos pocos momentos fluidos ($n$, $\mathbf u$, $T$) y la descripción de Navier-Stokes es exacta hasta correcciones pequeñas. En el viento solar, en cambio, la VDF medida **retiene estructura cinética** que una descripción fluida promedia y descarta:

- **Anisotropía de temperatura** $T_\perp \neq T_\parallel$ respecto al campo magnético, sostenida porque el campo $\mathbf B$ rompe la isotropía y la giromoción termaliza solo las dos direcciones perpendiculares. Los protones del viento rápido muestran de forma persistente $T_{\perp p} > T_{\parallel p}$, una signatura que la expansión adiabática por sí sola no puede producir y que exige calentamiento perpendicular activo.
- **Haces (beams) alineados al campo**: una segunda población de protones que fluye más rápido que el núcleo a lo largo de $\mathbf B$ y casi siempre en sentido anti-solar.
- **Colas supratérmicas** (no Maxwellianas, de tipo ley de potencias) que delatan la naturaleza colisionless del plasma y se modelan con distribuciones $\kappa$.

### 1.2 La VDF como densidad en el espacio de fases

La VDF $f(\mathbf v)$ es la **densidad de partículas en el espacio de fases de velocidades**: $f(\mathbf v)\,d^3v$ es el número de partículas por unidad de volumen físico cuya velocidad cae en el elemento $d^3v$ alrededor de $\mathbf v$. Todas las cantidades fluidas son **momentos** de $f$, es decir integrales en el espacio de velocidades ponderadas por potencias de $\mathbf v$:

$$
n=\int f\,d^3v,\qquad
\mathbf u=\frac1n\int \mathbf v\,f\,d^3v,\qquad
\mathsf P=m\int(\mathbf v-\mathbf u)(\mathbf v-\mathbf u)\,f\,d^3v,
$$

con el tensor de temperatura $\mathsf T=\mathsf P/(n k_B)$ y el flujo de calor $\mathbf q=\tfrac{m}{2}\int(\mathbf v-\mathbf u)|\mathbf v-\mathbf u|^2 f\,d^3v$. La densidad $n$ es el momento de orden 0, la velocidad de flujo $\mathbf u$ el de orden 1, la presión/temperatura el de orden 2 y el flujo de calor el de orden 3.

La lógica del proyecto se desprende directamente de aquí: puesto que toda la información cinética vive en $f(\mathbf v)$ y los momentos solo son proyecciones parciales, **recuperar y descomponer la $f(\mathbf v)$ medida es la manera natural de estudiar la física cinética directamente**, sin colapsarla prematuramente a un puñado de momentos. El resto del documento desarrolla, paso a paso, cómo (a) se reconstruye $f(\mathbf v)$ desde las cuentas del instrumento, (b) se la expresa en coordenadas físicamente significativas alineadas al campo, y (c) se la descompone sobre una base ortonormal completa que cuantifica, modo a modo, cuánto se aparta del equilibrio Maxwelliano.

---

## 2. El instrumento SWA-PAS y el producto L2 `swa-pas-vdf`

### 2.1 La suite SWA y el sensor PAS

La suite **Solar Wind Analyser (SWA)** a bordo de **Solar Orbiter** consta de tres analizadores electrostáticos de tipo *top-hat*: **EAS** (electrones), **PAS** (protones y partículas alfa) y **HIS** (iones pesados) (Owen et al. 2020, §1, §3.3). El **Proton-Alpha Sensor (PAS)** es un analizador electrostático *top-hat* de sección parcialmente esférica (cuarto de esfera), optimizado para el haz de iones del viento solar, que es estrecho y fuertemente radial. PAS mide la VDF 3D completa de los iones del viento solar **sin separación de masa ni de carga** (en la práctica protones $\mathrm{H^+}$ y alfas $\mathrm{He^{2+}}$) en el rango **200 eV/q – 20 keV/q**, montado tras una abertura dedicada en el escudo térmico de la nave de modo que apunta aproximadamente hacia el Sol (Owen et al. 2020, §3.3.1–3.3.2).

### 2.2 Analizador electrostático *top-hat* y los tres ejes de medida

PAS forma su matriz de cuentas 3D a partir de **tres ejes de medida** físicamente distintos:

**1. Energía por carga ($E/q$) — el analizador electrostático en sí.** Un analizador de sección esférica solo transmite los iones cuya $E/q$ cae en una banda estrecha fijada por el voltaje del hemisferio interior, $E = k\,V_{\rm hemi}$, con **constante del analizador $k\approx 13$–$14$ eV/V** (Owen, Tabla 8). PAS barre **96 niveles de energía espaciados logarítmicamente** (cuasi-geométricos), de los cuales se usan típicamente 92 por muestra, cubriendo 200 eV/q – 20 keV/q, con **resolución relativa $\Delta E/E\approx 5.5\%$** de diseño (3.0–9.3% medida). Por ser logarítmico el espaciado, cada canal tiene **semianchos asimétricos** `delta_p_Energy` y `delta_m_Energy` (no iguales entre sí).

**2. Elevación (ángulo polar) — los deflectores de entrada.** Dos placas deflectoras curvas más un electrodo de tapa superior (*top-cap*) dirigen electrostáticamente los iones hacia la entrada del analizador; la razón de voltajes deflector-a-analizador $U_{\rm def}/U_{\rm an}$ selecciona el índice de elevación. Se reparte en **9 *bins* de elevación** que cubren $\approx -22.5^\circ$ a $+22.5^\circ$ (resolución $\sim 5^\circ$).

**3. Azimut — el arreglo de 11 CEM.** La posición azimutal del impacto sobre el detector codifica el índice de azimut. PAS usa un arreglo de **11 multiplicadores de electrones de canal (CEM, *channeltrons*)**, leídos **simultáneamente**, que cubren **$-24^\circ$ a $+42^\circ$** (un abanico de $66^\circ$ desplazado $+9^\circ$ respecto a la dirección del Sol para captar la aberración del viento solar), con resolución $\sim 5^\circ$.

**Esquema de barrido.** A energía fija, PAS barre la elevación mientras acumula los 11 *bins* de azimut a la vez; al completar el barrido de elevación, escalona a la siguiente energía. Una matriz completa $(96\times 9\times 11)$ se adquiere en $\sim 1$ s, con acumulación por elemento de $\approx 1$ ms. Un algoritmo de **seguimiento de pico** (*peak-tracking*) recentra una ventana reducida de energía-elevación sobre el haz para habilitar los modos rápidos.

**Cadencia y modos** (Fedorov §1):
- **Normal**: una VDF completa cada **~4 s** (el modo de este proyecto).
- **Snapshot**: ráfagas de ~9 s cada 300 s a ~4 Hz.
- **Burst**: 300 s continuos hasta ~15–20 Hz, con espacio de fases reducido.

### 2.3 El producto CDF `swa-pas-vdf`

El producto 3D es el conjunto de datos CDF **`solo_L2_swa-pas-vdf`** (productos hermanos: `swa-pas-eflux`, flujo de energía diferencial omnidireccional 1D; `swa-pas-grnd-mom`, momentos calculados en tierra). La referencia autoritativa para nombres de variables, unidades y convenciones de marco es la **Guía de Usuario de Datos L2 de SWA-PAS** (Fedorov 2020, V02).

**Variable de datos principal.**

| Variable | Dimensión | Unidades | Notas |
|---|---|---|---|
| `vdf` | **[tiempo, 11, 9, 96]** | **$\mathrm{s^3\,m^{-6}}$** | Densidad en espacio de fases $f$, con `DEPEND_1=Azimuth`, `DEPEND_2=Elevation`, `DEPEND_3=Energy`. |

**Orden de las dimensiones (crítico).** La matriz por registro se almacena como **[azimut = 11, elevación = 9, energía = 96]**. Owen et al. (2020) la escriben conceptualmente como $(96,9,11)=$ energía×elevación×azimut, que es **el orden inverso** del de los ejes CDF. El orden autoritativo es el de `DEPEND_1/2/3`. En la implementación, `match_axes()` en `src/run_analysis.py` resuelve los ejes **en tiempo de ejecución** comparando los tamaños de la matriz `vdf` con `(n_E, n_el, n_az)`, lo cual evita por completo cualquier ambigüedad de orden.

**Unidades.** $\mathrm{s^3\,m^{-6}}$ es la densidad en espacio de fases en SI, tal que $n=\int f\,d^3v$ con $d^3v$ en $(\mathrm{m/s})^3$ da $\mathrm{m^{-3}}$. Para convertir a la unidad CGS habitual del plasma $\mathrm{s^3\,cm^{-6}}$ se multiplica por $10^{-6}$.

**Variables de coordenadas y de soporte** (Fedorov 2020, Tabla 1):

| Variable | Dimensión | Unidades | Significado |
|---|---|---|---|
| `Epoch` | [tiempo] | ns (TT2000) | Tiempo del registro. |
| `Info` | [tiempo] | — | Categoría de muestreo: 0 Tierra, 1 Normal, 2 Snapshot, 3 Burst, 4 Ingeniería, 5 Calibración. **Usar solo Info = 1/2/3.** |
| `Energy` | [96] | eV | Centros de *bin*; semianchos asimétricos `delta_p_Energy`, `delta_m_Energy`. |
| `Azimuth` | [11] | grados | Centros de los *bins* CEM; semiancho `delta_Azimuth` (simétrico). |
| `Elevation` | [9] | grados | Centros de los *bins* de elevación; semiancho `delta_Elevation` (simétrico). |
| `Full_azimuth` | [11, 9] | grados | Tabla completa $az(i_{az},i_{el})$ — el azimut depende de la elevación. |
| `Full_elevation` | [11, 9] | grados | Tabla completa $el(i_{az},i_{el})$. |
| `Elevation_correction` | [96] | grados | Desplazamiento "diente de sierra" por energía; el signo alterna con la paridad de la energía. |
| `PAS_to_RTN` | **[tiempo, 3, 3]** | adim. | Matriz de rotación instrumento → RTN (§3.4). |
| `start_*`, `nb_*` | [tiempo] | — | Índice de inicio y conteo realmente barridos (subconjunto del *peak-tracking*). |

Como el *peak-tracking* mueve la ventana activa, un registro puede poblar solo un subbloque del arreglo $11\times 9\times 96$; las variables `start_*`/`nb_*` indican qué índices son válidos. Se proveen tanto `DELTA_PLUS_VAR` como `DELTA_MINUS_VAR` para la energía (asimétricos, por el espaciado logarítmico) y para los ángulos (simétricos): son exactamente lo que se necesita para el elemento de volumen del espacio de velocidades (§3, §4). En la implementación, `read_vdf_timestep()` de `src/cdf_io.py` lee el registro más cercano a `TARGET_TIME` y extrae `vdf`, `Energy`, `Azimuth`, `Elevation`, `PAS_to_RTN`, los semianchos `delta_*` y el indicador `Info`.

---

## 3. Reconstrucción del espacio de velocidades

El objetivo de esta sección es convertir cada *bin* instrumental $(i,j,k)$ — etiquetado por $(E,\,\text{elevación},\,\text{azimut})$ — en un punto del espacio de velocidades 3D $(\mathbf v,\,f)$, primero en el marco del instrumento, luego en el marco heliocéntrico RTN. Esta es la "Tarea 1" del proyecto (el diagrama de dispersión de la VDF).

### 3.1 Rapidez a partir de la energía

El analizador selecciona $E/q$; con la variable `Energy` en eV, el módulo de la velocidad es

$$
|\mathbf v|=\sqrt{\frac{2qE}{m}}=\sqrt{\frac{2E_{\rm eV}\,e}{m}} .
$$

Para **protones** ($m=m_p$, $q=e$) esto se reduce a la constante de la Guía de Usuario **`E2V = 13.85`** (Fedorov 2020):

$$
|\mathbf v|\,[\mathrm{km/s}] = 13.85\,\sqrt{E_{\rm eV}},\qquad
\sqrt{2e/m_p} = 1.385\times 10^{4}\ \mathrm{m\,s^{-1}\,eV^{-1/2}} .
$$

(Verificación numérica con las constantes CODATA del proyecto: $\sqrt{2e/m_p}=1.3841\times10^4\ \mathrm{m\,s^{-1}\,eV^{-1/2}}$, es decir $13.84\ \mathrm{km\,s^{-1}\,eV^{-1/2}}$.) En el código, `speed_from_energy(E_eV)` de `src/physics.py` implementa exactamente $|\mathbf v|=\sqrt{2qE/m}$ en SI.

Los límites inferior/superior de rapidez por *bin* usan los semianchos de energía:

$$
v_{\min}[i_e]=13.85\sqrt{E[i_e]-\texttt{delta\_m\_Energy}},\qquad
v_{\max}[i_e]=13.85\sqrt{E[i_e]+\texttt{delta\_p\_Energy}} .
$$

**Nota sobre alfas.** Para una partícula alfa, el mismo escalón $E/q$ implica $q=2e$, $m=4m_p$, de modo que $|\mathbf v|_\alpha=|\mathbf v|_p/\sqrt2$. Como PAS no separa por masa, la conversión por defecto supone protones; este desplazamiento $m/q=2$ es justamente lo que más tarde permite separar espectralmente las alfas (§5.3).

### 3.2 De esférico a cartesiano: dirección de llegada frente a velocidad de la partícula

Una partícula contada en el *bin* $(i,j,k)$ **llega desde** la dirección de mirada $\hat{\mathbf n}_{jk}$; por lo tanto su velocidad apunta en sentido **opuesto**. Con la elevación medida desde el plano $x$–$y$ del instrumento y el azimut en ese plano, la dirección de mirada (todo positivo, convención esférica estándar) es

$$
\hat L_x=\cos(el)\cos(az),\quad \hat L_y=\cos(el)\sin(az),\quad \hat L_z=\sin(el),
$$

que es la dirección **de la cual** llegó el ion, **no** su velocidad. La velocidad física es $\mathbf v=-|\mathbf v|\,\hat L$.

**La convención SWA (con los signos menos).** La Guía de Usuario de PAS (Fedorov 2020, §2.1) absorbe el signo menos directamente en los versores. En el marco del analizador PAS:

$$
\hat V^{\rm PAS}_X=-\cos(El)\cos(Az),\qquad
\hat V^{\rm PAS}_Y=-\cos(El)\sin(Az),\qquad
\hat V^{\rm PAS}_Z=+\sin(El).
$$

Esto corresponde a los signos $(-1,-1,+1)$, que es exactamente lo que usa `instrument_velocity(E_eV, azim_deg, elev_deg, signs=(-1.0,-1.0,+1.0))` en `src/physics.py`:

```
vx = sx * |v| * cos(el) * cos(az)
vy = sy * |v| * cos(el) * sin(az)
vz = sz * |v| * sin(el)        con (sx, sy, sz) = (-1, -1, +1)
```

**Por qué los signos menos son obligatorios.** El instrumento mide la dirección de **llegada**; la velocidad de la partícula es la opuesta. Para un haz nominal ($Az\approx 0$, $El\approx 0$) la convención SWA da $\hat V^{\rm PAS}\approx(-1,0,0)$, es decir velocidad iónica apuntando **anti-solar** (hacia $-X$), como exige físicamente un viento solar que fluye hacia afuera. Combinado con la matriz `PAS_to_RTN` exacta por registro, la componente radial sale $V_R\approx+|v|>0$ (p. ej. $+300$–$400$ km/s). La fórmula "todo positivo" devolvería en cambio la dirección de llegada (hacia el Sol) y un $V_R$ **negativo espurio**. El proyecto valida esto explícitamente: `run_analysis.py` calcula un $u_R$ de prueba con los signos $(-1,-1,+1)$ y, solo si saliera $u_R<0$, recurre a la negación global $(+1,+1,-1)$; en este caso de estudio los signos confirmados son $(-1,-1,+1)$ (véase `sign` en `data/processed/moments.json`).

**Construcción por esquinas del *bin* (para momentos precisos).** La Guía de Usuario construye cada celda del espacio de velocidades a partir de los ángulos de esquina `Full_azimuth`/`Full_elevation` (mín, centro, máx por $(i_{az},i_{el})$) y la corrección por energía `Elevation_correction`:

```
x_i = -cos(elArr_i) * cos(azArr_i)
y_i =  cos(elArr_i) * sin(azArr_i)
z_i = -sin(elArr_i)
si ((ie - se) % 2 == 0):  z_i += Elevation_correction[ie]
si no:                    z_i -= Elevation_correction[ie]
P_i = [x_i, y_i, z_i]   # luego escalar por v_min / v_max para las esquinas de la celda
```

La `Elevation_correction` dependiente de paridad produce las características fronteras de *bin* en "diente de sierra" en el plano $V_x$–$V_z$. (La implementación de este proyecto usa centros de *bin* con los semianchos `delta_*` para los volúmenes, suficiente para los momentos de orden $\le 2$; la construcción por esquinas es la versión de máxima fidelidad.)

### 3.3 El jacobiano del espacio de velocidades

El mapeo esférico → cartesiano $(v,\theta,\phi)\to(v_x,v_y,v_z)$ arrastra el elemento de volumen

$$
\boxed{\,d^3v=v^2\,dv\,d\Omega\,},\qquad d\Omega=\cos(el)\,d(el)\,d(az)
$$

(convención de elevación; con un ángulo polar $\vartheta$ desde el eje $z$ sería $d\Omega=\sin\vartheta\,d\vartheta\,d\phi$). El factor $v^2$ es el jacobiano radial y $\cos(el)$ el jacobiano angular; **ambos son esenciales** y constituyen una fuente frecuente de error si se omiten (Wilson et al. 2025). El módulo de velocidad por *bin* sigue de $E=\tfrac12 mv^2 \Rightarrow dE=mv\,dv$, de donde

$$
dv=\frac{dE}{m\,v}=\sqrt{\frac{q}{2mE}}\;dE .
$$

En el código, `velocity_volume(E_eV, dE_eV, elev_deg, d_elev_deg, d_azim_deg)` de `src/physics.py` calcula precisamente

$$
dV = v^2\;\underbrace{\sqrt{\tfrac{q}{2mE}}\,dE}_{dv}\;\underbrace{\cos(el)\,\Delta(el)\,\Delta(az)}_{d\Omega},
$$

con los anchos completos de *bin* ($\Delta E = $ `delta_p_E` + `delta_m_E`, $\Delta el = 2\,$`delta_El`, $\Delta az = 2\,$`delta_Az`) tomados directamente del CDF.

### 3.4 Rotación PAS → RTN

Cada velocidad cartesiana se rota del marco del instrumento al marco heliocéntrico **RTN** (R radial hacia afuera, T aproximadamente a lo largo del movimiento orbital, N completando la tríada dextrógira). Para Solar Orbiter esto se encapsula en la matriz `PAS_to_RTN` por registro, $M$, de dimensión [tiempo, 3, 3], con metadatos `COORDINATE_SYSTEM = SOLO_SWA_PAS`, `TARGET_SYSTEM = SOLO_SUN_RTN`:

$$
\begin{pmatrix}v_R\\v_T\\v_N\end{pmatrix}
= M\begin{pmatrix}v_X^{\rm PAS}\\v_Y^{\rm PAS}\\v_Z^{\rm PAS}\end{pmatrix}.
$$

$M$ es una rotación propia (ortonormal, $\det=+1$), de modo que RTN→PAS $=M^{\mathsf T}$. Es dependiente del tiempo porque la actitud de la nave respecto a RTN varía a lo largo de la órbita; la relación aproximada de signos $X_{\rm RTN}\approx -X_{\rm SRF}$, $Y_{\rm RTN}\approx -Y_{\rm SRF}$ solo se cumple "la mayoría del tiempo", así que se usa `PAS_to_RTN` para la transformación exacta. La función `rotate_pas_to_rtn(vx, vy, vz, M)` de `src/physics.py` aplica $v_a^{\rm RTN}=\sum_b M_{ab}\,v_b^{\rm PAS}$. Como las rotaciones preservan $d^3v$ (ortogonales, $\det=1$), **los momentos son independientes del marco**. Cada *bin* medido $(i,j,k)$ se convierte así en un punto $(\mathbf v^{\rm RTN}_{ijk},\,f_{ijk})$: la nube de dispersión 3D del espacio de velocidades $(v_R,v_T,v_N)$ que se grafica en la Tarea 1.

---

## 4. Marco alineado al campo y momentos del plasma

### 4.1 El campo magnético y el versor $\hat{\mathbf b}$

El campo magnético se toma del producto MAG L2 **`mag-rtn-normal`** de Solar Orbiter. En el código, `read_mag_at(path, target_iso, window_s=4.0)` de `src/cdf_io.py` promedia la variable `B_RTN` en una ventana de $\pm$2 s alrededor del instante objetivo (descartando relleno), de modo que $\mathbf B$ se promedia sobre la misma cadencia de ~4 s de una VDF normal. Se define el versor

$$
\hat{\mathbf b}=\frac{\mathbf B}{|\mathbf B|}.
$$

### 4.2 Construcción de la tríada alineada al campo (FA)

La velocidad **paralela** es la proyección sobre el campo,

$$
v_\parallel=\mathbf v\cdot\hat{\mathbf b}.
$$

El par perpendicular **no es único**: cualquier rotación rígida por un ángulo de girofase $\psi$ alrededor de $\hat{\mathbf b}$ es una elección igualmente válida. Se fija la libertad con un vector de referencia $\hat{\mathbf r}$. En `field_aligned_basis(B)` de `src/physics.py` se elige $\hat{\mathbf r}=\hat R$ (la dirección radial RTN; si $\hat{\mathbf b}$ fuera casi paralelo a $\hat R$, se conmuta a $\hat T$), y se construye

$$
\hat{\mathbf e}_{\perp 1}=\frac{\hat{\mathbf b}\times\hat{\mathbf r}}{|\hat{\mathbf b}\times\hat{\mathbf r}|},\qquad
\hat{\mathbf e}_{\perp 2}=\frac{\hat{\mathbf b}\times\hat{\mathbf e}_{\perp 1}}{|\hat{\mathbf b}\times\hat{\mathbf e}_{\perp 1}|},
$$

de forma que $(\hat{\mathbf e}_{\perp 1},\hat{\mathbf e}_{\perp 2},\hat{\mathbf b})$ es dextrógira. Las proyecciones son $v_{\perp 1}=\mathbf v\cdot\hat{\mathbf e}_{\perp 1}$, $v_{\perp 2}=\mathbf v\cdot\hat{\mathbf e}_{\perp 2}$, y la girofase es $\psi=\arctan2(v_{\perp 2},v_{\perp 1})$. La función `project_field_aligned()` realiza estas tres proyecciones. Una sutileza de implementación: `run_analysis.py` fija el signo de $\hat{\mathbf b}$ "a lo largo del flujo" (exigiendo $u_\parallel>0$) antes de construir la tríada, de modo que el haz y el núcleo se separen consistentemente.

### 4.3 Girotropía: la reducción 3D → 2D

En un plasma magnetizado de $\beta$ bajo, el campo de fondo intenso es un **eje de simetría**: $f$ es invariante bajo rotación alrededor de $\hat{\mathbf b}$ (independiente de la girofase $\psi$). Esta **girotropía** permite definir la magnitud perpendicular y reducir la distribución de 3D a 2D:

$$
\boxed{\,v_\perp=\sqrt{v_{\perp 1}^2+v_{\perp 2}^2}\,},\qquad
f(v_\parallel,v_\perp)=\frac{1}{2\pi}\int_0^{2\pi}f(v_\parallel,v_\perp,\psi)\,d\psi .
$$

Esta **función de distribución girotrópica (GDF)** 2D, $f(v_\parallel,v_\perp)$, es el objeto que el resto del análisis descompone. Es también el marco natural para resolver haces de protones alineados al campo, que aparecen como una segunda población desplazada en $v_\parallel$.

**Cuándo falla la girotropía (agirotropía).** Cerca de choques, sitios de reconexión y regiones de difusión, la giración se interrumpe en una giroperiodo y $f$ depende de $\psi$; el promedio de girofase descarta entonces estructura real. La desviación se cuantifica con escalares construidos a partir del tensor de presión en el marco FA, el más usado de los cuales es la $Q$ de Swisdak:

$$
\sqrt{Q}=\sqrt{\frac{P_{12}^2+P_{13}^2+P_{23}^2}{P_\perp^2+2P_\perp P_\parallel}}\in[0,1],
$$

con el eje "1" a lo largo de $\hat{\mathbf b}$ (Swisdak 2016). Diagnósticos relacionados son la $A\!\varnothing_e$ de Scudder & Daughton (autovalores del subbloque perpendicular $2\times2$) y la $D_{\rm ng}$ de Aunai et al. (RMS de los $P_{ij}$ fuera de la diagonal, normalizado a la energía térmica); todas miden la desviación de $\mathsf P$ respecto de la forma diagonal girotrópica $\mathrm{diag}(P_\parallel,P_\perp,P_\perp)$.

### 4.4 Momentos del plasma a partir de la VDF discreta

Reemplazando $d^3v\to v_i^2\,\Delta v_i\,\Delta\Omega_{jk}$ en cada *bin* medido, los momentos discretos son (la "receta" SWA-PAS de Fedorov §3):

$$
\boxed{\,n=\sum_{ijk} f_{ijk}\,v_i^2\,\Delta v_i\,\cos(el_j)\,\Delta el_j\,\Delta az_k\,},
$$
$$
n\,\mathbf u=\sum_{ijk}\mathbf v_{ijk}\,f_{ijk}\,v_i^2\,\Delta v_i\,\cos(el_j)\,\Delta el_j\,\Delta az_k,
$$
$$
\mathsf P=m\sum_{ijk}(\mathbf v_{ijk}-\mathbf u)(\mathbf v_{ijk}-\mathbf u)\,f_{ijk}\,v_i^2\,\Delta v_i\,\cos(el_j)\,\Delta el_j\,\Delta az_k .
$$

En el código, `compute_moments(vR, vT, vN, f, dV, b)` de `src/physics.py` introduce un **peso por *bin*** $w_{ijk}=f_{ijk}\,dV_{ijk}$ (con $dV$ ya conteniendo $v^2\,dv\,d\Omega$) y calcula:

$$
n=\sum w_{ijk},\qquad
\mathbf u=\frac1n\sum \mathbf v_{ijk}\,w_{ijk},\qquad
P_{ab}=m\sum (\mathbf v-\mathbf u)_a(\mathbf v-\mathbf u)_b\,w_{ijk},
$$

con $\mathbf c=\mathbf v-\mathbf u$ la **velocidad peculiar**. El tensor de temperatura es $\mathsf T=\mathsf P/(n k_B)$ en kelvin.

**Temperaturas paralela y perpendicular.** Proyectando el tensor sobre el campo:

$$
\boxed{\,T_\parallel=\hat{\mathbf b}\cdot\mathsf T\cdot\hat{\mathbf b}\,},\qquad
\boxed{\,T_\perp=\frac{\operatorname{Tr}\mathsf T-T_\parallel}{2}\,},
$$

con temperatura escalar $T=\tfrac13\operatorname{Tr}\mathsf T=\tfrac13(T_\parallel+2T_\perp)$, anisotropía $A=T_\perp/T_\parallel$ y beta del plasma $\beta_s=8\pi n_s k_B T_s/B^2$. En el código, `T_par = b @ T @ b` y `T_perp = 0.5*(trace(T) - T_par)`.

**Velocidades térmicas.** Con la convención de **rapidez más probable** estándar en cinética del viento solar,

$$
\boxed{\,w=\sqrt{\frac{2k_B T}{m}}\,},\qquad
w_\parallel=\sqrt{\frac{2k_B T_\parallel}{m}},\qquad
w_\perp=\sqrt{\frac{2k_B T_\perp}{m}} .
$$

(Coexisten convenciones rivales $\sqrt{k_BT/m}$, $\sqrt{3k_BT/m}$, $\sqrt{8k_BT/\pi m}$; la forma $\sqrt{2k_BT/m}$ es la que hace que el exponente Maxwelliano sea exactamente $-c^2/w^2$.) En el código, `w_par`, `w_perp` siguen esta definición.

**Advertencias que sesgan las sumas** (Fedorov §3; Wilson et al. 2025): campo de visión truncado (cobertura parcial del cielo), cortes finitos de energía baja/alta, desplazamiento de $E$ por el potencial de la nave, ruido de "una cuenta". La VDF es ruidosa entre 300–400 eV (el factor geométrico decae) y **poco fiable por debajo de 300 eV**; aparecen valores "fantasma" esporádicos en *bins* angulares extremos. En el producto `swa-pas-grnd-mom` conviene usar `validity` $\ge 2$.

### 4.5 Coordenadas normalizadas y velocidad peculiar

El argumento de la base de Hermite (§6) es la **velocidad peculiar normalizada**:

$$
\boxed{\;\xi=\frac{v-u}{w}\;},\qquad
\xi_\parallel=\frac{v_\parallel-u_\parallel}{w_\parallel},\qquad
\xi_\perp=\frac{v_\perp-u_\perp}{w_\perp}.
$$

La forma más limpia y robusta, que evita tener que definir $u_\perp$ o una girofase, es la **velocidad peculiar en el marco de reposo del flujo**: se resta el flujo masivo **antes** de proyectar,

$$
\mathbf c=\mathbf v-\mathbf u,\qquad
c_\parallel=\mathbf c\cdot\hat{\mathbf b},\qquad
\boxed{\,c_\perp=\big|\mathbf c-c_\parallel\hat{\mathbf b}\big|=\sqrt{|\mathbf c|^2-c_\parallel^2}\,}.
$$

Aquí $c_\perp\ge 0$ por construcción y es invariante de rotación alrededor de $\hat{\mathbf b}$ — la cantidad correcta para histogramar una $f(c_\parallel,c_\perp)$ girotrópica. Las coordenadas normalizadas son entonces $\xi_\parallel=c_\parallel/w_\parallel$ y $\xi_\perp=c_\perp/w_\perp$. Exactamente esto hace `normalized_coords(vR, vT, vN, u, b, w_par, w_perp)` de `src/physics.py`. Como la base de Hermite está **centrada en la Maxwelliana local** (corrida al marco de flujo y escalada por la velocidad térmica), las fluctuaciones a gran escala de densidad, flujo y temperatura **no contaminan** el espectro de orden alto.

---

## 5. Modelos paramétricos y sesgo humano

Antes de la descomposición *sin modelo* (§6–§9), conviene entender los modelos paramétricos clásicos, sus aciertos, y por qué su sesgo intrínseco motiva un método libre de modelo.

### 5.1 Maxwelliana y bi-Maxwelliana

La **Maxwelliana** con deriva es la solución de máxima entropía / equilibrio térmico:

$$
f_M(\mathbf v)=\frac{n}{\pi^{3/2}\,w^3}\exp\!\Big[-\frac{(\mathbf v-\mathbf u)^2}{w^2}\Big],\qquad w=\sqrt{\tfrac{2k_BT}{m}} .
$$

Su única "estructura" es densidad, deriva y una temperatura. En un plasma magnetizado débilmente colisional, la giromoción termaliza las dos direcciones perpendiculares pero $T_\parallel\neq T_\perp$ en general; la generalización mínima fuera del equilibrio es la **bi-Maxwelliana**:

$$
f_{bM}(v_\parallel,v_\perp)=\frac{n}{\pi^{3/2}\,w_\parallel w_\perp^{2}}
\exp\!\Big[-\frac{(v_\parallel-u_\parallel)^2}{w_\parallel^{2}}\Big]
\exp\!\Big[-\frac{(v_\perp-u_\perp)^2}{w_\perp^{2}}\Big],
$$

con $w_{\parallel,\perp}=\sqrt{2k_BT_{\parallel,\perp}/m}$ y anisotropía $A=T_\perp/T_\parallel$ (para un plasma girotrópico $u_\perp=0$). La bi-Maxwelliana introduce anisotropía de temperatura respecto a $\mathbf B$ pero **sigue siendo Maxwelliana en cualquier corte 1D**. En coordenadas normalizadas se escribe de forma compacta $f=\tfrac{n}{\pi^{3/2}w_\parallel w_\perp^2}\exp(-\xi_\parallel^2-\xi_\perp^2)$.

### 5.2 Las tres poblaciones a escala de protones

Las VDFs de protones del viento solar están bien organizadas respecto a $\mathbf B$ y se modelan típicamente como **suma de bi-Maxwellianas**:

- **Núcleo de protones (core)** — el pico global, $\sim 90$–95% de los protones, casi Maxwelliano (an)isótropo; define la velocidad masiva; es el más denso y "frío".
- **Haz de protones (beam)** — "una segunda componente que fluye más rápido que el núcleo a lo largo de $\mathbf B$ con una rapidez relativa $\gtrsim v_{Ap}$", deriva $\Delta u\sim(1$–$2)\,v_{Ap}$, "casi siempre alejándose del Sol". Se modela con una bi-Maxwelliana de haz propia ($n_b,u_b,w_{\parallel b},w_{\perp b}$).
- **Partículas alfa ($\mathrm{He^{2+}}$)** — especie distinta ($m_\alpha=4m_p$, $q=2e$, $m/q=2$), $\lesssim 20\%$ de la densidad de masa, derivando a lo largo de $\mathbf B$ alejándose del Sol a $\lesssim v_{Ap}$. Como los analizadores ordenan por $E/q=\tfrac12 m v^2/q$, una alfa a la *misma rapidez masiva* que el núcleo aparece al **doble del $E/q$ del protón**; si se (mal)convierte suponiendo masa de protón, aparece a $\sqrt2\times$ la rapidez del protón. Este desplazamiento $m/q=2$ es exactamente lo que permite separarlas espectralmente.

### 5.3 Distribución kappa y colas supratérmicas

Los plasmas colisionless muestran casi universalmente **colas de ley de potencias supratérmicas** que una Maxwelliana no puede ajustar. Vasyliūnas (1968) y Olbert (1968) introdujeron la **distribución kappa**:

$$
f_\kappa(\mathbf v)=\frac{n}{[\pi(2\kappa-3)]^{3/2}\,w^{3}}\,
\frac{\Gamma(\kappa+1)}{\Gamma(\kappa-\tfrac12)}\,
\Big[1+\frac{2(\mathbf v-\mathbf u)^2}{(2\kappa-3)\,w^{2}}\Big]^{-(\kappa+1)},
$$

donde el factor $(2\kappa-3)$ mantiene a $w$ como velocidad térmica con $T=mw^2/2k_B$ **independiente de $\kappa$** (para $\kappa>3/2$). El índice $\kappa$ controla las colas: más pronunciadas para $\kappa$ menor, con pendiente de alta energía $\propto v^{-2(\kappa+1)}$, y $f_\kappa\to f_M$ cuando $\kappa\to\infty$. Un $\kappa$ menor correlaciona con colisiones más débiles / mayor distancia heliocéntrica. Para **electrones**, la descomposición estándar es **núcleo Maxwelliano + halo kappa** (más un *strahl* alineado al campo como tercera componente). La kappa estándar tiene momentos divergentes para $\kappa$ bajo; la **kappa regularizada (RKD)** de Scherer, Fichtner & Lazar (2017) multiplica por un corte gaussiano $\exp(-\alpha^2 v^2/\theta^2)$, $0<\alpha\ll1$, manteniendo todos los momentos finitos.

### 5.4 Ajuste por $\chi^2$, GMM, y el sesgo humano

Los ajustes paramétricos minimizan un $\chi^2$ ponderado entre el modelo y la $f$ medida sobre los *bins* del espacio de velocidades,

$$
\chi^2(\boldsymbol\theta)=\sum_i\frac{[f_{\rm obs}(\mathbf v_i)-f(\mathbf v_i;\boldsymbol\theta)]^2}{\sigma_i^2},
$$

usualmente por **Levenberg-Marquardt (LM)**, que interpola entre Gauss-Newton y descenso de gradiente vía un parámetro de amortiguamiento $\lambda$:

$$
\big(\mathbf J^{\mathsf T}\mathbf J+\lambda\,\mathrm{diag}(\mathbf J^{\mathsf T}\mathbf J)\big)\,\delta\boldsymbol\theta=\mathbf J^{\mathsf T}\mathbf r,
$$

(con $\mathbf J$ el jacobiano de los residuos $\mathbf r$; implementación MPFIT, Markwardt 2009). Una alternativa no supervisada es un **Modelo de Mezcla Gaussiana (GMM)** ajustado por **Expectación-Maximización (EM)** (De Marco et al. 2023, sobre Solar Orbiter/PAS), que modela cada punto de velocidad como mezcla de $K$ gaussianas,

$$
P(\mathbf x_i\mid\boldsymbol\theta)=\sum_{k=1}^{K}w_k\,\mathcal N(\mathbf x_i\mid\boldsymbol\mu_k,\boldsymbol\Sigma_k),\qquad \sum_k w_k=1,
$$

asignando a cada punto una membresía probabilística (*responsibility*) en lugar de una categoría forzada, lo que elimina el sesgo de inicialización manual y enmascaramiento de regiones de LM (aunque sigue siendo de base gaussiana).

**El sesgo humano.** Todo ajuste paramétrico incrusta *a priori* que sesgan el resultado **antes de tocar los datos**: el número de componentes (1 vs 2 vs 3 Maxwellianas, Maxwelliana vs kappa) predetermina si existe un "haz" o una "cola"; las bi-Maxwellianas y kappas girotrópicas asumen girotropía y simetría de reflexión, descartando asimetría (*skewness*), agirotropía y distorsiones de flujo de calor que están físicamente presentes (el trabajo de la *normal inverse-Gaussian*, A&A 2024, muestra que las VDFs reales llevan *skewness*/*kurtosis* que ningún modelo simétrico captura). Por tanto los $n,\mathbf u,T_{\parallel,\perp},\kappa$ recuperados dependen del prior del modelador.

**Por qué es deseable un método sin modelo.** Esto motiva una **expansión en base** completa y físicamente natural que no imponga conteo de poblaciones, ni hipótesis de equilibrio, ni simetría — dejando que los datos pueblen los coeficientes. La descomposición de Hermite/Laguerre de las siguientes secciones es justo eso: (i) no está sesgada por conteo de poblaciones ni hipótesis de equilibrio, (ii) conserva *toda* la información, (iii) produce una "distancia al equilibrio Maxwelliano" cuantitativa, y (iv) conecta directamente con la teoría cinética (mezcla de fases, amortiguamiento de Landau, cascada en el espacio de velocidades).

---

## 6. Descomposición polinomial: la analogía de Fourier

### 6.1 La analogía de Fourier

Así como una transformada de Fourier escribe una señal como suma de ondas planas, la **transformada de Hermite escribe la VDF como una Maxwelliana multiplicada por una serie de polinomios de Hermite.** Larosa et al. (2025) lo enuncian explícitamente: *"un solo modo de Fourier describe una onda plana, y análogamente un solo modo de Hermite describe una Maxwelliana; cuantos más modos de Hermite sean no nulos, más lejos está la VDF de una Maxwelliana."* La expansión es exacta e invertible (una base espectral) y **se cumple Parseval-Plancherel**, de modo que la "energía" del espacio de velocidades (la enstrofía) se conserva entre el espacio físico y el espectral.

### 6.2 Las funciones de Hermite ortonormales

En coordenadas alineadas al campo $(\xi_\parallel,\xi_\perp)$ la VDF girotrópica se expande como

$$
f(\xi_\parallel,\xi_\perp)=\sum_{m,n}c_{mn}\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp),
$$

con las **funciones de Hermite ortonormales**

$$
\boxed{\;\psi_m(\xi)=\frac{H_m(\xi)}{\sqrt{2^m\,m!\,\sqrt\pi}}\;e^{-\xi^2/2}\;}
$$

donde $H_m$ son los polinomios de Hermite *de los físicos*, $H_m(\xi)=(-1)^m e^{\xi^2}\frac{d^m}{d\xi^m}e^{-\xi^2}$. Satisfacen la ortonormalidad "tipo Fourier"

$$
\int_{-\infty}^{\infty}\psi_m(\xi)\,\psi_n(\xi)\,d\xi=\delta_{mn} .
$$

Obsérvese que el peso es **1**, no $e^{-\xi^2}$: el envolvente $e^{-\xi^2/2}$ hace ortonormales a las *funciones* (no solo a los polinomios). Son las funciones de Hermite "simétricamente ponderadas", esto es, los autoestados del oscilador armónico cuántico.

Elecciones de normalización clave:
- El argumento es la **velocidad peculiar normalizada** $\xi=(v-u)/w$ — la base está *corrida al marco de flujo local* y *escalada por la velocidad térmica* en cada dirección.
- Por tanto $m=0$ es la Maxwelliana con deriva *local*; la proyección "cuantifica las correcciones de orden alto a la DF, ya que la base está corrida al marco de velocidad de fluido local y normalizada a la densidad y temperatura ambientes" (Pezzi). Consecuentemente, las fluctuaciones a gran escala de densidad, flujo masivo y temperatura (momentos 1° y 2°) **no contaminan el espectro de orden alto**.

En el código, `hermite_functions(x, mmax)` de `src/hermite.py` evalúa las $\psi_m$ con la **recurrencia estable** de las funciones (acotadas, $|\psi_m|\sim\pi^{-1/4}$), evitando el desbordamiento que sufrirían los $H_m$ crudos:

$$
\psi_0=\pi^{-1/4}e^{-\xi^2/2},\quad
\psi_1=\sqrt2\,\xi\,\psi_0,\quad
\psi_{m+1}=\sqrt{\tfrac{2}{m+1}}\,\xi\,\psi_m-\sqrt{\tfrac{m}{m+1}}\,\psi_{m-1}.
$$

### 6.3 Los coeficientes $c_{mn}$ y su significado físico

El coeficiente espectral es la integral de proyección

$$
\boxed{\;c_{mn}=\iint f(\xi_\parallel,\xi_\perp)\,\psi_m(\xi_\parallel)\,\psi_n(\xi_\perp)\,d\xi_\parallel\,d\xi_\perp\;}
$$

Como la base está centrada en la Maxwelliana local, **toda la física vive en $m,n\ge 1$**, organizada en una jerarquía que refleja los momentos:

| Índice | Contenido físico |
|---|---|
| $c_{00}$ | Amplitud Maxwelliana $\sim$ **densidad** (Maxwelliana con deriva local); referencia de normalización. |
| $m=1$ | Fluctuaciones de **flujo masivo / deriva**. |
| $m=2$ | Deformaciones de **temperatura** (anisotropía de presión). |
| $m=3$ | Perturbaciones de **flujo de calor**. |
| $m,n$ bajos | Momentos fluidos → "el acoplamiento fluido-cinético es una característica intrínseca". |
| $m,n$ altos | **Estructura no térmica de escala fina** en el espacio de velocidades (haces, *hammerheads*, filamentación de fase). |

En palabras de Pezzi: *"el índice $m$ corresponde aproximadamente a un orden de los momentos de velocidad: $m=1$ a fluctuaciones de flujo masivo; $m=2$ a deformaciones de temperatura; $m=3$ a perturbaciones de flujo de calor, y así sucesivamente."*

**Enstrofía / energía libre.** La suma de los coeficientes al cuadrado da la enstrofía del espacio de velocidades, la cantidad que cascadea,

$$
\Omega=\int \delta f^2\,d^3v=\sum_{m>0}\big[c_m\big]^2,
$$

ligada a un indicador de no-Maxwellianidad vía $\Omega=\epsilon^2 n^2$, y que es "esencialmente la energía libre en girocinética" (Pezzi, citando a Schekochihin).

### 6.4 Dos codificaciones del plano perpendicular: HH frente a HL

La girotropía implica $f=f(v_\parallel,v_\perp)$ con $v_\perp=\sqrt{v_{\perp1}^2+v_{\perp2}^2}$. Hay dos codificaciones espectrales naturales del plano perpendicular — la elección metodológica central.

**(a) Hermite-Hermite (HH), girotropía "tipo cartesiano" — Larosa et al. 2025 (PSP).** Se trata $v_\perp$ como un eje cartesiano y se usa $\psi_n(\xi_\perp)$ también en perpendicular. La girotropía se impone por **reflexión (*mirroring*)**: *"imponemos $f(-v_\perp)=f(v_\perp)$ y extendemos la malla a valores negativos en la dirección perpendicular… tal procedimiento implica coeficientes de Hermite impares nulos en la dirección perpendicular"*. Así $c_{m,n}=0$ para $n$ impar — **solo los modos perpendiculares pares llevan potencia**, lo cual es la causa del patrón oscilatorio par/impar en el espectro perpendicular 1D. **Esta es la base usada en este proyecto** (HH, malla $60\times60$ vía `numpy.polynomial.hermite.hermgauss`). En `run_analysis.py`, la reflexión se realiza concatenando los puntos con $\xi_\perp$ y $-\xi_\perp$ antes de interpolar:

```
src_xi_par  = concat([xipar,  xipar])
src_xi_perp = concat([xiperp, -xiperp])   # refleja la VDF perpendicular
src_logf    = concat([logf,   logf])
```

**(b) Hermite-Laguerre (HL), girotropía cilíndrica — Coburn et al. 2024 (electrones de Solar Orbiter).** La base *natural* para una función azimutalmente simétrica usa **polinomios de Laguerre asociados en $\mu=v_\perp^2$**:

$$
F_e(v_\parallel,v_\perp)=\sum_{m}\sum_{l}c_{ml}\,\psi_m(\hat v_\parallel)\,\Gamma_l^{1}(\hat v_\perp),
\qquad
\Gamma_l^{n}(\hat v_\perp)=e^{-\hat v_\perp^2/2}\,L_l^{n}(\hat v_\perp^2),
$$

con índice de peso $n=1$. El peso $n=1$ absorbe el jacobiano cilíndrico $v_\perp\,dv_\perp$, de modo que las $\{\Gamma_l^1\}$ son ortonormales bajo la medida girotrópica. Esta es la base del espacio de velocidades de la teoría girocinética (formulación pseudo-espectral Laguerre-Hermite; Mandell et al. 2018).

**Por qué importa.** HL es la representación *físicamente fiel* de una distribución girotrópica: un anillo/cáscara o una Maxwelliana perpendicular es un único modo de Laguerre de orden bajo ($\Gamma_0^1$ exactamente), sin modos impares espurios. HH es más simple/uniforme (misma base en ambos ejes) pero paga el precio de modos impares-$n$ nulos y un espectro oscilatorio, y una Maxwelliana perpendicular requiere muchos modos de Hermite pares.

---

## 7. Cuadratura gaussiana

La integral de proyección $c_{mn}$ se evalúa numéricamente por **cuadratura de Gauss-Hermite**. Esta sección desarrolla por qué los nodos deben ser las raíces del polinomio ortogonal, cómo se calculan, y cómo se implementa.

### 7.1 Planteamiento y las dos familias relevantes

Se aproxima una integral ponderada por una regla de $N$ puntos,

$$
\int_a^b p(x)\,f(x)\,dx\;\approx\;\sum_{i=1}^N w_i\,f(x_i),
$$

con peso fijo $p(x)\ge 0$, nodos $\{x_i\}$ y pesos $\{w_i\}$. Las dos familias relevantes para una VDF son:

| Familia | Peso $p(x)$ | Intervalo | Polinomio ortogonal |
|---|---|---|---|
| Gauss-**Hermite** | $e^{-x^2}$ | $(-\infty,\infty)$ | $H_N$ |
| Gauss-**Laguerre** | $e^{-x}$ | $(0,\infty)$ | $L_N$ |

Una regla de $N$ puntos tiene $2N$ parámetros libres ($N$ nodos + $N$ pesos), de modo que se puede aspirar a que sea exacta para todo polinomio de grado $\le 2N-1$. La cuadratura gaussiana logra precisamente este orden máximo **cuando los nodos son las $N$ raíces de $P_N$.**

### 7.2 Por qué los nodos deben ser las raíces de $P_N$

Sea $\{P_k\}$ la familia ortogonal respecto a $p$: $\langle P_j,P_k\rangle=\int_a^b p\,P_jP_k\,dx=h_k\delta_{jk}$. Tómense los nodos $x_1,\dots,x_N$ como las $N$ raíces (reales, simples, interiores) de $P_N$. Para cualquier polinomio $f$ de grado $\le 2N-1$, dividir entre $P_N$:

$$
f(x)=q(x)P_N(x)+r(x),\qquad \deg q\le N-1,\ \deg r\le N-1.
$$

Integrando contra $p$, el primer término **se anula por ortogonalidad** ($q$ es combinación de $P_0,\dots,P_{N-1}$, todos ortogonales a $P_N$):

$$
\int_a^b p\,f\,dx=\underbrace{\int_a^b p\,q\,P_N\,dx}_{=0}+\int_a^b p\,r\,dx=\int_a^b p\,r\,dx.
$$

En cada nodo $P_N(x_i)=0\Rightarrow f(x_i)=r(x_i)$, de modo que $\sum_i w_i f(x_i)=\sum_i w_i r(x_i)$. Si los pesos integran exactamente todo polinomio de grado $\le N-1$ (§7.3), entonces $\sum_i w_i r(x_i)=\int_a^b p\,r\,dx$, y por tanto

$$
\sum_{i=1}^N w_i f(x_i)=\int_a^b p\,f\,dx .
$$

La elección $x_i=$ raíces de $P_N$ es exactamente lo que hace que el resto de alto grado $q\,P_N$ integre a cero; **cualquier otro conjunto de nodos pierde los $N$ órdenes adicionales**. Las raíces están garantizadas reales, simples e interiores porque $P_N$ es ortogonal respecto a un peso positivo.

### 7.3 La fórmula del peso

Exigiendo exactitud hasta grado $N-1$ vía la base de Lagrange $\ell_i(x)=\prod_{j\ne i}\frac{x-x_j}{x_i-x_j}$ (con $\ell_i(x_j)=\delta_{ij}$):

$$
\boxed{\,w_i=\int_a^b p(x)\,\ell_i(x)\,dx\,}.
$$

Usando las raíces de $P_N$, $\ell_i(x)=\frac{P_N(x)}{P_N'(x_i)(x-x_i)}$, y la identidad de Christoffel-Darboux colapsa la integral a la forma de los **números de Christoffel**:

$$
\boxed{\,w_i=\frac{1}{\sum_{k=0}^{N-1}\hat P_k(x_i)^2}\,}.
$$

Especializando (Abramowitz & Stegun, Cap. 25):

$$
\textbf{Hermite (25.4.46):}\quad w_i=\frac{2^{\,N-1}\,N!\,\sqrt\pi}{N^2\,[H_{N-1}(x_i)]^2},
\qquad
\textbf{Laguerre (25.4.45):}\quad w_i=\frac{x_i}{(N+1)^2\,[L_{N+1}(x_i)]^2}.
$$

### 7.4 Cálculo de los nodos: Newton-Raphson frente a Golub-Welsch

**Newton-Raphson** resuelve $P_N(x_i)=0$ directamente, iterando $x^{(k+1)}=x^{(k)}-P_N(x^{(k)})/P_N'(x^{(k)})$ vía la recurrencia de tres términos. Es **frágil**: las raíces se agrupan, $P_N$ tiene un rango dinámico enorme, y sin buenas semillas asintóticas Newton puede converger a una raíz ya hallada, **saltarse una raíz**, o desbordarse; requiere deflación y siembra cuidadosas.

**Golub-Welsch (1969)** es el método robusto moderno. Los polinomios ortonormales obedecen la recurrencia simétrica de tres términos

$$
x\,P_m(x)=a_{m+1}P_{m+1}(x)+b_m P_m(x)+a_m P_{m-1}(x),
$$

que se apila en $x\,\mathbf P(x)=J\,\mathbf P(x)+a_N P_N(x)\,\mathbf e_N$ con la **matriz de Jacobi tridiagonal simétrica**

$$
J=\begin{pmatrix}
b_0 & a_1 & & \\
a_1 & b_1 & a_2 & \\
& a_2 & b_2 & \ddots \\
& & \ddots & \ddots & a_{N-1}\\
& & & a_{N-1} & b_{N-1}
\end{pmatrix}.
$$

En un nodo $x_i$ el último término se anula, de modo que $J\mathbf P(x_i)=x_i\mathbf P(x_i)$:

- **Nodos** $x_i$ = **autovalores** de $J$ (reales, simples, calculados de forma estable por QR);
- **Pesos** $w_i=\mu_0\,(v_{0,i})^2$, con $v_{0,i}$ la primera componente del $i$-ésimo autovector normalizado y $\mu_0=\int_a^b p\,dx$ el momento de orden 0.

Los coeficientes de recurrencia son: **Hermite** $b_m=0$, $a_m=\sqrt{m/2}$, $\mu_0=\sqrt\pi$; **Laguerre** $b_m=2m+1$, $a_m=m$, $\mu_0=1$. La ventaja decisiva: **sin semillas, sin raíces faltantes** — el problema de raíces se convierte en un problema de autovalores simétrico bien condicionado.

### 7.5 En la práctica: `numpy.polynomial.hermite.hermgauss`

En el código, `gauss_hermite(n)` de `src/hermite.py` es un envoltorio de

```python
x, w = np.polynomial.hermite.hermgauss(n)   # peso e^{-x^2} en (-inf, inf)
```

`hermgauss(deg)` devuelve nodos/pesos que "integran correctamente polinomios de grado $2\cdot\text{deg}-1$ o menor" con peso $e^{-x^2}$: internamente construye la matriz *companion* (Jacobi), toma sus autovalores como nodos (refinados con un paso de Newton) y calcula los pesos con la forma cerrada de §7.3 — un Golub-Welsch modificado para $N\lesssim 150$. **Nota de convención:** `hermgauss` es la de los *físicos* ($e^{-x^2}$); `hermegauss` es la de los *probabilistas* ($e^{-x^2/2}$) — hay que casar la convención con la base. Para la malla $60\times60$ de este proyecto, $N=60$ en cada dirección (`N_NODES = 60` en `src/config.py`). El auto-test de `hermite.py` verifica que $\int e^{-x^2}x^2\,dx=\sqrt\pi/2$ y la ortonormalidad discreta $\sum_i (w_i e^{x_i^2})\psi_m(x_i)\psi_n(x_i)=\delta_{mn}$.

---

## 8. Interpolación sobre la malla de nodos

### 8.1 El problema: datos dispersos frente a una malla fija

La VDF medida vive en **puntos dispersos** $(\xi_\parallel^{(k)},\xi_\perp^{(k)})$ con valores $f_k$ (un punto por *bin* válido del instrumento), pero la cuadratura requiere $f$ **en la malla fija de nodos** $(\xi_{\parallel,i},\xi_{\perp,j})$ — las raíces de Hermite, que casi nunca coinciden con las posiciones de medición. Hay que **interpolar** la nube dispersa sobre la malla de nodos.

### 8.2 Peso inverso a la distancia (IDW / Shepard)

Se usa **interpolación por peso inverso a la distancia (IDW)**:

$$
f(\xi_{\parallel,i},\xi_{\perp,j})=\frac{\sum_k W_k\,f_k}{\sum_k W_k},
\qquad
W_k=\frac{1}{r_k^2+\varepsilon},
\qquad
r_k^2=(\xi_{\parallel,i}-\xi_\parallel^{(k)})^2+(\xi_{\perp,j}-\xi_\perp^{(k)})^2 .
$$

- La forma de **combinación convexa** ($\sum_k W_k f_k/\sum_k W_k$) garantiza que el interpolante quede dentro del rango de los datos y reproduzca constantes exactamente.
- El regularizador $\varepsilon>0$ (una longitud de suavizado) evita la singularidad $1/r^2$ cuando un nodo casi coincide con una medición y controla la suavidad: $\varepsilon$ pequeño → agudo, casi vecino más cercano; $\varepsilon$ mayor → más suave pero más difusivo. Conviene elegir $\varepsilon$ comparable al cuadrado del espaciado local de muestras.
- **Este error de interpolación es el piso de error dominante de toda la tubería** (§9.4): es independiente de cuántos nodos de cuadratura se añadan, y es lo que hace que el espectro se aplane más allá de $N\approx 60$.

En el código, `idw_interpolate(xs, ys, vals, grid_x, grid_y, eps=1e-3, power=2)` de `src/hermite.py` implementa exactamente esto, procesando por bloques de 4096 puntos para no agotar memoria (la matriz $N_{\rm malla}\times N_{\rm puntos}$). El parámetro $\varepsilon$ es `IDW_EPS = 1e-3` en `src/config.py`. Un detalle de implementación importante: `run_analysis.py` interpola el **logaritmo** de la VDF ($\log_{10} f$), no $f$ directamente, y luego exponencia de vuelta a lineal antes de proyectar; esto preserva mejor el rango dinámico de varios órdenes de magnitud de la PSD.

### 8.3 Reflexión perpendicular (HH) frente a dominio cilíndrico (HL)

La forma de imponer la girotropía en la interpolación distingue ambos métodos:

- **HH (estilo Larosa, este proyecto):** se refleja el perpendicular, $f(\xi_\perp)=f(-\xi_\perp)$, concatenando los puntos $(\xi_\parallel,\xi_\perp)$ y $(\xi_\parallel,-\xi_\perp)$ antes de interpolar (ver el fragmento en §6.4). La malla de nodos cubre entonces todo el eje perpendicular $(-\infty,\infty)$, y la reflexión garantiza coeficientes impares-$n$ nulos.
- **HL (estilo Coburn):** el dominio perpendicular es cilíndrico $[0,\infty)$ y se usa cuadratura de Gauss-Laguerre en $\mu=\xi_\perp^2$, sin reflexión.

---

## 9. El espectro 2D de Hermite

### 9.1 La fórmula de cuadratura para los coeficientes

La proyección $c_{mn}=\iint f\,\psi_m(\xi_\parallel)\psi_n(\xi_\perp)\,d\xi_\parallel d\xi_\perp$ es una integral **sin peso** (peso 1), pero Gauss-Hermite acarrea el peso $e^{-\xi^2}$. Se inserta $1=e^{-\xi^2}e^{+\xi^2}$ de modo que la cantidad entre corchetes sea lo que la regla integra:

$$
\int g(\xi)\,d\xi=\int e^{-\xi^2}\big[g(\xi)e^{+\xi^2}\big]d\xi\approx\sum_i w_H(\xi_i)\,e^{+\xi_i^2}\,g(\xi_i).
$$

Aplicado a ambas dimensiones:

$$
\boxed{\;
c_{mn}=\sum_{i}\sum_{j}
f(\xi_i,\xi_j)\;
e^{\xi_i^2}\,w_H(\xi_i)\;
e^{\xi_j^2}\,w_H(\xi_j)\;
\psi_m(\xi_i)\,\psi_n(\xi_j)
\;}
$$

Los factores $e^{\xi_i^2}$, $e^{\xi_j^2}$ **deshacen exactamente** el peso $e^{-\xi^2}$ horneado en los pesos de Gauss-Hermite, convirtiendo la cuadratura ponderada en la proyección sin peso. En el código, `hh_coefficients(f_grid, nodes, weights, mmax, nmax)` de `src/hermite.py` define el **peso modificado** $W_i=w_i\,e^{\xi_i^2}$ (estable en float64 hasta $n\sim$ cientos) y calcula el coeficiente como un producto de matrices:

```
W = weights * exp(nodes**2)
A = (W[:,None] * f_grid) * W[None,:]
c = psi_par.T @ A @ psi_per        # forma (mmax+1, nmax+1)
```

**Estabilidad numérica.** El producto $e^{\xi^2}\psi_m(\xi)$ podría desbordarse para $|\xi_i|$ grande; combinándolos como $e^{\xi^2}\psi_m=H_m(\xi)e^{-\xi^2/2}\cdot e^{\xi^2}/\sqrt{2^m m!\sqrt\pi}$ y trabajando con la forma acotada de las $\psi_m$ (la recurrencia estable de §6.2) se evita el desbordamiento — la "transformada de Hermite estable". El uso del peso modificado $W_i=w_i e^{\xi_i^2}$ con las $\psi_m$ acotadas logra precisamente esto.

**Convergencia $(M,N)$.** Dos enteros controlan el cálculo: $M$ = número de modos de Hermite retenidos ($m,n\le M$), $N$ = nodos de cuadratura por dimensión. Si $f$ fuera ella misma una serie de Hermite de orden $\le M$, el integrando $f\,e^{\xi^2}\psi_m\psi_n$ sería (polinomio de grado $\le 2M$)$\times e^{-\xi^2}$; como Gauss-Hermite es exacto hasta grado $2N-1$, **$N\ge M+1$ integra cada modo retenido exactamente** (por debajo de eso, los coeficientes de orden bajo se contaminan por *aliasing* del contenido de orden alto). Una VDF medida/interpolada no es un polinomio de Hermite finito, así que se empuja $N$ más alto para resolver las colas; los coeficientes cambian rápido y luego **se estancan más allá de $N\approx 60$**, porque el error residual pasa a estar dominado por el error de interpolación IDW (§8), que $N$ no puede reducir. Por eso $N\simeq 60$ se sitúa justo más allá de la meseta de convergencia, en el piso limitado por interpolación, satisfaciendo cómodamente $N\ge M+1$ para el $M$ modesto de interés.

### 9.2 El espectro 2D y sus reducciones

Graficar $\log_{10}(c_{mn}^2/c_{00}^2)$ sobre la retícula entera $(m,n)$ es el análogo del espectro de potencia 2D en el espacio de velocidades. Normalizar por $c_{00}^2$ (la potencia Maxwelliana) lo hace adimensional y mide el apartamiento *fraccional* del equilibrio Maxwelliano; el término "DC" $c_{00}$ se excluye de la suma de la cascada (es el equilibrio, no una fluctuación). En el código, `plot_hh_spectrum(c, ...)` de `src/plots.py` produce esta figura (Tarea 3), con la escala de color `SPEC_CLIM = (-8, 0)` de `src/config.py`.

Dos reducciones complementarias:

- **Cortes direccionales.** El corte paralelo $c_{m0}^2/c_{00}^2$ (índice perpendicular fijo en 0) mide estructura *a lo largo* de $\mathbf B$; el corte perpendicular $c_{0n}^2/c_{00}^2$ mide estructura *a través* de $\mathbf B$. Los grafica `plot_spectrum_cuts(c, ...)`.
- **Espectro promediado por cáscaras (isótropo):** suma sobre cáscaras concéntricas de ancho unidad,
$$
P(m)=\!\!\sum_{m-\frac12<|\mathbf m'|\le m+\frac12}\!\!c_{\mathbf m'}^2,\qquad
\mathbf m=(m_\perp,m_\parallel),\ \ m=\sqrt{m_\perp^2+m_\parallel^2},
$$
el análogo directo del espectro omnidireccional en $k$ de la turbulencia de fluidos.

### 9.3 El estriado par/impar por la simetría de espejo

En HH, la imposición $f(-v_\perp)=f(v_\perp)$ anula los coeficientes de Hermite **impares** en perpendicular: $c_{m,n}=0$ para $n$ impar. Solo los modos perpendiculares pares llevan potencia. Esto es exactamente lo que produce el **patrón oscilatorio (par/impar)** característico en el espectro perpendicular 1D (Larosa, Fig. 3): *"el comportamiento oscilatorio… se debe a la falta de potencia en los modos impares del espectro perpendicular causada por la hipótesis de girotropía."* Esta es una signatura diagnóstica del método HH, no física del plasma; es la razón por la que, al medir pendientes locales, conviene diferenciar saltando de $m$ a $m+2$ (misma paridad).

### 9.4 La cascada anisótropa en el espacio de velocidades

El resultado físico central (Pezzi 2018, confirmando a Servidio 2017): **un campo magnético de fondo hace anisótropa la cascada del espacio de velocidades — se extiende a órdenes más altos (pendiente más somera) a lo largo de $\mathbf B$ y se suprime (pendiente más abrupta) a través de $\mathbf B$.** En palabras de Pezzi: *"los espectros están estirados en la dirección paralela… los gradientes de velocidad son más fuertes a lo largo del campo medio y la cascada está inhibida a través del campo medio"* — un análogo del efecto Shebalin. Mecánicamente, la extensión paralela refleja **mezcla de fases paralela / resonancias de Landau** que generan estructura fina en $v_\parallel$. Cuantitativamente, Pezzi (2018, Fig. 3) midió en simulación:

$$
P(m_z)\sim m_z^{-2.01}\ (\text{paralelo}),\qquad
P(m_\perp)\sim m_\perp^{-3.5}\ (\text{perpendicular, mucho menor y más abrupto}).
$$

El exponente promediado por cáscaras es *diagnóstico del proceso dominante* (Servidio 2017; Larosa §1):
- $P(m)\sim m^{-3/2}$ — régimen **débilmente magnetizado / isótropo** (mezcla de fases / dominado por campo eléctrico);
- $P(m)\sim m^{-2}$ — régimen **magnetizado / dominado por campo magnético**; más abrupto indica dinámica dominada por no linealidad.

Larosa et al. (2025, PSP) miden $P(m)\approx m^{-2}$ en $4\le m<12$ tanto para una corriente de *ondas* (~$28\,R_\odot$, superalfvénica, haces *hammerhead*) como para una *turbulenta* (~$11\,R_\odot$, subalfvénica), "consistente con las expectativas para un plasma de $\beta$ bajo (Servidio et al. 2017)", advirtiendo que la pendiente es "muy sensible al rango de $m$ elegido".

**Bajo orden frente a alto orden — haces frente a cascada** (el contraste clave de Larosa): a bajo $m$ ($1,2,3$) la corriente de *ondas* tiene *más* potencia — el haz de protones denso + *hammerhead* proyecta solo sobre los primeros modos de Hermite; a alto $m$ ($8\lesssim m\lesssim15$) la corriente *turbulenta* tiene *más* potencia, "interpretado como efecto de la cascada del espacio de velocidades más intensa… distorsiones de escala fina de la VDF".

> **Convenciones de trabajo de este proyecto (no en la literatura citada).** El par de pendientes paralela/perpendicular "$m^{-1.3}$ (paralela), $n^{-3.1}$ (perpendicular)" y el estimador de pendiente local
> $$\alpha_m=\frac{c_{m+2,0}^2-c_{m,0}^2}{2\,c_{00}^2}$$
> son *convenciones de trabajo*, no citas textuales de Pezzi/Servidio/Larosa/Coburn. Codifican la misma desigualdad robusta y citable $|\text{pendiente}_\parallel|<|\text{pendiente}_\perp|$ (cascada paralela extendida, perpendicular suprimida). El estimador es una diferencia finita de la potencia del corte paralelo, espaciada dos modos y normalizada por $c_{00}$ — un sustituto discreto de $dP/dm$ a lo largo de $\mathbf B$, que señala dónde la cascada se hace más abrupta y transita a la meseta de ruido. El **espaciado "+2" es deliberado**: como en HH los modos impares perpendiculares se anulan y los modos pares/impares de Hermite llevan potencia sistemáticamente distinta (el espectro oscilatorio), diferenciar $m\to m{+}2$ salta la oscilación de paridad. Se contrasta con el $P(m)$ promediado por cáscaras.

### 9.5 Aplanamiento espectral, truncamiento y reconstrucción

A orden alto el espectro medido **se aplana a una meseta** — esto es **ruido instrumental / error de interpolación, no física**:

- **Larosa 2025:** *"deliberadamente no discutimos los espectros para $m\gtrsim20$ porque el aplanamiento se debe probablemente a errores de interpolación"*; validan que la señal está por encima de las estimaciones de ruido de una cuenta y bi-Maxwelliano.
- **Coburn 2024 (la aplicación explícita de eliminación de ruido):** los coeficientes HL "decrecen en intensidad relativa hasta $m,l\approx14$ donde el espectro se aplana", identificado como el piso de ruido. Entonces **truncan** (ponen a cero todos los coeficientes más allá del piso) y **reconstruyen** una "VDF filtrada paso-bajo": suave, analítica, diferenciable, libre de ruido de disparo instrumental. Como la base es analítica, la **VDF sin ruido y sus derivadas de velocidad** se calculan exactamente, habilitando análisis de interacción onda-partícula.

Este es el análogo directo de filtrar paso-bajo un espectro de Fourier: conservar los modos de gran escala resueltos, descartar el ruido de "número de onda" alto. Para el análisis HH $60\times60$ de este proyecto aplica la misma lógica: leer los modos de orden bajo convergidos como físicos (deriva, anisotropía, flujo de calor, haz/*hammerhead* a bajo $m$; estructura de cascada a $m$ intermedio) y tratar la meseta de orden alto como el piso de ruido limitado por interpolación.

### 9.6 Síntesis del flujo de trabajo

Una VDF de protones del viento solar medida por SWA-PAS se reconstruye de $(E,\text{elev},\text{azim})$ a $\mathbf v_{xyz}$ (con la convención obligatoria de signos SWA para que $V_R>0$), se rota a RTN vía `PAS_to_RTN`, luego a coordenadas alineadas al campo, y se normaliza a $\xi=(v-u)/w$ usando los momentos $n,\mathbf u,T_{\parallel,\perp}$ y velocidades térmicas $w=\sqrt{2k_BT/m}$. Bajo la hipótesis de girotropía se descompone sobre la base ortonormal de Hermite $\psi_m=H_m(\xi)e^{-\xi^2/2}/\sqrt{2^m m!\sqrt\pi}$ (Hermite-Hermite, con $f(-v_\perp)=f(v_\perp)$ matando los modos impares perpendiculares), el análogo en el espacio de velocidades de una transformada de Fourier centrada en la Maxwelliana local. Los coeficientes $c_{mn}=\iint f\psi_m\psi_n\,d\xi$ se evalúan por cuadratura de Gauss-Hermite cuyos nodos son las raíces de Hermite (Golub-Welsch / `hermgauss`, $N=60$) tras interpolar IDW los datos dispersos sobre la malla de nodos. El espectro cuadrado y normalizado $\log_{10}(c_{mn}^2/c_{00}^2)$ expone la jerarquía $c_{00}\sim n$, $m=1,2,3\sim$ deriva/temperatura/flujo de calor, $m,n$ altos $\sim$ estructura fina no térmica; sus cortes direccionales revelan una cascada anisótropa (somera/extendida paralela, abrupta/suprimida perpendicular), y su aplanamiento de orden alto marca un piso de ruido instrumental/interpolación que motiva el truncamiento espectral y la reconstrucción de una VDF analítica sin ruido.

---

## 10. Resultados de esta reproducción

### 10.1 El instante analizado

El caso de estudio es el instante **2022-03-08 14:45:22**, fijado por `TARGET_TIME = "2022-03-08T14:45:22"` en `src/config.py`. El registro CDF más cercano es el índice **idx 13203**, con tiempo exacto **$t=$ 14:45:22.535** (`Info = 1`, modo Normal). Los datos provienen de:
- `data/raw/solo_L2_swa-pas-vdf_20220308_V02.cdf` (la VDF),
- `data/raw/solo_L2_mag-rtn-normal_20220308_V02.cdf` (el campo magnético).

Los productos procesados se guardan en `data/processed/vdf_processed.npz` y `data/processed/moments.json`.

### 10.2 Campo magnético y momentos del protón

Los valores obtenidos (de `data/processed/moments.json`), promediando $\mathbf B$ en una ventana de $\pm$2 s:

| Cantidad | Valor | Notas |
|---|---|---|
| $\mathbf B$ (RTN) | $[-21.5,\,-8.1,\,16.3]$ nT | de MAG L2 `mag-rtn-normal` |
| $\lvert\mathbf B\rvert$ | $28.2$ nT | |
| $n$ | $39.3$ cm$^{-3}$ | densidad de protones |
| $\lvert\mathbf u\rvert$ | $309.7$ km/s | $\mathbf u_{\rm RTN}=[304.4,\,-56.9,\,6.2]$ km/s |
| $u_\parallel$ | $212.7$ km/s | proyección de $\mathbf u$ sobre $\hat{\mathbf b}$ |
| $T_\parallel$ | $2.32\times10^5$ K | $=\hat{\mathbf b}\cdot\mathsf T\cdot\hat{\mathbf b}$ |
| $T_\perp$ | $1.85\times10^5$ K | $=(\operatorname{Tr}\mathsf T-T_\parallel)/2$ |
| $T_\perp/T_\parallel$ | $0.80$ | anisotropía (protones más calientes a lo largo de $\mathbf B$) |
| $w_\parallel$ | $61.9$ km/s | $=\sqrt{2k_BT_\parallel/m_p}$ |
| $w_\perp$ | $55.3$ km/s | $=\sqrt{2k_BT_\perp/m_p}$ |
| Signos de reconstrucción | $(-1,-1,+1)$ | convención SWA confirmada ($u_R>0$) |

**El ángulo flujo-campo.** Con $u_\parallel=212.7$ km/s y $\lvert\mathbf u\rvert=309.7$ km/s, el ángulo entre el flujo masivo y $\mathbf B$ es

$$
\theta_{uB}=\arccos\!\left(\frac{u_\parallel}{\lvert\mathbf u\rvert}\right)=\arccos\!\left(\frac{212.7}{309.7}\right)\approx 47^\circ .
$$

Es decir, el flujo forma $\sim47^\circ$ con el campo; este gran ángulo es la razón del **gran arrastre perpendicular** del flujo (una fracción importante de $\mathbf u$ está en el plano perpendicular), lo que tiene consecuencias directas en cómo se distribuyen los puntos en el plano $(\xi_\parallel,\xi_\perp)$ tras centrar por la velocidad peculiar.

**Puntos válidos.** Del arreglo completo, **614 *bins* tienen PSD$>0$** y finita; estos son los puntos que entran a la nube de dispersión y a la interpolación. (El número total de elementos del arreglo es mayor, pero el *peak-tracking* y los cortes de validez dejan 614 puntos útiles.)

### 10.3 La malla y el espectro

La descomposición usa una malla **Hermite-Hermite de $H_{60}\times H_{60}$** (`N_NODES = M_ORDER = N_ORDER = 60` en `src/config.py`), con la VDF perpendicular reflejada ($f(\xi_\perp)=f(-\xi_\perp)$, estilo Larosa) e interpolada por IDW ($\varepsilon=10^{-3}$) sobre los $60\times60$ nodos de Gauss-Hermite. El espectro $\log_{10}(c_{mn}^2/c_{00}^2)$ resultante (figura `figures/task3_hh_spectrum.png`):

- **Reproduce el estriado par/impar** en el índice $n$ perpendicular, la signatura esperada de la simetría de espejo de la girotropía HH (los coeficientes impares-$n$ son nulos).
- **Concentra la potencia en los órdenes bajos** ($m,n$ pequeños), consistente con que la VDF es cercana a una bi-Maxwelliana con correcciones cinéticas de orden bajo (deriva, anisotropía, algo de flujo de calor), más que con una distribución fuertemente filamentada.

**Los cortes direccionales** (figura `figures/task3_spectrum_cuts.png`) dan una **pendiente paralela más somera** ($\sim-1.7$) que la **perpendicular** ($\sim-4.6$):

$$
|\text{pendiente}_\parallel|\approx 1.7 \;<\; |\text{pendiente}_\perp|\approx 4.6 ,
$$

lo que **confirma la cascada anisótropa**: la estructura del espacio de velocidades se extiende a órdenes más altos a lo largo de $\mathbf B$ (corte paralelo somero) y se suprime a través de $\mathbf B$ (corte perpendicular abrupto). Esto es exactamente la desigualdad robusta $|\text{pendiente}_\parallel|<|\text{pendiente}_\perp|$ de Pezzi/Servidio (cf. los valores de simulación $m_z^{-2.01}$ vs $m_\perp^{-3.5}$ de §9.4), reproducida aquí con datos reales de SWA-PAS. (Las pendientes numéricas exactas son convenciones de trabajo de este proyecto, según la corrección de atribución (iii); lo citable es la desigualdad.)

### 10.4 Figuras producidas

El orquestador `src/run_analysis.py` genera:
- **Tarea 1** — `figures/task1_vdf_scatter.png`: dispersión 3D de la VDF en coordenadas FA.
- **Tarea 2** — `figures/task2_hermite_grid.png`: VDF observada e interpolada sobre la malla de Hermite.
- **Tarea 3** — `figures/task3_hh_spectrum.png` (espectro 2D) y `figures/task3_spectrum_cuts.png` (cortes paralelo/perpendicular).

---

## 11. Referencias

### Instrumento y producto de datos

- **Owen, C. J., et al. 2020**, "The Solar Orbiter Solar Wind Analyser (SWA) suite," *A&A* **642**, A16. DOI: [10.1051/0004-6361/201937259](https://doi.org/10.1051/0004-6361/201937259). https://www.aanda.org/articles/aa/full_html/2020/10/aa37259-19/aa37259-19.html
- **Fedorov, A. 2020**, *SWA-PAS L2 Data User Guide*, V02 (16 nov 2020), ESA COSMOS. https://www.cosmos.esa.int/documents/3689933/11863901/PAS_L2_Data_User_Guide_20201116.pdf — variables, unidades, convenciones de marco, fórmulas de reconstrucción de velocidad, E2V=13.85, advertencias. Versión en línea: https://swa-pas-data-user-guide.readthedocs.io/en/latest/products.html
- **Louarn, P., et al. 2024**, "…solar-wind proton moments from SWA-PAS," *A&A* **682**, A44. DOI: [10.1051/0004-6361/202347874](https://doi.org/10.1051/0004-6361/202347874) (momentos de PAS / calibración cruzada en vuelo).
- **Livi, S., et al. 2023**, "First results from the Solar Orbiter Heavy Ion Sensor," *A&A* **676**, A36 (SWA-**HIS**, **no** PAS; citar solo para desambiguar la referencia "Livi 2023"). DOI: [10.1051/0004-6361/202346304](https://doi.org/10.1051/0004-6361/202346304).
- IRAP SWA/PAS: https://www.irap.omp.eu/en/project/solar-orbiter-pas-2/ ; CNES SWA: https://solar-orbiter.cnes.fr/en/SOLO/GP_swa.htm

### Cuentas → densidad de espacio de fases, reconstrucción, momentos, girotropía

- **Livi, R., et al. 2022**, "The Solar Probe ANalyzer—Ions on PSP," *ApJ* **938**, 138. https://iopscience.iop.org/article/10.3847/1538-4357/ac93f5 (marco instrumental SPAN-i, factor $k$, cuentas→flujo).
- **Verscharen, D., Klein, K. G. & Maruca, B. A. 2019**, "The multi-scale nature of the solar wind," *Living Rev. Solar Phys.* **16**:5. https://pmc.ncbi.nlm.nih.gov/articles/PMC6934245/ ; arXiv:1902.03448 (definiciones de momentos, $T_\parallel,T_\perp$, velocidades térmicas, $\beta$, bi-Maxwelliana Eq. 61, kappa Eq. 62, derivas núcleo/haz/alfa).
- **Broderick / Klein et al. 2025**, "Recovering Ion Distribution Functions: Slepian Reconstruction (MMS & Solar Orbiter)," I: arXiv:2501.17294, https://arxiv.org/html/2501.17294 ; II (girotrópica), *ApJ*, DOI: [10.3847/1538-4357/ae1d71](https://doi.org/10.3847/1538-4357/ae1d71) y DOI: [10.3847/1538-4357/adb6a0](https://doi.org/10.3847/1538-4357/adb6a0).
- **Wilson, L. B., et al. 2025**, "How limited resolution of plasma analyzers affects moment accuracy," arXiv:2505.09869. https://arxiv.org/pdf/2505.09869 (sumas discretas de momentos, discretización $v^2\,dv\,d\Omega$, sesgos de resolución).
- **Verniero, J. L., et al. 2020**, "PSP Observations of Proton Beams Simultaneous with Ion-scale Waves," *ApJS* **248**, 5. https://iopscience.iop.org/article/10.3847/1538-4365/ab86af (VDFs alineadas al campo y haces).
- **Swisdak, M. 2016**, "Quantifying gyrotropy in magnetic reconnection," *GRL*. DOI: [10.1002/2015GL066980](https://doi.org/10.1002/2015GL066980) (medida de agirotropía $Q$; refs. a Scudder & Daughton 2008 $A\!\varnothing$, Aunai et al. 2013 $D_{\rm ng}$).
- Wikipedia, "Thermal velocity." https://en.wikipedia.org/wiki/Thermal_velocity (convenciones de velocidad térmica).

### Modelos paramétricos y ajuste

- **Marsch, E. 2006**, "Kinetic Physics of the Solar Corona and Solar Wind," *Living Rev. Solar Phys.* https://link.springer.com/article/10.12942/lrsp-2006-1
- **De Marco, R., et al. 2023**, "Innovative technique for separating proton core, proton beam, and alpha particles" (GMM + EM, Solar Orbiter/PAS), *A&A*. https://www.aanda.org/articles/aa/full_html/2023/01/aa43719-22/aa43719-22.html
- **Dupuis, R., et al. 2020**, "Characterizing magnetic reconnection regions using Gaussian mixture models on particle VDFs," arXiv:1910.10012. https://arxiv.org/pdf/1910.10012
- **Štverák, Š., et al. 2022**, "Implications of Kappa Suprathermal Halo of the Solar Wind Electrons," *Front. Astron. Space Sci.* https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2022.892236/full
- **Lazar, M., et al. 2017**, "Dual Maxwellian–Kappa modeling of the solar wind electrons," *A&A*. https://www.aanda.org/articles/aa/full_html/2017/06/aa30194-16/aa30194-16.html
- **Pierrard, V. & Lazar, M. 2010**, "Kappa Distributions: Theory and Applications in Space Plasmas," *Solar Phys.*
- **Scherer, K., Fichtner, H. & Lazar, M. 2017** — distribuciones kappa regularizadas (RKD). https://www.mdpi.com/2571-6182/6/3/36 ; "Regularized kappa-halos with ALPS," arXiv:2504.15955.
- "Skewness and kurtosis of solar wind proton distribution functions" (normal inverse-Gaussian) 2024, *A&A*. https://www.aanda.org/articles/aa/full_html/2024/02/aa47874-23/aa47874-23.html
- **Tu, C.-Y., et al. 2004**, "Dependence of the proton beam drift velocity on proton core plasma beta," *JGR*. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2004JA010391

### Descomposición espectral Hermite / Hermite-Laguerre y cascada en el espacio de velocidades

- **Servidio, S., et al. 2017**, "MMS observation of plasma velocity-space cascade," *PRL* **119**, 205101. https://link.aps.org/doi/10.1103/PhysRevLett.119.205101 ; arXiv:1707.08180 (primer espectro $P(m)$ observado; $m^{-3/2}$ vs $m^{-2}$).
- **Pezzi, O., et al. 2018**, "Velocity-space cascade in magnetized plasmas," *Phys. Plasmas* **25**, 060704. https://pubs.aip.org/aip/pop/article/25/6/060704/320076 ; arXiv:1803.01633 (cascada anisótropa $m_z^{-2.01}$ vs $m_\perp^{-3.5}$; base de Hermite, enstrofía).
- **Coburn, J., et al. 2024**, *ApJ* **964**, 100 (**Solar Orbiter** electrones, Hermite-Laguerre, $60\times60$, eliminación de ruido por truncamiento). DOI: [10.3847/1538-4357/ad1329](https://doi.org/10.3847/1538-4357/ad1329) ; ADS: https://ui.adsabs.harvard.edu/abs/2024ApJ...964..100C
- **Larosa, A., et al. 2025**, "Velocity-space turbulent cascade in the near-Sun solar wind" (**Parker Solar Probe**, SPAN-i, Hermite-Hermite, $50\times50$). arXiv:2512.01492. https://arxiv.org/abs/2512.01492
- **Roytershteyn, V. & Delzanno, G. L. 2018**, *Front. Astron. Space Sci.* **5**, 27. https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2018.00027/full
- **Mandell, N. R., et al. 2018**, "Laguerre–Hermite pseudo-spectral velocity formulation of gyrokinetics," arXiv:1708.04029. https://arxiv.org/pdf/1708.04029
- **Vencels, J., et al. 2018**, "Spectral Approach to Plasma Kinetic Simulations Based on Hermite Decomposition," *Front. Astron. Space Sci.*

### Numérica de la cuadratura

- **Golub, G. H. & Welsch, J. H. 1969**, "Calculation of Gauss Quadrature Rules," *Math. Comp.* **23**(106), 221–230. https://bibbase.org/network/publication/golub-welsch-calculationofgaussquadraturerules-1969
- Reichel, L., "Computation of Gauss-type quadrature rules." https://www.math.kent.edu/~reichel/publications/DC.pdf
- Abramowitz, M. & Stegun, I. A., *Handbook of Mathematical Functions*, Cap. 25 (fórmulas 25.4.45, 25.4.46).
- Wolfram MathWorld, "Hermite–Gauss Quadrature." https://mathworld.wolfram.com/Hermite-GaussQuadrature.html
- NumPy: `numpy.polynomial.hermite.hermgauss`, https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermgauss.html ; `numpy.polynomial.laguerre.laggauss`.
- "Stable Hermite transforms via the Golub–Welsch algorithm," arXiv:2604.02041. https://arxiv.org/pdf/2604.02041

---

### Apéndice: mapa de la teoría a la implementación

| Sección | Concepto teórico | Función / archivo |
|---|---|---|
| §2.3 | Lectura del CDF, ejes en tiempo de ejecución | `read_vdf_timestep`, `match_axes` (`cdf_io.py`, `run_analysis.py`) |
| §3.1 | $\lvert v\rvert=\sqrt{2qE/m}$ | `speed_from_energy` (`physics.py`) |
| §3.2 | Convención de signos SWA $(-1,-1,+1)$ | `instrument_velocity` (`physics.py`) |
| §3.3 | Volumen $d^3v=v^2\,dv\,d\Omega$ | `velocity_volume` (`physics.py`) |
| §3.4 | Rotación PAS → RTN | `rotate_pas_to_rtn` (`physics.py`) |
| §4.1 | $\mathbf B$ de MAG L2 | `read_mag_at` (`cdf_io.py`) |
| §4.2 | Tríada alineada al campo | `field_aligned_basis`, `project_field_aligned` (`physics.py`) |
| §4.4 | Momentos $n,\mathbf u,\mathsf P,T_{\parallel,\perp},w$ | `compute_moments` (`physics.py`) |
| §4.5 | $\xi=(v-u)/w$, $c_\perp=\sqrt{\lvert c\rvert^2-c_\parallel^2}$ | `normalized_coords` (`physics.py`) |
| §6.2 | Funciones de Hermite ortonormales | `hermite_functions` (`hermite.py`) |
| §7.5 | Nodos/pesos Gauss-Hermite (Golub-Welsch) | `gauss_hermite` (`hermite.py`) |
| §8.2 | Interpolación IDW | `idw_interpolate` (`hermite.py`) |
| §9.1 | Coeficientes $c_{mn}$ | `hh_coefficients` (`hermite.py`) |
| §9.2–9.5 | Espectro 2D y cortes | `plot_hh_spectrum`, `plot_spectrum_cuts` (`plots.py`) |
| §10 | Caso 2022-03-08 14:45:22 | `run_analysis.py`, `data/processed/moments.json` |
