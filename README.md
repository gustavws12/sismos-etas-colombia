# ¿Dónde y con qué probabilidad va a temblar en Colombia? Pronóstico de sismos con un modelo ETAS

Pronóstico probabilístico de sismos para Colombia a 30 días, 90 días y 12 meses, con un modelo ETAS espacio-temporal (secuencias de réplicas de tipo epidémico) ajustado por máxima verosimilitud al catálogo del USGS de 1973 a 2026.

Gustavo A Romero, estadístico y científico de datos. Noviembre de 2026.

## Lo que se puede y lo que no

Un sismo no se puede predecir en fecha, lugar y magnitud. Lo que sí se puede es pronosticar probabilidades, como en meteorología: cuántos sismos de cierta magnitud se esperan en una zona y un plazo, y la probabilidad de que ocurra al menos uno. Ese es el alcance de este repositorio. No reemplaza los estudios de amenaza sísmica del Servicio Geológico Colombiano ni la norma sismorresistente.

## Datos

| Serie | Fuente | Periodo | Archivo |
|---|---|---|---|
| Sismos de magnitud 4 o más, 5° S a 14° N y 82° a 66° O (7.436 eventos) | USGS, catálogo ComCat (servicio FDSN) | 1973-01-05 a 2026-11-13 | `datos/sismos_usgs_colombia_1973_2026_M4.csv` |
| Geometría municipal simplificada (para departamentos y mapa) | DANE, Marco Geoestadístico Nacional 2025 | 2025 | `datos/mpios_simplificado.geojson` |

El modelo usa los 4.025 sismos de magnitud 4,5 o más, umbral desde el que el catálogo es razonablemente completo.

## Método

La tasa esperada de sismos en el instante $t$ y el punto $(x,y)$ es

$$\lambda(t,x,y) = \nu\,u(x,y) + \sum_{i:\,t_i<t} K\,e^{\alpha(m_i - M_c)}\; g(t-t_i)\; f(r_i;\,m_i)$$

donde $\nu\,u(x,y)$ es la sismicidad de fondo; $K\,e^{\alpha(m_i - M_c)}$ es el número esperado de réplicas directas de un sismo de magnitud $m_i$; $g(t) = (p-1)\,c^{p-1}(t+c)^{-p}$ es la ley de Omori-Utsu, con $t$ en días; y $f(r;m) = \frac{q-1}{\pi D^2}\left(1+\frac{r^2}{D^2}\right)^{-q}$, con $D^2 = d^2 e^{\gamma(m-M_c)}$ y $r$ en kilómetros, es el núcleo espacial. $M_c = 4{,}5$ es la magnitud de completitud.

El fondo $u(x,y)$ se estima con núcleos gaussianos de ancho variable ponderados por la probabilidad de que cada sismo sea de fondo y no réplica (desagrupamiento estocástico de Zhuang, Ogata y Vere-Jones, 2002). Los ocho parámetros se estiman por máxima verosimilitud, alternando con la actualización del fondo hasta converger.

Parámetros estimados: $K = 0{,}37$, $\alpha = 1{,}37$, $c = 0{,}013$ días, $p = 1{,}05$, $d = 8{,}1$ km, $\gamma = 0{,}96$, $q = 2{,}11$. El 63 % de los sismos son de fondo y el 37 % réplicas.

El pronóstico se obtiene simulando el proceso 4.000 veces: sismos de fondo, réplicas pendientes de todo el catálogo y sus cascadas. Las magnitudes se remuestrean de la distribución observada.

## Scripts

Se corren en este orden desde la raíz del repositorio. Dependencias: `pip install -r requirements.txt`.

1. `python etas_colombia.py` ajusta el modelo a todo el catálogo (unos 40 segundos). Produce `resultados/etas_parametros.json` y `resultados/catalogo_etas.csv`.
2. `python etas_colombia.py 2024-01-01 _val` ajusta el modelo solo hasta 2023, para la validación.
3. `python pronostico_etas.py - 365 4000` simula 12 meses; `python pronostico_etas.py _val 730 2000` simula 2024-2025 para la validación.
4. `python resumen_pronostico.py` calcula probabilidades por magnitud, ciudad, departamento y celda, valida y dibuja la figura.

## Resultados principales

Probabilidad de al menos un sismo en Colombia (territorio y hasta 50 km mar afuera):

| Magnitud | 30 días | 90 días | 12 meses | Promedio histórico a 12 meses |
|---|---|---|---|---|
| 5 o más | 52 % | 88 % | 100 % | 100 % |
| 5,5 o más | 20 % | 50 % | 92 % | 89 % |
| 6 o más | 9 % | 23 % | 63 % | 55 % |
| 6,5 o más | 3 % | 10 % | 34 % | 30 % |
| 7 o más | 2 % | 4 % | 17 % | 15 % |

Probabilidad de un sismo de magnitud 5 o más a menos de 100 km en 12 meses: Bucaramanga 94 %, Armenia 63 %, Pereira 60 %, Manizales 50 %, Buenaventura 43 %, Quibdó 38 %, Ibagué 37 %, Cali 35 %, Cúcuta 31 %, Medellín 18 %, Bogotá 14 %, costa Caribe menos de 7 %.

El Eje Cafetero y el norte del Valle están entre 15 y 17 puntos por encima de su promedio histórico por las réplicas pendientes del sismo de magnitud 7,4 de San José del Palmar (10 de agosto de 2026). El país en conjunto está apenas sobre su promedio: el exceso de sismicidad de 2026 viene sobre todo de Venezuela y Panamá.

Validación: con datos hasta 2023 el modelo esperaba 148 sismos de magnitud 4,5 o más para 2024-2025 (intervalo de 90 %: 121 a 181) y ocurrieron 165; para magnitud 5 o más esperaba 40 (29 a 53) y ocurrieron 29. AUC espacial de 0,90 frente a 0,89 del conteo histórico por celda.

## Figura

`fig1_pronostico_sismos.png`: mapa de probabilidad por celda, probabilidad por magnitud y horizonte, y probabilidad por ciudad frente al promedio histórico.

## Limitaciones

El modelo es bidimensional y no distingue la profundidad: el nido de Bucaramanga está a unos 160 km y un sismo a esa profundidad se siente pero rara vez daña. La probabilidad de sismo no es la probabilidad de daño. La magnitud se sortea igual en todas partes, lo que sobrestima los sismos grandes en el nido de Bucaramanga (un solo sismo de magnitud 6 o más en 54 años) y los subestima en la subducción del Pacífico. Se usa una sola ley de productividad, aunque los sismos profundos generan menos réplicas, de modo que el aumento del Eje Cafetero probablemente está sobrestimado. El catálogo solo es completo desde magnitud 4,5, y 53 años es poco frente al ciclo de los grandes terremotos.

## Notas

`notas/01_pronostico_sismos.md` contiene la discusión completa, la predicción previa y los contrapesos.

## Referencias

Ogata, Y. (1988). Statistical models for earthquake occurrences and residual analysis for point processes. Journal of the American Statistical Association, 83, 9-27. Ogata, Y. (1998). Space-time point-process models for earthquake occurrences. Annals of the Institute of Statistical Mathematics, 50, 379-402. Zhuang, J., Ogata, Y. y Vere-Jones, D. (2002). Stochastic declustering of space-time earthquake occurrences. Journal of the American Statistical Association, 97, 369-380.

## Licencia

MIT para el código. Los datos son del USGS (dominio público) y del DANE, y se redistribuyen tal como se descargaron de sus portales, con las derivaciones indicadas.
