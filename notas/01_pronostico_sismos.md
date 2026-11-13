# ¿Dónde y con qué probabilidad va a temblar en Colombia? Pronóstico con un modelo ETAS

Gustavo A Romero. 13 de noviembre de 2026.

## Pregunta

En 2026 se sintieron varios sismos fuertes en Colombia y su entorno. ¿Qué metodología existe para pronosticar sismos y qué dice aplicada al país para los próximos doce meses?

## Lo que se puede y lo que no

Un sismo no se puede predecir en fecha, lugar y magnitud; ningún método lo ha logrado de manera verificable. Lo que sí se puede es pronosticar probabilidades, como en meteorología: cuántos sismos de cierta magnitud se esperan en una zona y un plazo, y cuál es la probabilidad de que ocurra al menos uno. Ese es el alcance de este ejercicio.

## Predicción previa

Antes de ajustar el modelo escribí que el nido de Bucaramanga concentraría la mayor probabilidad de sismos de magnitud 4 a 5, que la costa Pacífica de Nariño y Chocó concentraría la de magnitud 6 o más, y que la zona del sismo de agosto seguiría con tasa elevada por réplicas. Acertó la primera y la tercera. La segunda salió a medias: Chocó sí queda segundo, pero Nariño queda por debajo de su promedio histórico, porque ese promedio lo inflan las réplicas del terremoto de 1979 y hoy la zona está quieta.

## Datos

Catálogo ComCat del Servicio Geológico de Estados Unidos (USGS): 7.436 sismos de magnitud 4 o más entre el 5 de enero de 1973 y el 13 de noviembre de 2026, en el rectángulo de 5° sur a 14° norte y 82° a 66° oeste, que cubre Colombia y las fuentes vecinas de Panamá, Venezuela, Ecuador y el norte de Perú. Se usan los 4.025 de magnitud 4,5 o más, que es el umbral desde el que el catálogo es razonablemente completo.

El año 2026 es anómalo: 375 sismos de magnitud 4 o más hasta el 13 de noviembre, el doble de un año típico (unos 185), con cuatro eventos mayores: Venezuela el 24 de junio (7,2 y 7,5), San José del Palmar, Chocó, el 10 de agosto (7,4 a 116 km de profundidad), y el sur de Panamá el 9 de octubre (7,7, con más de 150 réplicas en el mes).

## Método

El modelo ETAS (secuencias de réplicas de tipo epidémico, Ogata 1988 y 1998) es el estándar de los pronósticos operativos de Italia, Nueva Zelanda y California. Trata la sismicidad como un proceso puntual en el que cada sismo puede disparar otros. La tasa esperada de sismos en el instante $t$ y el punto $(x,y)$ es

$$\lambda(t,x,y) = \nu\,u(x,y) + \sum_{i:\,t_i<t} K\,e^{\alpha(m_i - M_c)}\; g(t-t_i)\; f(r_i;\,m_i)$$

donde $\nu\,u(x,y)$ es la sismicidad de fondo, que no depende del tiempo; $K\,e^{\alpha(m_i - M_c)}$ es el número esperado de réplicas directas de un sismo de magnitud $m_i$ (ley de productividad); $g(t) = (p-1)\,c^{p-1}(t+c)^{-p}$ es la ley de Omori-Utsu, con la que las réplicas decaen en el tiempo; y $f(r;m) = \frac{q-1}{\pi D^2}\left(1+\frac{r^2}{D^2}\right)^{-q}$, con $D^2 = d^2 e^{\gamma(m-M_c)}$, es un núcleo espacial de cola potencial cuyo alcance crece con la magnitud. $M_c = 4{,}5$ es la magnitud de completitud, $r_i$ la distancia en kilómetros al sismo $i$ y $t$ se mide en días.

El fondo $u(x,y)$ se estima con núcleos gaussianos de ancho variable ponderados por la probabilidad de que cada sismo sea de fondo y no réplica (desagrupamiento estocástico de Zhuang, Ogata y Vere-Jones, 2002). Los ocho parámetros se estiman por máxima verosimilitud, alternando con la actualización del fondo hasta converger. Resultado: $K = 0{,}37$, $\alpha = 1{,}37$, $c = 0{,}013$ días, $p = 1{,}05$, $d = 8{,}1$ km, $\gamma = 0{,}96$, $q = 2{,}11$; el 63 % de los sismos son de fondo y el 37 % réplicas.

El pronóstico se obtiene simulando el proceso 4.000 veces hacia adelante: sismos de fondo, réplicas pendientes de todo el catálogo y las cascadas que generen. Las magnitudes se remuestrean de la distribución observada, porque una ley de Gutenberg-Richter con un solo valor de $b$ (1,09) sobrestimaba los sismos de magnitud 5 y subestimaba los de 7 en este catálogo.

## Validación

Con el modelo ajustado solo hasta diciembre de 2023 se pronosticaron 2024 y 2025. Para sismos de magnitud 4,5 o más el modelo esperaba 148 (intervalo de 90 %: 121 a 181) y ocurrieron 165. Para magnitud 5 o más esperaba 40 (29 a 53) y ocurrieron 29, en el borde inferior del intervalo: 2024 fue un año inusualmente tranquilo en ese rango. En el espacio, el modelo ordena las celdas donde hubo sismos con un AUC de 0,90, frente a 0,89 de usar el conteo histórico de cada celda. La ganancia sobre la historia es pequeña en un periodo sin grandes secuencias; la ventaja del ETAS aparece después de un sismo grande, que es la situación actual.

## Resultados

Probabilidad de al menos un sismo en Colombia (territorio y hasta 50 km mar afuera):

| Magnitud | 30 días | 90 días | 12 meses | Promedio histórico a 12 meses |
|---|---|---|---|---|
| 5 o más | 52 % | 88 % | 100 % | 100 % |
| 5,5 o más | 20 % | 50 % | 92 % | 89 % |
| 6 o más | 9 % | 23 % | 63 % | 55 % |
| 6,5 o más | 3 % | 10 % | 34 % | 30 % |
| 7 o más | 2 % | 4 % | 17 % | 15 % |

Se esperan 37 sismos de magnitud 4,5 o más y 10 de magnitud 5 o más en doce meses dentro de Colombia, apenas por encima del promedio histórico (36 y 9). En toda la región de estudio se esperan 92, frente a unos 75 de un año promedio: el exceso está en Panamá, no en Colombia.

Por zona, a doce meses, probabilidad de un sismo de magnitud 5 o más a menos de 100 km: Bucaramanga 94 % (histórico 85 %), Armenia 63 % (47 %), Pereira 60 % (45 %), Manizales 50 % (44 %), Buenaventura 43 % (26 %), Quibdó 38 % (39 %), Ibagué 37 % (37 %), Cali 35 % (31 %), Cúcuta 31 % (37 %), Tunja 22 % (17 %), Tumaco 21 % (46 %), Medellín 18 % (12 %), Villavicencio 17 % (23 %), Popayán 15 %, Neiva 14 %, Bogotá 14 % (15 %), Pasto 12 %, y menos de 7 % en Santa Marta, Montería, Barranquilla y Cartagena.

Por departamento: Santander 94 %, Chocó 67 %, Valle del Cauca 40 %, Antioquia 40 %, Norte de Santander 18 %, Boyacá 17 %, Tolima 17 %, Meta 17 %, Casanare 12 %, Cauca 12 %, Cundinamarca 10 %.

## Lectura

Lo que está elevado respecto a la historia es el Eje Cafetero y el norte del Valle: Armenia, Pereira y Buenaventura suben entre 15 y 17 puntos sobre su promedio, por las réplicas pendientes del sismo de 7,4 de San José del Palmar del 10 de agosto. Con $p = 1{,}05$ las réplicas decaen despacio, así que el efecto dura años, no semanas. Bucaramanga está alta siempre: el nido sísmico produce por sí solo cerca de un tercio de los sismos del país.

Lo que no está elevado es el país en conjunto. La sensación de que "está temblando más" es correcta para 2026, pero viene sobre todo de Venezuela y Panamá. El sismo de Panamá explica el exceso regional y casi no cambia la probabilidad dentro de Colombia.

Para Bogotá, Medellín y la costa Caribe el pronóstico es prácticamente el promedio histórico.

## Contrapesos

El modelo es bidimensional y no distingue la profundidad. El nido de Bucaramanga está a unos 160 km y el sismo de agosto fue a 116 km: un sismo de magnitud 5 a esa profundidad se siente pero rara vez causa daño, mientras uno superficial de la misma magnitud bajo una ciudad sí. La probabilidad de sismo no es la probabilidad de daño.

La magnitud se sortea igual en todas partes. Eso sobrestima los sismos grandes donde la fuente no los produce (el nido de Bucaramanga tuvo un solo sismo de magnitud 6 o más en 54 años, con máximo de 6,2, así que el 30 % de magnitud 6 o más que el modelo le asigna a Santander está inflado) y los subestima donde sí, que es la zona de subducción del Pacífico.

Las réplicas de sismos profundos son mucho menos numerosas que las de sismos superficiales, y el modelo usa una sola ley de productividad. Por eso el aumento que asigna al Eje Cafetero por el sismo de San José del Palmar probablemente está sobrestimado.

El catálogo del USGS solo es completo desde magnitud 4,5. El del Servicio Geológico Colombiano baja hasta magnitud 2 y permitiría un modelo más fino, con profundidad. Y 53 años es poco frente al ciclo de los grandes terremotos: el modelo no sabe nada de fallas que no han roto en ese lapso, como el segmento del Pacífico que rompió en 1906, ni de la falla frontal de la Cordillera Oriental.

Este pronóstico no reemplaza los estudios de amenaza sísmica del Servicio Geológico Colombiano ni la norma sismorresistente. Sirve para responder con números a la pregunta de si está temblando más de lo normal y dónde.

## Archivos

`etas_colombia.py` ajusta el modelo; `pronostico_etas.py` simula; `resumen_pronostico.py` calcula probabilidades, valida y dibuja. Resultados en `resultados/`: `etas_parametros.json`, `pronostico_resumen.json`, `prob_ciudades.csv`, `prob_departamentos.csv`, `prob_rejilla_05.csv`, `catalogo_etas.csv`. Figura: `fig1_pronostico_sismos.png`. Datos: `datos/sismos_usgs_colombia_1973_2026_M4.csv`.
