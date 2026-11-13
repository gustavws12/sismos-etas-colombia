# -*- coding: utf-8 -*-
"""Pronóstico por simulación del ETAS ajustado: fondo + réplicas del catálogo + cascadas, magnitudes remuestreadas de la distribución empírica del catálogo.
Uso: python pronostico_etas.py [sufijo] [horizonte_dias] [n_sim]"""
import numpy as np, pandas as pd, geopandas as gpd, json, os, sys
from shapely.geometry import Point
RAIZ = os.path.dirname(os.path.abspath(__file__)); D = f"{RAIZ}/datos"; OUT = f"{RAIZ}/resultados"
SUF = sys.argv[1] if len(sys.argv) > 1 else ""; SUF = "" if SUF == "-" else SUF; H = float(sys.argv[2]) if len(sys.argv) > 2 else 365.0; NS = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
P = json.load(open(f"{OUT}/etas_parametros{SUF}.json")); c = pd.read_csv(f"{OUT}/catalogo_etas{SUF}.csv"); par = P["parametros"]; MC, MMAX, beta, T = P["Mc"], 8.2, P["beta"], P["T_dias"]
K, alfa, cc, p, d, gama, q, nu = par["K"], par["alfa"], par["c"], par["p"], par["d"], par["gama"], par["q"], par["nu"]
LAT0 = 4.5; kx = 111.32 * np.cos(np.deg2rad(LAT0)); ky = 110.57; rng = np.random.default_rng(20261113)
G = lambda s: 1 - (cc / (s + cc)) ** (p - 1)
MAGS = c.mag.values
def mags(n): return rng.choice(MAGS, n)   # distribución empírica de magnitudes del catálogo (la Gutenberg-Richter con b único sobrestima M5 y subestima M7 en este catálogo)
def hijos(tp, xp, yp, mp, desde):
    """réplicas directas de los padres (tp, xp, yp, mp) que caen en (desde, T+H]; devuelve arreglos t, x, y, m"""
    lo = G(np.maximum(desde - tp, 0)); hi = G(T + H - tp); n = rng.poisson(K * np.exp(alfa * (mp - MC)) * (hi - lo)); tot = n.sum()
    if tot == 0: return [np.empty(0)] * 4
    idx = np.repeat(np.arange(len(tp)), n); u = lo[idx] + rng.random(tot) * (hi[idx] - lo[idx]); s = cc * ((1 - u) ** (-1 / (p - 1)) - 1)
    D2 = d ** 2 * np.exp(gama * (mp[idx] - MC)); r = np.sqrt(D2 * ((1 - rng.random(tot)) ** (-1 / (q - 1)) - 1)); a = rng.random(tot) * 2 * np.pi
    return tp[idx] + s, xp[idx] + r * np.cos(a), yp[idx] + r * np.sin(a), mags(tot)
tc, xc, yc, mc_, phi, h = c.t.values, c.x.values, c.y.values, c.mag.values, c.phi.values, c.h.values; pb = phi / phi.sum(); tasa_fondo = nu * phi.sum() / T
sims = []
for s in range(NS):
    nb = rng.poisson(tasa_fondo * H); j = rng.choice(len(c), nb, p=pb)
    ev = [(T + rng.random(nb) * H, xc[j] + rng.normal(0, h[j]), yc[j] + rng.normal(0, h[j]), mags(nb), np.zeros(nb))]
    g1 = hijos(tc, xc, yc, mc_, T); ev.append((*g1, np.ones(len(g1[0]))))
    gen = (np.concatenate([ev[0][0], g1[0]]), np.concatenate([ev[0][1], g1[1]]), np.concatenate([ev[0][2], g1[2]]), np.concatenate([ev[0][3], g1[3]]))
    while len(gen[0]):
        g = hijos(*gen, T); g = tuple(v for v in g)
        if len(g[0]) == 0: break
        ev.append((*g, np.full(len(g[0]), 2.0))); gen = g
    A = np.column_stack([np.concatenate([e[k] for e in ev]) for k in range(5)]); sims.append(np.column_stack([np.full(len(A), s), A]))
S = pd.DataFrame(np.vstack(sims), columns=["sim", "t", "x", "y", "mag", "origen"]); S["dia"] = S.t - T; S["lon"] = S.x / kx - 74; S["lat"] = S.y / ky + LAT0
S.to_pickle(f"{OUT}/simulaciones{SUF}.pkl"); print("eventos simulados:", len(S), "por simulación:", round(len(S) / NS, 1), "| origen (0 fondo, 1 réplica del catálogo, 2 cascada):", S.origen.value_counts(normalize=True).round(3).to_dict())
