# -*- coding: utf-8 -*-
"""ETAS espacio-temporal para Colombia y su entorno (catálogo USGS, M >= 4,5, 1973-2026).
Intensidad condicional: lambda(t,x,y) = nu*u(x,y) + sum_{t_i<t} K exp(alfa (m_i-Mc)) g(t-t_i) f(r; m_i)
  g(t) = (p-1) c^(p-1) (t+c)^(-p)                       (Omori-Utsu normalizada, t en días)
  f(r; m) = (q-1)/(pi D^2) (1 + r^2/D^2)^(-q),  D^2 = d^2 exp(gama (m-Mc))   (núcleo espacial potencial, r en km)
  u(x,y): densidad de fondo por núcleos gaussianos de ancho variable ponderados por la probabilidad de ser evento de fondo
          (desagrupamiento estocástico de Zhuang, Ogata y Vere-Jones, 2002).
Estimación por máxima verosimilitud alternada con la actualización del fondo. Magnitudes: Gutenberg-Richter con b por Aki-Utsu.
Pronóstico: simulación del proceso (fondo + réplicas del catálogo + cascadas) a 30, 90 y 365 días.
Uso: python etas_colombia.py [fecha_corte AAAA-MM-DD] [sufijo]   (sin argumentos: todo el catálogo)"""
import numpy as np, pandas as pd, json, os, sys, time
from scipy.optimize import minimize
from scipy.spatial import cKDTree
RAIZ = os.path.dirname(os.path.abspath(__file__)); D = f"{RAIZ}/datos"; OUT = f"{RAIZ}/resultados"; os.makedirs(OUT, exist_ok=True)
CORTE = pd.Timestamp(sys.argv[1], tz="UTC") if len(sys.argv) > 1 else None; SUF = sys.argv[2] if len(sys.argv) > 2 else ""
MC, MMAX, RMAX, LAT0 = 4.5, 8.2, 400.0, 4.5
cat = pd.read_csv(f"{D}/sismos_usgs_colombia_1973_2026_M4.csv", parse_dates=["time"]); cat = cat[(cat.type == "earthquake") & (cat.mag >= MC)].sort_values("time").reset_index(drop=True)
FIN = CORTE if CORTE is not None else cat.time.max() + pd.Timedelta(hours=1); full = cat.copy(); cat = cat[cat.time < FIN].reset_index(drop=True)
t0 = cat.time.min(); t = ((cat.time - t0).dt.total_seconds() / 86400).values; T = (FIN - t0).total_seconds() / 86400
kx = 111.32 * np.cos(np.deg2rad(LAT0)); ky = 110.57
x = (cat.longitude.values + 74) * kx; y = (cat.latitude.values - LAT0) * ky; m = cat.mag.values; N = len(cat)
b = np.log10(np.e) / (m.mean() - (MC - 0.05)); beta = b * np.log(10)
print(f"N={N} eventos M>={MC}, {t0.date()} a {FIN.date()} (T={T:.0f} días), b={b:.3f}")
# pares (i hijo, j padre) con t_j < t_i y distancia < RMAX
tree = cKDTree(np.c_[x, y]); pr = tree.query_pairs(RMAX, output_type="ndarray"); a_, b_ = pr[:, 0], pr[:, 1]
late = np.where(t[a_] > t[b_], a_, b_); early = np.where(t[a_] > t[b_], b_, a_); ok = t[late] > t[early]; I = late[ok]; J = early[ok]
dt = t[I] - t[J]; r2 = (x[I] - x[J]) ** 2 + (y[I] - y[J]) ** 2; dm = m[J] - MC; print("pares:", len(I))
def fondo(phi):
    """densidad de fondo u en cada evento (por km2 y por día, antes de escalar por nu), núcleos de ancho variable"""
    dk, _ = tree.query(np.c_[x, y], k=6); h = np.clip(dk[:, -1], 15.0, 150.0)
    u = np.zeros(N); pr2 = tree.query_pairs(450.0, output_type="ndarray"); p, q = pr2[:, 0], pr2[:, 1]; d2 = (x[p] - x[q]) ** 2 + (y[p] - y[q]) ** 2
    np.add.at(u, p, phi[q] * np.exp(-d2 / (2 * h[q] ** 2)) / (2 * np.pi * h[q] ** 2)); np.add.at(u, q, phi[p] * np.exp(-d2 / (2 * h[p] ** 2)) / (2 * np.pi * h[p] ** 2))
    u += phi / (2 * np.pi * h ** 2); return u / T, h
def partes(th, u):
    nu, K, alfa, c, p, d, gama, q = np.exp(th[0]), np.exp(th[1]), np.exp(th[2]), np.exp(th[3]), 1 + np.exp(th[4]), np.exp(th[5]), np.exp(th[6]), 1 + np.exp(th[7])
    D2 = d ** 2 * np.exp(gama * dm); g = (p - 1) * c ** (p - 1) * (dt + c) ** (-p); f = (q - 1) / (np.pi * D2) * (1 + r2 / D2) ** (-q)
    trig = np.bincount(I, weights=K * np.exp(alfa * dm) * g * f, minlength=N); lam = nu * u + trig
    G = 1 - (c / (T - t + c)) ** (p - 1); integ = nu * PHI + (K * np.exp(alfa * (m - MC)) * G).sum()
    return lam, trig, integ, dict(nu=nu, K=K, alfa=alfa, c=c, p=p, d=d, gama=gama, q=q)
def nll(th, u):
    lam, _, integ, _ = partes(th, u); return -(np.log(lam + 1e-300).sum() - integ)
th = np.log([1.0, 0.3, 1.2, 0.02, 0.15, 12.0, 1.0, 0.7]); phi = np.full(N, 0.6); hist = []
for it in range(12):
    u, h = fondo(phi); PHI = phi.sum(); tic = time.time()
    res = minimize(nll, th, args=(u,), method="L-BFGS-B", bounds=[(-6, 4), (-6, 2), (-3, 1.2), (-9, 2), (-5, 1.5), (0, 5), (-4, 1.2), (-4, 2)], options=dict(maxiter=150))
    th = res.x; lam, trig, integ, par = partes(th, u); phi_n = par["nu"] * u / lam; cambio = np.abs(phi_n - phi).max(); phi = phi_n
    hist.append(dict(it=it, nll=float(res.fun), frac_fondo=float(phi.mean()), cambio=float(cambio), **{k: float(v) for k, v in par.items()}))
    print(f"it {it}: nll={res.fun:.1f} fondo={phi.mean():.3f} cambio={cambio:.4f} {time.time()-tic:.0f}s", {k: round(float(v), 4) for k, v in par.items()}, flush=True)
    if cambio < 0.01 and it >= 3: break
u, h = fondo(phi); PHI = phi.sum(); lam, trig, integ, par = partes(th, u)
ramif = float(par["K"] * beta / (beta - par["alfa"]) * (1 - np.exp(-(beta - par["alfa"]) * (MMAX - MC))) / (1 - np.exp(-beta * (MMAX - MC)))) if beta != par["alfa"] else None
out = dict(Mc=MC, N=int(N), inicio=str(t0), fin=str(FIN), T_dias=T, b=float(b), beta=float(beta), parametros={k: float(v) for k, v in par.items()}, fraccion_fondo=float(phi.mean()), tasa_fondo_dia=float(par["nu"] * PHI / T),
           razon_ramificacion=ramif, nll=float(nll(th, u)), historia=hist)
json.dump(out, open(f"{OUT}/etas_parametros{SUF}.json", "w"), indent=1); print(json.dumps({k: v for k, v in out.items() if k != "historia"}, indent=1))
cat["x"], cat["y"], cat["t"], cat["phi"], cat["h"], cat["lam"] = x, y, t, phi, h, lam
cat[["time", "latitude", "longitude", "depth", "mag", "place", "x", "y", "t", "phi", "h", "lam"]].to_csv(f"{OUT}/catalogo_etas{SUF}.csv", index=False)
