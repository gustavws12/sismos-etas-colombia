# -*- coding: utf-8 -*-
"""Resume las simulaciones del ETAS: probabilidades por horizonte y magnitud, por ciudad, por departamento y por celda; validación 2024-2025; figuras."""
import numpy as np, pandas as pd, geopandas as gpd, json, os, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from sklearn.metrics import roc_auc_score
RAIZ = os.path.dirname(os.path.abspath(__file__)); D = f"{RAIZ}/datos"; OUT = f"{RAIZ}/resultados"
AZUL, NARANJA, AQUA, AMAR, GRIS, GRISC, TXT, TXT2, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#9a9893", "#d9d8d4", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": GRISC, "axes.labelcolor": TXT2, "xtick.color": TXT2, "ytick.color": TXT2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.grid": True, "grid.color": "#e8e7e3"})
S = pd.read_pickle(f"{OUT}/simulaciones.pkl"); NS = int(S.sim.max()) + 1; P = json.load(open(f"{OUT}/etas_parametros.json")); cat = pd.read_csv(f"{OUT}/catalogo_etas.csv", parse_dates=["time"])
mp = gpd.read_file(f"{D}/mpios_simplificado.geojson"); mp["geometry"] = mp.buffer(0.012); dep = mp.dissolve("dpto_cnmbr").reset_index()[["dpto_cnmbr", "geometry"]]; dep["geometry"] = dep.buffer(-0.012); col = dep.dissolve().geometry.iloc[0]; colb = col.buffer(0.45)   # ~50 km mar afuera
pts = gpd.GeoDataFrame(S, geometry=gpd.points_from_xy(S.lon, S.lat), crs=4326); pts["en_col"] = pts.within(colb)
pts = gpd.sjoin(pts, dep, how="left", predicate="within").drop(columns="index_right")
def prob(mask, sub=None):
    z = pts[mask]; return z.sim.nunique() / NS
HZ = [30, 90, 365]; MG = [4.5, 5.0, 5.5, 6.0, 6.5, 7.0]; res = dict(n_sim=NS, fin_catalogo=P["fin"], parametros=P["parametros"], b=P["b"], fraccion_fondo=P["fraccion_fondo"], razon_ramificacion=P["razon_ramificacion"])
res["region_completa"] = {f"{h}d": {f"M{m}": dict(prob=round(prob((pts.dia <= h) & (pts.mag >= m)), 4), esperados=round(((pts.dia <= h) & (pts.mag >= m)).sum() / NS, 2)) for m in MG} for h in HZ}
res["colombia"] = {f"{h}d": {f"M{m}": dict(prob=round(prob((pts.dia <= h) & (pts.mag >= m) & pts.en_col), 4), esperados=round(((pts.dia <= h) & (pts.mag >= m) & pts.en_col).sum() / NS, 2)) for m in MG} for h in HZ}
# tasas históricas en Colombia para comparar (M>=5, M>=6, M>=7 por año)
cg = gpd.GeoDataFrame(cat, geometry=gpd.points_from_xy(cat.longitude, cat.latitude), crs=4326); cat["en_col"] = cg.within(colb).values; anos = P["T_dias"] / 365.25
res["historico_colombia_por_ano"] = {f"M{m}": round(float(((cat.mag >= m) & cat.en_col).sum() / anos), 3) for m in MG}; res["historico_prob_poisson_365d"] = {k: round(float(1 - np.exp(-v)), 4) for k, v in res["historico_colombia_por_ano"].items()}
# ciudades: al menos un sismo M>=5 (y M>=6) a menos de 100 km en 365 días
CI = {"Bogotá": (4.71, -74.07), "Medellín": (6.25, -75.56), "Cali": (3.45, -76.53), "Barranquilla": (10.96, -74.80), "Cartagena": (10.39, -75.51), "Bucaramanga": (7.12, -73.12), "Cúcuta": (7.89, -72.50), "Pereira": (4.81, -75.69), "Manizales": (5.07, -75.52), "Armenia": (4.53, -75.68),
      "Ibagué": (4.44, -75.24), "Popayán": (2.44, -76.61), "Pasto": (1.21, -77.28), "Neiva": (2.93, -75.28), "Villavicencio": (4.14, -73.63), "Quibdó": (5.69, -76.66), "Tumaco": (1.81, -78.76), "Buenaventura": (3.88, -77.03), "Santa Marta": (11.24, -74.20), "Montería": (8.75, -75.88), "Tunja": (5.53, -73.36), "Yopal": (5.34, -72.39)}
kx = 111.32 * np.cos(np.deg2rad(4.5)); ky = 110.57; rows = []
for n, (la, lo) in CI.items():
    dkm = np.hypot((pts.lon - lo) * kx, (pts.lat - la) * ky); near = dkm <= 100; hk = np.hypot((cat.longitude - lo) * kx, (cat.latitude - la) * ky) <= 100
    rows.append(dict(ciudad=n, prob_M5_365d=prob(near & (pts.mag >= 5)), prob_M5_90d=prob(near & (pts.mag >= 5) & (pts.dia <= 90)), prob_M6_365d=prob(near & (pts.mag >= 6)), esperados_M45_365d=round((near.sum()) / NS, 2),
                     historico_M5_por_ano=round(float((hk & (cat.mag >= 5)).sum() / anos), 3), prob_hist_M5_365d=round(float(1 - np.exp(-(hk & (cat.mag >= 5)).sum() / anos)), 3)))
ciu = pd.DataFrame(rows).sort_values("prob_M5_365d", ascending=False); ciu.round(4).to_csv(f"{OUT}/prob_ciudades.csv", index=False); print(ciu.round(3).to_string(index=False))
# departamentos
dp = pts.dropna(subset=["dpto_cnmbr"]); rows = []
for n, z in dp.groupby("dpto_cnmbr"):
    rows.append(dict(departamento=n, prob_M5_365d=z[z.mag >= 5].sim.nunique() / NS, prob_M6_365d=z[z.mag >= 6].sim.nunique() / NS, prob_M45_30d=z[z.dia <= 30].sim.nunique() / NS, esperados_M45_365d=len(z) / NS))
dpt = pd.DataFrame(rows).sort_values("prob_M5_365d", ascending=False); dpt.round(4).to_csv(f"{OUT}/prob_departamentos.csv", index=False); print(dpt.head(15).round(3).to_string(index=False))
# rejilla 0,5 grados
g = 0.5; pts["gx"] = np.floor(pts.lon / g) * g; pts["gy"] = np.floor(pts.lat / g) * g; m5 = pts[pts.mag >= 5]
grid = m5.groupby(["gx", "gy"]).sim.nunique().div(NS).rename("prob_M5_365d").reset_index().merge(pts.groupby(["gx", "gy"]).size().div(NS).rename("esperados_M45_365d").reset_index(), how="outer").fillna(0); grid.round(4).to_csv(f"{OUT}/prob_rejilla_05.csv", index=False)
# ---------- validación: ajuste hasta 2023, pronóstico 2024-2025 ----------
V = pd.read_pickle(f"{OUT}/simulaciones_val.pkl"); NV = int(V.sim.max()) + 1; full = pd.read_csv(f"{D}/sismos_usgs_colombia_1973_2026_M4.csv", parse_dates=["time"]); full = full[(full.type == "earthquake") & (full.mag >= 4.5)]
ob = full[(full.time >= "2024-01-01") & (full.time < "2025-12-31")]; cv = pd.read_csv(f"{OUT}/catalogo_etas_val.csv")
nsim = V.groupby("sim").size().reindex(range(NV)).fillna(0); n5 = V[V.mag >= 5].groupby("sim").size().reindex(range(NV)).fillna(0)
V["gx"] = np.floor(V.lon / g) * g; V["gy"] = np.floor(V.lat / g) * g; rate = V.groupby(["gx", "gy"]).size().div(NV)
gx = np.arange(-82, -66, g); gy = np.arange(-5, 14, g); cells = pd.MultiIndex.from_product([gx, gy], names=["gx", "gy"]); r_etas = rate.reindex(cells).fillna(0)
hist_rate = cv.assign(gx=np.floor(cv.longitude / g) * g, gy=np.floor(cv.latitude / g) * g).groupby(["gx", "gy"]).size().reindex(cells).fillna(0)
oc = ob.assign(gx=np.floor(ob.longitude / g) * g, gy=np.floor(ob.latitude / g) * g).groupby(["gx", "gy"]).size().reindex(cells).fillna(0); yb = (oc > 0).astype(int)
ll = lambda r: float((oc * np.log(r / r.sum() * oc.sum() + 1e-12) - r / r.sum() * oc.sum()).sum()); rs = r_etas + 1e-4; hs = hist_rate + 0.05
res["validacion_2024_2025"] = dict(observados_M45=int(len(ob)), simulados_M45_mediana=float(nsim.median()), simulados_M45_p05_p95=[float(nsim.quantile(.05)), float(nsim.quantile(.95))], observados_M5=int((ob.mag >= 5).sum()), simulados_M5_mediana=float(n5.median()), simulados_M5_p05_p95=[float(n5.quantile(.05)), float(n5.quantile(.95))],
    auc_espacial_etas=float(roc_auc_score(yb, r_etas)), auc_espacial_conteo_historico=float(roc_auc_score(yb, hist_rate)), ganancia_logverosimilitud_por_evento_vs_historico=float((ll(rs) - ll(hs)) / len(ob)), celdas_con_evento=int(yb.sum()), celdas=int(len(yb)))
json.dump(res, open(f"{OUT}/pronostico_resumen.json", "w"), indent=1, ensure_ascii=False); print(json.dumps({k: res[k] for k in ["colombia", "region_completa", "historico_colombia_por_ano", "historico_prob_poisson_365d", "validacion_2024_2025"]}, indent=1, ensure_ascii=False))
# ---------- Figura 1: mapa ----------
cmap = LinearSegmentedColormap.from_list("p", ["#f3f2ee", "#fbd9a8", AMAR, NARANJA, "#b8321a", "#5c0f0a"]); bounds = [0, .01, .03, .06, .10, .20, .35, 1.0001]; norm = BoundaryNorm(bounds, cmap.N)
fig = plt.figure(figsize=(19.2, 9.6), dpi=100); ax = fig.add_axes([0.0, 0.17, 0.52, 0.69]); ax.set_axis_off()
G = grid.pivot(index="gy", columns="gx", values="prob_M5_365d").reindex(index=gy, columns=gx).fillna(0)
ax.pcolormesh(np.append(gx, gx[-1] + g), np.append(gy, gy[-1] + g), G.values, cmap=cmap, norm=norm); dep.boundary.plot(ax=ax, color="#5a5955", linewidth=0.5)
rec = cat[cat.time >= "2026-05-01"]; ax.scatter(rec.longitude, rec.latitude, s=(rec.mag - 4) ** 2.6 * 6, facecolors="none", edgecolors=AZUL, linewidths=0.9, label="Sismos M ≥ 4,5 desde mayo de 2026")
for n in ["Bogotá", "Medellín", "Cali", "Bucaramanga", "Pasto", "Barranquilla", "Cúcuta", "Quibdó", "Villavicencio"]:
    la, lo = CI[n]; ax.plot(lo, la, "o", ms=3.5, color=TXT); ax.annotate(n, (lo, la), xytext=(5, 3), textcoords="offset points", fontsize=9.5, color=TXT, path_effects=[pe.withStroke(linewidth=2.5, foreground=SURF)])
ax.set_xlim(-82, -66); ax.set_ylim(-4.6, 13.2); ax.set_aspect(1.0); ax.legend(frameon=False, fontsize=10, loc="lower right")
cax = fig.add_axes([0.07, 0.155, 0.38, 0.018]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal", ticks=bounds[:-1] + [1]); cb.set_ticklabels(["0", "1 %", "3 %", "6 %", "10 %", "20 %", "35 %", "100 %"]); cb.outline.set_visible(False)
cb.set_label("Probabilidad de al menos un sismo de magnitud 5 o más en la celda (0,5° × 0,5°) en los próximos 12 meses", color=TXT2, fontsize=10)
ax2 = fig.add_axes([0.60, 0.54, 0.38, 0.28]); c12 = res["colombia"]; lab = ["M ≥ 5", "M ≥ 5,5", "M ≥ 6", "M ≥ 6,5", "M ≥ 7"]; ks = ["M5.0", "M5.5", "M6.0", "M6.5", "M7.0"]; xx = np.arange(5); w = .26
for i, (hh, cc, ln) in enumerate([("30d", GRIS, "30 días"), ("90d", AMAR, "90 días"), ("365d", NARANJA, "12 meses")]):
    v = [c12[hh][k]["prob"] * 100 for k in ks]; ax2.bar(xx + (i - 1) * w, v, w, color=cc, label=ln)
    for xi, a in zip(xx, v): ax2.text(xi + (i - 1) * w, a + 1.5, f"{a:.0f}", ha="center", fontsize=9.5, color=TXT)
ax2.set_xticks(xx); ax2.set_xticklabels(lab); ax2.set_ylim(0, 112); ax2.set_yticks([0, 25, 50, 75, 100]); ax2.set_yticklabels(["0", "25 %", "50 %", "75 %", "100 %"]); ax2.grid(axis="x", visible=False); ax2.legend(frameon=False, fontsize=10, ncol=3, loc="upper right")
ax2.set_title("Probabilidad de al menos un sismo en Colombia (territorio y hasta 50 km mar afuera)", loc="left", fontsize=12, color=TXT, pad=8)
ax3 = fig.add_axes([0.66, 0.17, 0.32, 0.27]); t = ciu.head(10).iloc[::-1]; yy = np.arange(len(t)); ax3.barh(yy + .2, t.prob_hist_M5_365d * 100, .38, color=GRISC, label="Promedio histórico 1973-2026"); ax3.barh(yy - .2, t.prob_M5_365d * 100, .38, color=NARANJA, label="Pronóstico ETAS")
for yi, a in zip(yy, t.prob_M5_365d * 100): ax3.text(a + 1, yi - .2, f"{a:.0f} %", va="center", fontsize=9.5, color=TXT2)
ax3.set_yticks(yy); ax3.set_yticklabels(t.ciudad); ax3.grid(axis="y", visible=False); ax3.legend(frameon=False, fontsize=9.5, loc="lower right"); ax3.set_xlabel("Probabilidad en 12 meses (%)")
ax3.set_title("Sismo de magnitud 5 o más a menos de 100 km de la ciudad", loc="left", fontsize=12, color=TXT, pad=8)
fig.text(0.012, 0.955, "Colombia tembló más en 2026 y el modelo espera que siga temblando: dónde y con qué probabilidad", fontsize=21, color=TXT, va="center")
v = res["validacion_2024_2025"]
fig.text(0.012, 0.905, f"Pronóstico probabilístico de sismos para los próximos 12 meses con un modelo ETAS (secuencias de réplicas de tipo epidémico) ajustado a {P['N']:,} sismos de magnitud 4;5 o más desde 1973.\nNo predice fecha ni lugar exactos de un sismo: estima la probabilidad de que ocurra al menos uno en cada zona.".replace(",", ".").replace(";", ","), fontsize=12, color=TXT2, va="center")
fig.text(0.012, 0.045, f"Fuente: Elaboración propia con el catálogo del Servicio Geológico de Estados Unidos (USGS, ComCat) a 13 de noviembre de 2026 y el Marco Geoestadístico Nacional del DANE. {NS:,} simulaciones.".replace(",", ".") + f"\nValidación: con datos hasta 2023 el modelo pronosticó {v['simulados_M45_mediana']:.0f} sismos M ≥ 4,5 para 2024-2025 (intervalo 90 %: {v['simulados_M45_p05_p95'][0]:.0f} a {v['simulados_M45_p05_p95'][1]:.0f}) y ocurrieron {v['observados_M45']}.\nEl modelo es bidimensional: no distingue la profundidad, y un sismo profundo (nido de Bucaramanga, 150 km) se siente menos que uno superficial de la misma magnitud. Gustavo A Romero.", fontsize=9.5, color=TXT2)
fig.savefig(f"{RAIZ}/fig1_pronostico_sismos.png"); print("ok")
