# -*- coding: utf-8 -*-
"""
NHANES 2009–2010 睡眠障碍—慢性腰痛 重建验证脚本(终版, 与稿件口径一致)
Python 3.12 | 用 statsmodels GLM(var_weights) + 手动 PSU 级得分聚合
复现 R survey 包 Taylor 线性化设计方差。运行: python3 nhanes_sleep_clbp_verification.py [csv目录]
"""
import sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, patsy, statsmodels.api as sm
from scipy import stats as st_

DIR = sys.argv[1] if len(sys.argv) > 1 else "."
def ld(f): return pd.read_csv(os.path.join(DIR, f), encoding="latin1", low_memory=False)

# ---------- 1. 合并 ----------
demo = ld("DEMO_F.csv")[["SEQN","RIDAGEYR","RIAGENDR","RIDRETH1","DMDEDUC2","INDFMPIR","WTMEC2YR","SDMVSTRA","SDMVPSU"]]
d = (demo
     .merge(ld("SLQ_F.csv")[["SEQN","SLD010H","SLQ050","SLQ060"]], on="SEQN", how="left")
     .merge(ld("ARQ_F.csv")[["SEQN","ARQ010","ARQ020D","ARQ022AD","ARQ024D"]], on="SEQN", how="left")
     .merge(ld("BMX_F.csv")[["SEQN","BMXBMI"]], on="SEQN", how="left")
     .merge(ld("BPX_F.csv")[["SEQN","BPXSY1","BPXSY2","BPXSY3"]], on="SEQN", how="left")
     .merge(ld("SMQ_F.csv")[["SEQN","SMQ020","SMQ040"]], on="SEQN", how="left")
     .merge(ld("ALQ_F.csv")[["SEQN","ALQ130"]], on="SEQN", how="left")
     .merge(ld("MCQ_F.csv")[["SEQN","MCQ160A"]], on="SEQN", how="left")
     .merge(ld("CBC_F.csv")[["SEQN","LBDNENO","LBDLYMNO","LBXPLTSI"]], on="SEQN", how="left")
     .merge(ld("CRP_F.csv")[["SEQN","LBXCRP"]], on="SEQN", how="left")
     .merge(ld("DPQ_F.csv")[["SEQN"]+[f"DPQ0{i}0" for i in range(1,10)]], on="SEQN", how="left"))

# ---------- 2. 队列 (n=4,998, CLBP=600) ----------
d = d[(d.ARQ010.notna()) & (d.WTMEC2YR > 0)]
d = d[((d.ARQ022AD.isin([1,2])) & (d.ARQ024D.isin([1,2]))) | (d.ARQ020D.isna()) | (d.ARQ010 == 2)].copy()
d["clbp"] = ((d.ARQ022AD == 1) & (d.ARQ024D == 1)).astype(int)
assert len(d) == 4998 and d.clbp.sum() == 600

# ---------- 3. 变量 ----------
d["short7"]  = (d.SLD010H < 7).astype(int)
d["long9"]   = (d.SLD010H > 9).astype(int)
d["trouble"] = (d.SLQ050 == 1).astype(int)
d["slpdx"]   = (d.SLQ060 == 1).astype(int)
d["any_slp"] = ((d.short7 | d.long9 | d.trouble | d.slpdx) > 0).astype(int)
d["age"] = d.RIDAGEYR; d["female"] = (d.RIAGENDR == 2).astype(int)
d["race"] = d.RIDRETH1; d["edu"] = d.DMDEDUC2; d["pir"] = d.INDFMPIR; d["bmi"] = d.BMXBMI
d["sbp"] = d.BPXSY1.fillna(d.BPXSY2).fillna(d.BPXSY3)          # 首个可用读数
d["smoke"] = np.where(d.SMQ020 == 2, 0, np.where(d.SMQ040 == 3, 1,
             np.where(d.SMQ040.isin([1,2]), 2, np.nan)))        # 0从不/1既往/2现在
d["alq130"] = d.ALQ130
d["arth"] = (d.MCQ160A == 1).astype(float).where(d.MCQ160A.isin([1,2]))
items = [f"DPQ0{i}0" for i in range(1,10)]
d["phq9"] = d[items].where(d[items] <= 3).sum(axis=1, min_count=9)
d["dep"] = (d.phq9 >= 10).astype(float).where(d.phq9.notna())
d["crp"] = d.LBXCRP; d["nlr"] = d.LBDNENO / d.LBDLYMNO
d["sii"] = d.LBXPLTSI * d.LBDNENO / d.LBDLYMNO

# ---------- 4. 复杂抽样 logistic (survey 复现) ----------
def svy_logit(data, formula):
    y, X = patsy.dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    dd = data.loc[y.index]; yv = y.iloc[:,0].values.astype(float); w = dd.WTMEC2YR.values
    res = sm.GLM(yv, X, family=sm.families.Binomial(), var_weights=w).fit()
    p = res.fittedvalues.values
    u = (yv - p)[:, None] * X.values * w[:, None]
    t = pd.DataFrame(u, columns=X.columns); t["st"] = dd.SDMVSTRA.values; t["ps"] = dd.SDMVPSU.values
    D = np.zeros((X.shape[1],) * 2)
    for _, g in t.groupby("st"):
        S = g.groupby("ps")[X.columns].sum().values; Sc = S - S.mean(0)
        D += S.shape[0] / (S.shape[0] - 1) * Sc.T @ Sc
    I = (X.values * (w * p * (1 - p))[:, None]).T @ X.values
    V = np.linalg.inv(I) @ D @ np.linalg.inv(I)
    return res, V, X.columns.tolist()

def or_ci(res, V, names, term):
    i = names.index(term); b = res.params[i]; se = np.sqrt(V[i, i])
    return np.exp(b), np.exp(b - 1.959964 * se), np.exp(b + 1.959964 * se), 2 * st_.norm.sf(abs(b / se))

def svy_ols(data, formula):
    y, X = patsy.dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    dd = data.loc[y.index]; yv = y.iloc[:,0].values.astype(float); w = dd.WTMEC2YR.values
    res = sm.GLM(yv, X, family=sm.families.Gaussian(), var_weights=w).fit()
    p = res.fittedvalues.values
    u = (yv - p)[:, None] * X.values * w[:, None]
    t = pd.DataFrame(u, columns=X.columns); t["st"] = dd.SDMVSTRA.values; t["ps"] = dd.SDMVPSU.values
    D = np.zeros((X.shape[1],) * 2)
    for _, g in t.groupby("st"):
        S = g.groupby("ps")[X.columns].sum().values; Sc = S - S.mean(0)
        D += S.shape[0] / (S.shape[0] - 1) * Sc.T @ Sc
    I = (X.values * w[:, None]).T @ X.values
    V = np.linalg.inv(I) @ D @ np.linalg.inv(I)
    return res, V, X.columns.tolist()

M2 = "age + female + C(race) + C(edu) + pir"
M3 = M2 + " + bmi + sbp + smoke + alq130 + phq9 + arth"
cc2 = d.dropna(subset=["age","female","race","edu","pir"])
fin = d.dropna(subset=["age","female","race","edu","pir","bmi","sbp","smoke","alq130","phq9","arth","crp"])
assert len(fin) == 2789 and int(fin.dep.sum()) == 262
f = lambda a: f"{a[0]:.2f} ({a[1]:.2f}\u2013{a[2]:.2f})"

print(f"队列 n={len(d)}  CLBP={d.clbp.sum()}  未加权={d.clbp.mean()*100:.1f}%  "
      f"加权={(d.clbp*d.WTMEC2YR).sum()/d.WTMEC2YR.sum()*100:.1f}%")
print(f"完全协变量 n={len(fin)}  抑郁 n={int(fin.dep.sum())}  "
      f"加权比例={(fin.dep*fin.WTMEC2YR).sum()/fin.WTMEC2YR.sum()*100:.2f}%\n")

print("== Table 2 (M1 n=4,998 / M2 n=4,513 / M3 n=2,789) ==")
for e, lab in [("any_slp","Any"),("short7","Short"),("long9","Long"),("trouble","Trouble"),("slpdx","Dx")]:
    a1 = or_ci(*svy_logit(d,  f"clbp ~ {e}"), e)
    a2 = or_ci(*svy_logit(cc2, f"clbp ~ {e} + {M2}"), e)
    a3 = or_ci(*svy_logit(fin, f"clbp ~ {e} + {M3}"), e)
    print(f"{lab:8s} {f(a1):18s} {f(a2):18s} {f(a3):18s} P={'<0.001' if a3[3]<.001 else f'{a3[3]:.3f}'}")

print("\n== 抑郁分层与交互 ==")
nd = or_ci(*svy_logit(fin[fin.dep == 0], f"clbp ~ any_slp + {M3}"), "any_slp")
dp = or_ci(*svy_logit(fin[fin.dep == 1], f"clbp ~ any_slp + {M3}"), "any_slp")
ix = or_ci(*svy_logit(fin, f"clbp ~ any_slp*dep + {M3}"), "any_slp:dep")
print(f"非抑郁 {f(nd)} | 抑郁 {f(dp)} | 交互 {f(ix)} P={ix[3]:.3f}")

print("\n== RCS (3 knots 10/50/90pct, 均值剖面预测) ==")
def rcs(x, k):
    # 3-knot 限制性立方样条第二基 (Harrell): (x-k1)³₊ - (x-k2)³₊
    return np.maximum(x - k[0], 0)**3 - np.maximum(x - k[1], 0)**3
S = fin.dropna(subset=["SLD010H"]).copy(); S = S[S.SLD010H.between(2, 14)]
kn = np.percentile(S.SLD010H, [10, 50, 90])
S["slin"] = S.SLD010H; S["snl"] = rcs(S.SLD010H.values, kn)
res2, V2, nm2 = svy_logit(S, f"clbp ~ slin + snl + {M3}")
i = nm2.index("snl"); p_nl = st_.chi2.sf(res2.params[i]**2 / V2[i, i], 1)
Xd = patsy.dmatrix(f"slin + snl + {M3}", S, return_type="dataframe")[nm2]
Xp = Xd.copy()
for c in Xd.columns:
    if c not in ("Intercept", "slin", "snl"): Xp[c] = Xd[c].mean()
pm = lambda h: float((1/(1+np.exp(-(Xp.assign(slin=h, snl=rcs(np.array([float(h)]),kn)[0]).values @ res2.params.values)))).mean()*100)
print(f"knots={kn.round(1)}  非线性 Wald P={p_nl:.3f}  5/7/10h 预测概率 = {pm(5):.1f}%/{pm(7):.1f}%/{pm(10):.1f}%")

print("\n== 中介 (a×b, delta法 CI) ==")
med_ = fin.dropna(subset=["nlr","sii"])
for med in ["nlr","sii","crp"]:
    ra, Va, nma = svy_ols(med_, f"{med} ~ any_slp + {M3}")
    ia = nma.index("any_slp"); a_, sa = ra.params[ia], np.sqrt(Va[ia, ia])
    rb, Vb, nmb = svy_logit(med_, f"clbp ~ any_slp + {med} + {M3}")
    ib = nmb.index(med); b_, sb = rb.params[ib], np.sqrt(Vb[ib, ib])
    ind = a_ * b_; sei = np.sqrt(b_**2*sa**2 + a_**2*sb**2)
    print(f"{med:4s} ab={ind:+.4f}  95%CI=({ind-1.96*sei:+.4f}, {ind+1.96*sei:+.4f})  跨零={abs(ind)<1.96*sei}")

d.to_csv("NHANES_睡眠-慢性腰痛_分析数据集.csv", index=False)
print("\n数据集已写出: NHANES_睡眠-慢性腰痛_分析数据集.csv (n=4,998)")
