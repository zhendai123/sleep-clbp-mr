# -*- coding: utf-8 -*-
"""OSA -> CLBP 重跑(500kb距离窗prune版) —— 放 D:\gwas 运行: python osa_pruned_rerun.py
输入(均已在你硬盘): 
  暴露: finngen_R11_G6_SLEEPAPNO_filtered_mlogp5.csv
  结局: LBP_extracted_OSA.csv (CLBP在OSA工具位点的效应)
输出: 独立工具数 / IVW / MR-Egger / FTO区删除敏感性 / 加权中位数
"""
import pandas as pd, numpy as np
from scipy import stats
import statsmodels.api as sm

W = 500_000  # 距离窗 500kb

# ---------- 1. 暴露工具: GWS + 距离窗prune ----------
osa = pd.read_csv("finngen_R11_G6_SLEEPAPNO_filtered_mlogp5.csv")
osa = osa[osa.pval < 5e-8].copy().sort_values("pval")
print(f"GWS候选: {len(osa)}")

kept = []
for chrom, grp in osa.groupby("#chrom"):
    taken = []          # 该染色体已保留位点
    for _, r in grp.sort_values("pval").iterrows():
        p = int(r["pos"])
        if all(abs(p - q) > W for q in taken):
            taken.append(p); kept.append(r)
iv = pd.DataFrame(kept)
print(f"500kb距离窗后独立工具: {len(iv)}")

iv = iv.rename(columns={"rsids": "snp", "beta": "bx", "sebeta": "sex", "af_alt": "eaf", "alt": "ea", "ref": "oa", "#chrom": "chrom"})
iv = iv[["snp", "chrom", "bx", "sex", "pval", "eaf", "ea", "oa"]]

# ---------- 2. 结局效应匹配 ----------
out = pd.read_csv("LBP_extracted_OSA.csv")
out["snp"] = out["rsids"].astype(str)
m = iv.merge(out, on="snp", suffixes=("", "_o"))
print(f"结局匹配上: {len(m)}")

# 等位基因对齐: 暴露EA=alt; 结局EA同为alt -> 同向; 暴露EA=结局ref -> 反向; 否则丢弃
def align(r):
    if r.ea == r.alt and r.oa == r.ref: return 1
    if r.ea == r.ref and r.oa == r.alt: return -1
    return 0
m["flip"] = m.apply(align, axis=1)
m = m[m.flip != 0].copy()
m["pal"] = m.apply(lambda r: {r.ea, r.oa} in [{"A","T"},{"C","G"}], axis=1)
m = m[~(m.pal & m.eaf.astype(float).between(0.42, 0.58))]
m["by"] = m.beta.astype(float) * m.flip
m["sy"] = m.sebeta.astype(float)
m = m.reset_index(drop=True)
print(f"等位基因对齐+去模糊回文后进入MR: {len(m)}")

# ---------- 3. MR估计 ----------
def ivw(bx, by, sy):
    w = bx**2/sy**2
    b = (w*by/bx).sum()/w.sum()
    Q = (w*(by/bx - b)**2).sum()
    se = np.sqrt(1/w.sum())
    se_re = se*np.sqrt(max(1, Q/(len(bx)-1)))
    p = 2*stats.norm.sf(abs(b/se_re))
    return b, se_re, p, Q, stats.chi2.sf(Q, len(bx)-1)

def egger(bx, by, sy):
    X = sm.add_constant(bx); Wt = 1/sy**2
    Xm = X*np.sqrt(Wt)[:,None]; ym = by*np.sqrt(Wt)
    c,_,_,_ = np.linalg.lstsq(Xm, ym, rcond=None)
    r = ym - Xm@c; s2 = (r**2).sum()/(len(bx)-2)
    cov = s2*np.linalg.lstsq(Xm.T@Xm, np.eye(2), rcond=None)[0]
    se = np.sqrt(np.diag(cov))
    return c[1], se[1], 2*stats.t.sf(abs(c[1]/se[1]), len(bx)-2), c[0], 2*stats.t.sf(abs(c[0]/se[0]), len(bx)-2)

def wmedian(bx, by, sy):
    bj = by/bx; wj = (bx/sy)**2
    o = np.argsort(bj); bj, wj = bj[o], wj[o]
    cw = np.cumsum(wj) - 0.5*wj; cw /= wj.sum()
    b = np.interp(0.5, cw, bj)
    # bootstrap SE
    rng = np.random.default_rng(42); bs = []
    for _ in range(2000):
        idx = rng.integers(0, len(bj), len(bj))
        bo, wo = bj[idx], wj[idx]; oo = np.argsort(bo); bo, wo = bo[oo], wo[oo]
        cwo = np.cumsum(wo) - 0.5*wo; cwo /= wo.sum()
        bs.append(np.interp(0.5, cwo, bo))
    se = np.std(bs)
    return b, se, 2*stats.norm.sf(abs(b/se))

bx, by, sy = m.bx.astype(float).values, m.by.values, m.sy.values
b, se, p, Q, Qp = ivw(bx, by, sy)
print(f"\nIVW(随机效应): OR={np.exp(b):.3f} 95%CI {np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f} P={p:.2e}  Q P={Qp:.3f}")
eb, ese, ep, ic, icp = egger(bx, by, sy)
print(f"MR-Egger: OR={np.exp(eb):.3f} P={ep:.2e}  截距P={icp:.3f}")
mb, mse, mp = wmedian(bx, by, sy)
print(f"加权中位数: OR={np.exp(mb):.3f} 95%CI {np.exp(mb-1.96*mse):.3f}-{np.exp(mb+1.96*mse):.3f} P={mp:.2e}")

# FTO区敏感性 (chr16 52.5-54.5Mb 删除)
fto = m[~((m.chrom==16) & m.pos.between(52_500_000, 54_500_000))]
if len(fto) >= 10:
    b2, se2, p2, _, _ = ivw(fto.bx.astype(float).values, fto.by.values, fto.sy.values)
    print(f"删FTO区后 IVW: OR={np.exp(b2):.3f} 95%CI {np.exp(b2-1.96*se2):.3f}-{np.exp(b2+1.96*se2):.3f} P={p2:.2e} (n={len(fto)})")
else:
    print(f"删FTO区后工具数不足: {len(fto)}")

input("\n按回车退出...")
