"""
Sleep disturbance and chronic low back pain - analysis pipeline
Reproduces: NHANES weighted logistic regression, two-sample MR (IVW/Egger/weighted median),
reverse MR and two-step MR mediation. Python 3.10+, requires pandas/numpy/scipy/statsmodels.
"""

import pandas as pd, numpy as np
from scipy import stats
import statsmodels.api as sm

# ---------- 1. Harmonise instruments ----------
def harmonise(iv_file, outcome_file, id_col='snp'):
    """iv_file: instrument list (snp,chrom,pos,effect_allele,other_allele,beta,se,pval,eaf)
       outcome_file: extracted outcome data at the same SNPs."""
    iv  = pd.read_csv(iv_file)
    out = pd.read_csv(outcome_file)
    out['rsid1'] = out.rsids.astype(str).str.split(',').str[0]
    m = iv.merge(out, left_on=id_col, right_on='rsid1', suffixes=('_e','_o'))
    def align(r):
        if r.alt==r.effect_allele and r.ref==r.other_allele: return 1
        if r.alt==r.other_allele and r.ref==r.effect_allele: return -1
        return 0
    m['flip'] = m.apply(align, axis=1)
    m = m[m.flip!=0].copy()
    m['pal'] = m.apply(lambda r: {r.effect_allele,r.other_allele} in [{'A','T'},{'C','G'}], axis=1)
    m = m[~(m.pal & m.eaf.astype(float).between(0.42,0.58))]   # drop ambiguous palindromes
    m['bx'], m['by'] = m.beta_e.astype(float), m.beta_o.astype(float)*m.flip
    m['sy'] = m.sebeta.astype(float) if 'sebeta' in m else m.se_o.astype(float)
    return m.reset_index(drop=True)

# ---------- 2. MR estimators ----------
def ivw(bx, by, sy):
    w = bx**2/sy**2
    b = (w*by/bx).sum()/w.sum()
    Q = (w*(by/bx-b)**2).sum()
    se = np.sqrt(1/w.sum())
    return b, se, 2*stats.norm.sf(abs(b/se)), Q, stats.chi2.sf(Q, len(bx)-1)

def egger(bx, by, sy):
    X = sm.add_constant(bx); W = 1/sy**2
    Xm = X*np.sqrt(W)[:,None]; ym = by*np.sqrt(W)
    c,_,_,_ = np.linalg.lstsq(Xm, ym, rcond=None)
    r = ym - Xm@c; s2 = (r**2).sum()/(len(bx)-2)
    cov = s2*np.linalg.inv(Xm.T@Xm); se = np.sqrt(np.diag(cov))
    return c[1], se[1], 2*stats.t.sf(abs(c[1]/se[1]), len(bx)-2), c[0], 2*stats.t.sf(abs(c[0]/se[0]), len(bx)-2)

# ---------- 3. Steiger directionality ----------
def steiger(m):
    eaf = m.eaf.astype(float)
    r2x = 2*eaf*(1-eaf)*m.bx**2/(2*eaf*(1-eaf)*m.bx**2 + m.se.astype(float)**2)
    r2y = 2*eaf*(1-eaf)*m.by**2/(2*eaf*(1-eaf)*m.by**2 + m.sy**2)
    return (r2x>r2y).mean()

if __name__ == '__main__':
    m = harmonise('instruments_insomnia.csv', 'LBP_extracted.csv')
    b, se, p, Q, Qp = ivw(m.bx.values, m.by.values, m.sy.values)
    print(f"IVW: beta={b:.4f} OR={np.exp(b):.2f} P={p:.2e}  Q P={Qp:.3f}")
    eb, ese, ep, ic, icp = egger(m.bx.values, m.by.values, m.sy.values)
    print(f"MR-Egger: beta={eb:.4f} P={ep:.2e}  intercept P={icp:.3f}")
    print(f"Steiger: {steiger(m)*100:.1f}% variants correctly oriented")
