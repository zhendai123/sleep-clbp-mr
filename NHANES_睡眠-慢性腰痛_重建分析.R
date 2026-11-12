# =============================================================================
# NHANES 2009–2010 睡眠障碍与慢性腰痛(CLBP) —— 重建分析脚本(终版, 与稿件口径一致)
# R 版本: 4.4.1
#
# 口径(已与稿件 Methods 逐项核对):
#   结局 CLBP(双条目验证定义): ARQ022AD==1(仍有低背痛) 且 ARQ024D==1(连续≥3个月)
#   暴露 复合表型 any_slp: 短睡(SLD010H<7) | 长睡(SLD010H>9)
#                          | 曾告知医生睡眠问题(SLQ050==1) | 曾诊断睡眠障碍(SLQ060==1)
#   模型 M1 未调整; M2 = 年龄+性别+种族+教育+PIR;
#        M3 = M2 + BMI + SBP(首个可用读数) + 吸烟 + 饮酒 + PHQ-9 + 关节炎
#   抽样设计: WTMEC2YR + SDMVSTRA + SDMVPSU (survey 包, Taylor 线性化)
# =============================================================================

library(dplyr)
library(survey)

dir0 <- "."   # NHANES CSV 所在目录

# ---------- 1. 读取与合并 ----------
read_xpt_csv <- function(f) read.csv(file.path(dir0, f), stringsAsFactors = FALSE)

demo <- read_xpt_csv("DEMO_F.csv") %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, DMDEDUC2, INDFMPIR,
         WTMEC2YR, SDMVSTRA, SDMVPSU)
slq <- read_xpt_csv("SLQ_F.csv") %>% select(SEQN, SLD010H, SLQ050, SLQ060)
arq <- read_xpt_csv("ARQ_F.csv") %>% select(SEQN, ARQ010, ARQ020D, ARQ022AD, ARQ024D)
bmx <- read_xpt_csv("BMX_F.csv") %>% select(SEQN, BMXBMI)
bpx <- read_xpt_csv("BPX_F.csv") %>% select(SEQN, BPXSY1, BPXSY2, BPXSY3)
smq <- read_xpt_csv("SMQ_F.csv") %>% select(SEQN, SMQ020, SMQ040)
alq <- read_xpt_csv("ALQ_F.csv") %>% select(SEQN, ALQ130)
mcq <- read_xpt_csv("MCQ_F.csv") %>% select(SEQN, MCQ160A)
cbc <- read_xpt_csv("CBC_F.csv") %>% select(SEQN, LBDNENO, LBDLYMNO, LBXPLTSI)
crp <- read_xpt_csv("CRP_F.csv") %>% select(SEQN, LBXCRP)
dpq <- read_xpt_csv("DPQ_F.csv") %>%
  select(SEQN, DPQ010, DPQ020, DPQ030, DPQ040, DPQ050,
         DPQ060, DPQ070, DPQ080, DPQ090)

dat <- demo %>%
  left_join(slq, by = "SEQN") %>% left_join(arq, by = "SEQN") %>%
  left_join(bmx, by = "SEQN") %>% left_join(bpx, by = "SEQN") %>%
  left_join(smq, by = "SEQN") %>% left_join(alq, by = "SEQN") %>%
  left_join(mcq, by = "SEQN") %>% left_join(cbc, by = "SEQN") %>%
  left_join(crp, by = "SEQN") %>% left_join(dpq, by = "SEQN")

# ---------- 2. 队列构建 (n = 4,998) ----------
dat <- dat %>%
  filter(!is.na(ARQ010), WTMEC2YR > 0) %>%
  filter((ARQ022AD %in% c(1, 2) & ARQ024D %in% c(1, 2)) | is.na(ARQ020D) | ARQ010 == 2)

# ---------- 3. 变量定义 ----------
phq_items <- paste0("DPQ0", 1:9, "0")
dat <- dat %>%
  mutate(
    # 结局: 双条目 CLBP (600 例)
    clbp   = as.integer(ARQ022AD == 1 & ARQ024D == 1),
    # 暴露
    short7 = as.integer(SLD010H < 7),
    long9  = as.integer(SLD010H > 9),
    trouble = as.integer(SLQ050 == 1),
    slpdx  = as.integer(SLQ060 == 1),
    any_slp = as.integer(short7 == 1 | long9 == 1 | trouble == 1 | slpdx == 1),
    # 协变量
    age    = RIDAGEYR,
    female = as.integer(RIAGENDR == 2),
    race   = factor(RIDRETH1),
    edu    = factor(DMDEDUC2),
    pir    = INDFMPIR,
    bmi    = BMXBMI,
    # SBP: 首个可用读数 (BPXSY1 -> 2 -> 3)
    sbp    = case_when(!is.na(BPXSY1) ~ BPXSY1,
                       !is.na(BPXSY2) ~ BPXSY2,
                       TRUE ~ BPXSY3),
    # 吸烟: 0=从不 1=既往 2=现在
    smoke  = case_when(SMQ020 == 2 ~ 0,
                       SMQ040 == 3 ~ 1,
                       SMQ040 %in% c(1, 2) ~ 2,
                       TRUE ~ NA_real_),
    alq130 = ALQ130,                    # 过去一年饮酒日平均杯数
    arth   = ifelse(MCQ160A %in% c(1, 2), as.numeric(MCQ160A == 1), NA_real_),
    # PHQ-9: 9 条目 (7/9 视为缺失), 全部有效才求和; dep = PHQ-9 >= 10
    phq9   = rowSums(across(all_of(phq_items),
                            ~ ifelse(.x <= 3, .x, NA_real_))),
    dep    = as.integer(phq9 >= 10),
    # 炎症指标
    crp    = LBXCRP,
    nlr    = LBDNENO / LBDLYMNO,
    sii    = LBXPLTSI * LBDNENO / LBDLYMNO
  )
stopifnot(nrow(dat) == 4998, sum(dat$clbp) == 600)

# ---------- 4. 抽样设计对象 ----------
des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                     weights = ~WTMEC2YR, nest = TRUE, data = dat)

m3_vars <- c("age", "female", "race", "edu", "pir", "bmi", "sbp",
             "smoke", "alq130", "phq9", "arth", "crp")
fin <- dat[complete.cases(dat[, m3_vars]), ]   # n = 2,789 (抑郁亚组 n = 262)
des_fin <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                     weights = ~WTMEC2YR, nest = TRUE, data = fin)
stopifnot(nrow(fin) == 2789, sum(fin$dep) == 262)

# ---------- 5. Table 2: 5 暴露 x 3 模型 ----------
F2 <- "clbp ~ %s + age + female + race + edu + pir"
F3 <- paste("clbp ~ %s +", "age + female + race + edu + pir +",
            "bmi + sbp + smoke + alq130 + phq9 + arth")
expos <- c("any_slp", "short7", "long9", "trouble", "slpdx")
orci <- function(fit, term) {
  b <- coef(fit)[term]; ci <- confint(fit)[term, ]
  sprintf("%.2f (%.2f\u2013%.2f)", exp(b), exp(ci[1]), exp(ci[2]))
}
tab2 <- lapply(expos, function(e) {
  m1 <- svyglm(as.formula(sprintf("clbp ~ %s", e)), des_all, family = quasibinomial)
  m2 <- svyglm(as.formula(sprintf(F2, e)), des_all, family = quasibinomial)
  m3 <- svyglm(as.formula(sprintf(F3, e)), des_fin, family = quasibinomial)
  c(orci(m1, e), orci(m2, e), orci(m3, e),
    sprintf("%.3f", summary(m3)$coefficients[e, 4]))
})
tab2 <- do.call(rbind, tab2)
rownames(tab2) <- expos
colnames(tab2) <- c("Model1", "Model2", "Model3", "P_M3")
print(tab2)

# ---------- 6. 抑郁分层与交互 ----------
des_nd <- subset(des_fin, dep == 0)   # n = 2,527
des_dp <- subset(des_fin, dep == 1)   # n = 262
f_nd <- svyglm(as.formula(gsub("%s", "any_slp", F3)), des_nd, family = quasibinomial)
f_dp <- svyglm(as.formula(gsub("%s", "any_slp", F3)), des_dp, family = quasibinomial)
f_ix <- svyglm(clbp ~ any_slp * dep + age + female + race + edu + pir +
                 bmi + sbp + smoke + alq130 + phq9 + arth,
               des_fin, family = quasibinomial)
cat("非抑郁:", orci(f_nd, "any_slp"), "\n")
cat("抑郁:  ", orci(f_dp, "any_slp"), "\n")
cat("交互:  ", orci(f_ix, "any_slp:dep"),
    " P =", sprintf("%.3f", summary(f_ix)$coefficients["any_slp:dep", 4]), "\n")

# ---------- 7. RCS 剂量反应 (3 knots: 10/50/90 分位) ----------
library(rms)
spl <- subset(fin, SLD010H >= 2 & SLD010H <= 14)
dd_ <- datadist(spl); options(datadist = "dd_")
des_spl <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                     weights = ~WTMEC2YR, nest = TRUE, data = spl)
f_rcs <- svyglm(clbp ~ rcs(SLD010H, 3) + age + female + race + edu + pir +
                  bmi + sbp + smoke + alq130 + phq9 + arth,
                des_spl, family = quasibinomial)
f_lin <- svyglm(clbp ~ SLD010H + age + female + race + edu + pir +
                  bmi + sbp + smoke + alq130 + phq9 + arth,
                des_spl, family = quasibinomial)
regTermTest(f_rcs, "rcs(SLD010H, 3)SLD010H'")   # 非线性 Wald 检验 (P = 0.11)

# ---------- 8. 炎症中介 (NLR / SII / CRP; 500 次 bootstrap) ----------
set.seed(2024)
B <- 500
med_boot <- function(med) {
  medf <- fin[complete.cases(fin[, med]), ]
  point <- {
    fa <- glm(as.formula(paste(med, "~ any_slp + age + female + race + edu + pir +",
                               "bmi + sbp + smoke + alq130 + phq9 + arth")),
              data = medf, weights = WTMEC2YR)
    fb <- glm(as.formula(paste("clbp ~ any_slp +", med, "+ age + female + race + edu +",
                               "pir + bmi + sbp + smoke + alq130 + phq9 + arth")),
              data = medf, family = binomial, weights = WTMEC2YR)
    coef(fa)["any_slp"] * coef(fb)[med]
  }
  boots <- replicate(B, {
    i <- sample.int(nrow(medf), replace = TRUE); dd <- medf[i, ]
    fa <- glm(as.formula(paste(med, "~ any_slp + age + female + race + edu + pir +",
                               "bmi + sbp + smoke + alq130 + phq9 + arth")),
              data = dd, weights = WTMEC2YR)
    fb <- glm(as.formula(paste("clbp ~ any_slp +", med, "+ age + female + race + edu +",
                               "pir + bmi + sbp + smoke + alq130 + phq9 + arth")),
              data = dd, family = binomial, weights = WTMEC2YR)
    coef(fa)["any_slp"] * coef(fb)[med]
  })
  c(indirect = point, quantile(boots, c(.025, .975), na.rm = TRUE))
}
for (m in c("nlr", "sii", "crp")) { cat(m, ": "); print(med_boot(m)) }
# 预期: 三个间接效应 95% CI 均跨零 (重建验证: 与稿件结论一致)

write.csv(dat, "NHANES_睡眠-慢性腰痛_分析数据集.csv", row.names = FALSE)
