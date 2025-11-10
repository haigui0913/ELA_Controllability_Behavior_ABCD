# ==== 加载包 ====
library(lavaan)
library(dplyr)
library(purrr)
library(tidyr)
library(glue)

# ==== 载入数据 ====
load("C:/future/ABCD_data/20250825/7_lmm_results_longitudinal/my_all_data_with_cbcl_control_year2_5172.RData")
data_raw <- merged_df2
data_raw <- data_raw %>% filter(!is.na(genetic_pc_1))
data_raw <- data_raw %>% filter(!is.na(cbcl_scr_syn_internal_t_Year2))

# ==== 设置变量 ====
brain_vars      <- colnames(data_raw)[4:13]       # Year1 脑网络变量
brain_vars_y2   <- colnames(data_raw)[97:106]      # Year2 脑网络变量   nihtb [78:87]

cog_vars        <- c("cbcl_scr_syn_internal_t","cbcl_scr_syn_external_t","cbcl_scr_syn_totprob_t")        # Year1 认知变量（可扩展）c("nihtbx_cryst_fc")  
cog_vars_y2     <- paste0(cog_vars, "_Year2")     # Year2 认知变量

covariates <- c("age", "sex", "race", "handness", "site", "pubertyZ", 
                "genetic_pc_1", "genetic_pc_2", "genetic_pc_3", 
                "genetic_pc_4", "genetic_pc_5")

categorical_vars <- c("sex", "race", "handness", "site")
data_raw[categorical_vars] <- lapply(data_raw[categorical_vars], as.factor)

# ==== 残差化函数 ====
residualize <- function(varname, data, covars) {
  fml <- reformulate(covars, response = varname)
  model <- lm(fml, data = data)
  return(residuals(model))
}

# ==== 执行残差化 ====
all_vars <- c(brain_vars, brain_vars_y2, cog_vars, cog_vars_y2)
data_resid <- data.frame(subject_id = data_raw$subject_id)

for (v in all_vars) {
  data_resid[[v]] <- residualize(v, data_raw, covariates)
}

# ==== 可选：对残差标准化（z-score） ====
data_resid[all_vars] <- lapply(data_resid[all_vars], scale)


# ==== 建模函数 ====
run_clpm_resid <- function(brain, brain_y2, cog, cog_y2, data) {
  model_str <- glue('
    {brain_y2} ~ a1*{brain} + c1*{cog}
    {cog_y2} ~ a2*{cog} + c2*{brain}
    
    {brain} ~~ {cog}
    {brain_y2} ~~ {cog_y2}
  ')
  
  fit <- tryCatch({
    sem(model_str, data = data, missing = "ML", estimator = "MLR")
  }, error = function(e) return(NULL))
  
  if (is.null(fit)) return(NULL)
  
  result <- parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "~") %>%
    mutate(brain = brain, cog = cog)
  return(result)
}

# ==== 批量建模 ====
all_results <- list()

for (i in seq_along(brain_vars)) {
  for (j in seq_along(cog_vars)) {
    res <- run_clpm_resid(
      brain     = brain_vars[i],
      brain_y2  = brain_vars_y2[i],
      cog       = cog_vars[j],
      cog_y2    = cog_vars_y2[j],
      data      = data_resid
    )
    if (!is.null(res)) {
      key <- paste0(brain_vars[i], "_", cog_vars[j])
      all_results[[key]] <- res
    }
  }
}

# ==== 合并 & FDR 校正 ====
df_clpm_all <- bind_rows(purrr::compact(all_results), .id = "id") %>%
  select(id, lhs, rhs, est, se, z, pvalue, std.all, brain, cog) %>%
  mutate(p_fdr = p.adjust(pvalue, method = "fdr"))

df_clpm_sig <- df_clpm_all %>% filter(p_fdr < 0.05)

# ==== 保存结果 ====
outdir <- "C:/future/ABCD_data/20250825/7_lmm_results_longitudinal/"
write.csv(df_clpm_all, file.path(outdir, "cbcl_clpm_all_residualized_scaled.csv"), row.names = FALSE)
write.csv(df_clpm_sig, file.path(outdir, "cbcl_clpm_fdr_sig_residualized_scaled.csv"), row.names = FALSE)

