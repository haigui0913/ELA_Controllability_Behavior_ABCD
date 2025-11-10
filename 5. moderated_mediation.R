library(lme4)
library(mediation)   # 简单/经典中介
library(processR)    # Hayes 风格的调节中介
library(dplyr)
library(purrr)
library(stringr)
library(lavaan)
library(readr)
library(openxlsx)
library(glue)



load ("C:\\future\\ABCD_data\\20250818\\euro_prsice_7970_factors_3967subjs\\euro_PRSice_norm_3976_with_testprs.RData")

# 确保分类变量为因子
data_eur$sex      <- as.factor(data_eur$sex)
data_eur$race <- as.factor(data_eur$race)
data_eur$handness <- as.factor(data_eur$handness)
data_eur$site     <- as.factor(data_eur$site)

# 协变量（与你之前一致）
ctrl_vars <- c("age_z","sex","pubertyZ_z","handness",
               "genetic_pc_1_z","genetic_pc_2_z","genetic_pc_3_z",
               "genetic_pc_4_z","genetic_pc_5_z")

# 你定义的变量集
behavior_vars <- c("cbcl_scr_syn_totprob_t_z","nihtbx_totalcomp_fc_z") #colnames(data_eur)[c(121:123, 140:142)]
brain_vars    <- c("MF") #,"FP","Mot","Vas","SAL","CBL"
env_vars      <- c("F4")      
prs_vars      <- c("PRS_adhd_z","PRS_asd_z")

# 1. site 层面残差化函数
residualize_by_site <- function(var){
  fml  <- as.formula(paste0(var, " ~ 1 + (1|site)"))
  fit  <- lmer(fml, data = data_eur, REML = TRUE)
  rcol <- paste0(var, "_resid_site")
  data_eur[[rcol]] <<- resid(fit)
  return(rcol)
}

residualize_by_site <- function(var_name, data) {
  fml <- as.formula(paste0(var_name, " ~ 1 + (1|site)"))
  fit <- lmer(fml, data = data, REML = TRUE)
  resid(fit)
}

# -------------------------
# 建议文件夹名
result_dir <- "C:\\future\\ABCD_data\\20250818\\euro_prsice_7970_factors_3967subjs\\moderated_mediation_results_test"
dir.create(result_dir, showWarnings = FALSE)

for (y_var in behavior_vars) {
  for (mi_var in brain_vars) {
    for (x_var in env_vars) {
      for (w_var in prs_vars) {
        
        # 构造文件名作为唯一标识
        fname <- paste0(y_var, "__", mi_var, "__", x_var, "__", w_var, ".xlsx")
        fpath <- file.path(result_dir, fname)
        
        # ✅ 跳过已存在文件（中断恢复核心）
        if (file.exists(fpath)) {
          cat("✅ 跳过已完成：", fname, "\n")
          next
        }
        
        cat("\n➡️ 正在分析：", fname, "\n")
        
        # 残差化脑和行为
        data_model <- data_eur
        data_model$Mi_resid <- residualize_by_site(mi_var, data_model)
        data_model$Y_resid  <- residualize_by_site(y_var,  data_model)
        
        df <- data_model[, c(x_var, w_var, "Mi_resid", "Y_resid")]
        colnames(df) <- c("X", "W", "Mi", "Y")
        df <- na.omit(df)
        if (nrow(df) < 30) next
        
        # ✅ 残差后的脑 & 行为变量需再次标准化
        df$Mi <- scale(df$Mi)
        df$Y  <- scale(df$Y)
        
        # 构造交互项
        df$XW  <- df$X * df$W
        df$MiW <- df$Mi * df$W
        
        # W 的估值（mean, ±1SD）
        w_mean <- mean(df$W)
        w_sd   <- sd(df$W)
        w_low  <- w_mean - w_sd
        w_high <- w_mean + w_sd
        
        # 构造lavaan模型（PROCESS Model 59）
        # === 构建 lavaan 模型语法 ===
        model_text <- glue("
          Mi ~ a1*X + a2*W + a3*X:W
          Y  ~ c1*X + c2*W + c3*X:W + b1*Mi + b2*Mi:W
          W ~ wmean*1
          W ~~ wvar*W

          indirect_mean := (a1 + a3*{w_mean}) * (b1 + b2*{w_mean})
          direct_mean   := c1 + c3*{w_mean}
          total_mean    := indirect_mean + direct_mean
          prop_mean     := indirect_mean / total_mean

          indirect_low  := (a1 + a3*{w_low}) * (b1 + b2*{w_low})
          direct_low    := c1 + c3*{w_low}
          total_low     := indirect_low + direct_low
          prop_low      := indirect_low / total_low

          indirect_high := (a1 + a3*{w_high}) * (b1 + b2*{w_high})
          direct_high   := c1 + c3*{w_high}
          total_high    := indirect_high + direct_high
          prop_high     := indirect_high / total_high
        ")

          # 拟合模型
        fit <- tryCatch({
          sem(model_text, data = df, se = "bootstrap", bootstrap = 5000)
        }, error = function(e) NULL)
        
        if (!is.null(fit) && lavInspect(fit, "converged")) {
          
          model_id <- gsub(".xlsx", "", fname) 
          
          est <- parameterEstimates(fit, boot.ci.type = "bca.simple", level = 0.95, ci = TRUE)
          est$model_id <- model_id
          
          # 添加拟合指标
          fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "srmr", "tli"))
          fit_info <- as.data.frame(as.list(fit_measures)) 
          fit_info$model_id <- model_id
          
          # 保存结果和拟合度
          write.xlsx(est, file = fpath, rowNames = FALSE)
          write.xlsx(fit_info, file = sub(".xlsx", "_fit.xlsx", fpath), rowNames = FALSE)
          
          cat("✅ 模型完成：", model_id, "\n")
        }
        
          
          }
        }
      }
    }
  }
}

library(openxlsx)
library(dplyr)

# 设置保存你输出文件的文件夹路径
result_dir <- "C:/future/ABCD_data/20250818/euro_prsice_7970_factors_3967subjs/moderated_mediation_results_test"

# 获取所有模型参数结果（不包含 _fit 后缀）
param_files <- list.files(result_dir, pattern = "^[^_]+.*\\.xlsx$", full.names = TRUE)
param_files <- param_files[!grepl("_fit\\.xlsx$", param_files)]  # 只保留主表

# 读取合并
all_params <- do.call(rbind, lapply(param_files, function(f) {
  df <- read.xlsx(f)
  df$File <- basename(f)
  return(df)
}))

fit_files <- list.files(result_dir, pattern = "_fit\\.xlsx$", full.names = TRUE)

all_fit <- do.call(rbind, lapply(fit_files, function(f) {
  df <- read.xlsx(f)
  df$File <- gsub("_fit\\.xlsx", ".xlsx", basename(f))  # 与主结果对齐
  return(df)
}))

all_params <- all_params %>%
  filter(!is.na(label) & label != "") %>%
  mutate(FDR_p = p.adjust(pvalue, method = "fdr"))

all_results <- left_join(all_params, all_fit, by = "File")

# 含全部参数估计和拟合指标
write.xlsx(all_results, file = file.path(result_dir, "summary_all_models_with_fit.xlsx"), rowNames = FALSE)

# 单独保存 indirect 路径 + FDR
indirect_results <- all_results %>%
  filter(grepl("^indirect_", label)) %>%
  arrange(FDR_p)

write.xlsx(indirect_results, file = file.path(result_dir, "summary_indirect_effects_FDR.xlsx"), rowNames = FALSE)

