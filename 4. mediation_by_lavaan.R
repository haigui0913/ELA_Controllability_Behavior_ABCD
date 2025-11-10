# ========== 0. 加载必要包 ==========
library(lme4)
library(lavaan)
library(dplyr)
library(openxlsx)

setwd("C:\\future\\ABCD_data\\20250818\\euro_prsice_7970_factors_3967subjs")
# ========== 1. 读取数据并设置分类变量 ==========
load("C:/future/ABCD_data/20250818/euro_prsice_7970_factors_3967subjs/my_all_data_7970.RData")

# ========== 2. 定义变量 ==========
# 变量组
behavior_vars <- colnames(data)[c(150:155)]   # 行为变量
brain_vars    <- c("MF", "FP","Mot", "Vas", "SAL") # 脑网络变量
ela_vars      <- c("F1", "F4", "F5")                        # ELA变量

# 协变量
ctrl_vars <- c("age_z","sex","race","pubertyZ_z","handness",
               "genetic_pc_1_z","genetic_pc_2_z","genetic_pc_3_z",
               "genetic_pc_4_z","genetic_pc_5_z")

data$sex      <- as.factor(data$sex)
data$race     <- as.factor(data$race)
data$handness <- as.factor(data$handness)
data$site     <- as.factor(data$site)

# 替换 NA 为 0（特别是 genetic_pc）
genetic_pcs <- grep("genetic_pc", ctrl_vars, value = TRUE)
data[genetic_pcs] <- lapply(data[genetic_pcs], function(x) { x[is.na(x)] <- 0; x })

# 创建输出文件夹
if (!dir.exists("model_results")) dir.create("model_results")

# 获取已完成模型（用于断点续跑）
done_files <- list.files("model_results", pattern = "\\.xlsx$")
done_keys <- gsub(".xlsx", "", done_files)

# ========== 1. 组合循环 ==========
for (xvar in ela_vars) {
  for (mvar in brain_vars) {
    for (yvar in behavior_vars) {
      
      model_id <- paste(xvar, mvar, yvar, sep = "_")
      if (model_id %in% done_keys) next  # 跳过已完成
      
      cat("Running:", model_id, "\n")
      
      # ------- 1. 残差化 FP_resid 和 cbcl_resid ---------
      # 构建公式
      m_formula <- as.formula(paste0(mvar, " ~ ", paste(ctrl_vars, collapse = " + "), " + (1|site)"))
      y_formula <- as.formula(paste0(yvar, " ~ ", paste(ctrl_vars, collapse = " + "), " + (1|site)"))
      
      # 残差化：中介变量 M
      m_complete <- complete.cases(data[, c(mvar, ctrl_vars, "site")])
      m_fit <- lmer(m_formula, data = data[m_complete, ])
      m_resid <- rep(NA, nrow(data))
      m_resid[m_complete] <- resid(m_fit)
      data[[paste0(mvar, "_resid")]] <- m_resid
      m_resid_name <- paste0(mvar, "_resid")
      
      # 残差化：行为变量 Y
      y_complete <- complete.cases(data[, c(yvar, ctrl_vars, "site")])
      y_fit <- lmer(y_formula, data = data[y_complete, ])
      y_resid <- rep(NA, nrow(data))
      y_resid[y_complete] <- resid(y_fit)
      data[[paste0(yvar, "_resid")]] <- y_resid
      y_resid_name <- paste0(yvar, "_resid")
      
      # ------- 2. lavaan 模型拟合 ---------
      model_text <- paste0('
        ', m_resid_name, ' ~ a * ', xvar, '
        ', y_resid_name, ' ~ b * ', m_resid_name, ' + cp * ', xvar, '
        indirect := a * b
        total := cp + (a * b)
        PM := (a * b) / (cp + a * b)
      ')
      
      fit <- tryCatch({
        sem(model_text, data = data, se = "bootstrap", bootstrap = 5000)
      }, error = function(e) {
        cat("Model failed for", model_id, "\n")
        return(NULL)
      })
      
      if (is.null(fit)) next
      
      # ------- 3. 提取结果 ---------
      pe <- parameterEstimates(fit, boot.ci.type = "perc", standardized = TRUE) %>%
        filter(label %in% c("a", "b", "cp", "indirect", "total", "PM"))
      
      r2_vals <- tryCatch(inspect(fit, "r2"), error = function(e) rep(NA, 2))
      r2_m <- r2_vals[[m_resid_name]] %||% NA
      r2_y <- r2_vals[[y_resid_name]] %||% NA
      
      result_df <- data.frame(
        ELA = xvar,
        Brain = mvar,
        Behavior = yvar,
        Path = pe$label,
        Est = pe$est,
        SE = pe$se,
        pvalue = pe$pvalue,
        CI.lower = pe$ci.lower,
        CI.upper = pe$ci.upper,
        R2_M = r2_m,
        R2_Y = r2_y,
        stringsAsFactors = FALSE
      )
      
      # ------- 4. 保存结果 ---------
      out_path <- file.path("model_results", paste0(model_id, ".xlsx"))
      write.xlsx(result_df, file = out_path, rowNames = FALSE)
      
    }
  }
}


library(openxlsx)
library(dplyr)

# ========== 1. 批量读取所有 .xlsx 结果 ==========
folder <- "C:\\future\\ABCD_data\\20250818\\euro_prsice_7970_factors_3967subjs\\model_results"
files <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)

# 合并所有文件
all_results <- do.call(rbind, lapply(files, function(f) {
  df <- read.xlsx(f)
  df$File <- basename(f)
  return(df)
}))

# ========== 2. 对 indirect 路径进行多重校正 ==========
indirect_df <- all_results %>%
  filter(Path == "indirect") %>%
  mutate(
    FDR_p = p.adjust(pvalue, method = "fdr"))
  

# ========== 3. 合并回主表 ==========
all_results <- left_join(
  all_results,
  indirect_df[, c("File", "FDR_p")],
  by = "File"
)

# ========== 4. 保存结果 ==========
write.xlsx(all_results, "C:\\future\\ABCD_data\\20250818\\euro_prsice_7970_factors_3967subjs\\model_results\\all_mediation_results_with_fdr2.xlsx", rowNames = FALSE)




