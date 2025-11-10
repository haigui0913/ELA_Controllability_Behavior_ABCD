library(lmerTest)
library(broom.mixed)

load ("C:/future/ABCD_data/20250825/my_all_data_7970.RData")
#load ("C:/future/ABCD_data/20250825/my_euro_data_3967.RData")

# 提取非欧洲人
non_eur <- subset(final_merged, race != 1)
eur <- subset(data, race == 1)

data_eur <- final_merged

vars_to_scale <- c(2:24, 27, 30, 31:75)
data_eur[vars_to_scale] <- lapply(data_eur[vars_to_scale], function(x) as.numeric(scale(x)))

# 1. env/prs to brain
# ==== 设置变量 ====
brain_vars <- colnames(data_eur)[4:13] # modal 
#brain_vars <- colnames(data_eur)[14:23] # ave 
env_vars <- colnames(data_eur)[61:65] 
#prs_vars <- colnames(data_eur)[76:78]  # colnames(data_eur)[153:157]

ctrl_vars <- c("age", "sex", "race", "pubertyZ", "handness","genetic_pc_1", "genetic_pc_2", "genetic_pc_3", "genetic_pc_4", "genetic_pc_5") #

categorical_vars <- c("sex", "site", "race", "handness") #race
data_eur[categorical_vars] <- lapply(data_eur[categorical_vars], as.factor)


# ==== 建模函数 ====
run_lmm <- function(y_var, x_var, data, interaction = FALSE, prs = NULL) {
  if (interaction && !is.null(prs)) {
    formula_str <- paste0(y_var, " ~ ", x_var, "*", prs, " + ", 
                          paste(ctrl_vars, collapse = " + "), " + (1|site)")
  } else {
    formula_str <- paste0(y_var, " ~ ", x_var, " + ", 
                          paste(ctrl_vars, collapse = " + "), " + (1|site)")
  }
  model <- lmer(as.formula(formula_str), data = data, REML = FALSE)
  result <- tidy(model)
  result$formula <- formula_str
  return(result)
}

# ==== 主效应分析：环境因子 ====
results_env_main <- list()
for (brain in brain_vars) {
  for (env in env_vars) {
    res <- run_lmm(brain, env, data_eur, interaction = FALSE)
    res$model <- paste0(brain, "~", env)
    results_env_main[[paste0(brain, "_", env)]] <- res
  }
}

# ==== 主效应分析：PRS ====
results_prs_main <- list()
for (brain in brain_vars) {
  for (prs in prs_vars) {
    res <- run_lmm(brain, prs, data_eur, interaction = FALSE)
    res$model <- paste0(brain, "~", prs)
    results_prs_main[[paste0(brain, "_", prs)]] <- res
  }
}

# ==== G×E交互模型 ====
results_interact <- list()
for (brain in brain_vars) {
  for (env in env_vars) {
    for (prs in prs_vars) {
      res <- run_lmm(brain, env, data_eur, interaction = TRUE, prs = prs)
      res$model <- paste0(brain, "~", env, "*", prs)
      results_interact[[paste0(brain, "_", env, "_x_", prs)]] <- res
    }
  }
}

# ==== 合并结果 ====
main_env_df <- bind_rows(results_env_main, .id = "id")
main_prs_df <- bind_rows(results_prs_main, .id = "id")
interact_df <- bind_rows(results_interact, .id = "id")

# ==== FDR 校正 ====
main_env_df <- main_env_df %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_bonf = p.adjust(p.value, method = "bonferroni")
  ) %>%
  ungroup()


main_prs_df <- main_prs_df %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_bonf = p.adjust(p.value, method = "bonferroni")
  ) %>%
  ungroup()

interact_df <- interact_df %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_bonf = p.adjust(p.value, method = "bonferroni")
  ) %>%
  ungroup()

# ==== 保存结果 ====
write.csv(main_env_df, "C:\\future\\ABCD_data\\20250825\\main_effects_env_results_non_euro_3742_modal.csv", row.names = FALSE)
write.csv(main_prs_df, "C:\\future\\ABCD_data\\20250825\\main_effects_prs_results_euro_3967_ave.csv", row.names = FALSE)
write.csv(interact_df, "C:\\future\\ABCD_data\\20250825\\interaction_effects_results_euro_3967_ave.csv", row.names = FALSE)

##热图
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)


# 添加显著性标记（* 表示 FDR < 0.05）
df_env_brain <- main_env_df %>%
  filter(term %in% paste0("F", 1:5)) %>%  # 仅保留 F1–F5
  mutate(
    brain_net = str_extract(id, "^[^_]+"),
    Signif = ifelse(p.value < 0.05, "*", ""),
    term = factor(term, levels = paste0("F", 1:5)),
    brain_net = factor(brain_net, levels = unique(brain_net))
  )

# 绘图
p_env_brain <- ggplot(df_env_brain, aes(x = term, y = brain_net, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, name = "Beta"
  ) +
  labs(
    title = "Heatmap: ELA Factors (F1–F5) → Brain Network Controllability",
    x = "ELA Factor",
    y = "Network Controllability"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

# 保存图像
ggsave("C:\\future\\ABCD_data\\20250825\\Env_to_Brain_Heatmap_fdr_non_euro_3742_modal.tiff",
       plot = p_env_brain, width = 150, height = 90, units = "mm", dpi = 600)




## 2. brain to behavior
behavior_vars <- colnames(data_eur)[c(39:41,58:60)] # [c(121:129,140:142)]
modal_vars <- colnames(data_eur)[14:23]           # modal指标（前10列）##[14:23]

ctrl_vars <- c("age", "sex", "race", "pubertyZ", "handness","genetic_pc_1", "genetic_pc_2", "genetic_pc_3", "genetic_pc_4", "genetic_pc_5") #

categorical_vars <- c("sex", "site", "race", "handness") #race

data_eur[categorical_vars] <- lapply(data_eur[categorical_vars], as.factor)

# ==== 建模函数 ====
run_lmm_behavior <- function(y_var, x_var, data) {
  formula_str <- paste0(y_var, " ~ ", x_var, " + ", 
                        paste(ctrl_vars, collapse = " + "), " + (1|site)")
  model <- lmer(as.formula(formula_str), data = data, REML = FALSE)
  result <- tidy(model)
  result$behavior <- y_var
  result$brain_network <- x_var
  result$formula <- formula_str
  return(result)
}

# ==== 批量建模 ====
results_behavior <- list()
for (y in behavior_vars) {
  for (x in modal_vars) {
    key <- paste0(y, "_~_", x)
    results_behavior[[key]] <- run_lmm_behavior(y, x, data_eur)
  }
}

# ==== 合并结果 ====
df_behavior <- bind_rows(results_behavior, .id = "id")

# ==== 添加FDR和Bonf校正 ====
df_behavior <- df_behavior %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "fdr"),
    p_bonf = p.adjust(p.value, method = "bonferroni")
  ) %>%
  ungroup()

# ==== 导出 ====
write.csv(df_behavior, "C:\\future\\ABCD_data\\20250825\\Control_to_6_behaviors_lmm_results_non_euro_3742_ave.csv", row.names = FALSE)

##热图
library(ggplot2)
library(dplyr)
library(tidyr)

# 添加显著性标记（* 表示 FDR < 0.05）
heatmap_data <- df_behavior %>%
  filter(term == brain_network) %>%  # 仅保留脑网络本体主效应行
  mutate(
    Signif = ifelse(p_fdr < 0.05, "*", ""),
    behavior = factor(behavior, levels = unique(behavior)),
    brain_network = factor(brain_network, levels = unique(brain_network))
  )

# 绘图
p <- ggplot(heatmap_data, aes(x = brain_network, y = behavior, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Signif), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, name = "Beta"
  ) +
  labs(
    title = "Heatmap: Brain Controllability Predicting Behavior",
    x = "Network Controllability",
    y = "Behavioral Outcome"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

# 保存图像
ggsave("C:\\future\\ABCD_data\\20250825\\Control_6_Behavior_heatmap_fdr_non_euro_3742_ave.tiff",
       plot = p, width = 180, height = 100, units = "mm", dpi = 600)
