# 安装必要包
library(lavaan)
library(semTools)
library(psych)
library(GPArotation)
library(polycor)
library(openxlsx)

#读取数据
data <- read.csv("C:/future/ABCD_data/controllability_ELAs_behavior_7970.csv", stringsAsFactors = FALSE)

#1. 因子分析
# 提前定义变量类型
ela_data <- data[, 31:97]

# 自动识别变量
n_unique <- sapply(ela_data, function(x) length(unique(na.omit(x))))
binary_vars <- names(ela_data)[n_unique == 2]
ordinal_vars <- names(ela_data)[n_unique >= 3 & n_unique <= 5]
cont_vars <- names(ela_data)[n_unique > 5]

# 将变量统一为 ordered factor（lavaan对有序数据更友好）
ela_for_lavaan <- ela_data
ela_for_lavaan[binary_vars] <- lapply(ela_for_lavaan[binary_vars], function(x) as.ordered(as.factor(x)))
ela_for_lavaan[ordinal_vars] <- lapply(ela_for_lavaan[ordinal_vars], function(x) as.ordered(as.numeric(x)))
ela_for_lavaan[cont_vars] <- lapply(ela_for_lavaan[cont_vars], function(x) as.numeric(x))

# 检查类型
str(ela_for_lavaan[1:5])


# 将 ela_for_lavaan 中所有变量统一转为 numeric（用于检验）
ela_for_test <- as.data.frame(lapply(ela_for_lavaan, function(x) as.numeric(as.character(x))))

# 计算 spearman 相关矩阵（适用于有序数据）
ela_cor_test <- cor(ela_for_test, method = "spearman", use = "pairwise.complete.obs")

# ------------------ KMO 检验 ------------------
kmo_result <- KMO(ela_cor_test)
print(kmo_result$MSA)  # 输出整体KMO值（>0.6可接受）
print(kmo_result$MSAi) # 各变量的KMO值（>0.5为合格）

# ------------------ Bartlett's 球形检验 ------------------
bartlett_result <- cortest.bartlett(ela_cor_test, n = nrow(ela_for_test))
print(bartlett_result)



# 创建一个容错版本的 polychoric 相关计算函数
safe_hetcor <- function(data) {
  good_vars <- names(data)
  repeat {
    result <- tryCatch({
      hetcor(data[good_vars])$correlations
    }, error = function(e) {
      message("Error in hetcor: ", e$message)
      return(NULL)
    })
    
    if (!is.null(result) && !any(is.na(result))) break
    
    # 如果有 NA，找出有问题的变量，剔除后重试
    na_vars <- which(colSums(is.na(result)) > 0)
    if (length(na_vars) == 0) stop("仍存在未能解决的相关性问题")
    
    var_to_remove <- names(na_vars)[1]
    message("Removing problematic variable: ", var_to_remove)
    good_vars <- setdiff(good_vars, var_to_remove)
  }
  return(result)
}

# 将 ela_for_lavaan 转换为 numeric 类型用于 polychoric 计算
ela_numeric_check <- as.data.frame(lapply(ela_for_lavaan, function(x) as.numeric(as.character(x))))

# 然后运行安全的 polychoric 相关函数
polycor_mat <- safe_hetcor(ela_numeric_check)

# ------------------ KMO 检验 ------------------
kmo_result <- KMO(polycor_mat)
print(kmo_result$MSA)  # 输出整体KMO值（>0.6可接受）
print(kmo_result$MSAi) # 各变量的KMO值（>0.5为合格）

# ------------------ Bartlett's 球形检验 ------------------
bartlett_result <- cortest.bartlett(polycor_mat, n = 7970)

print(bartlett_result)


# 然后继续平行分析
fa.parallel(polycor_mat, fa = "fa", fm = "pa", n.obs = nrow(ela_numeric_check))


# ------------------ 4. VSS（Very Simple Structure） + MAP ------------------
vss_result <- VSS(polycor_mat, n = 10, fm = "pa", rotate = "oblimin", plot = TRUE)

# 查看最优因子数推荐结果
cat("建议因子数（基于 VSS complexity 1）:", which.max(vss_result$cfit.1), "\n")
cat("建议因子数（基于 VSS complexity 2）:", which.max(vss_result$cfit.2), "\n")
cat("建议因子数（基于 MAP）:", which.min(vss_result$map), "\n")

# 也可以用这个函数查看所有指标
print(vss_result)

# 保存结果为表格
vss_df <- data.frame(
  nfactors = 1:length(vss_result$cfit.1),
  VSS1 = round(vss_result$cfit.1, 3),
  VSS2 = round(vss_result$cfit.2, 3),
  MAP  = round(vss_result$map, 3)
)
write.csv(vss_df, "C:/future/ABCD_data/VSS_MAP_results.csv", row.names = FALSE)


# ------------------ 5. EFA + CFA ------------------
# 执行 EFA，指定 5 因子、WLSMV、迭代次数最多10000
efa_model <- efa('efa', data = ela_for_lavaan,
                 nfactors = 5,                      # 可替换为 6 或 7
                 estimator = "WLSMV",
                 rotation = "oblimin",
                 ordered = c(binary_vars, ordinal_vars),
                 control = list(iter.max = 10000))


# 查看载荷矩阵
print(efa_model$loadings, cutoff = 0.45)

best_k <- 5

# 提取并保存负载矩阵（加载 factor loading）
loadings_mat <- as.matrix(efa_model$loadings)
loadings_df <- as.data.frame(loadings_mat)
loadings_df$Variable <- rownames(loadings_mat)
loadings_df <- loadings_df[, c("Variable", paste0("f", 1:best_k))]

write.xlsx(loadings_df, file = paste0("C:/future/ABCD_data/EFA_", best_k, "f_loadings.xlsx"), rowNames = FALSE)

# 【自动生成 CFA model_text】
model_list <- list()
n_factor <- ncol(loadings_mat)
for (i in 1:n_factor) {
  vars <- rownames(loadings_mat)[abs(loadings_mat[, i]) > 0.45]
  if (length(vars) > 0) {
    factor_name <- paste0("F", i)
    eq <- paste0(factor_name, " =~ ", paste(vars, collapse = " + "))
    model_list[[i]] <- eq
  }
}
model_text <- paste(model_list, collapse = "\n")
cat("【自动生成 CFA model_text】\n", model_text, "\n")

# 拟合模型
fit_cfa <- cfa(model_text, data = ela_for_lavaan,
               estimator = "WLSMV",
               ordered = c(binary_vars, ordinal_vars))

# 提取得分
factor_scores <- lavPredict(fit_cfa)
head(factor_scores)

factor_scores_z <- scale(factor_scores) 

# 合并得分
data_with_scores <- cbind(data,factor_scores)

data_all <- data_with_scores[c(1:30,98:132)]

write.csv(data_all,"C:/future/ABCD_data/45_test_factor_scores_data_all_7970.csv",row.names = FALSE)

set.seed(123)  # 保证可重复性

# 1. 创建 10 折分组（按行号分组）
folds <- createFolds(1:nrow(ela_for_lavaan), k = 10, list = TRUE, returnTrain = FALSE)

# 2. 存储每一折的模型拟合结果（CFI、RMSEA等）
cv_results <- data.frame(
  Fold = integer(),
  CFI = numeric(),
  TLI = numeric(),
  RMSEA = numeric(),
  SRMR = numeric(),
  stringsAsFactors = FALSE
)

# 3. 开始 10 折交叉验证
for (i in seq_along(folds)) {
  cat(paste0("正在处理第 ", i, " 折...\n"))
  
  test_idx <- folds[[i]]
  test_data <- ela_for_lavaan[test_idx, ]
  
  # 运行 CFA 模型（固定结构）
  fit_fold <- tryCatch({
    cfa(model_text,
        data = test_data,
        estimator = "WLSMV",
        ordered = c(binary_vars, ordinal_vars))
  }, error = function(e) {
    message("Fold ", i, " CFA 模型报错: ", e$message)
    return(NULL)
  })
  
  # 如果模型正常运行，提取拟合指标
  if (!is.null(fit_fold)) {
    fit_measures <- fitMeasures(fit_fold, c("cfi", "tli", "rmsea", "srmr"))
    cv_results <- rbind(cv_results, data.frame(
      Fold = i,
      CFI = fit_measures["cfi"],
      TLI = fit_measures["tli"],
      RMSEA = fit_measures["rmsea"],
      SRMR = fit_measures["srmr"]
    ))
  }
}

# 4. 查看结果
print(cv_results)

# 5. 平均拟合指标
cat("=== 10折交叉验证平均模型拟合指标 ===\n")
print(colMeans(cv_results[, -1], na.rm = TRUE))

library(lavaan)
library(dplyr)

# Step 1: 全样本模型的因子得分
gold_scores <- lavPredict(fit_cfa)
gold_scores_df <- as.data.frame(gold_scores)

# 添加 row ID
gold_scores_df$ID <- 1:nrow(gold_scores_df)

# 初始化列表存储每一折的相关系数
cor_results_list <- list()

# Step 2: 针对每一折提取得分并计算与gold的相关
for (i in seq_along(folds)) {
  cat("处理第", i, "折...\n")
  test_idx <- folds[[i]]
  test_data <- ela_for_lavaan[test_idx, ]
  
  # 折内拟合模型
  fit_fold <- tryCatch({
    cfa(model_text,
        data = test_data,
        estimator = "WLSMV",
        ordered = c(binary_vars, ordinal_vars))
  }, error = function(e) {
    message("第", i, "折模型出错: ", e$message)
    return(NULL)
  })
  
  # 提取得分 & 计算相关
  if (!is.null(fit_fold)) {
    fold_scores <- lavPredict(fit_fold)
    fold_scores_df <- as.data.frame(fold_scores)
    
    # 匹配折内的 gold 标准得分（同样的 index）
    gold_sub <- gold_scores_df[test_idx, 1:ncol(fold_scores_df)]
    
    # 计算 Pearson 相关
    cor_mat <- cor(fold_scores_df, gold_sub, method = "pearson", use = "pairwise.complete.obs")
    
    # 提取对角线（每个因子对应自己的相关系数）
    diag_cor <- diag(cor_mat)
    
    cor_results_list[[i]] <- diag_cor
  }
}

# 合并结果
cor_df <- do.call(rbind, cor_results_list)
rownames(cor_df) <- paste0("Fold_", 1:length(cor_results_list))

# 输出
print(round(cor_df, 3))

# 平均一致性
cat("=== 各因子平均 Pearson 相关系数 ===\n")
print(round(colMeans(cor_df, na.rm = TRUE), 3))

