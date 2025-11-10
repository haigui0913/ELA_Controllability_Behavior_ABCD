# =============================
# 读取ELAs环境及精神健康、认知数据
# =============================# 设置路径
df <- read.csv("C:/future/ABCD_data/ABCD_behavior_envs/ELAs_ave_imputed_9796.csv", stringsAsFactors = FALSE)

# =============================
# 读取全脑controllability数据
# =============================

data_ave <- read.csv("C:/future/ABCD_data/ABCD_DTI/mean_ave_control_baseline.csv", stringsAsFactors = FALSE)
data_modal <- read.csv("C:/future/ABCD_data/ABCD_DTI/mean_modal_control_baseline.csv", stringsAsFactors = FALSE)

data_ave_subnet <- read.csv("C:/future/ABCD_data/ABCD_DTI/subnet_ave_ctrl.csv", stringsAsFactors = FALSE)
data_modal_subnet <- read.csv("C:/future/ABCD_data/ABCD_DTI/subnet_modal_ctrl.csv", stringsAsFactors = FALSE)

# Step 1：处理 data_ave 的 ID，提取核心部分并加上下划线
data_ave_subnet$participant_id_clean <- gsub("sub-", "", data_ave_subnet$participant_id)
data_ave_subnet$participant_id_clean <- gsub("NDARINV", "NDAR_INV", data_ave_subnet$participant_id_clean)


# Step 2：将 final_merged 的 ID 列重命名为统一名称
colnames(data_ave_subnet)[colnames(data_ave_subnet) == "participant_id_clean"] <- "subject_id"
data_ave_subnet$participant_id <- NULL

# Step 3：按统一的 ID 合并两个数据框
merged_df <- merge(data_ave_subnet, df, by = "subject_id")

# Step 1：处理 data_ave 的 ID，提取核心部分并加上下划线
data_modal_subnet$participant_id_clean <- gsub("sub-", "", data_modal_subnet$participant_id)
data_modal_subnet$participant_id_clean <- gsub("NDARINV", "NDAR_INV", data_modal_subnet$participant_id_clean)


# Step 2：将 final_merged 的 ID 列重命名为统一名称
colnames(data_modal_subnet)[colnames(data_modal_subnet) == "participant_id_clean"] <- "subject_id"
data_modal_subnet$participant_id <- NULL

# Step 3：按统一的 ID 合并两个数据框
merged_df2 <- merge(data_modal_subnet,merged_df, by = "subject_id")

# Step 1：处理 data_ave 的 ID，提取核心部分并加上下划线
data_ave$participant_id_clean <- gsub("sub-", "", data_ave$participant_id)
data_ave$participant_id_clean <- gsub("NDARINV", "NDAR_INV", data_ave$participant_id_clean)


# Step 2：将 final_merged 的 ID 列重命名为统一名称
colnames(data_ave)[colnames(data_ave) == "participant_id_clean"] <- "subject_id"
data_ave$participant_id <- NULL

# Step 3：按统一的 ID 合并两个数据框
merged_df3 <- merge(data_ave, merged_df2, by = "subject_id")

# Step 1：处理 data_modal 的 ID，提取核心部分并加上下划线
data_modal$participant_id_clean <- gsub("sub-", "", data_modal$participant_id)
data_modal$participant_id_clean <- gsub("NDARINV", "NDAR_INV", data_modal$participant_id_clean)

# Step 2：将 final_merged 的 ID 列重命名为统一名称
colnames(data_modal)[colnames(data_modal) == "participant_id_clean"] <- "subject_id"
data_modal$participant_id <- NULL

# Step 3：按统一的 ID 合并两个数据框
merged_df4 <- merge(data_modal, merged_df3, by = "subject_id")

write.csv(merged_df4, "C:/future/ABCD_data/controllability_ELAs_behavior_7970.csv", row.names = FALSE)

##merge_euro
## 1. 读取数据
data <- read.csv("C:/future/ABCD_data/45_test_factor_scores_data_all_7970.csv", stringsAsFactors = FALSE)
pc_file<-read.csv("C:/future/ABCD_data/ABCD_PRS/gen_y_pihat.csv", stringsAsFactors = FALSE)
colnames(pc_file)[colnames(pc_file) == "src_subject_id"] <- "subject_id"
pc_file2 <- pc_file[, c(1,21:30)] 
final_merged <- merge(data, pc_file2, by = "subject_id")
save(final_merged, file = "C:/future/ABCD_data/20250825/my_all_data_7970.RData")

load("C:/future/ABCD_data/20250825/my_all_data_7970.RData")

prs_file <- read.csv("C:/future/ABCD_data/ABCD_PRS/test/mean_model_control/model_B2_euro_prs.csv", stringsAsFactors = FALSE)

# 改列名、删除列
colnames(prs_file)[colnames(prs_file) == "FID"] <- "subject_id"
prs_df <- prs_file[, -c(4,5)]  # 删除第2列和第4列

data_eur <- merge(final_merged, prs_df, by = "subject_id")

save(data_eur, file = "C:/future/ABCD_data/20250825/my_euro_data_3967.RData")


