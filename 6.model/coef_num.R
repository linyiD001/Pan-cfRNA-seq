coef_matrix <- coef(glmnet_model)

best_lambda <- glmnet_model$lambda.min

# 提取最佳 λ 值下的系数
best_coef <- coef(glmnet_model, s = "lambda.min")

# 统计非零变量数量（不包括截距）
nonzero_count <- sum(best_coef != 0) - 1  # 去掉截距项

# 输出结果
cat("最佳 λ =", best_lambda, "\n")
cat("最佳 λ 下的非零变量数量为：", nonzero_count, "\n")

# 提取非零系数的索引（排除截距项）
nonzero_index <- which(best_coef != 0)[-1]  # 去掉截距项（通常是第一个）

# 获取非零变量的名称
nonzero_names <- rownames(best_coef)[nonzero_index]

# 输出结果
cat("最佳 λ 下的非零变量名称：\n")
print(nonzero_names)
