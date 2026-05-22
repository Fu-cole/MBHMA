# 加载成员3的完整结果
if (file.exists("hill_results.RData")) {
  load("hill_results.RData")  # 应包含 paramsHill, fitHill 等
  if (!exists("paramsHill")) stop("hill_results.RData 中缺少 paramsHill")
} else if (file.exists("hill_params.csv")) {
  paramsHill <- read.csv("hill_params.csv", stringsAsFactors = FALSE)
} else {
  stop("找不到 hill_results.RData 或 hill_params.csv，请先运行成员3的代码。")
}

# 检查必要的列：term, estimate, std.error
if (!all(c("term", "estimate", "std.error") %in% colnames(paramsHill))) {
  stop("paramsHill 缺少 term, estimate 或 std.error 列")
}

# 如果不存在 Biochar 和 Metal 列，则从 model 或 Biochar.Metal 列拆分
if (!all(c("Biochar", "Metal") %in% colnames(paramsHill))) {
  # 优先使用 model 列（成员3代码中创建的）
  if ("model" %in% colnames(paramsHill)) {
    paramsHill <- paramsHill %>%
      separate(model, into = c("Biochar", "Metal"), sep = "\\.", remove = FALSE)
  } else if ("Biochar.Metal" %in% colnames(paramsHill)) {
    paramsHill <- paramsHill %>%
      separate(Biochar.Metal, into = c("Biochar", "Metal"), sep = "\\.", remove = FALSE)
  } else {
    stop("paramsHill 中既没有 Biochar/Metal 列，也没有 model/Biochar.Metal 列，无法继续。")
  }
}

# 筛选参数 a (最大吸附量)
params_a <- paramsHill[paramsHill$term == "a", ]
if (nrow(params_a) == 0) {
  stop("未找到 term 为 'a' 的行，请检查 paramsHill 中的 term 列。")
}

# 确保 estimate 和 std.error 是数值型
params_a$estimate <- as.numeric(params_a$estimate)
params_a$std.error <- as.numeric(params_a$std.error)

# 去除缺失值
params_a <- params_a[!is.na(params_a$estimate) & !is.na(params_a$std.error), ]

# 转换为因子
params_a$Biochar <- as.factor(params_a$Biochar)
params_a$Metal <- as.factor(params_a$Metal)
params_a <- droplevels(params_a)

# 查看数据结构
print(head(params_a))

# ========= 1. 绘制 Qmax 比较图 =========
p <- ggplot(params_a, aes(x = Biochar, y = estimate, color = Metal)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(0.3)) +
  facet_wrap(~ Metal, scales = "free_y") +
  theme_bw() +
  labs(y = bquote('Qmax ('*mg.g^-1*')'), 
       title = "不同生物炭对三种重金属的最大吸附量 (Hill 模型参数 a)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Qmax_comparison.png", p, width = 8, height = 6)
print(p)

# ========= 2. 两两比较（基于标准误的 z 检验） =========
pairwise_z_test <- function(df) {
  biochars <- unique(df$Biochar)
  if (length(biochars) < 2) return(data.frame())
  results <- data.frame()
  for (i in 1:(length(biochars)-1)) {
    for (j in (i+1):length(biochars)) {
      sub_i <- df[df$Biochar == biochars[i], ]
      sub_j <- df[df$Biochar == biochars[j], ]
      if (nrow(sub_i) == 0 || nrow(sub_j) == 0) next
      est1 <- sub_i$estimate[1]
      se1 <- sub_i$std.error[1]
      est2 <- sub_j$estimate[1]
      se2 <- sub_j$std.error[1]
      diff <- est1 - est2
      se_diff <- sqrt(se1^2 + se2^2)
      z <- diff / se_diff
      p_val <- 2 * (1 - pnorm(abs(z)))
      results <- rbind(results, data.frame(
        Biochar1 = as.character(biochars[i]),
        Biochar2 = as.character(biochars[j]),
        diff = diff,
        z = z,
        p_value = p_val
      ))
    }
  }
  return(results)
}

# 对每种金属分别进行两两比较
metal_list <- unique(params_a$Metal)
all_results <- list()
for (m in metal_list) {
  subdf <- params_a[params_a$Metal == m, ]
  if (nrow(subdf) < 2) next
  res <- pairwise_z_test(subdf)
  if (nrow(res) > 0) {
    res$Metal <- m
    all_results[[m]] <- res
  }
}
final_results <- bind_rows(all_results)
write.csv(final_results, "pairwise_ztest_results.csv", row.names = FALSE)

# 输出显著差异的对
sig_results <- final_results[final_results$p_value < 0.05, ]
cat("\n=== 显著差异的生物炭对 (p < 0.05) ===\n")
if (nrow(sig_results) > 0) {
  print(sig_results)
} else {
  cat("无显著差异。\n")
}

# 保存成员4的输出
save(params_a, final_results, sig_results, file = "statistical_comparison.RData")