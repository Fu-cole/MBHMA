# 加载成员1清洗后的数据
load("data_cleaned.RData")   # 应包含 selmetals

# 1. 定义 Hill 函数
myFuncHill <- function(Ce_mg_l, a, b, c) {
  (a * Ce_mg_l^b) / (c^b + Ce_mg_l^b)
}

myFitHill <- function(df, a, b, c) {
  tryCatch(
    nls(Qe_mg_g ~ myFuncHill(Ce_mg_l, a, b, c),
        start = list(a = a, b = b, c = c),
        data = df),
    error = function(e) NULL   # 拟合失败返回 NULL，不中断流程
  )
}

# 2. 准备数据：按 Biochar 和 Metal 拆分
selmetals$Biochar.Metal <- interaction(selmetals$Biochar, selmetals$Metal)
listmetals <- split(selmetals, selmetals$Biochar.Metal)

# 获取各组合名称
all_names <- names(listmetals)

# 根据重金属类型分别设定起始值（参考原文经验）
# Cd: a=20, b=4, c=0.8
# Pb: a=45, b=9, c=0.1
# Zn: a=8,  b=8, c=0.05

fitHill <- list()

for (nm in all_names) {
  # 判断金属类型
  if (grepl("Cd", nm)) {
    fit <- myFitHill(listmetals[[nm]], a=20, b=4, c=0.8)
  } else if (grepl("Pb", nm)) {
    fit <- myFitHill(listmetals[[nm]], a=45, b=9, c=0.1)
  } else if (grepl("Zn", nm)) {
    fit <- myFitHill(listmetals[[nm]], a=8, b=8, c=0.05)
  } else {
    fit <- NULL
  }
  fitHill[[nm]] <- fit
}

# 移除拟合失败的组合（如果有）
fitHill <- fitHill[!sapply(fitHill, is.null)]
cat("成功拟合的 Hill 模型数量：", length(fitHill), "\n")

# 3. 提取参数（系数、标准误、p值、AIC）
extract_hill_params <- function(fit, model_name) {
  if (is.null(fit)) return(NULL)
  tidy_fit <- tidy(fit)   # 包含 term, estimate, std.error, statistic, p.value
  aic_val <- AIC(fit)
  tidy_fit <- tidy_fit %>%
    mutate(model = model_name,
           AIC = aic_val)
  return(tidy_fit)
}

paramsHill <- bind_rows(lapply(names(fitHill), function(nm) {
  extract_hill_params(fitHill[[nm]], nm)
}))
write.csv(paramsHill, "hill_params.csv", row.names = FALSE)

# 4. 可选：与 Langmuir 模型的 AIC 比较（如果成员2的结果可用）
if (file.exists("langmuir_results.RData")) {
  load("langmuir_results.RData")   # 期望得到 listfitLang
  if (exists("listfitLang")) {
    # 只比较两者都存在的组合
    common_names <- intersect(names(fitHill), names(listfitLang))
    aic_compare <- data.frame(
      Biochar.Metal = common_names,
      AIC_Langmuir = sapply(listfitLang[common_names], AIC),
      AIC_Hill     = sapply(fitHill[common_names], AIC)
    ) %>%
      mutate(Delta_AIC = AIC_Hill - AIC_Langmuir,
             Better = ifelse(Delta_AIC < 0, "Hill", "Langmuir"))
    write.csv(aic_compare, "AIC_comparison.csv", row.names = FALSE)
  }
}

# 5. 生成拟合曲线预测数据（平滑曲线）
newdata <- data.frame(Ce_mg_l = c(seq(0, 1.9, 0.02), seq(2, 30, 1)))
predHill <- lapply(fitHill, function(fit) {
  if (is.null(fit)) return(NULL)
  data.frame(Ce_mg_l = newdata,
             fitQe = predict(fit, newdata))
})
predHill_df <- bind_rows(predHill, .id = "Biochar.Metal")

# 添加元数据（Biochar 和 Metal 列）
metadata <- data.frame(Biochar.Metal = names(fitHill)) %>%
  mutate(
    Biochar = str_split_fixed(Biochar.Metal, "\\.", 2)[,1],
    Metal   = str_split_fixed(Biochar.Metal, "\\.", 2)[,2]
  )
predHill_df <- left_join(predHill_df, metadata, by = "Biochar.Metal")

# 6. 残差诊断（使用 fitted() 和 residuals()，不依赖 nlstools）
resid_list <- lapply(fitHill, function(fit) {
  if (is.null(fit)) return(NULL)
  data.frame(
    fitted = fitted(fit),
    residuals = residuals(fit)
  )
})
resid_df <- bind_rows(resid_list, .id = "Biochar.Metal")

# 绘制残差图（所有组合分面）
p_resid <- ggplot(resid_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Biochar.Metal, scales = "free") +
  theme_bw() +
  labs(title = "Hill 模型残差", x = "拟合值", y = "残差")
ggsave("hill_residuals.png", p_resid, width = 10, height = 8)

# 7. 绘制 Hill 模型拟合曲线总图（全范围）
p_hill_curves <- ggplot() +
  geom_point(data = selmetals, aes(x = Ce_mg_l, y = Qe_mg_g, color = Metal)) +
  geom_line(data = predHill_df, aes(x = Ce_mg_l, y = fitQe, group = Biochar.Metal)) +
  facet_grid(Metal ~ Biochar, scales = "free") +
  theme_bw() +
  labs(title = "Hill 模型拟合曲线", 
       x = bquote('Ce ('*mg.L^-1*')'), 
       y = bquote('Qe ('*mg.g^-1*')'))
ggsave("hill_fit_curves_full.png", p_hill_curves, width = 10, height = 6)

# 8. 低浓度范围（0-2 mg/L）放大图，检查 S 型起始部分
p_hill_low <- ggplot() +
  geom_point(data = selmetals, aes(x = Ce_mg_l, y = Qe_mg_g, color = Metal)) +
  geom_line(data = predHill_df, aes(x = Ce_mg_l, y = fitQe, group = Biochar.Metal)) +
  facet_grid(Metal ~ Biochar, scales = "free") +
  xlim(c(0, 2)) +
  theme_bw() +
  labs(title = "Hill 模型拟合（低浓度区域放大）", 
       x = bquote('Ce ('*mg.L^-1*')'), 
       y = bquote('Qe ('*mg.g^-1*')'))
ggsave("hill_fit_curves_low.png", p_hill_low, width = 10, height = 6)


# 10. 保存所有结果供成员5使用
save(fitHill, paramsHill, predHill_df, resid_df, file = "hill_results.RData")