# 加载成员1清洗后的数据
load("data_cleaned.RData")   # 包含 selmetals

# 1. 定义 Langmuir 函数
myFuncLang <- function(qm, Ce_mg_l, kl) {
  (qm * kl * Ce_mg_l) / (1 + (kl * Ce_mg_l))
}

myFitLang <- function(df, kl, qm) {
  tryCatch(
    nls(Qe_mg_g ~ myFuncLang(qm, Ce_mg_l, kl),
        start = list(kl = kl, qm = qm),
        data = df),
    error = function(e) NULL   # 拟合失败返回 NULL，避免中断
  )
}

# 2. 按 Biochar 和 Metal 拆分数据
selmetals$Biochar.Metal <- interaction(selmetals$Biochar, selmetals$Metal)
listmetals <- split(selmetals, selmetals$Biochar.Metal)

# 3. 批量拟合 Langmuir 模型（含错误处理）
listfitLang <- lapply(listmetals, function(x) myFitLang(x, kl=1, qm=70))

# 移除拟合失败的组合（如果有）
listfitLang <- listfitLang[!sapply(listfitLang, is.null)]
cat("成功拟合的模型数量：", length(listfitLang), "\n")

# 4. 提取参数与 AIC
extract_params_aic <- function(fit, model_name) {
  if (is.null(fit)) return(NULL)
  tidy_fit <- tidy(fit)
  aic_val <- AIC(fit)
  tidy_fit <- tidy_fit %>%
    mutate(model = model_name,
           AIC = aic_val)
  return(tidy_fit)
}

paramsLang <- bind_rows(lapply(names(listfitLang), function(nm) {
  extract_params_aic(listfitLang[[nm]], nm)
}))
write.csv(paramsLang, "langmuir_params.csv", row.names=FALSE)

# 5. 残差诊断（使用 fitted 和 residuals，避免 nlsResiduals 的坑）
resid_list <- lapply(listfitLang, function(fit) {
  data.frame(
    fitted = fitted(fit),
    residuals = residuals(fit)
  )
})
resid_df <- bind_rows(resid_list, .id = "Biochar.Metal")

# 绘制残差图
p_resid <- ggplot(resid_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Biochar.Metal, scales = "free") +
  theme_bw() +
  labs(title = "Langmuir 模型残差", x = "拟合值", y = "残差")
ggsave("langmuir_residuals.png", p_resid, width = 10, height = 8)

# 6. 生成 Langmuir 拟合曲线数据（用于绘制平滑曲线）
newdata <- data.frame(Ce_mg_l = c(seq(0, 1.9, 0.02), seq(2, 30, 1)))
predLang <- lapply(listfitLang, function(fit) {
  data.frame(Ce_mg_l = newdata,
             fitQe = predict(fit, newdata))
})
predLang_df <- bind_rows(predLang, .id = "Biochar.Metal")

# 分离 Biochar 和 Metal 信息
metadata <- data.frame(Biochar.Metal = names(listfitLang)) %>%
  mutate(
    Biochar = str_split_fixed(Biochar.Metal, "\\.", 2)[,1],
    Metal   = str_split_fixed(Biochar.Metal, "\\.", 2)[,2]
  )
predLang_df <- left_join(predLang_df, metadata, by = "Biochar.Metal")

# 7.绘制 Langmuir 拟合曲线（快速检查）
p_lang_curves <- ggplot() +
  geom_point(data = selmetals, aes(x = Ce_mg_l, y = Qe_mg_g, color = Metal)) +
  geom_line(data = predLang_df, aes(x = Ce_mg_l, y = fitQe, group = Biochar.Metal)) +
  facet_grid(Metal ~ Biochar, scales = "free") +
  theme_bw() +
  labs(title = "Langmuir 模型拟合", x = "Ce (mg/L)", y = "Qe (mg/g)")
ggsave("langmuir_fit_curves.png", p_lang_curves, width = 10, height = 6)

# 8. 保存中间结果供其他成员使用
save(listfitLang, paramsLang, predLang_df, file = "langmuir_results.RData")