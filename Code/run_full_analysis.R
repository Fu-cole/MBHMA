# ============================================================
# 一键复现脚本：植物基生物炭对 Cd, Pb, Zn 的吸附等温线分析
# 整合成员1-5 的全部功能，输出最终报告
# 运行要求：R (>= 4.0) 以及以下包：tidyverse, minpack.lm, rmarkdown
# 数据文件：allmetals2_isotherm.csv (请放置于当前工作目录)
# ============================================================

# 设置工作目录（可修改为实际路径，或保持当前）
# setwd("your_path")

# 安装缺失的包
required_pkgs <- c("tidyverse", "minpack.lm", "rmarkdown", "knitr", "broom")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(minpack.lm)
library(broom)
library(knitr)

# 创建输出目录
if(!dir.exists("output")) dir.create("output")

# ========================= 1. 数据清洗（成员1） =========================
cat("步骤1：加载与清洗数据...\n")
if(!file.exists("allmetals2_isotherm.csv")) {
  stop("未找到 allmetals2_isotherm.csv，请将数据文件放在当前目录下。")
}
allmetals <- read.csv("allmetals2_isotherm.csv", stringsAsFactors = FALSE)

# 剔除无效生物炭和高浓度异常点
selmetals <- allmetals %>%
  filter(!Biochar %in% c("SWP550", "SWP550AW"), Ce_mg_l < 50) %>%
  mutate(Biochar = as.factor(Biochar),
         Metal = as.factor(Metal))

# 描述性统计
desc_stats <- selmetals %>%
  group_by(Biochar, Metal) %>%
  summarise(n = n(),
            Ce_min = min(Ce_mg_l), Ce_max = max(Ce_mg_l),
            Qe_mean = mean(Qe_mg_g), Qe_sd = sd(Qe_mg_g),
            .groups = "drop")
write.csv(desc_stats, "output/descriptive_stats.csv", row.names = FALSE)

# 保存清洗后数据供后续使用
saveRDS(selmetals, "output/selmetals.rds")

# 清洗后散点图
p1 <- ggplot(selmetals, aes(x = Ce_mg_l, y = Qe_mg_g, color = Metal)) +
  geom_point() + facet_wrap(~Biochar, scales = "free") + theme_bw() +
  labs(title = "清洗后的吸附等温线数据")
ggsave("output/cleaned_scatter.png", p1, width = 10, height = 6)

# ========================= 2. Langmuir 模型拟合（成员2） =========================
cat("步骤2：拟合 Langmuir 模型...\n")

Langmuir_func <- function(Ce, Qmax, KL) { (Qmax * KL * Ce) / (1 + KL * Ce) }

safe_fit_lang <- function(df) {
  tryCatch(
    nlsLM(Qe_mg_g ~ Langmuir_func(Ce_mg_l, Qmax, KL),
          data = df,
          start = list(Qmax = 70, KL = 1),
          control = nls.lm.control(maxiter = 500)),
    error = function(e) NULL
  )
}

# 按生物炭和金属分组
df_list <- selmetals %>% group_by(Biochar, Metal) %>% group_split()
names_df <- selmetals %>% group_by(Biochar, Metal) %>% group_keys() %>%
  mutate(name = paste(Biochar, Metal, sep = "."))

lang_fits <- list()
for(i in seq_along(df_list)) {
  fit <- safe_fit_lang(df_list[[i]])
  lang_fits[[names_df$name[i]]] <- fit
}
lang_fits <- lang_fits[!sapply(lang_fits, is.null)]

# 提取参数
lang_params <- bind_rows(lapply(names(lang_fits), function(nm) {
  fit <- lang_fits[[nm]]
  tid <- tidy(fit) %>% mutate(model = nm, AIC = AIC(fit))
  tid$Biochar <- str_split(nm, "\\.")[[1]][1]
  tid$Metal <- str_split(nm, "\\.")[[1]][2]
  tid
}))
write.csv(lang_params, "output/langmuir_params.csv", row.names = FALSE)

# 残差图（简单绘制）
resid_lang <- bind_rows(lapply(names(lang_fits), function(nm) {
  fit <- lang_fits[[nm]]
  data.frame(model = nm,
             fitted = fitted(fit),
             resid = residuals(fit))
}))
p_resid_lang <- ggplot(resid_lang, aes(fitted, resid)) + geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~model, scales = "free") + theme_bw()
ggsave("output/langmuir_residuals.png", p_resid_lang, width = 12, height = 8)

# 预测曲线
newdata <- data.frame(Ce_mg_l = c(seq(0,1.9,0.02), seq(2,30,1)))
pred_lang <- bind_rows(lapply(names(lang_fits), function(nm) {
  fit <- lang_fits[[nm]]
  pred <- predict(fit, newdata)
  data.frame(model = nm, Ce_mg_l = newdata$Ce_mg_l, Qe_pred = pred)
})) %>% separate(model, into = c("Biochar","Metal"), sep = "\\.")

p_lang_curve <- ggplot() +
  geom_point(data = selmetals, aes(Ce_mg_l, Qe_mg_g, color = Metal)) +
  geom_line(data = pred_lang, aes(Ce_mg_l, Qe_pred, group = Metal), color = "black") +
  facet_grid(Metal ~ Biochar, scales = "free") + theme_bw()
ggsave("output/langmuir_fit_curves.png", p_lang_curve, width = 10, height = 6)

# ========================= 3. Hill 模型拟合（成员3） =========================
cat("步骤3：拟合 Hill 模型...\n")

Hill_func <- function(Ce, a, b, c) { (a * Ce^b) / (c^b + Ce^b) }

safe_fit_hill <- function(df, a0, b0, c0) {
  tryCatch(
    nlsLM(Qe_mg_g ~ Hill_func(Ce_mg_l, a, b, c),
          data = df,
          start = list(a = a0, b = b0, c = c0),
          control = nls.lm.control(maxiter = 500)),
    error = function(e) NULL
  )
}

# 根据金属分别设定初值
hill_fits <- list()
for(i in seq_along(df_list)) {
  df <- df_list[[i]]
  metal <- unique(df$Metal)
  if(metal == "Cd") { a0 <- 20; b0 <- 4; c0 <- 0.8 }
  else if(metal == "Pb") { a0 <- 45; b0 <- 9; c0 <- 0.1 }
  else { a0 <- 8; b0 <- 8; c0 <- 0.05 }
  fit <- safe_fit_hill(df, a0, b0, c0)
  hill_fits[[names_df$name[i]]] <- fit
}
hill_fits <- hill_fits[!sapply(hill_fits, is.null)]

# 提取参数
hill_params <- bind_rows(lapply(names(hill_fits), function(nm) {
  fit <- hill_fits[[nm]]
  tid <- tidy(fit) %>% mutate(model = nm, AIC = AIC(fit))
  tid$Biochar <- str_split(nm, "\\.")[[1]][1]
  tid$Metal <- str_split(nm, "\\.")[[1]][2]
  tid
}))
write.csv(hill_params, "output/hill_params.csv", row.names = FALSE)

# 残差图
resid_hill <- bind_rows(lapply(names(hill_fits), function(nm) {
  fit <- hill_fits[[nm]]
  data.frame(model = nm,
             fitted = fitted(fit),
             resid = residuals(fit))
}))
p_resid_hill <- ggplot(resid_hill, aes(fitted, resid)) + geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~model, scales = "free") + theme_bw()
ggsave("output/hill_residuals.png", p_resid_hill, width = 12, height = 8)

# Hill 预测曲线
pred_hill <- bind_rows(lapply(names(hill_fits), function(nm) {
  fit <- hill_fits[[nm]]
  pred <- predict(fit, newdata)
  data.frame(model = nm, Ce_mg_l = newdata$Ce_mg_l, Qe_pred = pred)
})) %>% separate(model, into = c("Biochar","Metal"), sep = "\\.")

p_hill_curve <- ggplot() +
  geom_point(data = selmetals, aes(Ce_mg_l, Qe_mg_g, color = Metal)) +
  geom_line(data = pred_hill, aes(Ce_mg_l, Qe_pred, group = Metal), color = "black") +
  facet_grid(Metal ~ Biochar, scales = "free") + theme_bw()
ggsave("output/hill_fit_curves_full.png", p_hill_curve, width = 10, height = 6)

# 低浓度放大图
p_hill_low <- p_hill_curve + xlim(c(0,2))
ggsave("output/hill_fit_curves_low.png", p_hill_low, width = 10, height = 6)

# 参数 a (Qmax) 图
params_a <- hill_params %>% filter(term == "a")
p_a <- ggplot(params_a, aes(x = Biochar, y = estimate, fill = Metal)) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Metal, scales = "free_y") + theme_bw() +
  labs(y = bquote('Qmax ('*mg.g^-1*')')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("output/hill_param_a.png", p_a, width = 8, height = 5)

# 参数 c 和 b 类似（可选）
params_c <- hill_params %>% filter(term == "c")
if(nrow(params_c) > 0){
  p_c <- ggplot(params_c, aes(x = Biochar, y = estimate, fill = Metal)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  width = 0.2, position = position_dodge(0.9)) +
    facet_wrap(~Metal, scales = "free_y") + theme_bw()
  ggsave("output/hill_param_c.png", p_c, width = 8, height = 5)
}

params_b <- hill_params %>% filter(term == "b")
if(nrow(params_b) > 0){
  p_b <- ggplot(params_b, aes(x = Biochar, y = estimate, fill = Metal)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  width = 0.2, position = position_dodge(0.9)) +
    facet_wrap(~Metal, scales = "free_y") + theme_bw()
  ggsave("output/hill_param_b.png", p_b, width = 8, height = 5)
}

# ========================= 4. 统计比较（成员4简化版） =========================
cat("步骤4：统计比较（基于 Hill 参数的 z 检验）...\n")

# 两两 z 检验函数
pairwise_z <- function(df_param) {
  biochars <- unique(df_param$Biochar)
  if(length(biochars) < 2) return(NULL)
  res <- data.frame()
  for(i in 1:(length(biochars)-1)) {
    for(j in (i+1):length(biochars)) {
      est1 <- df_param$estimate[df_param$Biochar == biochars[i]]
      se1 <- df_param$std.error[df_param$Biochar == biochars[i]]
      est2 <- df_param$estimate[df_param$Biochar == biochars[j]]
      se2 <- df_param$std.error[df_param$Biochar == biochars[j]]
      if(length(est1)==0 || length(est2)==0) next
      diff <- est1 - est2
      se_diff <- sqrt(se1^2 + se2^2)
      z <- diff / se_diff
      p_val <- 2 * (1 - pnorm(abs(z)))
      res <- rbind(res, data.frame(Biochar1 = biochars[i], Biochar2 = biochars[j],
                                   diff = diff, z = z, p_value = p_val))
    }
  }
  return(res)
}

# 按金属分组进行两两比较
metal_list <- unique(params_a$Metal)
ztest_list <- list()
for(m in metal_list) {
  sub <- params_a[params_a$Metal == m, ]
  if(nrow(sub) >= 2) {
    tmp <- pairwise_z(sub)
    if(!is.null(tmp) && nrow(tmp)>0) {
      tmp$Metal <- m
      ztest_list[[m]] <- tmp
    }
  }
}
ztest_results <- bind_rows(ztest_list)
write.csv(ztest_results, "output/pairwise_ztest_results.csv", row.names = FALSE)

# 显著差异输出
sig_res <- ztest_results %>% filter(p_value < 0.05)
if(nrow(sig_res) > 0) {
  write.csv(sig_res, "output/significant_pairs.csv", row.names = FALSE)
}

# ========================= 5. 生成最终报告（HTML） =========================
cat("步骤5：生成最终报告...\n")

# 创建临时 Rmd 文件
rmd_content <- '
---
title: "复现论文：植物基生物炭对Cd、Pb、Zn吸附等温线的建模与比较"
author: "小组作业（一键复现版本）"
date: "`r Sys.Date()`"
output: html_document
---