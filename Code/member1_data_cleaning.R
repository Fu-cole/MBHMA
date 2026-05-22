# 1. 读取原始数据
allmetals <- read.csv('allmetals2_isotherm.csv')
cat("原始数据维度:", dim(allmetals), "\n")
str(allmetals)

# 3. 初步可视化（全数据）
p1 <- ggplot(data=allmetals, aes(x=Ce_mg_l, y=Qe_mg_g, group=interaction(Biochar,Metal))) +
  theme_bw() +
  geom_point() +
  facet_grid(Metal ~ Biochar, scales="free")
ggsave("initial_scatter_all.png", p1, width=10, height=6)

# 4. 剔除无效生物炭（SWP550, SWP550AW）
selmetals <- allmetals[!allmetals$Biochar %in% c("SWP550","SWP550AW"), ]
cat("剔除SWP550系列后维度:", dim(selmetals), "\n")

# 5. 剔除高浓度异常点（Ce > 50 mg/L）
selmetals <- selmetals[selmetals$Ce_mg_l < 50, ]
cat("剔除Ce>50后维度:", dim(selmetals), "\n")

# 6. 将Biochar和Metal列转为因子（去除多余水平）
selmetals <- droplevels(selmetals)

# 7. 描述性统计
desc_stats <- selmetals %>%
  group_by(Biochar, Metal) %>%
  summarise(n = n(),
            Ce_mean = mean(Ce_mg_l),
            Ce_sd = sd(Ce_mg_l),
            Qe_mean = mean(Qe_mg_g),
            Qe_sd = sd(Qe_mg_g),
            .groups = "drop")
write.csv(desc_stats, "descriptive_stats.csv", row.names=FALSE)

# 8. 清洗后可视化（最终用于建模的数据）
p2 <- ggplot(selmetals, aes(x=Ce_mg_l, y=Qe_mg_g, color=Metal)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~ Biochar, scales="free") +
  labs(title = "Cleaned isotherm data")
ggsave("cleaned_scatter.png", p2, width=8, height=6)

# 9. 保存清洗后的数据供其他成员使用
save(selmetals, file = "data_cleaned.RData")