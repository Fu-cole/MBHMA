# 复现论文：植物基生物炭对 Cd(II)、Pb(II)、Zn(II) 的吸附等温线分析
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue)](https://www.r-project.org/)
## 项目简介
本仓库为 **D2RS-2026spring** 课程小组作业，复现论文 **Essibu et al. (2025)** 中关于四种植物基生物炭（OSR550、MSP550、MSP700、PKS550）对 Cd(II)、Pb(II)、Zn(II) 吸附等温线的数据分析部分。使用 Langmuir 和 Hill 模型进行非线性拟合，比较模型优劣，并基于 Hill 模型参数进行统计检验。

**原始论文**：Essibu, A. K., et al. (2025). *Contrasting effects of Corn cob and Cocoa pod husk biochars on Heavy metal Bioavailability, Speciation, and Uptake by Maize in a Mining-Contaminated soil*. West African Journal of Applied Ecology.（原文链接：https://doi.org/10.1016/j.heliyon.2020.e05388 ）  
**数据来源**：Mendeley Data （数据链接：https://data.mendeley.com/datasets/wk8m4t64dh/2 ）

## 仓库结构
├── Code/ # 所有 R 脚本（按成员分工组织）

│ ├── member1_data_cleaning.R # 数据清洗与探索性分析

│ ├── member2_langmuir.R # Langmuir 模型拟合

│ ├── member3_hill.R # Hill 模型拟合

│ ├── member4_statistical_tests.R# 统计比较（z 检验）

│ └── run_full_analysis.R # 一键运行所有脚本

├── Data/ # 原始数据（只读，不修改）

│ └── allmetals2_isotherm.csv # 吸附等温线数据

├── output/ # 所有生成的结果

├── renv.lock # R 包版本锁定文件（用于完全复现环境）

└── README.md # 项目说明文档（本文件）

## 小组分工

| 成员   | GitHub 用户名 | 负责模块                 | 贡献文件                           |
|--------|---------------|--------------------------|------------------------------------|
| 成员1  | @member1      | 数据清洗与探索性分析     | `Code/member1_data_cleaning.R`     |
| 成员2  | @member2      | Langmuir 模型拟合        | `Code/member2_langmuir.R`          |
| 成员3  | @member3      | Hill 模型拟合            | `Code/member3_hill.R`              |
| 成员4  | @Hackysaw     | 统计比较（z 检验）       | `Code/member4_statistical_tests.R` |
| 成员5  | @Fu-cole      | 可视化与报告整合         | `Code/run_full_analysis.R`（以及报告文档） |


## 复现步骤

### 1. 克隆仓库
```bash
git clone https://github.com/D2RS-2026spring/MBHMA.git
cd MBHMA
```
### 2. 准备R环境
```bash
install.packages("renv")
renv::activate()   # 初始化项目环境（创建 renv 文件夹）
renv::restore()    # 根据 renv.lock 安装所有包
```
### 3. 运行完整分析
```bash
Rscript scripts/run_full_analysis.R
```
