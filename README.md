# 复现论文：植物基生物炭对 Cd(II)、Pb(II)、Zn(II) 的吸附等温线分析
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue)](https://www.r-project.org/)
## 项目简介
本仓库为 **D2RS-2026spring** 课程小组作业，复现论文 **Essibu et al. (2025)** 中关于四种植物基生物炭（OSR550、MSP550、MSP700、PKS550）对 Cd(II)、Pb(II)、Zn(II) 吸附等温线的数据分析部分。使用 Langmuir 和 Hill 模型进行非线性拟合，比较模型优劣，并基于 Hill 模型参数进行统计检验。

**原始论文**：Essibu, A. K., et al. (2025). *Contrasting effects of Corn cob and Cocoa pod husk biochars on Heavy metal Bioavailability, Speciation, and Uptake by Maize in a Mining-Contaminated soil*. West African Journal of Applied Ecology.（原文链接：https://doi.org/10.1016/j.heliyon.2020.e05388 ）  
**数据来源**：Mendeley Data （数据链接：https://data.mendeley.com/datasets/wk8m4t64dh/2 ）

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
