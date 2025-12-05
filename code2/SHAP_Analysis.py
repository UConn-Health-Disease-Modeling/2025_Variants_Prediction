import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import pandas as pd
import numpy as np

# 让 rpy2 自动把 R data.frame 转成 pandas DataFrame
pandas2ri.activate()

# 绑定 R 的 readRDS 函数
read_rds = ro.r["readRDS"]

# ----------- 读取 RDS -----------
res_peak_r = read_rds("code2/peak_results_1117.rds")
res_dur_r  = read_rds("code2/duration_results_1117.rds")
all_data_r = read_rds("code2/all_data.rds")

# ----------- 提取内部元素 -----------
# R: res_peak <- res_peak_r$feat_35_ensemble$final_df

res_peak = pandas2ri.rpy2py(
    res_peak_r.rx2("feat_35_ensemble").rx2("final_df")
)

res_dur = pandas2ri.rpy2py(
    res_dur_r.rx2("feat_35_ensemble").rx2("final_df")
)

# all_data 是一个 list, measurements 是其中的一个元素
all_data = pandas2ri.rpy2py(all_data_r)
measurements = all_data["measurements"]


# # ----------- peak_share 分类 -----------

# bins = [-np.inf, 0.05, 0.20, np.inf]
# labels = ["<0.10", "0.10–0.20", ">0.20"]

# measurements["peak_share_cat"] = pd.cut(
#     measurements["peak_share"],
#     bins=bins,
#     labels=labels,
#     right=True
# )

# print(measurements.head())
