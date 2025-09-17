import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_file_path = os.path.expanduser("./cellphonedb_v500_NatProtocol/v5.0.0/cellphonedb.zip")
meta_file_path = "./cellphonedb_meta_Normal_P3.txt"
counts_file_path = "./Normal_P3.h5ad"
out_path = "./cellphonedb_output_Normal"
os.makedirs(out_path, exist_ok=True)







cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path=cpdb_file_path,
    meta_file_path=meta_file_path,
    counts_file_path=counts_file_path,
    counts_data='hgnc_symbol',
    iterations=500,
    threshold=0.1,
    threads=32,
    pvalue=0.05,
    # 可选的高级参数
    # active_tfs_file_path=active_tf_path,    # 如果使用转录因子分析
    # microenvs_file_path=microenvs_file_path, # 如果使用微环境定义
    # subsampling=True,                     # 大数据集时可启用子抽样
    output_path=out_path,
    output_suffix='Normal')

print("CellPhoneDB 分析已完成！")
print(f"结果保存在: {out_path}")
