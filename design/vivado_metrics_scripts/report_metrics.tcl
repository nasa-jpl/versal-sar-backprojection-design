# By: Austin Owens
# Date: 7/9/2024
# Desc: Runs vivado on built XSA to get power and resource utilization

open_run impl_1

# Utilization
report_utilization -file utilization_hierarchical.txt -hierarchical_depth 2 -hierarchical
report_utilization -file utilization_report.txt

# Power
report_power -file power_hierarchical.txt -hierarchical_depth 3 
report_power -file {power.txt} -xpe {power.xpe} -rpx {power.rpx}
