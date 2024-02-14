library(data.table)

# Feedback from Hanqing
## {"AMY-1": "AMY-2", "AMY-2": "AMY-1"}, {"ACB-2": "CP-1", "CP-1": "ACB-2"}
## 1. AMY那个我不确定是我们谁反了，因为cell-type 比较像，但我们肯定是不一样的。
## 2. CP那个应该是ATAC 反了

# For current sample data we have
## ACB-2: 4E, CEMBA180110_4E and CEMBA180111_4E
## CP-1: 4D, CEMBA171214_4D and CEMBA171219_4D
## AMY-1: 7H, CEMBA200820_7H and CEMBA200827_7H
## AMY-2: 8H, CEMBA200903_8H and CEMBA200910_8H


# ATAC-seq experiments record:
# https://docs.google.com/spreadsheets/d/1HbYP0tLpv4rPwkJn6uZPjR_M7dnIRglzeTcY9Os-CVA/edit#gid=1893884661

# ATAC-seq LIMS spreadsheet
# https://docs.google.com/spreadsheets/d/1UPkKv3potJtNEbYxkpMY5X_V1xgkRbkvkW4Qi14o1x4/edit#gid=1307716493

# Comments from Yang
## brain dissection编号(4D: CP-1, etc) 肯定是一致的，毫无疑问。
## 上面的描述里，没懂CP-1跟那个region反了？可以查实验室记录吧，看反了的两个sample是不是同一天做的实验，如果不是，就好分辨了？
## 假设Marga没有把sample标记错，tissue被同时发给两个实验室后，某一个实验室的sample标记错了。
## 这样就核对marga，hanqing和我们三处的日期

## 一般是4个sample为一组做实验，比如1号做了CP-1+其他三个sample，5号做了ACB-2+其他3个sample。
## 假设：
##      Marga,   hanqing,   renlab
## CP-1    1         5         1
## others  1         5         1
## ACB-1   5         1         5
## 你就应该知道是哪里出错了 
