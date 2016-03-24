# Breast MRI datasets descriptive stats:
Cristina Gallego  
March 20, 2016  

Read in patient information from biomatrix and collects statistics of datasets used in study:
- Patient info:
  * #, Age,
  


```r
options(width = 500)
setwd("C:/Users/windows/Documents/repoCode-local/T2wR/lop_3Dtex_T2w_addedvalue")
library(R.utils)
library("RSQLite")
sqlite <- dbDriver("SQLite")
conn <- dbConnect(sqlite, "textureUpdatedFeatures.db")

# 2) all T1W features
lesionsQuery <- dbGetQuery(conn, "SELECT *
                           FROM  lesion 
                           INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                           INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                           INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                           INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)
                           INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)
                           INNER JOIN radiologyInfo ON (lesion.lesion_id = radiologyInfo.lesion_id)")

# prune entries and extract feature subsets
# corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
lesioninfo = lesionsQuery[c(1:26)]
dynfeatures = lesionsQuery[c(29:62)]
morphofeatures = lesionsQuery[c(65:83)]
texfeatures = lesionsQuery[c(86:129)]
T2info = lesionsQuery[c(259:270)]
T2features = lesionsQuery[c(259,267:291,232:251)]
stage1features = lesionsQuery[c(132:231)]
imagingdetails = lesionsQuery[c(293:318)]

##### set data splits
# combine all features and exclude foci lesions at this point
namest1w = names(cbind(dynfeatures, morphofeatures, texfeatures, stage1features))
namest2w = names(T2features)

# all lesions at the lesion id
allfeatures = cbind(lesioninfo[c("lesion_label")], dynfeatures, morphofeatures, texfeatures, stage1features, T2features)   
# select non foci
lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )

allfeatures = allfeatures[rownames(lesioninfo),]
allfeatures$origlesion_label = factor(allfeatures$lesion_label)
print(summary(allfeatures$origlesion_label))
```

```
##    massB    massM nonmassB nonmassM 
##      240      168      142       77
```

```r
########## Patient 
patientage = data.frame()
for (k in 1:length(lesioninfo$lesion_id)) {
    # find the age of patient at the time of imaging
    days = difftime(lesioninfo$exam_dt_datetime[k], lesioninfo$anony_dob_datetime[k], 
        units = "days")
    age = days[[1]]/365
    patientage = rbind(patientage, c(age))
}

# combine with other data
patientinfo = cbind(lesioninfo$cad_pt_no_txt, patientage, lesioninfo$lesion_label, lesioninfo$BIRADS)
colnames(patientinfo) <- c("CADid", "age", "type", "BIRADS")

# print summary statistics
print(summary(patientinfo))
```

```
##      CADid          age              type     BIRADS 
##  0880   :  5   Min.   :22.37   massB   :240   2: 22  
##  3055   :  5   1st Qu.:40.64   massM   :168   3: 56  
##  6054   :  5   Median :46.86   nonmassB:142   4:343  
##  0252   :  4   Mean   :48.75   nonmassM: 77   5:116  
##  0728   :  4   3rd Qu.:55.29                  6: 90  
##  3005   :  4   Max.   :85.59                         
##  (Other):600
```

```r
# print All cad ids
print(unique(patientinfo$CADid))
```

```
##   [1] 0002 0016 0025 0027 0066 0093 0102 0111 0114 0121 0122 0123 0127 0129 0130 0132 0133 0135 0168 0171 0172 0173 0177 0180 0186 0189 0103 0190 0196 0197 0198 0199 0205 0207 0212 0220 0229 0232 0246 0252 0259 0266 0276 0277 0280 0293 0311 0325 0331 0352 0357 0376 0388 0409 0420 0426 0442 0456 0462 0463 0465 0473 0503 0510 0513 0519 0536 0551 0552 0553 0561 0571 0572 0573 0576 0578 0580 0586 0595 0603 0606 0608 0613 0616 0619 0624 0635 0651 0657 0663 0664 0666 0667 0668 0672 0673 0679 0681 0682
## [100] 0683 0684 0685 0687 0689 0690 0691 0692 0696 0700 0705 0707 0710 0713 0714 0718 0720 0721 0722 0723 0724 0726 0727 0728 0729 0730 0731 0734 0735 0736 0737 0740 0742 0743 0744 0745 0748 0752 0755 0757 0760 0764 0765 0767 0771 0775 0778 0781 0782 0783 0758 0776 0779 0789 0790 0791 0792 0793 0795 0796 0799 0802 0803 0805 0807 0809 0810 0812 0813 0814 0815 0817 0818 0827 0828 0829 0830 0831 0834 0837 0839 0843 0844 0845 0846 0847 0850 0851 0853 0855 0856 0857 0861 0862 0863 0865 0867 0870 0871
## [199] 0873 0875 0876 0877 0880 0883 0884 0885 0887 0888 0896 0898 0900 0904 0913 0918 0920 0921 0924 0934 0937 0943 0944 0950 0952 0954 0956 0962 0965 0967 0978 0985 0993 0995 0997 0999 1003 1004 1006 1008 1012 1018 1021 1024 1025 1027 1044 1045 1050 1026 1053 1062 1065 1071 1072 1077 1078 1079 1081 1086 1087 1090 1092 1095 1099 2007 2016 2017 2023 2024 2027 2028 2029 2033 2042 2049 2050 2051 2053 2055 2059 2065 2068 2069 2071 2072 2073 2075 2078 2079 3004 3005 3010 3011 3017 3018 3020 3021 3023
## [298] 3025 3026 3028 3030 3031 3033 3035 3039 3045 3046 3052 3053 3054 3055 3057 3063 3065 3070 3072 3073 3075 3076 3077 3078 3080 3081 3082 3083 3086 3092 3093 3097 4002 4003 4008 4012 4017 4018 4019 4020 4021 4023 4026 4029 4039 4040 4041 4043 4044 4045 4047 4049 4055 6001 6004 6005 6008 6014 6015 6017 6018 6019 6020 6021 6022 6023 6024 6025 6026 6027 6029 6032 6034 6035 6037 6038 6039 6040 6041 6042 6043 6044 6045 6046 6047 6048 6050 6051 6052 6054 6069 6100 6101 6105 6114 6117 6141 6148 6150
## [397] 6164 6174 6223 6224 6226 6233 7008 7011 7018 7024 7029 7030 7043 7045 7053 7066 7076 7077 7085 7086 7088 7094 7096 7097 7104 7105 7127 7151 7159 7165 7178 7183 7186 7189 7190 7192 7193 7201 7220
## 435 Levels: 0002 0016 0025 0027 0066 0093 0102 0103 0111 0114 0121 0122 0123 0127 0129 0130 0132 0133 0135 0168 0171 0172 0173 0177 0180 0186 0189 0190 0196 0197 0198 0199 0205 0207 0212 0220 0229 0232 0246 0252 0259 0266 0276 0277 0280 0293 0311 0325 0331 0352 0357 0376 0388 0409 0420 0426 0442 0456 0462 0463 0465 0473 0503 0510 0513 0519 0536 0551 0552 0553 0561 0571 0572 0573 0576 0578 0580 0586 0595 0603 0606 0608 0613 0616 0619 0624 0635 0651 0657 0663 0664 0666 0667 0668 0672 0673 ... 7220
```

```r
print(length(unique(patientinfo$CADid)))
```

```
## [1] 435
```

```r
# print Age stats
print(summary(patientinfo$age))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   22.37   40.64   46.86   48.75   55.29   85.59
```

```r
print(sd(patientinfo$age))
```

```
## [1] 10.58338
```

```r
range(patientinfo$age)
```

```
## [1] 22.36712 85.59167
```

```r
# Range of imaging studies
range(lesioninfo$exam_dt_datetime)
```

```
## [1] "2007-03-15 00:00:00.000000" "2014-11-15 00:00:00.000000"
```

```r
# time between imaging and procedure for biopsies only
print(summary(factor(lesioninfo$proc_proc_source_int)))
```

```
##                                                            NA                                                     Radiology Surgical/Operating Rm (includes 'Sentinel Lymph Node Biopsy') 
##                                                             2                                                           566                                                            59
```

```r
print(summary(factor(lesioninfo$proc_proc_tp_int)))
```

```
##     Core Needle Biopsy Fine Needle Aspiration                    N/A                     NA Vacuum Assisted Biopsy 
##                    313                     14                     60                      5                    235
```

```r
biopsies = subset(lesioninfo, proc_proc_source_int!="Surgical/Operating Rm (includes 'Sentinel Lymph Node Biopsy')")
print(summary(factor(biopsies$proc_proc_tp_int)))
```

```
##     Core Needle Biopsy Fine Needle Aspiration                    N/A                     NA Vacuum Assisted Biopsy 
##                    313                     14                      4                      2                    235
```

```r
# for biopsies
procinfo = c()
for (k in 1:length(biopsies$proc_pt_procedure_id)) {
    # find the age of patient at the time of imaging
    days = difftime(biopsies$proc_proc_dt_datetime[k], biopsies$exam_dt_datetime[k], 
        units = "days")
    print(days[[1]])
    procinfo = c(procinfo, days[[1]])
}
```

```
## [1] 32
## [1] 19
## [1] 13
## [1] 68
## [1] 10
## [1] 25.95833
## [1] 24.04167
## [1] 43.04167
## [1] 43.04167
## [1] 28
## [1] 26
## [1] 23.95833
## [1] 26
## [1] 52
## [1] 38
## [1] 38
## [1] 21
## [1] 10
## [1] -461.9583
## [1] -461.9583
## [1] 23
## [1] 462
## [1] 462
## [1] 16.95833
## [1] 23
## [1] 25
## [1] 15
## [1] 17.95833
## [1] 54
## [1] 36
## [1] 48.04167
## [1] 28
## [1] 29
## [1] 19
## [1] 16
## [1] 16
## [1] 36
## [1] 36
## [1] 36
## [1] 34
## [1] 34
## [1] 34
## [1] 18
## [1] 30
## [1] 17
## [1] 38
## [1] 54.95833
## [1] 54.95833
## [1] 12
## [1] 9
## [1] 9
## [1] 42
## [1] 14
## [1] 35
## [1] 386
## [1] 386
## [1] 0
## [1] 18
## [1] 19
## [1] 19
## [1] 9
## [1] 0
## [1] 25
## [1] 39
## [1] 12
## [1] 13
## [1] 37.95833
## [1] 11
## [1] 30
## [1] 32
## [1] 10
## [1] 15
## [1] 28
## [1] 31
## [1] 24
## [1] 31
## [1] 23.95833
## [1] 29.95833
## [1] 10
## [1] 233.0417
## [1] 39
## [1] 32
## [1] 13
## [1] 13
## [1] 33
## [1] 32
## [1] 28
## [1] 35
## [1] 15
## [1] -163.9583
## [1] 27
## [1] 18
## [1] 32
## [1] 52
## [1] 53.04167
## [1] 21.04167
## [1] -12
## [1] 61.04167
## [1] 16
## [1] 10
## [1] 25.95833
## [1] 57.04167
## [1] 9
## [1] 30
## [1] 49
## [1] 6
## [1] 8
## [1] 8
## [1] 25
## [1] 13
## [1] 24
## [1] 2
## [1] 19
## [1] 27
## [1] 807.9583
## [1] 27
## [1] 25
## [1] 70.95833
## [1] 29
## [1] 79
## [1] 4
## [1] 12
## [1] 27
## [1] 30
## [1] 30
## [1] 30
## [1] 66
## [1] 268.0417
## [1] 268.0417
## [1] 13
## [1] 10
## [1] 10
## [1] 8
## [1] 27
## [1] 9
## [1] 11
## [1] 10
## [1] 10
## [1] 17
## [1] 17
## [1] 212.0417
## [1] 4
## [1] 4
## [1] 6
## [1] 6
## [1] 6
## [1] 14
## [1] 20
## [1] 20
## [1] 2
## [1] 5
## [1] 3
## [1] 45
## [1] 13
## [1] 13
## [1] 13
## [1] 13
## [1] 11
## [1] 25
## [1] 26
## [1] 2
## [1] 23
## [1] 24
## [1] 14
## [1] 25
## [1] 25
## [1] 17
## [1] 26
## [1] 15
## [1] -23
## [1] 1
## [1] 1
## [1] 1
## [1] 23
## [1] 23
## [1] 12
## [1] 15.04167
## [1] 15.04167
## [1] 35.04167
## [1] 29
## [1] -15
## [1] -15
## [1] 42
## [1] 27
## [1] 13
## [1] 136
## [1] 28.04167
## [1] 28.04167
## [1] 28
## [1] 28
## [1] 2
## [1] 2
## [1] 45
## [1] 25
## [1] 25
## [1] 48
## [1] 48
## [1] 18
## [1] -540
## [1] -11
## [1] -25
## [1] 12
## [1] 29.95833
## [1] 2
## [1] 12
## [1] 8
## [1] 32
## [1] 42
## [1] 7
## [1] 17
## [1] 17
## [1] 12
## [1] 25
## [1] 5
## [1] 15
## [1] -28
## [1] 38
## [1] 8
## [1] 26
## [1] 23
## [1] 23
## [1] 5
## [1] 12
## [1] 25.95833
## [1] 19.95833
## [1] 19.95833
## [1] 12
## [1] 63
## [1] 14.04167
## [1] 6
## [1] 3
## [1] 29
## [1] 29
## [1] 29
## [1] 670.0417
## [1] 670.0417
## [1] 945
## [1] 26
## [1] 23.04167
## [1] 23.04167
## [1] 19
## [1] 61
## [1] 412
## [1] 1
## [1] 37
## [1] 37
## [1] 17
## [1] 31
## [1] 2
## [1] 13
## [1] -29.04167
## [1] 32
## [1] 35.04167
## [1] 18
## [1] 24
## [1] 62
## [1] 20
## [1] 20
## [1] 20
## [1] 62
## [1] 4
## [1] 5
## [1] 5
## [1] 21
## [1] 32
## [1] 5
## [1] 17.95833
## [1] 23.04167
## [1] 37.95833
## [1] 61
## [1] 51
## [1] 24
## [1] 17
## [1] 38
## [1] 44
## [1] 44
## [1] 5
## [1] 13
## [1] 10
## [1] 19
## [1] 21
## [1] 28
## [1] 2
## [1] 2
## [1] 31
## [1] 37
## [1] 10
## [1] 23
## [1] 46
## [1] 28
## [1] 22
## [1] 22
## [1] 57
## [1] 30
## [1] 28
## [1] 52
## [1] 52
## [1] -298
## [1] 39.04167
## [1] 35
## [1] 18
## [1] 21
## [1] 4
## [1] 36.04167
## [1] 54.04167
## [1] 54.04167
## [1] 26
## [1] 24
## [1] 2
## [1] 2
## [1] 11
## [1] 30
## [1] 8.958333
## [1] 12
## [1] 38
## [1] 13.04167
## [1] 25.95833
## [1] 14
## [1] 14
## [1] 23
## [1] 11
## [1] 15
## [1] 15
## [1] 34
## [1] 37
## [1] 37
## [1] 18.95833
## [1] 6
## [1] 13
## [1] 7
## [1] 7
## [1] 42
## [1] 14
## [1] 14
## [1] 8
## [1] 25.04167
## [1] 25.04167
## [1] 15
## [1] 14
## [1] 22
## [1] -26
## [1] 25
## [1] 1
## [1] 1
## [1] 5
## [1] -13
## [1] 13
## [1] 8
## [1] -720
## [1] 6
## [1] 6
## [1] -17
## [1] 23
## [1] 23
## [1] 24
## [1] 43
## [1] 10
## [1] 16
## [1] 18
## [1] 17
## [1] 17
## [1] 21
## [1] 22
## [1] 22
## [1] 41.95833
## [1] 9
## [1] 9
## [1] -49
## [1] 22
## [1] 832.9583
## [1] 22
## [1] -40
## [1] 48
## [1] 54.95833
## [1] 47
## [1] 27
## [1] 17
## [1] 2
## [1] 40
## [1] 36
## [1] 41
## [1] 1
## [1] -7
## [1] 38
## [1] 80.04167
## [1] 29
## [1] 47.04167
## [1] 13
## [1] 31
## [1] 13
## [1] 87
## [1] 87
## [1] 17
## [1] 16
## [1] 24
## [1] 10
## [1] 8.958333
## [1] 10
## [1] 10
## [1] 8.958333
## [1] 20
## [1] 13
## [1] 29.95833
## [1] 5
## [1] 5
## [1] 26
## [1] 11
## [1] 6
## [1] 6
## [1] 16
## [1] 16
## [1] 19
## [1] 19
## [1] 36
## [1] 36
## [1] 13
## [1] 95.04167
## [1] 15
## [1] 11
## [1] 9
## [1] 13
## [1] 30
## [1] 53.95833
## [1] 53.95833
## [1] 21
## [1] 21
## [1] 21
## [1] 41.95833
## [1] 41.95833
## [1] 749
## [1] 70
## [1] 36
## [1] 11
## [1] 25
## [1] 28
## [1] 234.9583
## [1] 35
## [1] 17
## [1] 17
## [1] 15
## [1] -8
## [1] 10
## [1] 20
## [1] 42.95833
## [1] 12
## [1] 12
## [1] 12
## [1] 24
## [1] 34
## [1] 34
## [1] 9
## [1] -17
## [1] -17
## [1] 14
## [1] 14
## [1] 14
## [1] 6
## [1] 6
## [1] -42
## [1] -42
## [1] 17
## [1] -9
## [1] 6
## [1] 40
## [1] 8
## [1] 38
## [1] 38
## [1] 6
## [1] 6
## [1] 11
## [1] 11
## [1] 11
## [1] 22.95833
## [1] 310
## [1] 310
## [1] 620
## [1] 620
## [1] 51
## [1] 45
## [1] 1
## [1] 1
## [1] 11
## [1] 11
## [1] 11
## [1] 5
## [1] 5
## [1] 1340
## [1] 27
## [1] 27
## [1] 27
## [1] 16.04167
## [1] 16.04167
## [1] 283.9583
## [1] 12
## [1] 13
## [1] 13
## [1] 14
## [1] 14
## [1] 0
## [1] 0
## [1] 0
## [1] 0
## [1] -6
## [1] 4
## [1] 37
## [1] 3.041667
## [1] 49
## [1] 49
## [1] 49
## [1] 49
## [1] 49
## [1] 31.04167
## [1] 31.04167
## [1] 7.958333
## [1] 17.04167
## [1] -26
## [1] 10
## [1] 33
## [1] 73
## [1] 20
## [1] 23
## [1] 28
## [1] 98
## [1] 38
## [1] 53
## [1] 34
## [1] 29
## [1] 14
## [1] 27
## [1] 338
## [1] -345
## [1] 28
## [1] 23
## [1] 13
## [1] 13
## [1] -688.0417
## [1] 49
## [1] 22.04167
## [1] 32
## [1] 20
## [1] 24
## [1] 24
## [1] 16
## [1] 19
## [1] 31
## [1] 52
## [1] 5
## [1] 5
## [1] 24
## [1] 24
## [1] 180.0417
## [1] 1
## [1] 16.04167
## [1] 17
## [1] 39.95833
## [1] 19
## [1] 13
## [1] -31
## [1] 37
## [1] -9
## [1] -678.0417
## [1] 27.95833
## [1] 15
## [1] 8.041667
## [1] 22.95833
## [1] 22.95833
## [1] 16
## [1] 21
```

```r
# time after imgaing stats
summary(procinfo)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -720.00   11.00   20.00   33.82   32.00 1340.00
```

```r
# plot By default, geom_bar uses stat='bin'. This makes the height of each
# bar equal to the number of cases in each group, and it is incompatible
# with mapping values to the y aesthetic.
require(ggplot2)
IQR = c(summary(patientinfo$age)[[2]], summary(patientinfo$age)[[5]])
histages <- ggplot(patientinfo, aes(age))
histages + geom_bar(width=1.5) + 
  geom_vline(xintercept = range(patientinfo$age), linetype = "dotted") + 
  geom_vline(xintercept = mean(patientinfo$age), linetype = "longdash", colour = "red") + 
  geom_vline(xintercept = median(patientinfo$age), linetype = "longdash", colour = "blue") + 
  geom_vline(xintercept = IQR, linetype = "longdash", colour = "green") + 
  labs(x = "age in years", y = "# patients",  title = "Patients age at the time of imaging")
```

![](descriptive_Datasets_files/figure-html/unnamed-chunk-1-1.png) 

```r
bartypes <- ggplot(patientinfo, aes(type, fill = BIRADS))
bartypes + geom_bar() + labs(x = "type of lesion", y = "# patients", title = "Type of lesion and BIRADS category")
```

![](descriptive_Datasets_files/figure-html/unnamed-chunk-1-2.png) 

```r
library(reshape)
# subset by pathology
patho = summary(as.factor(lesioninfo$lesion_diagnosis))
for (k in 1:length(patho)) {
    print(names(patho[k]))
    print(cast(lesioninfo, ~lesion_label, subset = lesion_diagnosis == 
        names(patho[k]), value = "lesion_id", length))
}
```

```
## [1] "Adenocarcinoma"
##   value nonmassM
## 1 (all)        1
## [1] "ADENOSIS"
##   value massB nonmassB
## 1 (all)     9        2
## [1] "ADENOSIS, COLUMNAR CELL CHANGES"
##   value massB
## 1 (all)     1
## [1] "ATYPICAL DUCTAL HYPERPLASIA"
##   value massB massM nonmassB
## 1 (all)    10     1       16
## [1] "ATYPICAL LOBULAR HYPERPLASIA"
##   value massB nonmassB
## 1 (all)     9        3
## [1] "ATYPICAL PAPILLARY LESION"
##   value massB
## 1 (all)     1
## [1] "AtypicalCells"
##   value massB
## 1 (all)     1
## [1] "AtypicalPapilloma"
##   value massB
## 1 (all)     1
## [1] "BENIGN BREAST TISSUE"
##   value massB nonmassB
## 1 (all)    40       32
## [1] "BENIGN HAMARTOMA"
##   value massB
## 1 (all)     1
## [1] "BENIGN INTRAMAMMARY LYMPH NODE"
##   value massB
## 1 (all)     1
## [1] "benign lymphoid tissue"
##   value massB
## 1 (all)     1
## [1] "capillary hemangioma"
##   value massB
## 1 (all)     1
## [1] "COLUMNAR CELL CHANGES"
##   value massB nonmassB
## 1 (all)     5       11
## [1] "ColumnarAlterationwoAtypia"
##   value nonmassB
## 1 (all)        2
## [1] "COMPLEX FIBROEPITHELIAL LESION"
##   value massB
## 1 (all)     1
## [1] "COMPLEX PAPILLARY LESION"
##   value nonmassB
## 1 (all)        2
## [1] "Cyst"
##   value massB
## 1 (all)     3
## [1] "DENSE FIBROSIS"
##   value nonmassB
## 1 (all)        5
## [1] "DENSE FIBROSIS AND FIBROADENOMATOID CHANGE"
##   value massB
## 1 (all)     1
## [1] "DUCT PAPILLOMA"
##   value massB nonmassB
## 1 (all)    10        5
## [1] "DUCT PAPILLOMA WITH ATYPIA"
##   value massB nonmassB
## 1 (all)     1        1
## [1] "DUCTAL HYPERPLASIA"
##   value massB
## 1 (all)     1
## [1] "DUCTAL HYPERPLASIA WITHOUT ATYPIA"
##   value massB
## 1 (all)     2
## [1] "DYSTROPHICCALCIFICATION"
##   value massB
## 1 (all)     1
## [1] "FAT NECROSIS"
##   value massB
## 1 (all)     1
## [1] "FIBROADENOMA"
##   value massB nonmassB
## 1 (all)    59       10
## [1] "FIBROADENOMATOID"
##   value massB
## 1 (all)     1
## [1] "FIBROCYSTIC"
##   value massB nonmassB
## 1 (all)    38       32
## [1] "FIBROEPITHELIAL"
##   value massB
## 1 (all)     3
## [1] "FIBROSIS"
##   value massB nonmassB
## 1 (all)     2        3
## [1] "FIBROTIC STROMA"
##   value nonmassB
## 1 (all)        1
## [1] "FLAT EPITHELIAL ATYPIA"
##   value massB
## 1 (all)     1
## [1] "FLORID DUCT HYPERPLASIA"
##   value massB
## 1 (all)     1
## [1] "FLORID DUCTAL HYPERPLASIA"
##   value nonmassB
## 1 (all)        2
## [1] "FLORID HYPERPLASIA"
##   value massB
## 1 (all)     1
## [1] "FOCAL CELLULAR STROMA"
##   value massB
## 1 (all)     1
## [1] "FOCAL HYPERPLASIA"
##   value nonmassB
## 1 (all)        1
## [1] "FOCAL USUAL DUCTAL HYPERPLASIA"
##   value massB nonmassB
## 1 (all)     1        1
## [1] "GRANULOMATOUS LOBULAR MASTITIS"
##   value massB
## 1 (all)     2
## [1] "HYPERPLASIA"
##   value massB
## 1 (all)     4
## [1] "IN SITU PAPILLARY CARCINOMA"
##   value nonmassM
## 1 (all)        1
## [1] "INFLAMED CYST WALL"
##   value massB
## 1 (all)     1
## [1] "InsituDuctal"
##   value massM nonmassM
## 1 (all)    40       40
## [1] "InsituLobular"
##   value massB nonmassB
## 1 (all)     4        3
## [1] "INSITUPAPILLARYCARCINOMA"
##   value massM nonmassM
## 1 (all)     1        1
## [1] "InvasiveDuctal"
##   value massM nonmassM
## 1 (all)   107       27
## [1] "InvasiveDuctal micropapillary"
##   value massM
## 1 (all)     1
## [1] "InvasiveLobular"
##   value massM nonmassM
## 1 (all)    15        6
## [1] "LARGE DUCT PAPILLOMA"
##   value nonmassB
## 1 (all)        1
## [1] "LobularHyperplasia"
##   value massB nonmassB
## 1 (all)     2        2
## [1] "LYMPHOVASCULAR INVASION"
##   value nonmassM
## 1 (all)        1
## [1] "MetaplasticCarcinoma"
##   value massM
## 1 (all)     1
## [1] "PAPILLARY LESION"
##   value nonmassB
## 1 (all)        1
## [1] "Papillary(focalAtypia)"
##   value massB
## 1 (all)     1
## [1] "PHYLLODES TUMOR"
##   value massM
## 1 (all)     2
## [1] "PSEUDOANGIOMATOUS STROMAL HYPERPLASIA"
##   value massB
## 1 (all)     1
## [1] "RADIAL SCLEROSING LESION"
##   value massB
## 1 (all)     1
## [1] "RadialScar"
##   value massB
## 1 (all)     1
## [1] "SCLEROSING ADENOSIS"
##   value massB nonmassB
## 1 (all)     4        2
## [1] "SCLEROSING INTRADUCTAL PAPILLOMA"
##   value nonmassB
## 1 (all)        1
## [1] "SCLEROSING PAPILLARY LESION"
##   value massB
## 1 (all)     1
## [1] "STROMAL FIBROSIS"
##   value massB nonmassB
## 1 (all)     1        3
## [1] "STROMAL HYPERPLASIA"
##   value massB
## 1 (all)     2
## [1] "TUBULAR ADENOMA"
##   value massB
## 1 (all)     1
## [1] "USUAL DUCTAL HYPERPLASIA"
##   value massB
## 1 (all)     4
```

```r
benigns = c(72, 34, 16, 7, 14, 4, 19, 7, 70, 20, 70, 19, 18, 12) # n=382
maligns = c(80, 3, 136, 21, 1, 1, 2, 1)
print(sum(benigns))
```

```
## [1] 382
```

```r
print(sum(maligns))
```

```
## [1] 245
```

```r
# percentages
print(sum(benigns))/print(sum(benigns)+sum(maligns))
```

```
## [1] 382
## [1] 627
```

```
## [1] 0.6092504
```

```r
print(sum(maligns))/print(sum(benigns)+sum(maligns))
```

```
## [1] 245
## [1] 627
```

```
## [1] 0.3907496
```

```r
print(benigns/sum(benigns)*100)
```

```
##  [1] 18.848168  8.900524  4.188482  1.832461  3.664921  1.047120  4.973822  1.832461 18.324607  5.235602 18.324607  4.973822  4.712042  3.141361
```

```r
print(maligns/sum(maligns)*100)
```

```
## [1] 32.6530612  1.2244898 55.5102041  8.5714286  0.4081633  0.4081633  0.8163265  0.4081633
```

Analize type of scans
==========

```r
seqdf = data.frame()
caseids = lesioninfo$lesion_id
for(k in 1:length(caseids)){
  df = data.frame(id=caseids[k])
  radreport = imagingdetails[imagingdetails$lesion_id==caseids[k],"original_report_txt"]
  #cat(radreport)
  if(grepl("VIBRANT", radreport)){df$T1w="VIBRANT"}else{
    df$T1w="REVIEW"; 
    }
  if(grepl("T2 weighted FSE", radreport)){df$T2w="T2 weighted FSE"}else{
    df$T2w="REVIEW"; print(caseids[k])
  }
  
  # so far print update
  print(df)
  print(lesioninfo[lesioninfo$lesion_id==caseids[k],c(1,3,6:7,13,24:26)])
  ## find prior histories
  lesioncomts = lesioninfo[lesioninfo$lesion_id==caseids[k],"proc_lesion_comments_txt"]
  cat(lesioncomts)
  surgery=c(" ")
  if(grepl("mammoplasties", lesioncomts)){surgery = list(c(surgery, "mammoplasties, "))}
  if(grepl("lumpectomy", lesioncomts)){surgery = list(c(surgery, "lumpectomy, "))}
  
  df$surgery = list(c(surgery,lesioncomts))
  
  # asppend
  seqdf = rbind(seqdf, df)
  
   ## pause  
   #   cat ("Press [enter] to continue")
   #   line <- readline()
}
```

```
##   id     T1w             T2w
## 1  1 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 1         1          0002           6745896 2011-04-03 00:00:00.000000                 2067      4     nonmassM     InsituDuctal
## None  id     T1w             T2w
## 1  2 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label        lesion_diagnosis
## 2         2          0016           6920252 2011-09-24 00:00:00.000000                 2677      4        massB FLORID DUCT HYPERPLASIA
## None  id     T1w             T2w
## 1  3 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 3         3          0025           7002835 2011-12-28 00:00:00.000000                 2678      4     nonmassB   DENSE FIBROSIS
## None  id     T1w             T2w
## 1  4 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 4         4          0027           7171944 2012-08-15 00:00:00.000000                 3739      4     nonmassB STROMAL FIBROSIS
## CLINICAL INDICATION:  rt sided multifocal idc with large area of dcis. 1st post-op mri for re-evaluation  id     T1w             T2w
## 1  5 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 5         5          0027           6805483 2011-05-27 00:00:00.000000                 2680      4     nonmassM     InsituDuctal
## None  id     T1w             T2w
## 1  6 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 6         6          0066           4583735 2008-02-17 00:00:00.000000                 2685      3        massB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1  7 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 7         7          0066           7556910 2013-10-12 00:00:00.000000                 3637      4     nonmassB   DENSE FIBROSIS
## CLINICAL INDICATION:  OBSP HR >25% LT risk  id     T1w             T2w
## 1  8 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 8         8          0093           7156466 2012-10-16 00:00:00.000000                 2692      4     nonmassM   InvasiveDuctal
## None  id     T1w             T2w
## 1  9 VIBRANT T2 weighted FSE
##   lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 9         9          0093           7156466 2012-10-16 00:00:00.000000                 2692      4     nonmassM   InvasiveDuctal
## None  id     T1w             T2w
## 1 10 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 10        10          0102           4755778 2009-03-19 00:00:00.000000                 2708      4        massB      FIBROCYSTIC
## None[1] 11
##   id     T1w    T2w
## 1 11 VIBRANT REVIEW
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 11        11          0111           6907205 2011-09-15 00:00:00.000000                 2726      4     nonmassB   DUCT PAPILLOMA
## None  id     T1w             T2w
## 1 12 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 12        12          0114           6896014 2011-10-02 00:00:00.000000                 2267      4        massM   InvasiveDuctal
## High risk screening study. BRCA 1 mutation carrier. 6 month follow up probably benign enhancment. Reduction mammoplasties 1997. No HRT or supplements, has gained weight.  id     T1w             T2w
## 1 13 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 13        13          0121           6714524 2011-02-22 00:00:00.000000                 1988      4        massB         ADENOSIS
## None  id     T1w             T2w
## 1 14 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 14        14          0121           7091267 2012-08-29 00:00:00.000000                 2727      4        massB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 15 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 15        15          0122           5108281 2009-12-14 00:00:00.000000                  611      3        massB             Cyst
## None[1] 16
##   id     T1w    T2w
## 1 16 VIBRANT REVIEW
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 16        16          0123           6909758 2011-09-16 00:00:00.000000                 2728      4     nonmassB COLUMNAR CELL CHANGES
## None[1] 17
##   id     T1w    T2w
## 1 17 VIBRANT REVIEW
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 17        17          0123           6909758 2011-09-16 00:00:00.000000                 2728      4     nonmassB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 18 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 18        18          0127           4696964 2008-09-11 00:00:00.000000                 2729      4     nonmassB     FIBROADENOMA
## unclassified variant in BRCA 2  id     T1w             T2w
## 1 19 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 19        19          0129           5326737 2010-07-17 00:00:00.000000                  113      4        massB BENIGN BREAST TISSUE
## INDICATION: History of brown/bloody left nipple discharge.
## Failed ductogram June 25 2010. H/o brownish nipple discharge
## (left side). No palpable abnormality.  id     T1w             T2w
## 1 20 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 20        20          0130           5017534 2010-03-29 00:00:00.000000                 1489      2        massB ATYPICAL LOBULAR HYPERPLASIA
## None  id     T1w             T2w
## 1 21 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 21        21          0130           5017534 2010-03-29 00:00:00.000000                 1489      2        massB ATYPICAL LOBULAR HYPERPLASIA
## None  id     T1w             T2w
## 1 22 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 22        22          0130           7347205 2013-05-05 00:00:00.000000                 3822      4        massB ATYPICAL LOBULAR HYPERPLASIA
## CLINICAL INDICATION:  OBSP high risk screening.  id     T1w             T2w
## 1 23 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 23        23          0132           5154279 2010-04-08 00:00:00.000000                 3823      3        massB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 24 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 24        24          0132           5154279 2010-04-08 00:00:00.000000                 3823      3        massB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 25 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 25        25          0133           7072006 2012-03-10 00:00:00.000000                 3749      4        massB      FIBROCYSTIC
## CLINICAL INDICATION: High risk Mutation carrier.  TAH \T\ BSO in 2005.  On Vagifem  id     T1w             T2w
## 1 26 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 26        26          0135           7777131 2014-05-17 00:00:00.000000                 4111      4        massB      FIBROCYSTIC
## CLINICAL INDICATION:  4 mth f/u as per rad (mass for bx not seen). BRCA1 carrier. Benign left surgical biopsy. Benign right upper outer MRI guided biopsy.  id     T1w             T2w
## 1 27 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 27        27          0135           5083620 2010-07-19 00:00:00.000000                 2269      4     nonmassB      FIBROCYSTIC
## BRCA 1. Family history breast cancer. Prior benign surgical biopsy left 1999. Post menopausal (hysterectomy May
## 2009).  id     T1w             T2w
## 1 29 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 29        29          0168           5240535 2010-04-12 00:00:00.000000                 2798      4        massB      FIBROCYSTIC
## CLINICAL INDICATION: Intraductal papillary lesion on recent
## breast ultrasound 11 o'clock subareolar position left breast,
## recommended for surgical excision biopsy. Recommended for MRI
## evaluation. Bilateral reduction mammoplasty 1990.  id     T1w             T2w
## 1 30 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 30        30          0171           4751079 2009-02-22 00:00:00.000000                 2272      4        massM     InsituDuctal
## High risk screening detected linear enhancement left subareolar  id     T1w             T2w
## 1 31 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 31        31          0172           4703102 2008-08-24 00:00:00.000000                 2799      4        massB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 32 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 32        32          0173           5123923 2009-11-30 00:00:00.000000                 2800      4     nonmassB   DUCT PAPILLOMA
## None  id     T1w             T2w
## 1 33 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 33        33          0177           6996979 2011-12-09 00:00:00.000000                 2828      3        massM     InsituDuctal
## None  id     T1w             T2w
## 1 34 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 34        34          0177           6996979 2011-12-09 00:00:00.000000                 2828      3     nonmassM     InsituDuctal
## None  id     T1w             T2w
## 1 35 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 35        35          0180           4632561 2008-10-25 00:00:00.000000                 2276      4     nonmassB BENIGN BREAST TISSUE
## BRCA+ve on annual screening MRI  id     T1w             T2w
## 1 36 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 36        36          0180           5254957 2010-12-14 00:00:00.000000                 2275      4     nonmassB     FIBROADENOMA
## BRCA1 positive, previous h/o MR guided biopsy of left lower outer quadrant with benign results.  id     T1w             T2w
## 1 37 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 37        37          0186           6869828 2011-08-02 00:00:00.000000                 2805      4     nonmassM     InsituDuctal
## None  id     T1w             T2w
## 1 38 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 38        38          0189           5057674 2009-10-10 00:00:00.000000                 2832      4     nonmassB SCLEROSING ADENOSIS
## None  id     T1w             T2w
## 1 39 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 39        39          0103           6836585 2011-06-26 00:00:00.000000                 2711      5        massM  PHYLLODES TUMOR
## None  id     T1w             T2w
## 1 40 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 40        40          0190           6760690 2011-04-12 00:00:00.000000                 2092      4        massM   InvasiveDuctal
## None  id     T1w             T2w
## 1 41 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 41        41          0190           6760690 2011-04-12 00:00:00.000000                 2092      4     nonmassM   InvasiveDuctal
## None  id     T1w             T2w
## 1 42 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 42        42          0196           5289117 2010-11-29 00:00:00.000000                 1734      4     nonmassB SCLEROSING ADENOSIS
## None  id     T1w             T2w
## 1 43 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label           lesion_diagnosis
## 43        43          0196           5289117 2010-11-29 00:00:00.000000                 1733      4     nonmassB ColumnarAlterationwoAtypia
## None  id     T1w             T2w
## 1 44 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label           lesion_diagnosis
## 44        44          0196           5289117 2010-11-29 00:00:00.000000                 1733      4     nonmassB ColumnarAlterationwoAtypia
## None  id     T1w             T2w
## 1 45 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 45        45          0197           6667696 2011-05-10 00:00:00.000000                 2204      4     nonmassB LobularHyperplasia
## None  id     T1w             T2w
## 1 46 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 46        46          0197           6667696 2011-05-10 00:00:00.000000                 2204      4     nonmassB LobularHyperplasia
## None  id     T1w             T2w
## 1 47 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 47        47          0197           6667696 2011-05-10 00:00:00.000000                 2204      4        massB LobularHyperplasia
## None  id     T1w             T2w
## 1 48 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 48        48          0198           4809893 2009-05-03 00:00:00.000000                 2278      2        massB      FIBROCYSTIC
## 52 years old BRCA 2 positive. Prior surgical
## excision of right breast fibroadenoma in 2002. BSO 2007.  id     T1w             T2w
## 1 49 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 49        49          0198           4809893 2009-05-03 00:00:00.000000                 2278      2        massB     FIBROADENOMA
## 52 years old BRCA 2 positive. Prior surgical
## excision of right breast fibroadenoma in 2002. BSO 2007.  id     T1w             T2w
## 1 50 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 50        50          0199           4362726 2007-05-18 00:00:00.000000                 2835      4        massB ATYPICAL LOBULAR HYPERPLASIA
## None  id     T1w             T2w
## 1 51 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 51        51          0205           5085133 2010-04-17 00:00:00.000000                  435      4        massB  FIBROEPITHELIAL
## None  id     T1w             T2w
## 1 52 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 52        52          0207           4982884 2009-11-09 00:00:00.000000                 2562      4        massB         FIBROSIS
## 6 month follow-up right breast mass (rim
## enhancing). Prior right lumpectomy, axillary dissection and
## radiation in 1992. BRCA +.  id     T1w             T2w
## 1 53 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 53        53          0212           4734525 2008-09-08 00:00:00.000000                 2877      4        massB     FIBROADENOMA
## None  id     T1w             T2w
## 1 54 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 54        54          0220           6715021 2011-02-22 00:00:00.000000                 2057      5        massB       RadialScar
## None  id     T1w             T2w
## 1 55 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 55        55          0220           6715021 2011-02-22 00:00:00.000000                 2057      5        massB     FIBROADENOMA
## None  id     T1w             T2w
## 1 56 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 56        56          0229           6831376 2011-06-15 00:00:00.000000                 2881      5     nonmassB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 57 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 57        57          0232           6671713 2011-01-11 00:00:00.000000                 1759      5     nonmassB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 58 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 58        58          0232           6671713 2011-01-11 00:00:00.000000                 1760      5     nonmassB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 59 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 59        59          0246           7485590 2013-05-14 00:00:00.000000                 3824      4        massB BENIGN BREAST TISSUE
## HISTORY: 55 year old female. Strong family history of breast cancer. BRCA 1/2 negative. Nodularity palpable right breast, 11 o'clock 2 cm from nipple, negative mammography and ultrasound.  id     T1w             T2w
## 1 60 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 60        60          0252           5142106 2009-12-04 00:00:00.000000                  672      4        massB     FIBROADENOMA
## None  id     T1w             T2w
## 1 61 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 61        61          0252           5142106 2009-12-04 00:00:00.000000                 2883      4        massB     FIBROADENOMA
## None  id     T1w             T2w
## 1 62 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 62        62          0252           6700964 2011-07-18 00:00:00.000000                 2885      3     nonmassB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 63 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 63        63          0252           6700964 2011-07-18 00:00:00.000000                 2885      3        massB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 65 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 65        65          0259           7364573 2013-02-02 00:00:00.000000                   NA      2     nonmassB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 66 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 66        66          0266           5254958 2010-07-16 00:00:00.000000                 2886      4     nonmassB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 67 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 67        67          0276           6952525 2011-12-31 00:00:00.000000                 2570      4        massM   InvasiveDuctal
## OBSP High Risk Screen.  id     T1w             T2w
## 1 68 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 68        68          0276           6952525 2011-12-31 00:00:00.000000                 2570      4        massM   InvasiveDuctal
## OBSP High Risk Screen.  id     T1w             T2w
## 1 69 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 69        69          0277           5077098 2009-09-22 00:00:00.000000                 2571      5     nonmassM     InsituDuctal
## 51 years-old female. BRCA2 carrier. Left ductal carcinoma in situ at biopsy. Evaluate extent of disease.  id     T1w             T2w
## 1 70 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 70        70          0280           5091695 2009-12-07 00:00:00.000000                 2298      4        massB BENIGN BREAST TISSUE
## BRCA 2. LMP 8 years ago (TAH-BSO)  id     T1w             T2w
## 1 72 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 72        72          0293           7491268 2013-06-09 00:00:00.000000                   NA      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk.  Annual recall  id     T1w             T2w
## 1 73 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 73        73          0311           6677243 2011-01-10 00:00:00.000000                 1796      4        massB     FIBROADENOMA
## bilateral breast implants  id     T1w             T2w
## 1 74 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 74        74          0325           4696948 2008-12-01 00:00:00.000000                 2317      4        massB      FIBROCYSTIC
## CLINICAL INDICATION: MRI study  id     T1w             T2w
## 1 75 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 75        75          0331           4722659 2009-01-24 00:00:00.000000                 2318      2        massB      FIBROCYSTIC
## High risk screening study. Bilateral
## surgical biopsies in 1996 (right lower outer quadrant and left upper inner quadrant). Prior US showed bilateral cysts.  id     T1w             T2w
## 1 76 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 76        76          0331           7347095 2013-02-16 00:00:00.000000                 4033      4        massB      FIBROCYSTIC
## None  id     T1w             T2w
## 1 77 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 77        77          0331           7347095 2013-02-16 00:00:00.000000                 4269      4        massB capillary hemangioma
## None  id     T1w             T2w
## 1 78 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 78        78          0352           4785776 2009-01-19 00:00:00.000000                 2328      4        massB     FIBROADENOMA
## Screening. BRCA 2.  id     T1w             T2w
## 1 79 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 79        79          0357           5137030 2009-12-15 00:00:00.000000                 2332      4     nonmassB      FIBROCYSTIC
## BRCA2 carrier, right benign excision 2001, Follow up right breast biopsy July 2009 showing small intraductal
## papilloma and fibrocystic changes without atypia, prior left lumpectomy with radiation and axillary dissection 1998  id     T1w             T2w
## 1 80 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 80        80          0376           4609403 2008-04-04 00:00:00.000000                 2340      4        massB BENIGN HAMARTOMA
## High risk screening study. LMP March 28, 2008.  id     T1w             T2w
## 1 81 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 81        81          0388           7395410 2013-02-26 00:00:00.000000                 2974      5        massM   InvasiveDuctal
## LMP:hysterectomy  id     T1w             T2w
## 1 82 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 82        82          0409           5161803 2010-06-08 00:00:00.000000                 2353      4        massB BENIGN BREAST TISSUE
## Follow up probably benign findings from December 2009. Lumpectomies lower inner and upper outer right breast 1985 and 2000, axillary dissection and radiation (2000) for previous carcinoma. BRCA 2 carrier.  id     T1w             T2w
## 1 83 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 83        83          0420           6738142 2011-03-22 00:00:00.000000                 2065      3     nonmassM     InsituDuctal
## None  id     T1w             T2w
## 1 84 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 84        84          0426           7169326 2012-06-29 00:00:00.000000                 2889      4     nonmassB STROMAL FIBROSIS
## None  id     T1w             T2w
## 1 85 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 85        85          0426           7169326 2012-06-29 00:00:00.000000                 2769      4     nonmassB BENIGN BREAST TISSUE
## None  id     T1w             T2w
## 1 86 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 86        86          0426           7169326 2012-06-29 00:00:00.000000                 2889      4     nonmassB STROMAL FIBROSIS
## None  id     T1w             T2w
## 1 87 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 87        87          0442           4936886 2010-03-13 00:00:00.000000                 2370      4        massB BENIGN BREAST TISSUE
## 40 years-old female. Previous breast
## carcinoma. BRCA 2 gene carrier.  id     T1w             T2w
## 1 88 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 88        88          0456           6689214 2011-02-12 00:00:00.000000                 2706      4        massM   InvasiveDuctal
## None  id     T1w             T2w
## 1 89 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 89        89          0462           5466989 2010-12-23 00:00:00.000000                 2704      3     nonmassM   InvasiveDuctal
## BRCA 1 variant. Known multiple fibroadenomas.  id     T1w             T2w
## 1 90 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 90        90          0462           5466989 2010-12-23 00:00:00.000000                 2703      4        massB     FIBROADENOMA
## BRCA 1 variant. Known multiple fibroadenomas.  id     T1w             T2w
## 1 91 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 91        91          0463           7626269 2014-01-14 00:00:00.000000                 4145      4        massB FLORID HYPERPLASIA
## CLINICAL INDICATION:  OBSP high risk  id     T1w             T2w
## 1 92 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 92        92          0465           4885863 2009-06-07 00:00:00.000000                 2378      2     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## 34 years old, prior right lumpectomy for
## ADH. BRCA 1 positive. Stopped lactating Dec/2008. Family history
## of breast cancer. High risk screening MRI. LMP May/26/2009.  id     T1w             T2w
## 1 93 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 93        93          0473           7364625 2013-12-19 00:00:00.000000                 3326      4        massB      FIBROCYSTIC
## CLINICAL INDICATION: OBSP annual high risk screening.  id     T1w             T2w
## 1 95 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 95        95          0503           6697826 2011-03-19 00:00:00.000000                 2068      3        massM   InvasiveDuctal
## CO-REGISTRATION STUDY 
## 
## INDICATION:  Surveillance.  BRCA 1.  Post menopausal.  id     T1w             T2w
## 1 96 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 96        96          0510           7662547 2014-04-22 00:00:00.000000                 4155      4     nonmassB COLUMNAR CELL CHANGES
## CLINICAL INDICATION:  OBSP HR MRI - returning annual  id     T1w             T2w
## 1 97 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 97        97          0510           7662547 2014-04-22 00:00:00.000000                 4155      4     nonmassB COLUMNAR CELL CHANGES
## CLINICAL INDICATION:  OBSP HR MRI - returning annual  id     T1w             T2w
## 1 98 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 98        98          0513           5043867 2010-09-16 00:00:00.000000                 3693      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: High risk screening study. BRCA 2 carrier.
## LMP 2008. On HRT.  id     T1w             T2w
## 1 99 VIBRANT T2 weighted FSE
##    lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label       lesion_diagnosis
## 99        99          0519           4937737 2009-05-14 00:00:00.000000                 2403      4        massB FLAT EPITHELIAL ATYPIA
## High risk screening study. LMP May 6, 2009.   id     T1w             T2w
## 1 100 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 100       100          0536           7786869 2014-05-06 00:00:00.000000                 3668      4        massB     FIBROADENOMA
## CLINICAL INDICATION: 48 year old female for OBSP High Risk screening; BRCA2 carrier.   id     T1w             T2w
## 1 101 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 101       101          0551           4804820 2008-11-06 00:00:00.000000                 2413      4        massB  FIBROEPITHELIAL
## High risk screening. Postmenopausal.   id     T1w             T2w
## 1 102 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 102       102          0552           4663314 2008-05-15 00:00:00.000000                 2414      4        massB         ADENOSIS
## Routine screen. BRCA 2 carrier. Post
## pregnancy/lactation.   id     T1w             T2w
## 1 103 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 103       103          0553           6687000 2011-07-28 00:00:00.000000                 3831      2        massB BENIGN BREAST TISSUE
## HISTORY: high risk screening.  BRCA 1 carrier. Benign ultrasound guided core biopsy 8 mm hypoechoic mass 4 o'clock left breast February 2011   id     T1w             T2w
## 1 104 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 104       104          0561           4668611 2008-05-30 00:00:00.000000                 2668      4        massB     FIBROADENOMA
## High-risk screening.   id     T1w             T2w
## 1 106 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 106       106          0571           4902166 2009-08-30 00:00:00.000000                 2427      4        massM   InvasiveDuctal
## 40 year-old female. BRCA 1 carrier.
## Hysterectomy and bilateral salpingo-oophorectomy in 2006 .   id     T1w             T2w
## 1 107 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 107       107          0571           4902166 2009-08-30 00:00:00.000000                 2423      4     nonmassB   DUCT PAPILLOMA
## 40 year-old female. BRCA 1 carrier.
## Hysterectomy and bilateral salpingo-oophorectomy in 2006 .   id     T1w             T2w
## 1 108 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 108       108          0572           4681582 2008-07-07 00:00:00.000000                 2428      4     nonmassM   InvasiveDuctal
## High risk screening study. TAHBSO 2006   id     T1w             T2w
## 1 109 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 109       109          0573           5142109 2010-11-06 00:00:00.000000                 3666      4     nonmassB COLUMNAR CELL CHANGES
## CLINICAL INDICATION: 65 year old BRCA 2. TAH BSO.   id     T1w             T2w
## 1 110 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 110       110          0576           6905042 2011-10-25 00:00:00.000000                 3832      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk screening   id     T1w             T2w
## 1 111 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 111       111          0578           6765702 2011-04-13 00:00:00.000000                 2020      6        massM   InvasiveDuctal
## CLINICAL INDICATION:  Biopsy proven invasive ductal carcinoma right breast 7 o'clock, 2 cm and the nipple presenting as a palpable finding.  Family history of breast cancer.   id     T1w             T2w
## 1 112 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 112       112          0580           6855384 2011-09-07 00:00:00.000000                 3833      4        massB     FIBROADENOMA
## CLINICAL INDICATION: BRCA 1 mutation carrier.  PHx ovarian CA and melanoma.  Post TAH/BSO 2006.  On ERT since June 2010   id     T1w             T2w
## 1 113 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 113       113          0586           5332925 2010-08-24 00:00:00.000000                  657      4     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 114 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 114       114          0595           7441706 2013-03-30 00:00:00.000000                 3560      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk.  Previous right lumpectomy for DCIS.  BRCA 2 positive.   id     T1w             T2w
## 1 115 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 115       115          0603           4593568 2008-03-06 00:00:00.000000                 2447      4     nonmassB         FIBROSIS
## High risk screening study. LMP 2004. TAH 2004. BSO 2007.   id     T1w             T2w
## 1 116 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 116       116          0606           6781309 2011-09-21 00:00:00.000000                 3414      4        massB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION:  High risk screening study.  39-year-old woman is followed because of a BRCA-1 mutation.  Prophylactic bilateral salpingo-oophorectomy last year and is currently on estrogen 1.3 mg daily.  She has no new breast concerns.   id     T1w             T2w
## 1 117 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 117       117          0608           5094101 2010-06-27 00:00:00.000000                  267      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION: BRCA1 mutation carrier. Reduction
## mammoplasty 1999.   id     T1w             T2w
## 1 118 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 118       118          0613           4681594 2008-07-15 00:00:00.000000                 2452      4     nonmassM     InsituDuctal
## High risk screening study   id     T1w             T2w
## 1 119 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 119       119          0613           4681594 2008-07-15 00:00:00.000000                 2451      3     nonmassB BENIGN BREAST TISSUE
## High risk screening study   id     T1w             T2w
## 1 120 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 120       120          0616           7910718 2014-11-15 00:00:00.000000                 4146      2        massM     InsituDuctal
## None   id     T1w             T2w
## 1 121 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 121       121          0619           7250777 2013-04-09 00:00:00.000000                 2975      5        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 122 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 122       122          0619           7250777 2013-04-09 00:00:00.000000                 2975      5     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 123 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 123       123          0624           4894714 2009-04-26 00:00:00.000000                 2458      5        massB      FIBROCYSTIC
## 44 years old, BRCA 2 positive. Prior left
## MRI guided biopsy with subsequent left lumpectomy for focal radical scar. Follow up. The patient had TAH and BSO. On estrogen since June/2008.   id     T1w             T2w
## 1 124 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 124       124          0635           7092156 2012-04-03 00:00:00.000000                 2699      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 125 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 125       125          0651           4695822 2008-09-07 00:00:00.000000                 2465      4        massB     FIBROADENOMA
## 6 month follow-up of non-mass enhancement left breast   id     T1w             T2w
## 1 126 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 126       126          0657           6980780 2012-02-14 00:00:00.000000                 2702      4        massM   InvasiveDuctal
## None[1] 127
##    id    T1w    T2w
## 1 127 REVIEW REVIEW
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 127       127          0663           4804825 2008-11-07 00:00:00.000000                 2467      4        massB     FIBROADENOMA
## High risk patient. Part of high risk screening study. Positive family history of breast cancer.   id     T1w             T2w
## 1 128 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 128       128          0664           7081071 2012-05-01 00:00:00.000000                 2795      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 129 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 129       129          0666           5088826 2010-01-02 00:00:00.000000                 2697      3        massM     InsituDuctal
## 56 years-old female. 6 months follow-up of
## bilateral probably benign masses. High risk screening   id     T1w             T2w
## 1 130 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 130       130          0667           4864590 2009-01-29 00:00:00.000000                 2750      3        massB     FIBROADENOMA
## High risk screening study. Postpartum June, 2008. Did not nurse. LMP January 23, 2009.   id     T1w             T2w
## 1 131 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 131       131          0667           4864590 2009-01-29 00:00:00.000000                   43      4        massM     InsituDuctal
## High risk screening study. Postpartum June, 2008. Did not nurse. LMP January 23, 2009.   id     T1w             T2w
## 1 132 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 132       132          0667           6993980 2012-01-11 00:00:00.000000                 2542      4     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 133 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 133       133          0668           6989634 2012-02-09 00:00:00.000000                 2696      4     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 134 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 134       134          0672           4899757 2009-03-03 00:00:00.000000                 1264      5     nonmassM     InsituDuctal
## Suspicious enhancement left breast on
## outside MRI   id     T1w             T2w
## 1 135 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 135       135          0673           4585908 2008-04-03 00:00:00.000000                 2473      4        massB      FIBROCYSTIC
## family history of breast cancer   id     T1w             T2w
## 1 136 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 136       136          0679           4994641 2009-06-12 00:00:00.000000                 2673      6        massM   InvasiveDuctal
## Known malignancy at 9 o'clock right breast.
## A questionable lesion at 7 o'clock   id     T1w             T2w
## 1 137 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 137       137          0681           4999374 2009-06-23 00:00:00.000000                  864      3        massB      FIBROCYSTIC
## Persistent LUOQ mammographic asymmetry with
## distortion. Ultrasound shows cysts and elongated duct.   id     T1w             T2w
## 1 138 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 138       138          0682           5050826 2009-09-06 00:00:00.000000                 2972      6     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 139 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 139       139          0683           5226149 2010-03-25 00:00:00.000000                  509      5        massM   InvasiveDuctal
## Palpable mass right breast, sonographically
## suspicious for multicentric carcinoma   id     T1w             T2w
## 1 140 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 140       140          0684           5266209 2010-07-24 00:00:00.000000                  173      4     nonmassM     InsituDuctal
## Life time risk > 25%. Left lumpectomy for ADH in
## 2004. Family history of breast cancer. LMP: Early June 2010.
##    id     T1w             T2w
## 1 141 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 141       141          0685           5456684 2010-12-14 00:00:00.000000                 1750      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 142 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 142       142          0687              1201 2008-06-22 00:00:00.000000                 2655      5        massM   InvasiveDuctal
## Screening breast MRI, family history of
## breast cancer and prior bilateral excisional biopsies for ADH
## thyroidectomy for thyroid cancer and hysterectomy, life time >30%.
##    id     T1w             T2w
## 1 143 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 143       143          0687              1201 2008-06-22 00:00:00.000000                 2655      5     nonmassM   InvasiveDuctal
## Screening breast MRI, family history of
## breast cancer and prior bilateral excisional biopsies for ADH
## thyroidectomy for thyroid cancer and hysterectomy, life time >30%.
##    id     T1w             T2w
## 1 144 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 144       144          0687              1201 2008-06-22 00:00:00.000000                 2655      5        massM   InvasiveDuctal
## Screening breast MRI, family history of
## breast cancer and prior bilateral excisional biopsies for ADH
## thyroidectomy for thyroid cancer and hysterectomy, life time >30%.
##    id     T1w             T2w
## 1 145 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 145       145          0689           5205923 2010-08-10 00:00:00.000000                 1549      2        massB COLUMNAR CELL CHANGES
## INDICATION: 6 month follow up left enhancement. History of left
## LCIS.   id     T1w             T2w
## 1 146 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 146       146          0690           5180451 2010-05-09 00:00:00.000000                 1790      3     nonmassB      FIBROCYSTIC
## BRCA 2. High risk screening. LMP April 28, 2010.
##    id     T1w             T2w
## 1 147 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 147       147          0690           5180451 2010-05-09 00:00:00.000000                 1790      3        massB USUAL DUCTAL HYPERPLASIA
## BRCA 2. High risk screening. LMP April 28, 2010.
##    id     T1w             T2w
## 1 148 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 148       148          0690           6681276 2011-01-19 00:00:00.000000                 1790      4     nonmassB COLUMNAR CELL CHANGES
## 6 month follow-up of right breast lesion.  BRCA2 mutation carrier.   id     T1w             T2w
## 1 149 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 149       149          0691           5178056 2010-01-23 00:00:00.000000                  590      5        massM   InvasiveDuctal
## Right 6 o'clock palpable finding. Additional lesions seen on ultrasound. LMP January 4 2010 (day 19).   id     T1w             T2w
## 1 150 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 150       150          0691           5178056 2010-01-23 00:00:00.000000                  590      5        massM   InvasiveDuctal
## Right 6 o'clock palpable finding. Additional lesions seen on ultrasound. LMP January 4 2010 (day 19).   id     T1w             T2w
## 1 151 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 151       151          0692           5199366 2010-02-11 00:00:00.000000                  591      4        massB     FIBROADENOMA
## Right upper calcifications, for further evaluation.
## Mother and sister with breast cancer. Post menopausal.   id     T1w             T2w
## 1 152 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 152       152          0696           6983274 2011-11-23 00:00:00.000000                 2543      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 153 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 153       153          0700           4660805 2008-05-25 00:00:00.000000                 2474      5        massM   InvasiveDuctal
## Right upper outer breast suspicious mass per
## mammogram and ultrasound with abnormal right axillary nodes. MRI
## for extent of disease.   id     T1w             T2w
## 1 154 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 154       154          0705           4648471 2008-05-05 00:00:00.000000                 2659      5        massM   InvasiveDuctal
## New lump UOQ left breast 1.5 cm spiculated
## mass on mammo with enlarged axillary nodes dense breasts on mammo
## - MRI for extent of disease. LMP April 24 2008, second week of the
## menstrual cycle.   id     T1w             T2w
## 1 155 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 155       155          0707           5184832 2010-01-30 00:00:00.000000                  612      4     nonmassB COMPLEX PAPILLARY LESION
## For further evaluation. Right brown nipple
## discharge, bilateral calcifications and right ultrasound findings
## (please refer to the recent previous imaging work up). LMP
## October 2009.   id     T1w             T2w
## 1 156 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 156       156          0707           5184832 2010-01-30 00:00:00.000000                 1078      4     nonmassB COMPLEX PAPILLARY LESION
## For further evaluation. Right brown nipple
## discharge, bilateral calcifications and right ultrasound findings
## (please refer to the recent previous imaging work up). LMP
## October 2009.   id     T1w             T2w
## 1 157 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 157       157          0710           5282770 2010-05-29 00:00:00.000000                  233      4        massB     FIBROADENOMA
## Persistent intermittent bilateral clear
## discharge. Bilateral breast pain. Nodule right breast 8 o'clock
## shown to be duct ectasia.
##    id     T1w             T2w
## 1 158 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 158       158          0710           5282770 2010-05-29 00:00:00.000000                  234      5        massB   DUCT PAPILLOMA
## Persistent intermittent bilateral clear
## discharge. Bilateral breast pain. Nodule right breast 8 o'clock
## shown to be duct ectasia.
##    id     T1w             T2w
## 1 159 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 159       159          0710           6798490 2011-05-31 00:00:00.000000                 3566      2        massB DUCTAL HYPERPLASIA
## CLINICAL INDICATION:  Excisional biopsy right upper outer quadrant with a diagnosis of intraductal papilloma with ADH, radial scar in the colon or cell changes with a TPN.  2 excisional biopsies left breast, medially complex fibroadenoma and
## laterally benign intraductal papilloma.  Now blood the left nipple discharge   id     T1w             T2w
## 1 160 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 160       160          0713           5150291 2009-12-11 00:00:00.000000                 2649      5        massM   InvasiveDuctal
## Patient with highly suspicious
## microcalcifications left breast   id     T1w             T2w
## 1 161 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 161       161          0713           5150291 2009-12-11 00:00:00.000000                 2649      5     nonmassM   InvasiveDuctal
## Patient with highly suspicious
## microcalcifications left breast   id     T1w             T2w
## 1 162 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 162       162          0714           5324209 2010-07-09 00:00:00.000000                  118      5        massM   InvasiveDuctal
## Probable multifocal left breast carcinoma on
## imaging.
##    id     T1w             T2w
## 1 163 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 163       163          0714           5324209 2010-07-09 00:00:00.000000                  118      5        massM   InvasiveDuctal
## Probable multifocal left breast carcinoma on
## imaging.
##    id     T1w             T2w
## 1 164 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 164       164          0714           5324209 2010-07-09 00:00:00.000000                  118      5     nonmassM     InsituDuctal
## Probable multifocal left breast carcinoma on
## imaging.
##    id     T1w             T2w
## 1 165 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 165       165          0718           4962581 2009-06-05 00:00:00.000000                  967      4        massB      FIBROCYSTIC
## Suspicious right breast lesion, radial scar
## versus carcinoma on mammography, sonographically occult.   id     T1w             T2w
## 1 166 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 166       166          0720           4965525 2009-05-29 00:00:00.000000                 2479      4        massB      FIBROCYSTIC
## For evaluation of sonographic mass and cysts
## found on baseline screening mammogram. Family history of breast
## cancer.   id     T1w             T2w
## 1 167 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 167       167          0720           4965525 2009-05-29 00:00:00.000000                 2480      4     nonmassB     FIBROADENOMA
## For evaluation of sonographic mass and cysts
## found on baseline screening mammogram. Family history of breast
## cancer.   id     T1w             T2w
## 1 168 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 168       168          0721           4961869 2009-05-08 00:00:00.000000                 1182      6        massM   InvasiveDuctal
## 3.5 cm mass within the right upper outer
## quadrant with suspicious anteriorly located micocalcifications and
## prominent axillary lymph nodes. Assess extent of disease.
##    id     T1w             T2w
## 1 169 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 169       169          0722           5366177 2010-08-24 00:00:00.000000                  983      5        massM   InvasiveDuctal
## Right inferior palpable mass.   id     T1w             T2w
## 1 170 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 170       170          0723           4884108 2009-02-12 00:00:00.000000                 2976      6        massM     InsituDuctal
## None   id     T1w             T2w
## 1 171 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 171       171          0724           5141876 2009-12-03 00:00:00.000000                  988      5        massM   InvasiveDuctal
## Suspicious mass in the left breast with skin
## retraction. Called back from consultation for BIRADS 5 lesion in
## the left breast. LMP 1985.
##    id     T1w             T2w
## 1 172 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 172       172          0726           5304228 2010-06-14 00:00:00.000000                 1008      5        massM   InvasiveDuctal
## Locally advanced breast cancer   id     T1w             T2w
## 1 173 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 173       173          0727           4803733 2009-03-17 00:00:00.000000                 1148      4        massM     InsituDuctal
## CLINICAL INDICATION: high risk screening 25% risk. History
## lumpiness superior central left breast, benign core biopsy 2002.
##    id     T1w             T2w
## 1 174 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 174       174          0728           5304244 2010-06-17 00:00:00.000000                  190      6        massB     FIBROADENOMA
## Left breast Ca 1 o'clock. Outside MRI
## images describe three separate lesions around mass and left lower
## outer left breast enhancement significance unknown
##    id     T1w             T2w
## 1 175 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label           lesion_diagnosis
## 175       175          0728           5304244 2010-06-17 00:00:00.000000                  190      4        massB DUCT PAPILLOMA WITH ATYPIA
## Left breast Ca 1 o'clock. Outside MRI
## images describe three separate lesions around mass and left lower
## outer left breast enhancement significance unknown
##    id     T1w             T2w
## 1 176 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 176       176          0728           5304244 2010-06-17 00:00:00.000000                  190      4        massB     FIBROADENOMA
## Left breast Ca 1 o'clock. Outside MRI
## images describe three separate lesions around mass and left lower
## outer left breast enhancement significance unknown
##    id     T1w             T2w
## 1 177 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 177       177          0728           5304244 2010-06-17 00:00:00.000000                  192      6     nonmassM   InvasiveDuctal
## Left breast Ca 1 o'clock. Outside MRI
## images describe three separate lesions around mass and left lower
## outer left breast enhancement significance unknown
##    id     T1w             T2w
## 1 178 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 178       178          0729           4805710 2009-04-20 00:00:00.000000                 1852      4        massB ATYPICAL LOBULAR HYPERPLASIA
## Previous right breast surgical biopsy 1996,
## ALH. Family history of breast cancer. 2 previous biopsies,
## including ultrasound guided an MRI guided for an MRI detected
## right breast lesion (January and April, 2008), both benign. No
## longer on HRT (according to patient history).
##    id     T1w             T2w
## 1 179 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 179       179          0730           5009497 2009-07-02 00:00:00.000000                  907      5        massM   InvasiveDuctal
## Indeterminate right breast mass and left
## breast calcifications on recent routine screening.   id     T1w             T2w
## 1 180 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label        lesion_diagnosis
## 180       180          0731           5265417 2010-05-06 00:00:00.000000                 2641      4        massB DYSTROPHICCALCIFICATION
## Mammographic architectural distortion benign
## on stereotactic guided VAB Post hysterectomy   id     T1w             T2w
## 1 181 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 181       181          0734           4532660 2008-04-26 00:00:00.000000                 2485      4        massB ATYPICAL DUCTAL HYPERPLASIA
## 53 year old high risk patient considering
## breast reduction, mother and sister with breast cancer, 2 year
## follow up for small bilateral enhancing foci. LMP April 26, 2008   id     T1w             T2w
## 1 182 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 182       182          0735           5276000 2010-05-11 00:00:00.000000                  404      5        massM   InvasiveDuctal
## New segmental linear and pleomorphic
## calcifications in the upper outer quadrant of the right breast
## with a possible, associated obscured mass on screening mammograms.
## Patient complains of pain in both breasts, with a lump in the
## right breast and left nipple discharge (yeast infection). LMP
## May   id     T1w             T2w
## 1 183 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 183       183          0736           4963473 2009-05-19 00:00:00.000000                 2486      4        massB BENIGN BREAST TISSUE
## follow up probably benign mass 6 o'clock left
## breast
##    id     T1w             T2w
## 1 184 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 184       184          0737           4559808 2008-02-12 00:00:00.000000                 2487      3        massM     InsituDuctal
## ACRIN 6666 study   id     T1w             T2w
## 1 185 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                 lesion_diagnosis
## 185       185          0740           4842984 2008-12-22 00:00:00.000000                 1153      4     nonmassB SCLEROSING INTRADUCTAL PAPILLOMA
## Lumpectomy upper inner left breast April 2007
## for 5 mm focus DCIS arising within an intraductal papilloma. 0/9
## notes negative. No radiation. Stereotactic core biopsy with clip
## placement for calcifications upper outer left breast July 2007
## negative. New calcifications lower inner left breast but    id     T1w             T2w
## 1 186 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 186       186          0742           5329785 2010-07-23 00:00:00.000000                 2488      4        massB   DUCT PAPILLOMA
## Left sided nipple discharge with 2
## unsuccessful ductograms   id     T1w             T2w
## 1 187 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 187       187          0742           5329785 2010-07-23 00:00:00.000000                 2488      4     nonmassB   DUCT PAPILLOMA
## Left sided nipple discharge with 2
## unsuccessful ductograms   id     T1w             T2w
## 1 188 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 188       188          0743           4827839 2009-06-06 00:00:00.000000                 2489      4        massM   InvasiveDuctal
## Family history of breast cancer.Exteremly
## dense breast and fibrocystic changes. LMP May/28/2009.   id     T1w             T2w
## 1 189 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 189       189          0744           4848278 2008-12-20 00:00:00.000000                 1396      5        massM   InvasiveDuctal
## Left mammogram with suspicious density,
## distortion and calcifications in 12 o'clock position. Post
## menopausal.   id     T1w             T2w
## 1 190 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 190       190          0745           4881779 2009-02-19 00:00:00.000000                 1307      4        massM   InvasiveDuctal
## Left lower inner quadrant mass and
## calcifiations seen mammographically and sonographically. Left 2
## o'clock small mass seen sonographically. Post menopausal.   id     T1w             T2w
## 1 191 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 191       191          0748           4940559 2009-04-23 00:00:00.000000                 1211      6        massM ATYPICAL DUCTAL HYPERPLASIA
## Right breast 1 o'clock lesion found
## incidentally at ultrasound, biopsy shows at least ADH. MRI for
## extent of disease. Post menopausal on hormones.   id     T1w             T2w
## 1 192 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 192       192          0752           4940477 2009-04-13 00:00:00.000000                 1155      4     nonmassB BENIGN BREAST TISSUE
## Further evaluation of right upper outer
## quadrant distortion and calcifications, prior to recommended core
## biopsy.
##    id     T1w             T2w
## 1 193 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 193       193          0752           4940477 2009-04-13 00:00:00.000000                 1156      4     nonmassB      FIBROCYSTIC
## Further evaluation of right upper outer
## quadrant distortion and calcifications, prior to recommended core
## biopsy.
##    id     T1w             T2w
## 1 194 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 194       194          0752           4940477 2009-04-13 00:00:00.000000                 1155      4     nonmassB BENIGN BREAST TISSUE
## Further evaluation of right upper outer
## quadrant distortion and calcifications, prior to recommended core
## biopsy.
##    id     T1w             T2w
## 1 195 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 195       195          0755           5059877 2009-09-05 00:00:00.000000                  723      4     nonmassM     InsituDuctal
## Bilateral reduction mammoplasty. Bilateral
## calcifications and left upper inner thickening. Post menopausal.   id     T1w             T2w
## 1 196 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 196       196          0755           5059877 2009-09-05 00:00:00.000000                  723      4     nonmassB BENIGN BREAST TISSUE
## Bilateral reduction mammoplasty. Bilateral
## calcifications and left upper inner thickening. Post menopausal.   id     T1w             T2w
## 1 197 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 197       197          0757           4779344 2008-10-20 00:00:00.000000                 2635      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## 59yo, ADH on stereobiopsy of left lower outer
## quadrant microcalcification. To assess extent of the disease. LMP
## , 9 yrs ago.   id     T1w             T2w
## 1 198 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 198       198          0760           4750742 2008-09-07 00:00:00.000000                 1521      5     nonmassM     InsituDuctal
## Patient with highly suspicious segmental
## microcalcifications right breast   id     T1w             T2w
## 1 199 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 199       199          0764           5088503 2009-10-29 00:00:00.000000                 1022      5        massM   InvasiveDuctal
## Right 10 and 7 o'clock suspicious nodules. No
## personal or family history of breast cancer. LMP October 12 2009.   id     T1w             T2w
## 1 200 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 200       200          0764           5088503 2009-10-29 00:00:00.000000                 1022      5        massM   InvasiveDuctal
## Right 10 and 7 o'clock suspicious nodules. No
## personal or family history of breast cancer. LMP October 12 2009.   id     T1w             T2w
## 1 201 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 201       201          0765           5094113 2009-10-15 00:00:00.000000                 2491      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## Right upper outer calcifications for further
## evaluation. Left mastectomy 2007 (ILC), right MRI guided biopsy
## 2007. Postmenopausal.   id     T1w             T2w
## 1 202 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label          lesion_diagnosis
## 202       202          0767           5306672 2010-06-24 00:00:00.000000                  144      4        massB ATYPICAL PAPILLARY LESION
## Consultation of outside imaging identified
## asymmetric tissue in the upper outer left breast. No abnormality
## identified on repeat mammogram. Targeted ultrasound of a palpable
## finding in the left breast at 12 o'clock demonstrated a dilated
## duct with an intraductal mass. MRI for further evaluation.   id     T1w             T2w
## 1 203 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 203       203          0771           4680997 2008-06-10 00:00:00.000000                 3456      4        massB   DUCT PAPILLOMA
## Serous left nipple discharge. Ultrasound May
## 2008 demonstrated dilated branching duct system lower outer left
## breast. At 4 o'clock 3 cm from the nipple lobulated 11 x 3 x 8 mm
## hypoechoic mass which was biopsied and showed fragments of large
## duct papilloma with no atypia. For assessment extent of   id     T1w             T2w
## 1 204 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 204       204          0771           4680997 2008-06-10 00:00:00.000000                 3456      4        massB   DUCT PAPILLOMA
## Serous left nipple discharge. Ultrasound May
## 2008 demonstrated dilated branching duct system lower outer left
## breast. At 4 o'clock 3 cm from the nipple lobulated 11 x 3 x 8 mm
## hypoechoic mass which was biopsied and showed fragments of large
## duct papilloma with no atypia. For assessment extent of   id     T1w             T2w
## 1 205 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 205       205          0775           5437780 2010-11-25 00:00:00.000000                 3458      3        massB BENIGN BREAST TISSUE
## Left breast lesion, assess for other abnormalities.
## LMP 4 weeks ago.   id     T1w             T2w
## 1 206 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 206       206          0775           5437780 2010-11-25 00:00:00.000000                 3458      3     nonmassB BENIGN BREAST TISSUE
## Left breast lesion, assess for other abnormalities.
## LMP 4 weeks ago.   id     T1w             T2w
## 1 207 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 207       207          0775           5437780 2010-11-25 00:00:00.000000                 3458      3     nonmassB BENIGN BREAST TISSUE
## Left breast lesion, assess for other abnormalities.
## LMP 4 weeks ago.   id     T1w             T2w
## 1 208 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 208       208          0778           4794199 2008-11-03 00:00:00.000000                 2492      5        massB     FIBROADENOMA
## Left breast lump  MALIGNANT  PHYLLOIDES TUMOUR-Dec 05, 2007
## 
## 54 yo , prior left mastectomy for phyllodes
## tumor, for follow up of right breast focus of enhancement seen on
## an outside MRI in June/08. LMP Nov/01/2008   id     T1w             T2w
## 1 209 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 209       209          0781           4738440 2008-09-18 00:00:00.000000                 2630      5     nonmassM   InvasiveDuctal
## Suspicious microcalcifications medial aspect
## of right breast on mammography. Hysterectomy more than 15 years
## ago.   id     T1w             T2w
## 1 210 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 210       210          0782           4775699 2008-10-04 00:00:00.000000                 2631      5        massM   InvasiveDuctal
## Highly suspicious mass on mammogram   id     T1w             T2w
## 1 211 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 211       211          0783           4758418 2009-03-24 00:00:00.000000                  952      3        massB      FIBROCYSTIC
## follow up enhancing mass lower inner left
## breast first identified August 2008 in investigation progressive
## thickening superior left breast. No ultrasound correlate. Post
## menopausal.   id     T1w             T2w
## 1 212 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 212       212          0758           4796378 2008-10-30 00:00:00.000000                 2637      4        massB      FIBROCYSTIC
## Right breast biopsy August 2008 for
## microcalcifications. Pathology reveals atypical lobular
## hyperplasia. Biopsy site was right upper inner quadrant. LMP Oct
## 16/08   id     T1w             T2w
## 1 213 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 213       213          0758           4796378 2008-10-30 00:00:00.000000                 2637      4        massB      FIBROCYSTIC
## Right breast biopsy August 2008 for
## microcalcifications. Pathology reveals atypical lobular
## hyperplasia. Biopsy site was right upper inner quadrant. LMP Oct
## 16/08   id     T1w             T2w
## 1 214 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 214       214          0776           5352670 2010-08-16 00:00:00.000000                  317      5        massB    AtypicalCells
## 48 year old female with new palpable lump in
## the left retroareolar region. LMP 21 july 2010
##    id     T1w             T2w
## 1 215 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 215       215          0776           5352670 2010-08-16 00:00:00.000000                 1007      5        massM   InvasiveDuctal
## 48 year old female with new palpable lump in
## the left retroareolar region. LMP 21 july 2010
##    id     T1w             T2w
## 1 216 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 216       216          0779           4934249 2009-04-06 00:00:00.000000                 1190      5        massB SCLEROSING ADENOSIS
## 49 years old with right breast skin
## thickening and nipple inversion. Previous mammogram and ultrasound
## from March/2009 demonstrated large right central mass with
## spiculation and abnormal right axillary lymph nodes. MRI for
## extent of the disease.   id     T1w             T2w
## 1 217 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 217       217          0779           4934249 2009-04-06 00:00:00.000000                 1190      5        massM   InvasiveDuctal
## 49 years old with right breast skin
## thickening and nipple inversion. Previous mammogram and ultrasound
## from March/2009 demonstrated large right central mass with
## spiculation and abnormal right axillary lymph nodes. MRI for
## extent of the disease.   id     T1w             T2w
## 1 218 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 218       218          0789           4785741 2008-11-22 00:00:00.000000                 1403      3        massM     InsituDuctal
## 58 yo ,family history of breast cancer.
## Previous right lateral breast biopsies in 1997 and 2004, pathology
## was bening (papilloma). LMP in 1995. screening MRI   id     T1w             T2w
## 1 219 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 219       219          0790           4708057 2009-01-17 00:00:00.000000                 2494      4     nonmassB   DUCT PAPILLOMA
## Baseline study. Family history of breast
## cancer, risk 25%.   id     T1w             T2w
## 1 220 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 220       220          0791           5365218 2010-08-26 00:00:00.000000                 1049      5        massM  InvasiveLobular
## Previous left lumpectomy for LCIS. Suspicion
## of multicentric disease in the right breast. LMP Aug 9/10   id     T1w             T2w
## 1 221 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 221       221          0791           5365218 2010-08-26 00:00:00.000000                 1049      5     nonmassM  InvasiveLobular
## Previous left lumpectomy for LCIS. Suspicion
## of multicentric disease in the right breast. LMP Aug 9/10   id     T1w             T2w
## 1 222 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 222       222          0792           5264066 2010-06-26 00:00:00.000000                  207      3        massB   DUCT PAPILLOMA
## Bilateral breast pain. Family history of
## breast cancer. Radiologist recommended breast MRI.   id     T1w             T2w
## 1 223 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 223       223          0792           5264066 2010-06-26 00:00:00.000000                  207      3        massB   DUCT PAPILLOMA
## Bilateral breast pain. Family history of
## breast cancer. Radiologist recommended breast MRI.   id     T1w             T2w
## 1 224 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 224       224          0793           4988020 2009-08-28 00:00:00.000000                  726      4        massB     FIBROADENOMA
## Screening. Mother had bilateral DCIS at age
## 38. LMP August 19 2009.
##    id     T1w             T2w
## 1 225 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 225       225          0793           7135216 2012-09-11 00:00:00.000000                 1992      2        massB COMPLEX FIBROEPITHELIAL LESION
## None   id     T1w             T2w
## 1 226 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 226       226          0795           5188009 2010-03-12 00:00:00.000000                  562      6     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## Left spontaneous nipple discharge (greenish) 11
## o'clock duct. Strong family history of breast/ovarian cancer.
## LMP February 28 2010.
##    id     T1w             T2w
## 1 227 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 227       227          0796            860773 2009-10-10 00:00:00.000000                  728      4     nonmassB ATYPICAL LOBULAR HYPERPLASIA
## Recent left breast ultrasound guided core
## biopsy, ALH on pathology.   id     T1w             T2w
## 1 228 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 228       228          0799           5372294 2010-09-11 00:00:00.000000                  354      4        massB      FIBROCYSTIC
## Left upper outer calcifications, with core biopsy of
## atypia (performed elsewhere). For evaluation of adjacent
## calcifications and left 4 o'clock sonographic mass. LMP: 9/4/10.   id     T1w             T2w
## 1 229 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 229       229          0802           4600874 2008-03-09 00:00:00.000000                 2496      4        massM   InvasiveDuctal
## Mammographic and sonographic medial right
## breast nodule (probably benign). Additional left lateral breast
## asymmetry with no sonographic correlate. MRI for problem solving.   id     T1w             T2w
## 1 230 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 230       230          0803           5058195 2009-08-31 00:00:00.000000                  320      5        massM   InvasiveDuctal
## Suspicious mass left breast seen on
## ultrasound   id     T1w             T2w
## 1 231 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 231       231          0805           5059167 2009-09-13 00:00:00.000000                  754      4     nonmassB      FIBROCYSTIC
## For further evaluation of right breast
## mammographic calcifications, radiologist recommended. Family
## history of breast cancer.   id     T1w             T2w
## 1 232 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 232       232          0807           5235491 2010-04-01 00:00:00.000000                  522      5     nonmassM     InsituDuctal
## Suspicious left breast calcifications
## recommended for stereotactic biopsy. For extent of disease.   id     T1w             T2w
## 1 233 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 233       233          0809           5016014 2009-07-13 00:00:00.000000                  827      4        massB      FIBROCYSTIC
## Suspicious calcifications in left breast.
## Post menopausal.
##    id     T1w             T2w
## 1 234 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 234       234          0810           4622489 2008-06-26 00:00:00.000000                 2615      4        massM     InsituDuctal
## 6 month follow up   id     T1w             T2w
## 1 235 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 235       235          0812           4700538 2008-06-30 00:00:00.000000                 2617      5        massM   InvasiveDuctal
## Locally advanced breast cancer. Assess
## extent of disease   id     T1w             T2w
## 1 236 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 236       236          0813           5378164 2010-09-14 00:00:00.000000                 1578      5     nonmassB      FIBROCYSTIC
## bloody left nipple discharge. Mammograms and ultrasound
## findings suspicious for malignancy medial left breast 8 o'clock 5
## cm from nipple with suspicious low axillary lymph node. Other
## sonographic findings closer to nipple including 5-6 o'clock may
## indicate DCIS. Unsuccessful ductogram.   id     T1w             T2w
## 1 237 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 237       237          0813           5378164 2010-09-14 00:00:00.000000                 1578      5        massM  InvasiveLobular
## bloody left nipple discharge. Mammograms and ultrasound
## findings suspicious for malignancy medial left breast 8 o'clock 5
## cm from nipple with suspicious low axillary lymph node. Other
## sonographic findings closer to nipple including 5-6 o'clock may
## indicate DCIS. Unsuccessful ductogram.   id     T1w             T2w
## 1 238 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 238       238          0814           4704240 2008-08-07 00:00:00.000000                 2618      5        massM   InvasiveDuctal
## Left upper inner mass with calcifications,
## for extent of disease.   id     T1w             T2w
## 1 239 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 239       239          0814           6667547 2011-01-06 00:00:00.000000                 3294      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION:  Recent biopsy proven invasive ductal carcinoma right axillary tail.  Left breast conserving therapy October, 2008.   id     T1w             T2w
## 1 240 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 240       240          0815           4828432 2008-12-04 00:00:00.000000                 1433      5        massM   InvasiveDuctal
## Highly suspicious lesion visualized on
## mammography and ultrasound. Assess for extent of disease. Last
## menstrual period was on December 2, 2008.   id     T1w             T2w
## 1 241 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 241       241          0817           5363917 2010-08-24 00:00:00.000000                 2499      6        massM   InvasiveDuctal
## Investigation in Bangladesh of palpable abnormality 1
## o'clock right breast with FNAB positive for malignancy. FNAB of
## prominent node negative for malignancy.   id     T1w             T2w
## 1 242 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 242       242          0818           5021762 2009-07-28 00:00:00.000000                  884      4        massB   DUCT PAPILLOMA
## Known papillary lesion on fine needle
## aspiration done on a palpable lump at 530. Patient with bloody
## nipple discharge   id     T1w             T2w
## 1 243 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 243       243          0827           4985128 2009-06-13 00:00:00.000000                 2620      4        massM INSITUPAPILLARYCARCINOMA
## 72 years old, increasing right lateral
## breast asymmetry suspicious for malignancy. Had FNA at an outside
## institution, pathology is suggestive of angiolipoma. Questionable
## angiosarcoma. Menopausal.   id     T1w             T2w
## 1 244 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 244       244          0827           4985128 2009-06-13 00:00:00.000000                 2620      4     nonmassM INSITUPAPILLARYCARCINOMA
## 72 years old, increasing right lateral
## breast asymmetry suspicious for malignancy. Had FNA at an outside
## institution, pathology is suggestive of angiolipoma. Questionable
## angiosarcoma. Menopausal.   id     T1w             T2w
## 1 245 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 245       245          0827           4985128 2009-06-13 00:00:00.000000                  924      4        massB USUAL DUCTAL HYPERPLASIA
## 72 years old, increasing right lateral
## breast asymmetry suspicious for malignancy. Had FNA at an outside
## institution, pathology is suggestive of angiolipoma. Questionable
## angiosarcoma. Menopausal.   id     T1w             T2w
## 1 246 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 246       246          0828           4787730 2008-10-13 00:00:00.000000                 1407      6        massM     InsituDuctal
## Left breast calcifications and sonographic
## lesion. Outside biopsy of ADH.   id     T1w             T2w
## 1 247 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 247       247          0829           5264139 2010-05-02 00:00:00.000000                  384      5        massM   InvasiveDuctal
## Incidental mass in right subareolar region
## on MRI from St. Michael's hospital. No mammographic abnormality.
## Post menopausal.
##    id     T1w             T2w
## 1 248 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 248       248          0830           4863868 2009-01-13 00:00:00.000000                 2793      5        massB ATYPICAL DUCTAL HYPERPLASIA
## Suspicious masses in the right breast on
## Mammograms and ultrasound, MRI to determine extent of disease.
## Targeted left breast lower inner quadrant was normal   id     T1w             T2w
## 1 249 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 249       249          0830           4863868 2009-01-13 00:00:00.000000                 2794      5        massB             Cyst
## Suspicious masses in the right breast on
## Mammograms and ultrasound, MRI to determine extent of disease.
## Targeted left breast lower inner quadrant was normal   id     T1w             T2w
## 1 250 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 250       250          0831           4633368 2008-04-13 00:00:00.000000                 2500      6     nonmassM     InsituDuctal
## Known left DCIS. Assess extent of disease.   id     T1w             T2w
## 1 251 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 251       251          0834           4614262 2008-03-27 00:00:00.000000                 2621      5        massM   InvasiveDuctal
## none seen on report   id     T1w             T2w
## 1 252 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 252       252          0837           4559849 2008-02-28 00:00:00.000000                 1839      5        massM   InvasiveDuctal
## High risk screening study. LMP February 11,
## 2008 (week 3 of menstrual cycle).
## 
## DICOM images saved as pt id837   id     T1w             T2w
## 1 253 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 253       253          0839           4739257 2009-03-07 00:00:00.000000                  954      4        massB BENIGN BREAST TISSUE
## Family history of breast cancer .   id     T1w             T2w
## 1 254 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 254       254          0839           4739257 2009-03-07 00:00:00.000000                  954      4        massB BENIGN BREAST TISSUE
## Family history of breast cancer .   id     T1w             T2w
## 1 255 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 255       255          0843           4798594 2009-04-19 00:00:00.000000                 1198      4        massB      FIBROCYSTIC
## 41 years old, family history of breast
## cancer. High risk screening MRI. Life time risk is 26%. LMP
## April/01/2009.   id     T1w             T2w
## 1 256 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 256       256          0843           6792402 2011-05-10 00:00:00.000000                 2607      4        massB      FIBROCYSTIC
## 43 year old with family history of breast cancer.   id     T1w             T2w
## 1 257 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 257       257          0844           4795902 2008-10-27 00:00:00.000000                 2608      4     nonmassM   InvasiveDuctal
## 48 yo with history of left axially mass in
## Aug/2008 of which was biopsied and positive for malignancy from
## questionable primary breast cancer. An outside MG showed two
## cluster of microcalcification at 5-6 Oclock and outside MRI showed
## area of enhancement at 4-5 and 6 oclock position. LMP
## Oct/26   id     T1w             T2w
## 1 258 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 258       258          0845           5433683 2010-11-12 00:00:00.000000                 2502      5        massM  InvasiveLobular
## Mass right breast suspicious of malignancy
## Postmenopausal
##    id     T1w             T2w
## 1 259 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 259       259          0846           4800867 2008-11-10 00:00:00.000000                 2503      5        massM MetaplasticCarcinoma
## 62 yo with history of right breast
## suspicious mass since o8/08   id     T1w             T2w
## 1 260 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 260       260          0847           5064132 2009-09-24 00:00:00.000000                 2610      4        massM   InvasiveDuctal
## Known fibrocystic breasts. Right breast
## enlarging subareolar cyst with mural nodules,? Nodules due to
## inspissated debris versus papillary solid nodules. Recommended for
## MRI assessment. Reduction mammoplasty 2004.   id     T1w             T2w
## 1 261 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 261       261          0850           5380609 2010-09-23 00:00:00.000000                 2504      5        massB BENIGN BREAST TISSUE
## New focal asymmetry in upper outer right breast on
## mammogram (July, 2010), not seen on ultrasound. Stereotactic
## biopsy demonstrated benign breast tissue. This MRI study is for
## further assessment.
##    id     T1w             T2w
## 1 262 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 262       262          0850           5380609 2010-09-23 00:00:00.000000                 2504      5     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## New focal asymmetry in upper outer right breast on
## mammogram (July, 2010), not seen on ultrasound. Stereotactic
## biopsy demonstrated benign breast tissue. This MRI study is for
## further assessment.
##    id     T1w             T2w
## 1 263 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 263       263          0850           5380609 2010-09-23 00:00:00.000000                 2504      5        massB      FIBROCYSTIC
## New focal asymmetry in upper outer right breast on
## mammogram (July, 2010), not seen on ultrasound. Stereotactic
## biopsy demonstrated benign breast tissue. This MRI study is for
## further assessment.
##    id     T1w             T2w
## 1 264 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 264       264          0851           4593282 2008-04-10 00:00:00.000000                  573      4        massB    InsituLobular
## Recurrent mastitis right breast.
## Intraductal echogenic filling defect on the right noted on recent
## ultrasound, scheduled for excision. Exclusion of additional
## pathology.   id     T1w             T2w
## 1 265 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 265       265          0851           4593282 2008-04-10 00:00:00.000000                  573      4        massM   InvasiveDuctal
## Recurrent mastitis right breast.
## Intraductal echogenic filling defect on the right noted on recent
## ultrasound, scheduled for excision. Exclusion of additional
## pathology.   id     T1w             T2w
## 1 266 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 266       266          0853           4798586 2009-03-27 00:00:00.000000                 2507      2     nonmassB      FIBROCYSTIC
## Follow-up Rt non-mass enhancement. Post left
## surgical biopsy for papillomas. Post hysterectomy   id     T1w             T2w
## 1 267 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 267       267          0853           4745782 2008-09-08 00:00:00.000000                 2044      3     nonmassB      FIBROCYSTIC
## Patient with left nipple discharge   id     T1w             T2w
## 1 268 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 268       268          0853           6696534 2011-03-16 00:00:00.000000                 2044      4     nonmassM     InsituDuctal
## 1 year f/u of non mass enhancement right breast.  Significant family history of breast CA.  Previous contralateral excisional biopsies of papillomas 2008, 2009.   id     T1w             T2w
## 1 269 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 269       269          0855           4641315 2008-04-29 00:00:00.000000                 2613      6        massB         FIBROSIS
## October 2007 lumpectomy for flat epithelial atypia, no ADH. Positive family hx breast Ca. New nodularity
## adjacent to surgical scar upper outer quadrant right breast, positive for ADH under ultrasound core bx   id     T1w             T2w
## 1 270 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 270       270          0855           4641315 2008-04-29 00:00:00.000000                 2613      6     nonmassB         FIBROSIS
## October 2007 lumpectomy for flat epithelial atypia, no ADH. Positive family hx breast Ca. New nodularity
## adjacent to surgical scar upper outer quadrant right breast, positive for ADH under ultrasound core bx   id     T1w             T2w
## 1 271 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 271       271          0856           4986174 2009-10-18 00:00:00.000000                 2508      4        massB     FIBROADENOMA
## 39 year old. LMP October 10, 2009. Very
## strong family history of breast and ovarian cancer. Twin sister
## is BRCA1 variant. No mammographic evidence of malignancy on
## consult of outside imaging.
##    id     T1w             T2w
## 1 272 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 272       272          0856           4986174 2009-10-18 00:00:00.000000                 2508      4        massB     FIBROADENOMA
## 39 year old. LMP October 10, 2009. Very
## strong family history of breast and ovarian cancer. Twin sister
## is BRCA1 variant. No mammographic evidence of malignancy on
## consult of outside imaging.
##    id     T1w             T2w
## 1 273 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 273       273          0856           6871177 2012-01-04 00:00:00.000000                 2572      2        massB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 274 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 274       274          0857           4870283 2009-04-16 00:00:00.000000                 2509      4        massB      FIBROCYSTIC
## Screening. Right mastectomy. Left
## calcifications with a benign biopsy. Post menopausal.   id     T1w             T2w
## 1 275 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 275       275          0857           5013393 2009-12-01 00:00:00.000000                 1775      2        massM     InsituDuctal
## Follow-up post benign biopsy   id     T1w             T2w
## 1 276 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 276       276          0861           5053396 2009-08-24 00:00:00.000000                  858      5        massM   InvasiveDuctal
## 49 years-old female. High breast density.
## Probable left breast carcinoma detected on recent imaging   id     T1w             T2w
## 1 277 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 277       277          0862           5395314 2010-10-04 00:00:00.000000                 2596      4        massM     InsituDuctal
## Recent stereotactic biopsy of left breast
## upper outer quadrant microcalcifications, atypical ductal
## hyperplasia on pathology.   id     T1w             T2w
## 1 278 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 278       278          0863           4969136 2009-06-02 00:00:00.000000                  929      4        massB   DUCT PAPILLOMA
## For problem solving. Ultrasound showed
## bilateral breast masses. Family history (2 sisters with breast
## cancer in their 20's). Hysterectomy 30 years ago.   id     T1w             T2w
## 1 279 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 279       279          0863           4969136 2009-06-02 00:00:00.000000                  928      4        massM  InvasiveLobular
## For problem solving. Ultrasound showed
## bilateral breast masses. Family history (2 sisters with breast
## cancer in their 20's). Hysterectomy 30 years ago.   id     T1w             T2w
## 1 280 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 280       280          0865           5267535 2010-05-03 00:00:00.000000                  422      5        massM   InvasiveDuctal
## Palpable mass in the right breast. LMP
## April 13th.   id     T1w             T2w
## 1 281 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 281       281          0865           5267535 2010-05-03 00:00:00.000000                 2511      5     nonmassM   InvasiveDuctal
## Palpable mass in the right breast. LMP
## April 13th.   id     T1w             T2w
## 1 282 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 282       282          0867           5372277 2010-09-07 00:00:00.000000                  331      5     nonmassM   InvasiveDuctal
## Palpable right breast mass and axillary adenopathy   id     T1w             T2w
## 1 283 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 283       283          0870           5141888 2009-12-03 00:00:00.000000                  679      6     nonmassM  InvasiveLobular
## 68 year old with node positive breast cancer
## without obvious breast primary lesion. Has had prior left
## axillary lymph node dissection after presenting post fall with a
## left axillary mass.   id     T1w             T2w
## 1 284 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 284       284          0871           5130094 2009-11-28 00:00:00.000000                 2597      4        massM  InvasiveLobular
## History of bilateral lumpectomies (11
## previous excisional biopsies for benign disease). Recent core
## biopsy on the left demonstrated ADH. The patient is at high risk
## for breast cancer and is considering prophylactic mastectomies
## (mother (age 38), maternal aunts).   id     T1w             T2w
## 1 285 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 285       285          0871           5130094 2009-11-28 00:00:00.000000                  808      4        massM     InsituDuctal
## History of bilateral lumpectomies (11
## previous excisional biopsies for benign disease). Recent core
## biopsy on the left demonstrated ADH. The patient is at high risk
## for breast cancer and is considering prophylactic mastectomies
## (mother (age 38), maternal aunts).   id     T1w             T2w
## 1 286 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 286       286          0871           5130094 2009-11-28 00:00:00.000000                 2597      4        massM  InvasiveLobular
## History of bilateral lumpectomies (11
## previous excisional biopsies for benign disease). Recent core
## biopsy on the left demonstrated ADH. The patient is at high risk
## for breast cancer and is considering prophylactic mastectomies
## (mother (age 38), maternal aunts).   id     T1w             T2w
## 1 287 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 287       287          0873           4956191 2009-05-09 00:00:00.000000                 2598      4        massB      FIBROCYSTIC
## 48 years old, previous biopsy showed ADH in
## left breast 3 O'clock in one out of five cores, to evaluate extent
## of disease. LMP Feb/2009.   id     T1w             T2w
## 1 288 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 288       288          0873           4956191 2009-05-09 00:00:00.000000                 2598      4        massB      FIBROCYSTIC
## 48 years old, previous biopsy showed ADH in
## left breast 3 O'clock in one out of five cores, to evaluate extent
## of disease. LMP Feb/2009.   id     T1w             T2w
## 1 289 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 289       289          0875           7141879 2012-05-17 00:00:00.000000                 3746      4        massB   DUCT PAPILLOMA
## CLINICAL INDICATION: Bilateral nipple discharge with prior excision of papillomas bilaterally  LMP 2003   id     T1w             T2w
## 1 290 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 290       290          0875           5396107 2010-10-14 00:00:00.000000                 1651      4     nonmassB PAPILLARY LESION
## None   id     T1w             T2w
## 1 291 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 291       291          0876           4719378 2008-09-18 00:00:00.000000                 1487      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 292 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 292       292          0877           4724338 2008-08-18 00:00:00.000000                 1546      4     nonmassB ATYPICAL LOBULAR HYPERPLASIA
## High Risk screening. Baseline MRI   id     T1w             T2w
## 1 293 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label       lesion_diagnosis
## 293       293          0880           4809515 2009-03-21 00:00:00.000000                 2514      4        massB Papillary(focalAtypia)
## Right blood nipple discharge with resultant resection
## of a right nipple adenoma with atypical ductal hyperplasia
## incompletely excised. Rule out residual disease. History of
## previous right benign surgical biopsy and bilateral
## sonographically seen breast masses. LMP March 9 2009.   id     T1w             T2w
## 1 294 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 294       294          0880           6778829 2011-05-05 00:00:00.000000                 2599      3        massM     InsituDuctal
## 31 yo female. Bilateral stable masses.  Papilloma with atypia in the right breast excised in 2009. Treated for adenoma and ADH in 2000.   id     T1w             T2w
## 1 295 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 295       295          0880           6778829 2011-05-05 00:00:00.000000                 2599      3     nonmassM     InsituDuctal
## 31 yo female. Bilateral stable masses.  Papilloma with atypia in the right breast excised in 2009. Treated for adenoma and ADH in 2000.   id     T1w             T2w
## 1 296 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 296       296          0880           6778829 2011-05-05 00:00:00.000000                 2599      3     nonmassM     InsituDuctal
## 31 yo female. Bilateral stable masses.  Papilloma with atypia in the right breast excised in 2009. Treated for adenoma and ADH in 2000.   id     T1w             T2w
## 1 297 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label  lesion_diagnosis
## 297       297          0880           4809515 2009-03-21 00:00:00.000000                 2514      4        massB AtypicalPapilloma
## Right blood nipple discharge with resultant resection
## of a right nipple adenoma with atypical ductal hyperplasia
## incompletely excised. Rule out residual disease. History of
## previous right benign surgical biopsy and bilateral
## sonographically seen breast masses. LMP March 9 2009.   id     T1w             T2w
## 1 298 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 298       298          0883           5177385 2010-01-22 00:00:00.000000                  605      5     nonmassM     InsituDuctal
## Left calcifications, for biopsy. Extent of disease.
## Postmenopausal.
##    id     T1w             T2w
## 1 299 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 299       299          0884           6876318 2011-08-05 00:00:00.000000                 2595      6        massM   InvasiveDuctal
## Kmown malignancy RUOQ - for extent of disease \T\ nodal assessment Post menopausal   id     T1w             T2w
## 1 300 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 300       300          0884           6876318 2011-08-05 00:00:00.000000                 2595      6     nonmassM     InsituDuctal
## Kmown malignancy RUOQ - for extent of disease \T\ nodal assessment Post menopausal   id     T1w             T2w
## 1 301 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 301       301          0885           6747175 2011-03-31 00:00:00.000000                 2073      4     nonmassB      FIBROCYSTIC
## Right bloody nipple discharge.  Prior right surgical biopsy 1998 with ADH.  Prior right benign biopsy of subareolar calcifications with clip placement 2005.  Post menopausal on HRT.
##    id     T1w             T2w
## 1 302 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 302       302          0887           6794529 2011-05-13 00:00:00.000000                 2675      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 303 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 303       303          0888           6744887 2011-03-23 00:00:00.000000                 2676      5        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 304 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 304       304          0896           6895982 2012-03-02 00:00:00.000000                 2683      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 305 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                  lesion_diagnosis
## 305       305          0898           5224531 2010-11-07 00:00:00.000000                 2687      4        massB DUCTAL HYPERPLASIA WITHOUT ATYPIA
## None   id     T1w             T2w
## 1 306 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 306       306          0900           6699226 2011-02-04 00:00:00.000000                 2690      4        massB             Cyst
## None   id     T1w             T2w
## 1 307 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 307       307          0904           7133915 2012-05-12 00:00:00.000000                 2695      3        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 308 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 308       308          0913           7350757 2013-01-13 00:00:00.000000                 2982      4        massB         ADENOSIS
## None   id     T1w             T2w
## 1 309 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 309       309          0918           6976567 2011-11-19 00:00:00.000000                 2903      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 310 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 310       310          0920           7095635 2012-06-16 00:00:00.000000                 2904      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 311 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 311       311          0921           6997232 2012-05-01 00:00:00.000000                 2908      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 312 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 312       312          0924           7532614 2013-06-26 00:00:00.000000                 3542      4        massB         ADENOSIS
## INDICATION:  Surveillance (OBSP).  Prior benign biopsy in 2009 (ultrasound and MRI finding, performed under ultrasound guidance).  Family history of breast cancer.   id     T1w             T2w
## 1 313 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 313       313          0924           7532614 2013-06-26 00:00:00.000000                 3543      4        massB         ADENOSIS
## INDICATION:  Surveillance (OBSP).  Prior benign biopsy in 2009 (ultrasound and MRI finding, performed under ultrasound guidance).  Family history of breast cancer.   id     T1w             T2w
## 1 314 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 314       314          0934           5314924 2010-07-24 00:00:00.000000                  139      4     nonmassB     FIBROADENOMA
## None   id     T1w             T2w
## 1 315 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 315       315          0937           7144673 2012-09-04 00:00:00.000000                 2919      4        massB  FIBROEPITHELIAL
## None   id     T1w             T2w
## 1 316 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 316       316          0943           5395204 2010-09-25 00:00:00.000000                 1572      4     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 317 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 317       317          0943           5395204 2010-09-25 00:00:00.000000                 1573      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 318 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 318       318          0944           7742881 2014-05-23 00:00:00.000000                 3615      4        massM     InsituDuctal
## CLINICAL INDICATION: High risk screening.  Benign right MRI guided biopsy  LMP May 16/14   id     T1w             T2w
## 1 319 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 319       319          0944           7092128 2012-04-17 00:00:00.000000                 2921      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 320 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 320       320          0950           6931716 2011-10-04 00:00:00.000000                 2922      5        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 321 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 321       321          0950           6931716 2011-10-04 00:00:00.000000                 2922      5     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 322 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 322       322          0952           7105222 2012-05-06 00:00:00.000000                 2923      4     nonmassB FOCAL USUAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 323 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 323       323          0952           7105222 2012-05-06 00:00:00.000000                 2924      4     nonmassB     FIBROADENOMA
## None   id     T1w             T2w
## 1 324 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 324       324          0954           7962026 2014-11-11 00:00:00.000000                 4144      4        massB ATYPICAL LOBULAR HYPERPLASIA
## None   id     T1w             T2w
## 1 325 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 325       325          0956           5062341 2010-02-16 00:00:00.000000                  559      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 326 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 326       326          0962           4755483 2009-03-09 00:00:00.000000                 1868      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 327 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 327       327          0965           6676125 2011-01-18 00:00:00.000000                 1892      3        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 328 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 328       328          0967           6938015 2012-03-29 00:00:00.000000                 2986      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION: Follow-up right breast mass  LMP Mar 17/12   id     T1w             T2w
## 1 329 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 329       329          0967           6938015 2012-03-29 00:00:00.000000                 2986      4     nonmassM   InvasiveDuctal
## CLINICAL INDICATION: Follow-up right breast mass  LMP Mar 17/12   id     T1w             T2w
## 1 330 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 330       330          0978           4851428 2009-06-15 00:00:00.000000                 2989      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 331 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 331       331          0985           7050619 2012-04-25 00:00:00.000000                 2994      4     nonmassB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 332 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 332       332          0993           6979299 2012-01-27 00:00:00.000000                 2997      4        massM  PHYLLODES TUMOR
## None   id     T1w             T2w
## 1 333 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 333       333          0995           6816787 2011-08-09 00:00:00.000000                 2998      4     nonmassB LARGE DUCT PAPILLOMA
## None   id     T1w             T2w
## 1 334 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 334       334          0995           6816787 2011-08-09 00:00:00.000000                 2999      3        massB   DUCT PAPILLOMA
## None   id     T1w             T2w
## 1 335 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 335       335          0997           7279207 2012-11-13 00:00:00.000000                 4173      3        massB ATYPICAL LOBULAR HYPERPLASIA
## None   id     T1w             T2w
## 1 336 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 336       336          0999           6925971 2011-10-29 00:00:00.000000                 3000      3        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 337 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 337       337          1003           6682777 2011-01-24 00:00:00.000000                 3005      4     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 338 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 338       338          1004           6801264 2011-05-29 00:00:00.000000                 3006      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 339 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 339       339          1006           4443563 2007-09-06 00:00:00.000000                 3010      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 340 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 340       340          1008           6745959 2011-03-24 00:00:00.000000                 2006      5        massM     InsituDuctal
## None   id     T1w             T2w
## 1 341 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 341       341          1012           7629993 2013-10-06 00:00:00.000000                 4372      6        massM   InvasiveDuctal
## CLINICAL INDICATION:  49-year-old woman with new right breast IDC, extent of disease.   id     T1w             T2w
## 1 342 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 342       342          1012           6940724 2011-10-22 00:00:00.000000                 3013      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 343 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 343       343          1012           6940724 2011-10-22 00:00:00.000000                 3013      4     nonmassB     FIBROADENOMA
## None   id     T1w             T2w
## 1 344 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 344       344          1018           4773924 2009-09-17 00:00:00.000000                 3016      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 345 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 345       345          1021           6760795 2011-08-26 00:00:00.000000                 3017      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 346 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 346       346          1024           6980462 2011-11-20 00:00:00.000000                 3020      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 347 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label           lesion_diagnosis
## 347       347          1024           6980462 2011-11-20 00:00:00.000000                 3020      4     nonmassB DUCT PAPILLOMA WITH ATYPIA
## None   id     T1w             T2w
## 1 348 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 348       348          1025           6703528 2011-02-11 00:00:00.000000                 3022      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 349 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 349       349          1027           6930730 2011-10-05 00:00:00.000000                 3025      3     nonmassB     FIBROADENOMA
## None   id     T1w             T2w
## 1 350 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 350       350          1044           7366817 2013-03-03 00:00:00.000000                 4323      5        massM   InvasiveDuctal
## CLINICAL INDICATION:  OBSP high risk.  First study   id     T1w             T2w
## 1 351 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 351       351          1045           7231265 2012-08-18 00:00:00.000000                 4324      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  New focal asymmetry in the right upper outer breast on screening mammogram that did not resolve with additional spot compression views.  No correlating abnormality identified on ultrasound.  Family history of breast cancer in mother at age 59 in sister at age 50.   id     T1w             T2w
## 1 352 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 352       352          1050           7296806 2012-12-04 00:00:00.000000                 4329      3     nonmassB   DENSE FIBROSIS
## CLINICAL INDICATION: OBSP high risk. Recommended 6 month follow up.   id     T1w             T2w
## 1 353 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                           lesion_diagnosis
## 353       353          1026           6907382 2011-10-26 00:00:00.000000                 3024      4        massB DENSE FIBROSIS AND FIBROADENOMATOID CHANGE
## None   id     T1w             T2w
## 1 354 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 354       354          1053           7748055 2014-03-09 00:00:00.000000                 4330      4        massB INFLAMED CYST WALL
## CLINICAL INDICATION:  OBSP high risk screening   id     T1w             T2w
## 1 355 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 355       355          1062           7408296 2013-02-19 00:00:00.000000                 4336      4        massB ATYPICAL LOBULAR HYPERPLASIA
## CLINICAL INDICATION:  Recent left lumpectomy.  Assess for residual disease   id     T1w             T2w
## 1 356 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 356       356          1062           7408296 2013-02-19 00:00:00.000000                 4335      4     nonmassB      FIBROCYSTIC
## CLINICAL INDICATION:  Recent left lumpectomy.  Assess for residual disease   id     T1w             T2w
## 1 358 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 358       358          1065           7741665 2014-03-19 00:00:00.000000                 4339      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk MRI   id     T1w             T2w
## 1 359 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 359       359          1071           7382882 2013-01-18 00:00:00.000000                 4340      4     nonmassM     InsituDuctal
## CLINICAL INDICATION:  Calcifications in left breast.  Atypia, suspicious for DCIS on previous biopsy.  Very dense breasts on mammogram.
## 
## LMP:  Postmenopausal   id     T1w             T2w
## 1 360 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 360       360          1072           7554174 2013-07-16 00:00:00.000000                 4342      6        massM   InvasiveDuctal
## CLINICAL INDICATION:  45 year-old female for assessment of extensive right breast carcinoma. The patient is also for assessment of nodules seen within the left breast on outside ultrasound.   id     T1w             T2w
## 1 361 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 361       361          1072           7554174 2013-07-16 00:00:00.000000                 4343      4        massB SCLEROSING PAPILLARY LESION
## CLINICAL INDICATION:  45 year-old female for assessment of extensive right breast carcinoma. The patient is also for assessment of nodules seen within the left breast on outside ultrasound.   id     T1w             T2w
## 1 362 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 362       362          1077           6890028 2011-09-28 00:00:00.000000                 4345      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  Assessment of area of palpable concern superior right breast.  History of lumpectomy inferior right breast (invasive mucinous cancer) with sentinel node biopsy December 2010 followed by radiation which finished March 2011.  Also in December 2010 DCIS removed from lateral left breast. 
## 
## Patient indicates persistent palpable concern on the MRI questionnaire, right upper, slightly outer breast.
## 
## LMP:  3 years ago   id     T1w             T2w
## 1 363 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 363       363          1078           7105247 2012-05-05 00:00:00.000000                 4347      4     nonmassB   DENSE FIBROSIS
## CLINICAL INDICATION:  OBSP high risk screening
## 
## LMP:  20/04/2012   id     T1w             T2w
## 1 364 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 364       364          1078           7105247 2012-05-05 00:00:00.000000                 4347      4     nonmassB   DENSE FIBROSIS
## CLINICAL INDICATION:  OBSP high risk screening
## 
## LMP:  20/04/2012   id     T1w             T2w
## 1 365 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 365       365          1079           7417880 2013-02-24 00:00:00.000000                 4348      4        massM     InsituDuctal
## CLINICAL INDICATION:  Biopsy proven left breast DCIS at 2 o'clock.  Right solid mass at 11 o'clock biopsied under ultrasound (pathology: Benign).
## 
## LMP:  Postmenopausal   id     T1w             T2w
## 1 366 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 366       366          1081           7078151 2012-07-25 00:00:00.000000                 4349      4        massB     FIBROADENOMA
## CLINICAL INDICATION:  OBSP high risk screen   id     T1w             T2w
## 1 367 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 367       367          1086           7173349 2012-06-13 00:00:00.000000                 4350      6     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: Known right DCIS for extent of disease  Post menopausal   id     T1w             T2w
## 1 368 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 368       368          1087           5360576 2010-08-26 00:00:00.000000                 4351      4        massB GRANULOMATOUS LOBULAR MASTITIS
## CLINICAL INDICATION: RIGHT inferior breast mass responding poorly
## to antibiotics with pain and nipple discharge. Post menopausal.   id     T1w             T2w
## 1 369 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 369       369          1087           5360576 2010-08-26 00:00:00.000000                 4351      4        massB GRANULOMATOUS LOBULAR MASTITIS
## CLINICAL INDICATION: RIGHT inferior breast mass responding poorly
## to antibiotics with pain and nipple discharge. Post menopausal.   id     T1w             T2w
## 1 370 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 370       370          1090           4288694 2007-03-15 00:00:00.000000                 4354      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: ACRIN 6666 study   id     T1w             T2w
## 1 371 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 371       371          1092           4951061 2009-04-30 00:00:00.000000                 4356      6        massB    InsituLobular
## CLINICAL INDICATION: Biopsy proven invasive lobular carcinoma
## right breast (microcalcifications on mammography). Preoperative
## extent of disease. LMP December, 2008.   id     T1w             T2w
## 1 372 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 372       372          1092           4951061 2009-04-30 00:00:00.000000                 4356      4     nonmassB    InsituLobular
## CLINICAL INDICATION: Biopsy proven invasive lobular carcinoma
## right breast (microcalcifications on mammography). Preoperative
## extent of disease. LMP December, 2008.   id     T1w             T2w
## 1 373 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 373       373          1095           4378323 2007-06-07 00:00:00.000000                 4357      3        massB     FIBROADENOMA
## CLINICAL INDICATION: Patient had right lumpectomy in January 2007
## for invasive and in situ carcinoma. Residual microcalcifications
## right breast, which are highly suspicious for residual disease.
## Probably benign mass left breast 9 o'clock for which the patient
## would prefer to have a core biopsy. MRI for residual assessment.
## 
## Comparison is made to the previous ultrasound and mammograms May
## 2007.   id     T1w             T2w
## 1 374 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 374       374          1095           4378323 2007-06-07 00:00:00.000000                 4358      5     nonmassM     InsituDuctal
## CLINICAL INDICATION: Patient had right lumpectomy in January 2007
## for invasive and in situ carcinoma. Residual microcalcifications
## right breast, which are highly suspicious for residual disease.
## Probably benign mass left breast 9 o'clock for which the patient
## would prefer to have a core biopsy. MRI for residual assessment.
## 
## Comparison is made to the previous ultrasound and mammograms May
## 2007.   id     T1w             T2w
## 1 375 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 375       375          1099           7646705 2013-10-31 00:00:00.000000                 4359      5        massM     InsituDuctal
## CLINICAL INDICATION:  Family history of breast cancer. Personal history of ADH left breast. Dense breast   id     T1w             T2w
## 1 376 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 376       376          1099           7646705 2013-10-31 00:00:00.000000                 4359      4     nonmassM     InsituDuctal
## CLINICAL INDICATION:  Family history of breast cancer. Personal history of ADH left breast. Dense breast   id     T1w             T2w
## 1 377 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 377       377          2007           7366811 2013-02-20 00:00:00.000000                 4276      4     nonmassB     FIBROADENOMA
## INDICATION:  Surveillance, OBSP.  Family history of breast cancer.  Left fibroadenoma excised in 2005.   id     T1w             T2w
## 1 378 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 378       378          2016           7052211 2012-02-14 00:00:00.000000                 4213      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 379 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 379       379          2017           7397047 2013-01-31 00:00:00.000000                 4263      6        massB RADIAL SCLEROSING LESION
## CLINICAL INDICATION: Known right invasive lobular carcinoma  Postmenopausal; Rt breast Ca lobular. r/o multicentric vs multifocal disease   id     T1w             T2w
## 1 380 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 380       380          2023           5141524 2010-01-07 00:00:00.000000                 3478      6        massM   InvasiveDuctal
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION: Recurrence right breast, for disease extent   id     T1w             T2w
## 1 381 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 381       381          2024           5190122 2010-02-08 00:00:00.000000                  648      6     nonmassM     InsituDuctal
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION: Known diagnosis of infiltrating lobular
## carcinoma right breast at 12 o'clock. Assess extent of disease   id     T1w             T2w
## 1 382 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 382       382          2027           5465838 2010-12-23 00:00:00.000000                 1772      6        massM     InsituDuctal
## CO-REGISTRATION STUDY 
## 
## HISTORY:  biopsy proven cancer upper outer right breast. Mammograms/ultrasound show irregular mass at 10 o'clock 17 cm from nipple with maximum dimension 3 cm on ultrasound. Calcifications extend anteriorly for 7 cm.  For assessment extent of disease and contralateral breast. FNAB of prominent right low axillary node: reactive.   id     T1w             T2w
## 1 383 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 383       383          2028           6702914 2011-02-08 00:00:00.000000                 3490      6     nonmassB      FIBROCYSTIC
## CO-REGISTRATION STUDY
## 
## 
## INDICATION:  Right IDC and left DCIS.  LMP January 16 2011.   id     T1w             T2w
## 1 384 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 384       384          2028           6702914 2011-02-08 00:00:00.000000                 3489      6        massM   InvasiveDuctal
## CO-REGISTRATION STUDY
## 
## 
## INDICATION:  Right IDC and left DCIS.  LMP January 16 2011.   id     T1w             T2w
## 1 385 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 385       385          2029           6716423 2011-02-24 00:00:00.000000                 1934      6        massB         ADENOSIS
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION: Known right invasive ductal carcinoma  LMP Feb 6/11   id     T1w             T2w
## 1 386 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 386       386          2033           6849696 2011-07-13 00:00:00.000000                 2806      6        massM  InvasiveLobular
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION:  Preoperative staging, left breast Ca   id     T1w             T2w
## 1 387 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 387       387          2033           6849696 2011-07-13 00:00:00.000000                 2807      4     nonmassB ATYPICAL LOBULAR HYPERPLASIA
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION:  Preoperative staging, left breast Ca   id     T1w             T2w
## 1 389 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label   lesion_diagnosis
## 389       389          2042           4964619 2009-05-28 00:00:00.000000                 3504      6        massB LobularHyperplasia
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION: Locally advanced right breast cancer.
## Follow-up enhancement within the left breast.   id     T1w             T2w
## 1 390 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 390       390          2042           5186978 2010-02-11 00:00:00.000000                  598      4        massB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 391 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 391       391          2042           7050570 2012-02-09 00:00:00.000000                  598      4        massB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 392 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 392       392          2049           5458850 2010-12-09 00:00:00.000000                 3596      5     nonmassM   InvasiveDuctal
## INDICATION: Bilateral breast masses, ? locally advanced breast
## cancer and inflammatory breast cancer. For extent of disease.
## LMP November 17 2010.   id     T1w             T2w
## 1 393 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 393       393          2049           5458850 2010-12-09 00:00:00.000000                 3595      5        massM   InvasiveDuctal
## INDICATION: Bilateral breast masses, ? locally advanced breast
## cancer and inflammatory breast cancer. For extent of disease.
## LMP November 17 2010.   id     T1w             T2w
## 1 394 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 394       394          2050           6689745 2011-01-27 00:00:00.000000                 1751      6        massM   InvasiveDuctal
## CO-REGISTRATION STUDY 
## 
## INDICATION:  Locally advanced breast cancer, assess extent of disease.  Biopsy performed January 10 2011 showed invasive duct cancer and positive axillary lymph node.  Post menopausal.   id     T1w             T2w
## 1 395 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 395       395          2051           6712632 2011-02-14 00:00:00.000000                 1952      6        massM     InsituDuctal
## CO-REGISTRATION STUDY
## 
## CLINICAL INDICATION:  Biopsy proven invasive carcinoma upper inner quadrant left breast corresponding to a palpable mass.  Biopsy proven ipsilateral metastatic axillary lymph node.  Recent clip placement in left carcinoma prior to neoadjuvant
## chemotherapy.   id     T1w             T2w
## 1 396 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 396       396          2051           6712632 2011-02-14 00:00:00.000000                 1953      6     nonmassB      FIBROCYSTIC
## CO-REGISTRATION STUDY
## 
## CLINICAL INDICATION:  Biopsy proven invasive carcinoma upper inner quadrant left breast corresponding to a palpable mass.  Biopsy proven ipsilateral metastatic axillary lymph node.  Recent clip placement in left carcinoma prior to neoadjuvant
## chemotherapy.   id     T1w             T2w
## 1 397 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 397       397          2053           6776964 2011-05-05 00:00:00.000000                 3284      6        massM   InvasiveDuctal
## CO-REGISTRATION STUDY 
## 
## CLINICAL INDICATION:  61 yo patient with LABC with positive node on FNA. Palpable nodes in both axilla.   id     T1w             T2w
## 1 398 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 398       398          2055           7041426 2012-02-03 00:00:00.000000                 3280      6     nonmassM   InvasiveDuctal
## CO-REGISTRATION STUDY 
## 
## INDICATION:  Bilateral breast cancer.  For extent of disease and to evaluate chest wall invasion.  LMP late December 2011.   id     T1w             T2w
## 1 399 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 399       399          2059           7749617 2014-08-05 00:00:00.000000                 4262      4     nonmassB COLUMNAR CELL CHANGES
## REGISTRATION STUDY: NON-MOTION CASE 
## 
## CLINICAL INDICATION: High risk screening   id     T1w             T2w
## 1 400 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 400       400          2065           7604632 2013-11-24 00:00:00.000000                 4360      4     nonmassB    InsituLobular
## CLINICAL INDICATION:  Right lumpectomy for LCIS. Screening   id     T1w             T2w
## 1 401 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 401       401          2068           7559583 2013-07-29 00:00:00.000000                 4362      5        massM   InvasiveDuctal
## CLINICAL INDICATION:  Left DCIS 2005 + >50% breast density.   id     T1w             T2w
## 1 402 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 402       402          2069           4976319 2009-05-26 00:00:00.000000                   49      6        massM     InsituDuctal
## CLINICAL INDICATION: 38 years old, biopsy proven DCIS left 6
## O'clock. Ultrasound describes a stellite nodule lateral to the
## biopsied mass. MRI for extent of the disease. LMP May/14/2009.   id     T1w             T2w
## 1 403 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 403       403          2071           7594721 2013-09-08 00:00:00.000000                 4364      4        massM   InvasiveDuctal
## CLINICAL INDICATION: Previous right mastectomy for DCIS  Post menopausal   id     T1w             T2w
## 1 404 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 404       404          2072           7256932 2012-09-13 00:00:00.000000                 4365      4        massM   InvasiveDuctal
## CLINICAL INDICATION:   39 year old female ca breast. right mastectomy and reconstruction with bilateral implants. new nodule R ant chest wall subcutaneous nodule lateral to manubrium, biopsy proven recurrent invasive ductal carcinoma.   id     T1w             T2w
## 1 405 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label              lesion_diagnosis
## 405       405          2073           4745825 2008-08-30 00:00:00.000000                 4366      5        massM InvasiveDuctal micropapillary
## CLINICAL INDICATION: Left axillary carcinoma. Searching for the
## primary   id     T1w             T2w
## 1 406 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 406       406          2073           4745825 2008-08-30 00:00:00.000000                 4366      5     nonmassM  InvasiveLobular
## CLINICAL INDICATION: Left axillary carcinoma. Searching for the
## primary   id     T1w             T2w
## 1 407 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 407       407          2075           6985605 2011-11-28 00:00:00.000000                 4369      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  Breast cancer with muscle invasion. Post-op. Query residual disease.   id     T1w             T2w
## 1 408 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 408       408          2075           6985605 2011-11-28 00:00:00.000000                 4370      4        massB     FIBROADENOMA
## CLINICAL INDICATION:  Breast cancer with muscle invasion. Post-op. Query residual disease.   id     T1w             T2w
## 1 409 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 409       409          2078           5116776 2010-04-25 00:00:00.000000                  434      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION: Prior left lumpectomy, sentinel lymph node
## biopsy and radiation (2008). Family history of breast cancer. 6
## month follow-up of two, probably benign areas of enhancement in
## the left breast on MRI of October 2009. LMP December 2008 (on
## Tamoxifen).   id     T1w             T2w
## 1 410 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label       lesion_diagnosis
## 410       410          2079           4591198 2009-02-12 00:00:00.000000                 1306      3        massB benign lymphoid tissue
## CLINICAL INDICATION: Prior left upper inner quadrant lumpectomy,
## axillary dissection and radiation. Follow-up probably benign
## right breast fibroadenoma. LMP October, 2008   id     T1w             T2w
## 1 411 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 411       411          3004           7691918 2014-02-16 00:00:00.000000                 4032      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 412 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 412       412          3004           7691918 2014-02-16 00:00:00.000000                 4032      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 413 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 413       413          3005           4974097 2009-06-19 00:00:00.000000                 1227      3     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION: Indeterminate microcalcifications right
## breast   id     T1w             T2w
## 1 414 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 414       414          3005           6757337 2012-04-04 00:00:00.000000                 4030      3     nonmassM     InsituDuctal
## CLINICAL INDICATION:  DCIS L BRCA lump + xrt   id     T1w             T2w
## 1 415 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 415       415          3005           5057668 2010-01-14 00:00:00.000000                 4030      2     nonmassB      FIBROCYSTIC
## CLINICAL INDICATION: Left lumpectomy and radiation 2000. Right
## surgical biopsy 6 o'clock 1997. Six-month follow of probably
## benign left breast subareolar enhancement. Stereotactic biopsy of
## small cluster calcifications upper outer right breast May, 2009
## benign on pathology. Attempt at stereotactic biopsy lower outer
## right breast calcifications too faint to visualize.   id     T1w             T2w
## 1 416 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 416       416          3005           6757337 2012-04-04 00:00:00.000000                 4030      4     nonmassM     InsituDuctal
## CLINICAL INDICATION:  DCIS L BRCA lump + xrt   id     T1w             T2w
## 1 417 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 417       417          3010           6828446 2011-06-21 00:00:00.000000                 1684      6        massM   InvasiveDuctal
## None[1] 419
##    id     T1w    T2w
## 1 419 VIBRANT REVIEW
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 419       419          3011           6898308 2011-09-06 00:00:00.000000                 4037      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 420 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label  lesion_diagnosis
## 420       420          3017           7014437 2012-01-31 00:00:00.000000                 4042      4     nonmassB FOCAL HYPERPLASIA
## None   id     T1w             T2w
## 1 421 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 421       421          3018           6865137 2011-07-27 00:00:00.000000                 4044      3        massB     FAT NECROSIS
## None   id     T1w             T2w
## 1 422 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 422       422          3020           7395195 2013-04-16 00:00:00.000000                 4048      4        massB STROMAL HYPERPLASIA
## None   id     T1w             T2w
## 1 423 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 423       423          3021           7019819 2012-01-07 00:00:00.000000                 4050      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 424 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 424       424          3023           7106703 2012-04-12 00:00:00.000000                 3952      6        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 425 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 425       425          3025           7103914 2013-02-05 00:00:00.000000                 4051      4        massB     FIBROADENOMA
## Invasive ductal carcinoma treated with lumpectomy, radiation and sentinel node biopsy 2008 in left breast.   id     T1w             T2w
## 1 426 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 426       426          3026           6830523 2011-06-23 00:00:00.000000                 4053      4        massB      FIBROCYSTIC
## CLINICAL INDICATION:  Multiple findings on ultrasound   id     T1w             T2w
## 1 427 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 427       427          3028           6991592 2011-12-12 00:00:00.000000                 4055      3        massB      HYPERPLASIA
## None   id     T1w             T2w
## 1 428 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 428       428          3030           7642998 2014-01-18 00:00:00.000000                 4056      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 429 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 429       429          3031           7106716 2012-04-08 00:00:00.000000                 4057      3        massB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 430 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 430       430          3033           5016967 2009-07-10 00:00:00.000000                  971      5        massM   InvasiveDuctal
## CLINICAL INDICATION: Known right multifocal invasive ductal
## cancer. IUD in place; no menstrual cycle.   id     T1w             T2w
## 1 431 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 431       431          3035           7002031 2012-01-06 00:00:00.000000                 4063      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 432 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 432       432          3035           7145247 2012-08-28 00:00:00.000000                 4064      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 433 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 433       433          3039           6894870 2011-08-29 00:00:00.000000                 4068      4        massB      HYPERPLASIA
## None   id     T1w             T2w
## 1 434 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 434       434          3045           7149704 2012-10-24 00:00:00.000000                 4071      4     nonmassB         ADENOSIS
## None   id     T1w             T2w
## 1 435 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 435       435          3046           7682447 2013-11-27 00:00:00.000000                 4072      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 436 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 436       436          3046           7289130 2012-12-22 00:00:00.000000                 2691      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 437 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 437       437          3046           7682447 2013-11-27 00:00:00.000000                 4072      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 438 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 438       438          3052           7100200 2012-03-31 00:00:00.000000                 4073      4        massB STROMAL FIBROSIS
## None   id     T1w             T2w
## 1 439 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 439       439          3052           7100200 2012-03-31 00:00:00.000000                 4073      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 440 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 440       440          3053           7041483 2012-01-31 00:00:00.000000                 4074      6        massB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 441 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 441       441          3053           7449310 2013-04-07 00:00:00.000000                 4075      5     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 442 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 442       442          3054           6714946 2011-06-18 00:00:00.000000                 4076      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 443 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 443       443          3055           7742700 2014-03-25 00:00:00.000000                 4077      4        massB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 444 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 444       444          3055           7060620 2012-03-07 00:00:00.000000                 2234      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 445 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 445       445          3055           7742700 2014-03-25 00:00:00.000000                 4077      4     nonmassB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 446 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 446       446          3055           7742700 2014-03-25 00:00:00.000000                 4077      4        massB STROMAL HYPERPLASIA
## None   id     T1w             T2w
## 1 447 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 447       447          3055           7060620 2012-03-07 00:00:00.000000                 2234      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 448 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 448       448          3057           7098623 2012-04-11 00:00:00.000000                 4079      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 449 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 449       449          3057           7098623 2012-04-11 00:00:00.000000                 4078      4        massB  TUBULAR ADENOMA
## None   id     T1w             T2w
## 1 450 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label        lesion_diagnosis
## 450       450          3063           7053508 2012-02-18 00:00:00.000000                 4047      6     nonmassM LYMPHOVASCULAR INVASION
## None   id     T1w             T2w
## 1 451 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 451       451          3065           7037223 2012-01-29 00:00:00.000000                 4087      4        massB         ADENOSIS
## None   id     T1w             T2w
## 1 452 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 452       452          3065           7037223 2012-01-29 00:00:00.000000                 4087      4        massB         ADENOSIS
## None   id     T1w             T2w
## 1 453 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                  lesion_diagnosis
## 453       453          3070           7085188 2012-04-04 00:00:00.000000                 4095      4        massB DUCTAL HYPERPLASIA WITHOUT ATYPIA
## None   id     T1w             T2w
## 1 454 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 454       454          3072           7054863 2012-02-12 00:00:00.000000                 4109      6     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 456 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 456       456          3073           7043941 2012-02-03 00:00:00.000000                 4104      6     nonmassM IN SITU PAPILLARY CARCINOMA
## None   id     T1w             T2w
## 1 457 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 457       457          3075           7064471 2012-02-24 00:00:00.000000                 4105      6        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 458 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 458       458          3075           7064471 2012-02-24 00:00:00.000000                 4106      6     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 459 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 459       459          3076           7053450 2012-02-15 00:00:00.000000                 4103      6        massM     InsituDuctal
## None   id     T1w             T2w
## 1 460 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 460       460          3076           7053450 2012-02-15 00:00:00.000000                 4103      6        massM     InsituDuctal
## None   id     T1w             T2w
## 1 461 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 461       461          3077           7042083 2012-01-31 00:00:00.000000                 4097      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 462 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 462       462          3077           7042083 2012-01-31 00:00:00.000000                 4096      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 463 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 463       463          3078           4836946 2009-02-14 00:00:00.000000                 4099      5     nonmassB    InsituLobular
## None   id     T1w             T2w
## 1 464 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 464       464          3080           7033654 2012-01-23 00:00:00.000000                 4100      6     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 465 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 465       465          3081           7041435 2012-01-29 00:00:00.000000                 4102      5     nonmassB COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 466 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 466       466          3081           7041435 2012-01-29 00:00:00.000000                 4102      5     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 467 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 467       467          3082           5355166 2010-08-19 00:00:00.000000                  347      6        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 468 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 468       468          3082           5355166 2010-08-19 00:00:00.000000                  347      6        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 469 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 469       469          3082           7080675 2012-03-20 00:00:00.000000                 4143      4     nonmassB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 470 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 470       470          3083           5345062 2010-08-12 00:00:00.000000                 4120      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 471 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 471       471          3086           7715466 2014-03-23 00:00:00.000000                 4115      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 472 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 472       472          3092           4462310 2007-09-17 00:00:00.000000                 4151      3        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 473 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 473       473          3093           7438787 2014-03-15 00:00:00.000000                 4152      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 474 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 474       474          3097           6909883 2012-01-03 00:00:00.000000                 4156      4        massB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 475 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 475       475          4002           6993690 2012-03-14 00:00:00.000000                 4158      5        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 476 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label          lesion_diagnosis
## 476       476          4003           7056445 2012-02-08 00:00:00.000000                 4159      4     nonmassB FLORID DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 477 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label          lesion_diagnosis
## 477       477          4003           7056445 2012-02-08 00:00:00.000000                 4159      4     nonmassB FLORID DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 478 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 478       478          4008           7014565 2012-01-03 00:00:00.000000                 4162      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 479 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 479       479          4008           7014565 2012-01-03 00:00:00.000000                 4161      6        massB      HYPERPLASIA
## None   id     T1w             T2w
## 1 480 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 480       480          4012           7002008 2012-03-15 00:00:00.000000                 4164      4        massB FOCAL USUAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 481 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 481       481          4017           6979356 2012-02-21 00:00:00.000000                 4167      4        massB USUAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 482 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label         lesion_diagnosis
## 482       482          4017           6979356 2012-02-21 00:00:00.000000                 4167      4        massB USUAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 483 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 483       483          4018           6983262 2011-11-25 00:00:00.000000                 4183      6     nonmassM     InsituDuctal
## None   id     T1w             T2w
## 1 484 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 484       484          4019           7151338 2012-09-22 00:00:00.000000                 3089      4        massB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 485 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 485       485          4020           6988975 2011-11-29 00:00:00.000000                 4169      6        massM  InvasiveLobular
## None   id     T1w             T2w
## 1 486 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 486       486          4021           6992707 2011-12-04 00:00:00.000000                 4172      4     nonmassB         ADENOSIS
## None   id     T1w             T2w
## 1 487 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 487       487          4023           7037125 2012-01-27 00:00:00.000000                 3850      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 488 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                lesion_diagnosis
## 488       488          4023           7037125 2012-01-27 00:00:00.000000                 3852      4        massB ADENOSIS, COLUMNAR CELL CHANGES
## None   id     T1w             T2w
## 1 489 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label               lesion_diagnosis
## 489       489          4023           7152678 2012-08-21 00:00:00.000000                 3849      4        massB BENIGN INTRAMAMMARY LYMPH NODE
## None   id     T1w             T2w
## 1 490 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 490       490          4023           7037125 2012-01-27 00:00:00.000000                 3849      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 491 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 491       491          4026           6998219 2012-01-24 00:00:00.000000                 4176      4        massB ATYPICAL DUCTAL HYPERPLASIA
## None   id     T1w             T2w
## 1 492 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 492       492          4029           7633460 2014-01-11 00:00:00.000000                 4178      4        massB    InsituLobular
## None   id     T1w             T2w
## 1 493 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 493       493          4029           7633460 2014-01-11 00:00:00.000000                 4178      4        massB    InsituLobular
## None   id     T1w             T2w
## 1 494 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 494       494          4039           7041331 2012-01-29 00:00:00.000000                 4181      6     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 495 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 495       495          4040           7003416 2011-12-17 00:00:00.000000                 2932      6        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 496 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 496       496          4040           7085105 2012-07-03 00:00:00.000000                 2935      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 497 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 497       497          4041           7003893 2011-12-20 00:00:00.000000                 4185      4     nonmassB BENIGN BREAST TISSUE
## None   id     T1w             T2w
## 1 498 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 498       498          4043           7041465 2012-01-31 00:00:00.000000                 4186      6     nonmassM   InvasiveDuctal
## None   id     T1w             T2w
## 1 499 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 499       499          4044           7066571 2012-02-26 00:00:00.000000                 4189      4        massB FIBROADENOMATOID
## None   id     T1w             T2w
## 1 500 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 500       500          4044           7066571 2012-02-26 00:00:00.000000                 4189      4        massB FOCAL CELLULAR STROMA
## None   id     T1w             T2w
## 1 501 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 501       501          4044           7066571 2012-02-26 00:00:00.000000                 4190      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 502 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 502       502          4045           7092118 2012-03-30 00:00:00.000000                 4187      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 503 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 503       503          4047           7009608 2011-12-21 00:00:00.000000                 4188      4        massB      FIBROCYSTIC
## None   id     T1w             T2w
## 1 504 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 504       504          4049           7009602 2012-01-05 00:00:00.000000                 4062      6        massM   InvasiveDuctal
## newly dx right breast ca with nodal involvement...for extent of disease
## 
## CLINICAL INDICATION:  Biopsy proven invasive ductal carcinoma and DCIS right breast. Biopsy proven right axillary lymph node metastases.   id     T1w             T2w
## 1 505 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label                      lesion_diagnosis
## 505       505          4055           7439091 2013-03-19 00:00:00.000000                 4205      4        massB PSEUDOANGIOMATOUS STROMAL HYPERPLASIA
## None   id     T1w             T2w
## 1 506 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 506       506          6001           4574766 2008-02-01 00:00:00.000000                    2      6     nonmassM   InvasiveDuctal
## Biopsy proven left breast invasive cancer with focal in situ component (prior ork-up done at an outside institution). Pre-operative MRI to determine extent of disease.   id     T1w             T2w
## 1 507 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 507       507          6004         ACC108249 2009-09-10 00:00:00.000000                  918      6        massM  InvasiveLobular
## known diagnosis of invasive lobular carcinoma   id     T1w             T2w
## 1 508 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 508       508          6004         ACC108249 2009-09-10 00:00:00.000000                  918      5     nonmassM  InvasiveLobular
## known diagnosis of invasive lobular carcinoma   id     T1w             T2w
## 1 509 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 509       509          6005         ACC108250 2009-09-08 00:00:00.000000                  844      5        massM  InvasiveLobular
## Bilateral breast carcinomas - according to EPR right
## mass is a biopsy proven invasive ductal cancer while the left
## breast biopsy showed LCIS (I do not have the biopsy or pathology
## reports). Evaluate extent of disease. LMP: September 7, 2009   id     T1w             T2w
## 1 510 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 510       510          6005         ACC108250 2009-09-08 00:00:00.000000                  844      5        massM  InvasiveLobular
## Bilateral breast carcinomas - according to EPR right
## mass is a biopsy proven invasive ductal cancer while the left
## breast biopsy showed LCIS (I do not have the biopsy or pathology
## reports). Evaluate extent of disease. LMP: September 7, 2009   id     T1w             T2w
## 1 511 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 511       511          6005         ACC108250 2009-09-08 00:00:00.000000                  844      5        massM  InvasiveLobular
## Bilateral breast carcinomas - according to EPR right
## mass is a biopsy proven invasive ductal cancer while the left
## breast biopsy showed LCIS (I do not have the biopsy or pathology
## reports). Evaluate extent of disease. LMP: September 7, 2009   id     T1w             T2w
## 1 512 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 512       512          6008           4644038 2008-04-24 00:00:00.000000                    8      6        massM   InvasiveDuctal
## Known right breast cancer   id     T1w             T2w
## 1 513 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 513       513          6014           5101372 2009-10-20 00:00:00.000000                   69      6        massM   InvasiveDuctal
## Left breast cancer, biopsy proven elsewhere.
## For pre operative evaluation. LMP October 14-18 2009.   id     T1w             T2w
## 1 514 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 514       514          6015           5082265 2009-09-26 00:00:00.000000                 3030      6        massM   InvasiveDuctal
## 34 years-old female. Locally advanced breast
## cancer right breast. Evaluate extent of disease and left breast.   id     T1w             T2w
## 1 515 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 515       515          6015           5082265 2009-09-26 00:00:00.000000                 3030      6     nonmassM   InvasiveDuctal
## 34 years-old female. Locally advanced breast
## cancer right breast. Evaluate extent of disease and left breast.   id     T1w             T2w
## 1 516 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 516       516          6017           5086121 2009-10-06 00:00:00.000000                  833      6        massM  InvasiveLobular
## Core biopsy suspicious mammographic and
## ultrasound findings 12 o'clock left breast showed invasive lobular
## carcinoma, classic type. Lumpectomy superior central right breast
## 1988, axillary dissection and radiation for 2.4 cm IDC, 5/7
## positive nodes. For assessment extent of disease and contralate   id     T1w             T2w
## 1 517 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 517       517          6017           5086121 2009-10-06 00:00:00.000000                  833      2        massB     FIBROADENOMA
## Core biopsy suspicious mammographic and
## ultrasound findings 12 o'clock left breast showed invasive lobular
## carcinoma, classic type. Lumpectomy superior central right breast
## 1988, axillary dissection and radiation for 2.4 cm IDC, 5/7
## positive nodes. For assessment extent of disease and contralate   id     T1w             T2w
## 1 518 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 518       518          6018           5088825 2009-10-06 00:00:00.000000                  799      5     nonmassM  InvasiveLobular
## Mass under the left nipple. Nipple discharge. Prior
## history of left lumpectomy. Post menopause.   id     T1w             T2w
## 1 519 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 519       519          6019         ACC109175 2009-10-10 00:00:00.000000                 1024      4        massM  InvasiveLobular
## 46 year old with dense breast. Left breast
## mass, biopsied at an outside institution with a diagnosis of
## fibroepithelial lesion. Repeat biopsy at this institution with a diagnosis of ALH and LCIS.   id     T1w             T2w
## 1 520 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 520       520          6020         ACC109177 2009-10-09 00:00:00.000000                 2518      6        massM   InvasiveDuctal
## Invasive ductal carcinoma right breast   id     T1w             T2w
## 1 521 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label             lesion_diagnosis
## 521       521          6021           4798692 2009-03-14 00:00:00.000000                 1094      4        massB ATYPICAL LOBULAR HYPERPLASIA
##  High risk screening MRI.
##    id     T1w             T2w
## 1 522 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 522       522          6022           5046558 2009-08-16 00:00:00.000000                 1133      4     nonmassB     FIBROADENOMA
## Known left locally advanced breast cancer.
## LMP August 14/09   id     T1w             T2w
## 1 523 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 523       523          6022           5046558 2009-08-16 00:00:00.000000                 3166      6        massM   InvasiveDuctal
## Known left locally advanced breast cancer.
## LMP August 14/09   id     T1w             T2w
## 1 524 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 524       524          6023           4697014 2009-06-23 00:00:00.000000                  821      3        massB      FIBROCYSTIC
## High risk screening LMP May 20   id     T1w             T2w
## 1 525 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 525       525          6023           4697014 2009-06-23 00:00:00.000000                  821      3        massB SCLEROSING ADENOSIS
## High risk screening LMP May 20   id     T1w             T2w
## 1 526 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 526       526          6024           5008021 2009-06-30 00:00:00.000000                 1058      5        massM   InvasiveDuctal
## Suspicious right breast mass and
## calcifications. LMP June 30 2009.   id     T1w             T2w
## 1 527 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 527       527          6024           5008021 2009-06-30 00:00:00.000000                 1058      5     nonmassM   InvasiveDuctal
## Suspicious right breast mass and
## calcifications. LMP June 30 2009.   id     T1w             T2w
## 1 528 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 528       528          6025           5111910 2009-11-05 00:00:00.000000                  732      6        massM   InvasiveDuctal
## Left breast cancer with dense breasts. LMP roughly 14 days ago.   id     T1w             T2w
## 1 529 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 529       529          6025           5111910 2009-11-05 00:00:00.000000                  732      6        massM   InvasiveDuctal
## Left breast cancer with dense breasts. LMP roughly 14 days ago.   id     T1w             T2w
## 1 530 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 530       530          6025           5111910 2009-11-05 00:00:00.000000                  732      6     nonmassM   InvasiveDuctal
## Left breast cancer with dense breasts. LMP roughly 14 days ago.   id     T1w             T2w
## 1 531 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 531       531          6026           4888386 2009-02-24 00:00:00.000000                 2519      4     nonmassM     InsituDuctal
## dcis rt breast, pathology shows close margins assess for extent of disease   id     T1w             T2w
## 1 532 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 532       532          6027           4770166 2008-09-28 00:00:00.000000                 3097      4        massB BENIGN BREAST TISSUE
## Left breast mass at two o'clock, pathology-invasive ductal carcinoma. MRI for extent of disease.   id     T1w             T2w
## 1 533 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 533       533          6027           4770166 2008-09-28 00:00:00.000000                 3097      4     nonmassB BENIGN BREAST TISSUE
## Left breast mass at two o'clock, pathology-invasive ductal carcinoma. MRI for extent of disease.   id     T1w             T2w
## 1 534 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 534       534          6029           5083338 2009-10-16 00:00:00.000000                 3307      6        massB     FIBROADENOMA
## Known right breast cancer.   id     T1w             T2w
## 1 535 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 535       535          6029           5083338 2009-10-16 00:00:00.000000                 3307      6        massB     FIBROADENOMA
## Known right breast cancer.   id     T1w             T2w
## 1 536 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 536       536          6029           6772981 2011-05-08 00:00:00.000000                 3307      4        massB     FIBROADENOMA
## None   id     T1w             T2w
## 1 537 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 537       537          6032           4982490 2009-06-01 00:00:00.000000                 2523      4     nonmassB      FIBROCYSTIC
## 50 years old, left lower outer lumpectomy
## and sentinel node biopsy (negative) in Feb/2009, close margins.
## Faint calcifications medial and lateral to surgical bed on
## mammogram. MRI to rule out residual disease.   id     T1w             T2w
## 1 538 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 538       538          6034           4997881 2009-06-15 00:00:00.000000                   67      6        massM   InvasiveDuctal
## 37 years old with right breast palpable abnormality.
## Biopsy proven right breast carcinoma with positive node(biopsy
## performed in an outside institution). Biopsy proven fat epithelial
## tissue with atypia of a second right breast mass. Family history
## of breast cancer (mother at age 50). LMP June/09   id     T1w             T2w
## 1 539 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 539       539          6034           4997881 2009-06-15 00:00:00.000000                 2586      6        massM   InvasiveDuctal
## 37 years old with right breast palpable abnormality.
## Biopsy proven right breast carcinoma with positive node(biopsy
## performed in an outside institution). Biopsy proven fat epithelial
## tissue with atypia of a second right breast mass. Family history
## of breast cancer (mother at age 50). LMP June/09   id     T1w             T2w
## 1 540 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 540       540          6034           4997881 2009-06-15 00:00:00.000000                 2586      6     nonmassM   InvasiveDuctal
## 37 years old with right breast palpable abnormality.
## Biopsy proven right breast carcinoma with positive node(biopsy
## performed in an outside institution). Biopsy proven fat epithelial
## tissue with atypia of a second right breast mass. Family history
## of breast cancer (mother at age 50). LMP June/09   id     T1w             T2w
## 1 541 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 541       541          6035           5062962 2009-09-05 00:00:00.000000                  800      5        massM   InvasiveDuctal
## Baseline mammogram and subsequent ultrasound showed
## left calcifications and masses. LMP September 5 2009.   id     T1w             T2w
## 1 542 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 542       542          6035           5062962 2009-09-05 00:00:00.000000                  800      5        massM   InvasiveDuctal
## Baseline mammogram and subsequent ultrasound showed
## left calcifications and masses. LMP September 5 2009.   id     T1w             T2w
## 1 543 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 543       543          6035           5062962 2009-09-05 00:00:00.000000                  800      5        massM   InvasiveDuctal
## Baseline mammogram and subsequent ultrasound showed
## left calcifications and masses. LMP September 5 2009.   id     T1w             T2w
## 1 544 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 544       544          6037           5043444 2009-08-13 00:00:00.000000                  933      5        massM   InvasiveDuctal
## Clinically apparent locally advanced breast
## carcinoma on the left, pre-therapeutic workup   id     T1w             T2w
## 1 545 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 545       545          6037           5043444 2009-08-13 00:00:00.000000                  934      5        massM   InvasiveDuctal
## Clinically apparent locally advanced breast
## carcinoma on the left, pre-therapeutic workup   id     T1w             T2w
## 1 546 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 546       546          6038           5044471 2009-08-14 00:00:00.000000                   66      6        massM   InvasiveDuctal
## Biopsy proven left breast cancer for
## staging.   id     T1w             T2w
## 1 547 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 547       547          6038           5044471 2009-08-14 00:00:00.000000                   66      6     nonmassB BENIGN BREAST TISSUE
## Biopsy proven left breast cancer for
## staging.   id     T1w             T2w
## 1 548 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 548       548          6038           5044471 2009-08-14 00:00:00.000000                 3043      6     nonmassB BENIGN BREAST TISSUE
## Biopsy proven left breast cancer for
## staging.   id     T1w             T2w
## 1 549 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 549       549          6039         ACC109197 2009-09-29 00:00:00.000000                  803      5        massM     InsituDuctal
## 44 years-old female . Evaluate extent of
## disease. Abnormal mammogram and ultrasound.   id     T1w             T2w
## 1 550 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 550       550          6039         ACC109197 2009-09-29 00:00:00.000000                  803      5     nonmassM     InsituDuctal
## 44 years-old female . Evaluate extent of
## disease. Abnormal mammogram and ultrasound.   id     T1w             T2w
## 1 551 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 551       551          6039         ACC109197 2009-09-29 00:00:00.000000                 1026      5        massM     InsituDuctal
## 44 years-old female . Evaluate extent of
## disease. Abnormal mammogram and ultrasound.   id     T1w             T2w
## 1 552 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 552       552          6040           5075204 2009-09-21 00:00:00.000000                 2529      5        massM   InvasiveDuctal
## Left lumpectomy August 13, 2009 with chest wall
## invasion. Post menopausal.   id     T1w             T2w
## 1 553 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 553       553          6041           5104414 2009-10-27 00:00:00.000000                 2532      6        massM   InvasiveDuctal
## palpable abnormality lateral right breast with core
## biopsy proven poorly differentiated carcinoma. For assessment
## extent of disease and contralateral breast.   id     T1w             T2w
## 1 554 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 554       554          6041           5104414 2009-10-27 00:00:00.000000                 2532      6        massM   InvasiveDuctal
## palpable abnormality lateral right breast with core
## biopsy proven poorly differentiated carcinoma. For assessment
## extent of disease and contralateral breast.   id     T1w             T2w
## 1 555 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 555       555          6042           4504274 2008-01-19 00:00:00.000000                 2587      3        massM   InvasiveDuctal
## Strong family history of breast cancer; 25%
## lifetime risk of breast cancer. New palpable nodule right upper
## outer breast.   id     T1w             T2w
## 1 556 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 556       556          6042           4504274 2008-01-19 00:00:00.000000                   25      3        massM   InvasiveDuctal
## Strong family history of breast cancer; 25%
## lifetime risk of breast cancer. New palpable nodule right upper
## outer breast.   id     T1w             T2w
## 1 557 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 557       557          6043           5249778 2010-04-17 00:00:00.000000                 2588      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## 38 year old with prior left mastectomy
## (2008) for IDC with 2/9 positive nodes. Increasing calcifications
## in the right breast. LMP March 2009. On Tamoxifen.   id     T1w             T2w
## 1 558 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 558       558          6044           5078981 2009-09-25 00:00:00.000000                  823      5        massB      FIBROCYSTIC
## Suspion of multifocal cancer right breast.
## For extent of disease and contralateral breast assessment. Remote
## benign biopsy right upper breast. LMP Aug 26/09   id     T1w             T2w
## 1 559 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 559       559          6044           5078981 2009-09-25 00:00:00.000000                 1043      5        massM   InvasiveDuctal
## Suspion of multifocal cancer right breast.
## For extent of disease and contralateral breast assessment. Remote
## benign biopsy right upper breast. LMP Aug 26/09   id     T1w             T2w
## 1 560 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 560       560          6045           5208117 2010-02-25 00:00:00.000000                 2583      6        massM   InvasiveDuctal
## Right breast cancer. Extent of disease. Mother with
## breast cancer age 55. LMP February 8 2010.   id     T1w             T2w
## 1 561 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 561       561          6045           5208117 2010-02-25 00:00:00.000000                 2583      6        massM   InvasiveDuctal
## Right breast cancer. Extent of disease. Mother with
## breast cancer age 55. LMP February 8 2010.   id     T1w             T2w
## 1 562 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 562       562          6046         ACC108189 2010-03-01 00:00:00.000000                  567      5        massM     InsituDuctal
## 36 years-old female.Probable large right
## breast carcinoma. Evaluate extent of disease   id     T1w             T2w
## 1 563 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 563       563          6046         ACC108189 2010-03-01 00:00:00.000000                 1037      5        massM     InsituDuctal
## 36 years-old female.Probable large right
## breast carcinoma. Evaluate extent of disease   id     T1w             T2w
## 1 564 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 564       564          6046         ACC108189 2010-03-01 00:00:00.000000                 1037      5        massM     InsituDuctal
## 36 years-old female.Probable large right
## breast carcinoma. Evaluate extent of disease   id     T1w             T2w
## 1 565 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 565       565          6046         ACC108189 2010-03-01 00:00:00.000000                 1037      5        massM     InsituDuctal
## 36 years-old female.Probable large right
## breast carcinoma. Evaluate extent of disease   id     T1w             T2w
## 1 566 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 566       566          6047           5275305 2010-05-18 00:00:00.000000                 2584      6     nonmassM   Adenocarcinoma
## Right axillary node excised with metastatic
## breast adenocarcinoma LMP about 3 weeks ago   id     T1w             T2w
## 1 567 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 567       567          6048           5284266 2010-05-30 00:00:00.000000                  402      6        massM   InvasiveDuctal
## 31 year-old female with positive FNA for
## breast cancer at 2 o'clock right breast   id     T1w             T2w
## 1 568 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 568       568          6050           5225817 2010-05-15 00:00:00.000000                 2536      4     nonmassB         FIBROSIS
## Surveillance. Right lumpectomy April 2009 with
## axillary node dissection, radiation, chemotherapy, and tamoxifen.
## BRCA 2 variant. LMP February 20, 2010.   id     T1w             T2w
## 1 569 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 569       569          6051           5426079 2010-11-05 00:00:00.000000                 2585      6        massM   InvasiveDuctal
## Evaluation of disease - newly diagnosed left invasive breast carcinoma.   id     T1w             T2w
## 1 570 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 570       570          6052           5369136 2010-09-10 00:00:00.000000                 2231      6     nonmassM   InvasiveDuctal
## Known right DCIS   id     T1w             T2w
## 1 571 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 571       571          6054           5425486 2010-11-29 00:00:00.000000                 2540      5        massM     InsituDuctal
## Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.
## Left breast severe atypical intraductal proliferative lesion.   id     T1w             T2w
## 1 572 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 572       572          6054           5425486 2010-11-29 00:00:00.000000                 2540      5        massM     InsituDuctal
## Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.
## Left breast severe atypical intraductal proliferative lesion.   id     T1w             T2w
## 1 573 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 573       573          6054           5425486 2010-11-29 00:00:00.000000                 2539      5        massM     InsituDuctal
## Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.
## Left breast severe atypical intraductal proliferative lesion.   id     T1w             T2w
## 1 574 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 574       574          6054           5425486 2010-11-29 00:00:00.000000                 2540      5        massM     InsituDuctal
## Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.
## Left breast severe atypical intraductal proliferative lesion.   id     T1w             T2w
## 1 575 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 575       575          6054           5425486 2010-11-29 00:00:00.000000                 2539      5        massM     InsituDuctal
## Bilateral needle biopsies. The right breast with a diagnosis of atypical ductal epithelial hyperplasia with microcalcificaitons and adenosis.
## Left breast severe atypical intraductal proliferative lesion.   id     T1w             T2w
## 1 576 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 576       576          6069           7581124 2013-10-29 00:00:00.000000                 3594      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  2yr hx of bil multiductal discharge.  occassionaly seen frank or old blood in diff ducts at diff times. ? ductectasia - r/o other, MR guided biopsy March 2011 left breast mass showed benign pathology.   id     T1w             T2w
## 1 577 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 577       577          6069           7581124 2013-10-29 00:00:00.000000                 3594      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  2yr hx of bil multiductal discharge.  occassionaly seen frank or old blood in diff ducts at diff times. ? ductectasia - r/o other, MR guided biopsy March 2011 left breast mass showed benign pathology.   id     T1w             T2w
## 1 578 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 578       578          6100           6722170 2011-03-07 00:00:00.000000                 1975      5        massM   InvasiveDuctal
## CLINICAL INDICATION:  DCIS left breast.  Multiple new masses right breast.   id     T1w             T2w
## 1 579 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 579       579          6101           5087078 2009-10-26 00:00:00.000000                 3516      4     nonmassB      FIBROCYSTIC
## CLINICAL INDICATION: Follow-up of probably benign lesion left
## breast upper outer quadrant. Prior right lumpectomy, axillary
## dissection and radiation 2005.   id     T1w             T2w
## 1 580 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 580       580          6101           7709238 2014-01-07 00:00:00.000000                 3517      6        massM   InvasiveDuctal
## CLINICAL INDICATION:  65-year-old female with history of previous right lumpectomy and axillary lymph node dissection for invasive ductal carcinoma. Prior benign biopsy left breast November 2009. Recent area of palpable concern upper outer right breast corresponding to superficial 0.6 cm mass 10 o'clock position 5 cm from the nipple, biopsy proven invasive ductal carcinoma. MRI to assess for multifocality and guide between mastectomy or repeat lumpectomy.   id     T1w             T2w
## 1 581 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 581       581          6105           5069712 2009-09-12 00:00:00.000000                 4018      4     nonmassM     InsituDuctal
## CLINICAL INDICATION: Right lumpectomy, axillary dissection and
## radiation 1993/1994. Increasing calcifications right breast
## recommended for stereotactic biopsy and MRI assessment.   id     T1w             T2w
## 1 582 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 582       582          6114           5148523 2009-12-10 00:00:00.000000                 3544      6     nonmassB BENIGN BREAST TISSUE
## HISTORY: locally advanced right breast carcinoma, on chemotherapy.
## Outside core biopsy right subareolar mass consistent with
## carcinoma with biopsy proven metastasis to a right axillary node.
## Clinical finding left breast consisting of fullness laterally. 5
## mm nodule 3 o'clock left breast described on outside ultrasound.
## History of benign surgery medial left breast.   id     T1w             T2w
## 1 583 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 583       583          6117           5154282 2010-05-04 00:00:00.000000                  111      3        massB     FIBROADENOMA
## CLINICAL INDICATION: Left lumpectomy May 2009 for multifocal
## breast cancer, positive nodes. For screening. LMP 2008.   id     T1w             T2w
## 1 584 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 584       584          6141           7044114 2012-05-05 00:00:00.000000                 3564      2        massB     FIBROADENOMA
## CLINICAL INDICATION:  High risk screening. Hodgkins disease and mantle radiation.   id     T1w             T2w
## 1 585 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 585       585          6148           7446343 2013-05-25 00:00:00.000000                 3587      4        massB SCLEROSING ADENOSIS
## CLINICAL INDICATION:  Previous left breast cancer (2004). Extremely dense breast.   id     T1w             T2w
## 1 586 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 586       586          6150           7128025 2012-05-08 00:00:00.000000                 3588      4     nonmassB COLUMNAR CELL CHANGES
## CLINICAL INDICATION:  fx breast ca  pre. mri here recent mri recalled for 2nd look u/s for rpt mri   id     T1w             T2w
## 1 587 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 587       587          6164           6971531 2012-04-24 00:00:00.000000                 3605      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  37 year old fat necrosis breast   id     T1w             T2w
## 1 588 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 588       588          6174           7009629 2012-05-05 00:00:00.000000                 3618      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION:  BRCA1. 6 month f/u R breast NMLE enhancement. History of benign MR biopsy 12 o'clock right breast November 2010. Benign US core biopsy 12 o'clock right breast October 2010.  Lumpectomy upper outer left breast 2009, axillary
## dissection and radiation for previous carcinoma. Right axillary sentinel node biopsy.   id     T1w             T2w
## 1 589 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 589       589          6223           7043947 2012-02-02 00:00:00.000000                 3861      4        massM   InvasiveDuctal
## INDICATION:  6 month follow up of left enhancement.  Prior right mastectomy (DCIS) and reconstruction.  Prior left lumpectomy and reduction.  Patient indicates left pain and lump on history sheet (site not specified).  LMP August 2010.   id     T1w             T2w
## 1 590 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 590       590          6224           4559525 2008-06-09 00:00:00.000000                 3873      4     nonmassB      FIBROCYSTIC
## CLINICAL INDICATION: Hodgkin's lymphoma with mantle radiation at
## age 26. High risk screening. 6 months follow up for a short linear
## enhancement in the right retroareolar area.   id     T1w             T2w
## 1 591 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 591       591          6226           6718391 2011-08-09 00:00:00.000000                 4034      4        massB BENIGN BREAST TISSUE
## INDICATION:  38 year-old female, previous right mastectomy   id     T1w             T2w
## 1 592 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 592       592          6233           7047121 2012-02-05 00:00:00.000000                 4090      6        massB     FIBROADENOMA
## CLINICAL INDICATION: Known invasive lobular carcinoma left breast.  For extent of disease Postmenopausal   id     T1w             T2w
## 1 593 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 593       593          7008           6875110 2011-08-04 00:00:00.000000                 3506      6        massM  InvasiveLobular
## None   id     T1w             T2w
## 1 594 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 594       594          7011           6918051 2011-09-20 00:00:00.000000                 3459      4        massB BENIGN BREAST TISSUE
## INDICATION:  Surveillance.  Family history of breast/ovarian  cancer.  LMP June 2011.   id     T1w             T2w
## 1 595 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 595       595          7018           6803089 2011-05-24 00:00:00.000000                 3315      4        massM     InsituDuctal
## CLINICAL INDICATION:  Suspicious segmental microcalcifications on the mammogram   id     T1w             T2w
## 1 596 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 596       596          7018           7138226 2012-05-13 00:00:00.000000                 2178      2        massM     InsituDuctal
## CLINICAL INDICATION:  Highly suspicious for ductal carcinoma in situ right lower inner quadrant   id     T1w             T2w
## 1 597 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 597       597          7024           6805356 2011-05-31 00:00:00.000000                 3317      4        massB     FIBROADENOMA
## CLINICAL INDICATION:  28 year old patient with infiltrating ductal carcinoma left upper quadrant with multiple foci plus 5 cm of DCIS   id     T1w             T2w
## 1 598 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 598       598          7029           7014263 2012-05-27 00:00:00.000000                 3324      4     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION:  OBSP HIGH RISK   id     T1w             T2w
## 1 599 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 599       599          7030           7538617 2013-08-07 00:00:00.000000                 3325      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk, annual recall   id     T1w             T2w
## 1 600 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 600       600          7030           7538617 2013-08-07 00:00:00.000000                 3325      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP high risk, annual recall   id     T1w             T2w
## 1 601 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 601       601          7043           7119983 2013-02-02 00:00:00.000000                 1983      4        massB     FIBROADENOMA
## CLINICAL INDICATION: High risk  LMP Jan 24/13   id     T1w             T2w
## 1 602 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 602       602          7045           6760802 2011-09-05 00:00:00.000000                 3340      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION: High risk - follow-up of right subareolar mass  LMP Sept 16/11   id     T1w             T2w
## 1 603 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 603       603          7053           7956343 2014-10-13 00:00:00.000000                 4027      4     nonmassB  FIBROTIC STROMA
## None   id     T1w             T2w
## 1 604 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 604       604          7066           6715383 2011-11-11 00:00:00.000000                 3365      4        massB     FIBROADENOMA
## CLINICAL INDICATION: High risk  LMP Oct 30/11   id     T1w             T2w
## 1 605 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 605       605          7066           7395276 2013-02-12 00:00:00.000000                 3371      4     nonmassB COLUMNAR CELL CHANGES
## CLINICAL INDICATION:  OBSP high risk.  Benign MRI guided biopsy right breast a year ago   id     T1w             T2w
## 1 606 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 606       606          7076           7267446 2013-04-27 00:00:00.000000                 3401      3     nonmassM     InsituDuctal
## CLINICAL INDICATION: HIgh risk  LMP Apr 16/13  Previous excision of breast tissue from axillary tails bilaterally   id     T1w             T2w
## 1 607 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 607       607          7076           7267446 2013-04-27 00:00:00.000000                 3400      3     nonmassB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: HIgh risk  LMP Apr 16/13  Previous excision of breast tissue from axillary tails bilaterally   id     T1w             T2w
## 1 608 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 608       608          7077           5077480 2009-09-24 00:00:00.000000                 1812      5        massM     InsituDuctal
## rior history of DCIS in 2000 treated with
## lumpectomy and radiation. Suspected local recurrence and
## indeterminant left breast mammographic findings. Strong family
## history of breast cancer.   id     T1w             T2w
## 1 609 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 609       609          7077           5077480 2009-09-24 00:00:00.000000                 1810      5        massM     InsituDuctal
## rior history of DCIS in 2000 treated with
## lumpectomy and radiation. Suspected local recurrence and
## indeterminant left breast mammographic findings. Strong family
## history of breast cancer.   id     T1w             T2w
## 1 610 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 610       610          7077           5077480 2009-09-24 00:00:00.000000                 1810      5     nonmassM     InsituDuctal
## rior history of DCIS in 2000 treated with
## lumpectomy and radiation. Suspected local recurrence and
## indeterminant left breast mammographic findings. Strong family
## history of breast cancer.   id     T1w             T2w
## 1 611 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 611       611          7085           7616788 2014-01-22 00:00:00.000000                 3409      2        massB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: OBSP high risk screening.   id     T1w             T2w
## 1 612 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 612       612          7086           6938067 2012-01-26 00:00:00.000000                 3416      4        massM     InsituDuctal
## CLINICAL INDICATION:  Screening - Greater than 25% lifetime risk.   id     T1w             T2w
## 1 613 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 613       613          7088           7066921 2012-07-17 00:00:00.000000                 3420      3     nonmassB      FIBROCYSTIC
## HISTORY: 44 year old female, high risk screening. Follow up probably benign findings. Unable to biopsy area of enhancement lateral right breast first seen in October 2011 due to presence of implant.
## 
## Silicone implants. FNAB prominent left axillary node June 2010: consistent with reactive node.   id     T1w             T2w
## 1 614 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 614       614          7094           7171259 2012-06-16 00:00:00.000000                 3437      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  OBSP High risk screen. IBIS 30.88%.   id     T1w             T2w
## 1 615 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 615       615          7096           6869668 2011-07-28 00:00:00.000000                  452      3        massB     FIBROADENOMA
## HISTORY: 37 year old female recently investigated for vague palpable abnormality right breast 7 o'clock. 8 mm intraductal mass seen on ultrasound 7 o'clock right breast and biopsy is being arranged. Left retroareolar dilated ducts with echogenic
## material: query debris versus intraductal lesions. For further MRI assessment. US 2000 demonstrated stable 2 cm probably fibroadenoma 3 o'clock subareolar region right breast and this was palpable.   id     T1w             T2w
## 1 616 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 616       616          7096           6869668 2011-07-28 00:00:00.000000                  452      3        massB     FIBROADENOMA
## HISTORY: 37 year old female recently investigated for vague palpable abnormality right breast 7 o'clock. 8 mm intraductal mass seen on ultrasound 7 o'clock right breast and biopsy is being arranged. Left retroareolar dilated ducts with echogenic
## material: query debris versus intraductal lesions. For further MRI assessment. US 2000 demonstrated stable 2 cm probably fibroadenoma 3 o'clock subareolar region right breast and this was palpable.   id     T1w             T2w
## 1 617 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label    lesion_diagnosis
## 617       617          7097           6805449 2011-06-24 00:00:00.000000                 3439      4        massB SCLEROSING ADENOSIS
## CLINICAL INDICATION:  Family history of breast cancer.  Risk of more than 25%.   id     T1w             T2w
## 1 618 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 618       618          7097           6805449 2011-06-24 00:00:00.000000                 3440      4        massB     FIBROADENOMA
## CLINICAL INDICATION:  Family history of breast cancer.  Risk of more than 25%.   id     T1w             T2w
## 1 619 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label            lesion_diagnosis
## 619       619          7097           7388464 2013-06-12 00:00:00.000000                 3441      2        massB ATYPICAL DUCTAL HYPERPLASIA
## CLINICAL INDICATION: High risk screening.  Previous benign left MRI guided biopsies  LMP June 2/13   id     T1w             T2w
## 1 620 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 620       620          7104           6941351 2011-10-25 00:00:00.000000                 3442      5        massM   InvasiveDuctal
## CLINICAL INDICATION: Suspicious mass left subareolar with nipple changes  LMP Sept 10/11   id     T1w             T2w
## 1 621 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 621       621          7105           7837892 2014-10-25 00:00:00.000000                 4061      4        massM   InvasiveDuctal
## None   id     T1w             T2w
## 1 622 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 622       622          7127           6989740 2012-01-28 00:00:00.000000                 3446      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  High risk screening. BRCA2. Followup right breast mass.   id     T1w             T2w
## 1 623 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 623       623          7151           7557684 2014-02-01 00:00:00.000000                 3453      2        massB      HYPERPLASIA
## CLINICAL INDICATION:  OBSP High Risk Screening.   id     T1w             T2w
## 1 624 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 624       624          7159           5435020 2010-11-11 00:00:00.000000                 1668      4     nonmassM  InvasiveLobular
## CLINICAL INDICATION: LABC left breast on neoadjuvant tx, extent
## of disease   id     T1w             T2w
## 1 625 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 625       625          7165           5021830 2010-01-09 00:00:00.000000                 4022      3        massB         ADENOSIS
## CLINICAL INDICATION: 41 years-old female. Family history breast
## and ovarian carcinoma. Lifetime risk greater than 25%. Left breast
## mass follow-up.   id     T1w             T2w
## 1 626 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 626       626          7178           7074874 2012-03-03 00:00:00.000000                 3834      6        massM   InvasiveDuctal
## CLINICAL INDICATION:  Known poorly differentiated invasive ductal carcinoma left breast   id     T1w             T2w
## 1 627 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 627       627          7183           7404761 2013-12-21 00:00:00.000000                 4036      4        massB      FIBROCYSTIC
## CLINICAL INDICATION:  Breast MRI. Annual Recall OBSP High Risk   id     T1w             T2w
## 1 628 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 628       628          7186           5263507 2010-04-30 00:00:00.000000                  462      6        massM   InvasiveDuctal
## INDICATION: Left breast mass seen on mammogram and ultrasound 5
## o'clock position - ultrasound guided core biopsy of the left
## breast mass shows invasive duct carcinoma. Left breast
## calcifications 11-12 o'clock position - biopsy not yet performed.
## Post menopausal.   id     T1w             T2w
## 1 629 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 629       629          7189           7068978 2012-02-28 00:00:00.000000                  501      4        massB     FIBROADENOMA
## CLINICAL INDICATION:  Increased in size of the left breast since January and hardening between 3 and 4 o'clock left breast. Recent mammogram and ultrasound negative for malignancy.   id     T1w             T2w
## 1 630 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 630       630          7189           7068978 2012-02-28 00:00:00.000000                 4116      4        massB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  Increased in size of the left breast since January and hardening between 3 and 4 o'clock left breast. Recent mammogram and ultrasound negative for malignancy.   id     T1w             T2w
## 1 631 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label      lesion_diagnosis
## 631       631          7190           7013378 2012-01-10 00:00:00.000000                 4118      3        massB COLUMNAR CELL CHANGES
## CLINICAL INDICATION: Right br microcalcs between 9:30 and 10:00 BX=FEA...dense breast MRI recommended   id     T1w             T2w
## 1 632 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 632       632          7192           7974056 2014-11-02 00:00:00.000000                 4124      4     nonmassB BENIGN BREAST TISSUE
## CLINICAL INDICATION:  6 month followup a focal enhancement in the upper aspect of the right breast.   id     T1w             T2w
## 1 633 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 633       633          7193           7347138 2013-03-09 00:00:00.000000                 4126      4     nonmassB BENIGN BREAST TISSUE
## INDICATION:  Surveillance (OBSP) and follow up   id     T1w             T2w
## 1 634 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label     lesion_diagnosis
## 634       634          7193           7347138 2013-03-09 00:00:00.000000                 4126      4     nonmassB BENIGN BREAST TISSUE
## INDICATION:  Surveillance (OBSP) and follow up   id     T1w             T2w
## 1 635 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 635       635          7201           5041620 2009-08-11 00:00:00.000000                  922      4     nonmassB     FIBROADENOMA
## CLINICAL INDICATION: on going pain and fullness left breast upper
## outer quadrant. No mammographic or ultrasound findings June 2009.   id     T1w             T2w
## 1 636 VIBRANT T2 weighted FSE
##     lesion_id cad_pt_no_txt exam_a_number_txt           exam_dt_datetime proc_pt_procedure_id BIRADS lesion_label lesion_diagnosis
## 636       636          7220           7288789 2013-01-11 00:00:00.000000                 4147      4     nonmassB      FIBROCYSTIC
## CLINICAL INDICATION: High risk  LMP Dec 24/12  Examination is poorly timed to the menstrual cycle.
```

```r
print(summary(factor(seqdf$T1w)))
```

```
##  REVIEW VIBRANT 
##       1     626
```

```r
print(summary(factor(seqdf$T2w)))
```

```
##          REVIEW T2 weighted FSE 
##               5             622
```


