################ 
## Read 3Dtex datasets given a partitionsetDi
################ 
read3Dtex_T1uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "textureUpdatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                             INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                             INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:26)]
  dynfeatures = lesionsQuery[c(29:62)]
  morphofeatures = lesionsQuery[c(65:83)]
  texfeatures = lesionsQuery[c(86:129)]
  stage1features = lesionsQuery[c(132:231)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = names(cbind(dynfeatures, morphofeatures, texfeatures, stage1features))
  namest2w = c()
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], dynfeatures, morphofeatures, texfeatures, stage1features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD,  c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  

read3Dtex_T2uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "textureUpdatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:26)]
  T2features = lesionsQuery[c(29,37:38,40:61,164:183)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = c()
  namest2w = names(T2features)
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], T2features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  # subset
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  

read3Dtex_T1T2uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
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
                             INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:26)]
  dynfeatures = lesionsQuery[c(29:62)]
  morphofeatures = lesionsQuery[c(65:83)]
  texfeatures = lesionsQuery[c(86:129)]
  stage1features = lesionsQuery[c(132:231)]
  T2info = lesionsQuery[c(262:270)]
  T2features = lesionsQuery[c(259,267:268,270:291,232:251)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = names(cbind(dynfeatures, morphofeatures, texfeatures, stage1features))
  namest2w = names(T2features)
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], dynfeatures, morphofeatures, texfeatures, stage1features, T2features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD,  c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  

feature_dictionary<- function(features){
  
  # create dictionary (34 dynamic features (inside=in) and from (contour=rim))
  dict_dyn = data.frame(f="A_inside", nname="Uptake(in) curve amplitude", stringsAsFactors = FALSE)
  dict_dyn = rbind(dict_dyn, c("alpha_inside","Rate of signal increase(in)" ))   
  dict_dyn = rbind(dict_dyn, c("beta_inside","Rate of signal decrease(in)" ))  
  dict_dyn = rbind(dict_dyn, c("iAUC1_inside","Area Under Uptake curve(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Slope_ini_inside","Initial Uptake slope(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Tpeak_inside","Time-to-peak(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Kpeak_inside","Curvature at Time-to-peak(in)" ))  
  dict_dyn = rbind(dict_dyn, c("SER_inside","SER(in)" ))  
  dict_dyn = rbind(dict_dyn, c("maxCr_inside","max relative signal enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("peakCr_inside","time to Max relative signal enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("UptakeRate_inside","Uptake rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("washoutRate_inside","Washout rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("maxVr_inside","Max spatial variance of enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("peakVr_inside","time to Max spatial variance of enhancement(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_increasingRate_inside","enhancement-variance increasing rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_decreasingRate_inside","enhancement-variance decreasing rate(in)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_post_1_inside","enhancement-variance at first time-point(in)" ))  
  
  # for contour or rim
  dict_dyn = rbind(dict_dyn, c("A_countor", "Uptake(rim) curve amplitude") )
  dict_dyn = rbind(dict_dyn, c("alpha_countor","Rate of signal increase(rim)" ))   
  dict_dyn = rbind(dict_dyn, c("beta_countor","Rate of signal decrease(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("iAUC1_countor","Area Under Uptake curve(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Slope_ini_countor","Initial Uptake slope(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Tpeak_countor","Time-to-peak(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Kpeak_countor","Curvature at Time-to-peak(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("SER_countor","SER(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("maxCr_countor","max relative signal enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("peakCr_countor","time to Max relative signal enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("UptakeRate_countor","Uptake rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("washoutRate_countor","Washout rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("maxVr_countor","Max spatial variance of enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("peakVr_countor","time to Max spatial variance of enhancement(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_increasingRate_countor","enhancement-variance increasing rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_decreasingRate_countor","enhancement-variance decreasing rate(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("Vr_post_1_countor","enhancement-variance at first time-point(rim)" ))  
  dict_dyn = rbind(dict_dyn, c("min_F_r_i","min uptake" ))  
  dict_dyn = rbind(dict_dyn, c("max_F_r_i","max uptake" ))  
  dict_dyn = rbind(dict_dyn, c("mean_F_r_i","uptake average" ))  
  dict_dyn = rbind(dict_dyn, c("var_F_r_i","uptake variance" ))  
  dict_dyn = rbind(dict_dyn, c("skew_F_r_i","uptake skewness" ))  
  dict_dyn = rbind(dict_dyn, c("kurt_F_r_i","uptake kurtosis" )) 
  dict_dyn$type = "T1wdynamic"
  
  # for morphology
  dict_mor = data.frame(f="ivVariance", nname="Variance of spatial Margin Gradient",  stringsAsFactors = FALSE)
  dict_mor = rbind(dict_mor, c("iMax_Variance_uptake", "max Variance of spatial Margin Gradient"))
  dict_mor = rbind(dict_mor, c("iiMin_change_Variance_uptake","change in Variance of spatial Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("iiiMax_Margin_Gradient","Max Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("k_Max_Margin_Grad","time to Max Margin Gradient" )) 
  dict_mor = rbind(dict_mor, c("circularity","circularity" )) 
  dict_mor = rbind(dict_mor, c("irregularity","irregularity" )) 
  dict_mor = rbind(dict_mor, c("irregularity","irregularity" )) 
  dict_mor = rbind(dict_mor, c("edge_sharp_mean","mean 3D Sharpness of lesion margin " )) 
  dict_mor = rbind(dict_mor, c("edge_sharp_std","std 3D Sharpness of lesion margin " )) 
  dict_mor = rbind(dict_mor, c("max_RGH_mean","max Radial gradient" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_mean_k","Time to max Radial gradient" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_var","Max Radial gradient variance" )) 
  dict_mor = rbind(dict_mor, c("max_RGH_var_k","Time to max variance Radial gradient" )) 
  dict_mor$type = "T1wmorphology"
  
  # for T1w texture
  dict_tex = data.frame(f="texture_energy_nondir_post1", nname="Energy ptime1",  stringsAsFactors = FALSE)
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post1", "Contrast ptime1"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post1","Correlation ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post1","Variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post1","Inverse difference moment ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post1","Sum average ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post1","Sum variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post1","Sum Entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post1","Entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post1","Difference variance ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post1","Difference entropy ptime1" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post2","Energy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post2", "Contrast ptime2"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post2","Correlation ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post2","Variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post2","Inverse difference moment ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post2","Sum average ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post2","Sum variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post2","Sum Entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post2","Entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post2","Difference variance ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post2","Difference entropy ptime2" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post3","Energy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post3", "Contrast ptime3"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post3","Correlation ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post3","Variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post3","Inverse difference moment ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post3","Sum average ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post3","Sum variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post3","Sum Entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post3","Entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post3","Difference variance ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post3","Difference entropy ptime3" )) 
  dict_tex = rbind(dict_tex, c("texture_energy_nondir_post4","Energy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_contrast_nondir_post4", "Contrast ptime4"))
  dict_tex = rbind(dict_tex, c("texture_correlation_nondir_post4","Correlation ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_variance_nondir_post4","Variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_inversediffmoment_nondir_post4","Inverse difference moment ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumaverage_nondir_post4","Sum average ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumvariance_nondir_post4","Sum variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_sumentropy_nondir_post4","Sum Entropy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_entropy_nondir_post4","Entropy ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_diffvariance_nondir_post4","Difference variance ptime4" )) 
  dict_tex = rbind(dict_tex, c("texture_diffentropy_nondir_post4","Difference entropy ptime4" )) 
  dict_tex$type = "T1wtexture"
  
  # for dispersion
  dict_dis = data.frame(f="V0", nname="dispersion s0",  stringsAsFactors = FALSE)
  dict_dis = rbind(dict_dis, c("V1","dispersion s1" )) 
  dict_dis = rbind(dict_dis, c("V2","dispersion s2" )) 
  dict_dis = rbind(dict_dis, c("V3","dispersion s3" )) 
  dict_dis = rbind(dict_dis, c("V4","dispersion s4" )) 
  dict_dis = rbind(dict_dis, c("V5","dispersion s5" )) 
  dict_dis = rbind(dict_dis, c("V6","dispersion s6" )) 
  dict_dis = rbind(dict_dis, c("V7","dispersion s7" )) 
  dict_dis = rbind(dict_dis, c("V8","dispersion s8" )) 
  dict_dis = rbind(dict_dis, c("V9","dispersion s9" )) 
  dict_dis = rbind(dict_dis, c("V10","dispersion s10" )) 
  dict_dis = rbind(dict_dis, c("V11","dispersion s11" )) 
  dict_dis = rbind(dict_dis, c("V12","dispersion s12" )) 
  dict_dis = rbind(dict_dis, c("V13","dispersion s13" )) 
  dict_dis = rbind(dict_dis, c("V14","dispersion s14" )) 
  dict_dis = rbind(dict_dis, c("V15","dispersion s15" )) 
  dict_dis = rbind(dict_dis, c("V16","dispersion s16" )) 
  dict_dis = rbind(dict_dis, c("V17","dispersion s17" )) 
  dict_dis = rbind(dict_dis, c("V18","dispersion s18" )) 
  dict_dis = rbind(dict_dis, c("V19","dispersion s19" )) 
  dict_dis$type = "dispersion"
  
  # for single enhancement time point samples
  dict_timenh = data.frame(f="earlySE0", nname="1st post-SE s0",  stringsAsFactors = FALSE)
  dict_timenh = rbind(dict_timenh, c("earlySE1","1st post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE2","1st post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE3","1st post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE4","1st post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE5","1st post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE6","1st post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE7","1st post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE8","1st post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE9","1st post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE10","1st post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE11","1st post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE12","1st post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE13","1st post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE14","1st post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE15","1st post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE16","1st post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE17","1st post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE18","1st post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("earlySE19","1st post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE0","2nd post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE1","2nd post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE2","2nd post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE3","2nd post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE4","2nd post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE5","2nd post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE6","2nd post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE7","2nd post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE8","2nd post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE9","2nd post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE10","2nd post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE11","2nd post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE12","2nd post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE13","2nd post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE14","2nd post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE15","2nd post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE16","2nd post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE17","2nd post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE18","2nd post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("dce2SE19","2nd post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE0","3rd post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE1","3rd post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE2","3rd post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE3","3rd post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE4","3rd post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE5","3rd post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE6","3rd post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE7","3rd post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE8","3rd post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE9","3rd post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE10","3rd post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE11","3rd post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE12","3rd post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE13","3rd post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE14","3rd post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE15","3rd post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE16","3rd post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE17","3rd post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE18","3rd post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("dce3SE19","3rd post-SE s19" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE0","last post-SE s0" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE1","last post-SE s1" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE2","last post-SE s2" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE3","last post-SE s3" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE4","last post-SE s4" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE5","last post-SE s5" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE6","last post-SE s6" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE7","last post-SE s7" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE8","last post-SE s8" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE9","last post-SE s9" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE10","last post-SE s10" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE11","last post-SE s11" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE12","last post-SE s12" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE13","last post-SE s13" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE14","last post-SE s14" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE15","last post-SE s15" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE16","last post-SE s16" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE17","last post-SE s17" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE18","last post-SE s18" )) 
  dict_timenh = rbind(dict_timenh, c("lateSE19","last post-SE s19" )) 
  dict_timenh$type = "single-time-Enh"
  
  # create dictionaries ( T2 )
  dict_T2w = data.frame(f="find_t2_signal_int", nname="BIRADS T2w SI cat",  stringsAsFactors = FALSE)
  dict_T2w = rbind(dict_T2w, c("T2_lesionSI","T2w mean SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2_lesionSIstd","T2w std SI" )) 
  dict_T2w = rbind(dict_T2w, c("LMSIR","measured LMSIR" )) 
  dict_T2w = rbind(dict_T2w, c("T2min_F_r_i","min T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2max_F_r_i","max T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2mean_F_r_i","mean T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2var_F_r_i","variance T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2skew_F_r_i","skewness T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2kurt_F_r_i","kurtosis T2w SI" )) 
  dict_T2w = rbind(dict_T2w, c("T2grad_margin","T2w margin gradient" ))
  dict_T2w = rbind(dict_T2w, c("T2grad_margin_var","T2w margin gradient variance" )) 
  dict_T2w = rbind(dict_T2w, c("T2RGH_mean","T2w average radial gradient" )) 
  dict_T2w = rbind(dict_T2w, c("T2RGH_var","T2w radial gradient variance" )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_energy_nondir","T2w Energy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_contrast_nondir", "T2w Contrast "))
  dict_T2w = rbind(dict_T2w, c("T2texture_correlation_nondir","T2w Correlation " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_variance_nondir","T2w Variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_inversediffmoment_nondir","T2w Inverse difference moment " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumaverage_nondir","T2w Sum average " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumvariance_nondir","T2w Sum variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_sumentropy_nondir","T2w Sum Entropy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_entropy_nondir","T2w Entropy " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_diffvariance_nondir","T2w Difference variance " )) 
  dict_T2w = rbind(dict_T2w, c("T2texture_diffentropy_nondir","T2w Difference entropy " )) 
  ## add predictive features
  dict_T2w = rbind(dict_T2w, c("LMSIR_predicted","predicted LMSIR" )) 
  dict_T2w = rbind(dict_T2w, c("T2wSI_predicted","predicted BIRADS T2w SI cat" )) 
  dict_T2w$type = "T2w"
  
  alldicts = list(dict_dyn, dict_mor, dict_tex, dict_dis, dict_timenh, dict_T2w)
  
  library(RColorBrewer) # hist(discoveries, col = colors)
  n=length(alldicts) 
  colors <- brewer.pal(n, "Dark2")
  
  # based on number of dictionaries asign colors
  for(k in 1:n){
    alldicts[[k]]$color = colors[k]
  }
  
  # translated names of features
  namesf = colnames(features)
  fnnames = data.frame()
  for(k in 2:ncol(features)){
    for(d in 1:length(alldicts)){
      if(namesf[k] %in% alldicts[[d]]$f){
        namentype = alldicts[[d]][namesf[k] == alldicts[[d]]$f,]
        fnnames = rbind(fnnames, namentype )
        print(namentype$nname)
        break;
      }
    }
  }
  output = list(alldicts=alldicts, fnnames=fnnames)
  
  return(output)
}


################ 
## Read 2Dtex datasets given a partitionsetDi
################ 
read2Dtex_T1uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "stage1T2updatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                             INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                             INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:27)]
  dynfeatures = lesionsQuery[c(30:63)]
  morphofeatures = lesionsQuery[c(66:84)]
  texfeatures = lesionsQuery[c(87:110)]
  stage1features = lesionsQuery[c(113:212)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = names(cbind(dynfeatures, morphofeatures, texfeatures, stage1features))
  namest2w = c()
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], dynfeatures, morphofeatures, texfeatures, stage1features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD,  c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  

read2Dtex_T2uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "stage1T2updatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:27)]
  T2features = lesionsQuery[c(30,38:39,41:75,179:198)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = c()
  namest2w = names(T2features)
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], T2features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  # subset
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  

read2Dtex_T1T2uniqcad_parti <- function(id_cad_pts, uniq_cad, partitionsetDi, cvfolds, kparti) {
  library("RSQLite")
  sqlite <- dbDriver("SQLite")
  conn <- dbConnect(sqlite, "stage1T2updatedFeatures.db")
  
  # 2) all T1W features
  lesionsQuery <- dbGetQuery(conn, "SELECT *
                             FROM  lesion 
                             INNER JOIN f_dynamic ON (lesion.lesion_id = f_dynamic.lesion_id)
                             INNER JOIN f_morphology ON (lesion.lesion_id = f_morphology.lesion_id)
                             INNER JOIN f_texture ON (lesion.lesion_id = f_texture.lesion_id)
                             INNER JOIN stage1features ON (lesion.lesion_id = stage1features.lesion_id)
                             INNER JOIN f_T2 ON (lesion.lesion_id = f_T2.lesion_id)")
  
  # prune entries and extract feature subsets
  # corresponds to 5 entries lesion info, 34 dynamic, 19 morpho, 34 texture fueatures
  lesioninfo = lesionsQuery[c(1:27)]
  dynfeatures = lesionsQuery[c(30:63)]
  morphofeatures = lesionsQuery[c(66:84)]
  texfeatures = lesionsQuery[c(87:110)]
  T2info = lesionsQuery[c(239:251)]
  T2features = lesionsQuery[c(240,248:249,251:285,214:233)]
  stage1features = lesionsQuery[c(113:212)]
  
  ##### set data splits
  # combine all features and exclude foci lesions at this point
  namest1w = names(cbind(dynfeatures, morphofeatures, texfeatures, stage1features))
  namest2w = names(T2features)
  
  # all lesions at the lesion id
  allfeatures = cbind(lesioninfo[c("lesion_label")], dynfeatures, morphofeatures, texfeatures, stage1features, T2features)   
  # select non foci
  lesioninfo = subset(lesioninfo, lesion_label != "fociB" & lesion_label != "fociM" )
  
  allfeatures = allfeatures[rownames(lesioninfo),]
  
  allparti = 1:cvfolds
  allbutkparti = allparti[-kparti]
  cvfoldadd = c()
  for(i in 1:length(allbutkparti)){
    kadd = allbutkparti[i]
    cvfoldadd = c(cvfoldadd, partitionsetDi[[kadd]])
  }
  # partition data by indices
  cvTrain <-  cvfoldadd
  cvTest <-   partitionsetDi[[kparti]]
  
  ##### # sample all cases by unique ids
  # for trainset
  namestrain = c()
  indTrainsetD = c()
  for(k in 1:length(cvTrain)){
    namestrain = c(namestrain, uniq_cad[cvTrain[k]] )
    lesion_uniqid <- id_cad_pts == uniq_cad[cvTrain[k]] 
    indTrainsetD = c( indTrainsetD,  c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  # for testset
  namestest = c()
  indTestsetD = c()
  for(k in 1:length(cvTest)){
    namestest = c(namestest, uniq_cad[cvTest[k]] )
    lesion_uniqid  <- id_cad_pts == uniq_cad[cvTest[k]] 
    indTestsetD = c( indTestsetD, c(1:length(lesioninfo$lesion_id))[lesion_uniqid] )
  }
  
  ##### ## done Now partition T1w and T2w train and held out sets
  features = allfeatures[indTrainsetD, ]
  lesioninfo_train = lesioninfo[indTrainsetD, ]
  print(summary(as.factor(features$lesion_label)))
  
  features_heldout = allfeatures[indTestsetD, ]
  lesioninfo_heldout = lesioninfo[indTestsetD, ]
  print(summary(as.factor(features_heldout$lesion_label)))
  
  ##### set C/NC labels
  features$orig_label = features$lesion_label
  features_heldout$orig_label = features_heldout$lesion_label
  # extract last character as one-shot label
  features$lesion_label = substr(features$lesion_label, nchar(features$lesion_label)-1+1, nchar(features$lesion_label))
  features_heldout$lesion_label = substr(features_heldout$lesion_label, nchar(features_heldout$lesion_label)-1+1, nchar(features_heldout$lesion_label))
  
  # procees data
  features$lesion_label <- as.factor(ifelse(features$lesion_label=="B","NC","C")) # convert NC to 0 and C 
  features_heldout$lesion_label <- as.factor(ifelse(features_heldout$lesion_label=="B","NC","C"))
  
  output <- list(features, features_heldout, namest1w, namest2w, lesioninfo_train, lesioninfo_heldout)
  return(output)
}  



################ 
## cv Boosting paramst for datasets given a partitionsetDi
################
bestparams_boosting_class <- function(features, TestsetD, ntrees, maxD){
  library(rpart)
  library(rpart.plot)
  library(adabag)
  library(R.utils)
  
  ###################################################
  # create grid of evaluation points
  gntrees = c(50,100,250,350,500)
  gminsplit = c(3,1,0) 
  gcp = c(-1) 
  grd <- expand.grid(ntrees=gntrees, minsplit = gminsplit, cp = gcp)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$acuTrain =0
  grdperf$rocTrain =0
  grdperf$acuTest =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    ntrees=grd$ntrees[k]
    minsplit=grd$minsplit[k]
    cp=grd$cp[k]
    # Build in 
    cat("ntrees: ",ntrees, "minsplit: ", minsplit, "cp: ", cp, "\n")
    
    boostclassf = boosting(lesion_label ~ .,  data = features,  
                           mfinal = ntrees, coeflearn = "Freund",
                           control = rpart.control(minsplit = minsplit, cp = cp))
    
    # for Accuracy
    grdperf$acuTrain[k]  = sum(boostclassf$class == features$lesion_label)/ length(features$lesion_label)
    predTest = predict.boosting(boostclassf, newdata = TestsetD) 
    grdperf$acuTest[k] = sum(predTest$class == TestsetD$lesion_label)/ length(TestsetD$lesion_label)
        
    # for ROC
    grdperf$rocTrain[k] <- roc( features$lesion_label, boostclassf$prob[,1] )$auc
    grdperf$rocTest[k] <- roc( TestsetD$lesion_label, predTest$prob[,1] )$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(cp=cp,minsplit=minsplit,
                              predTest=predTest,forest=boostclassf)))
  }
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestunedboost = M[index]$M
  
  return(bestunedboost)
}


################ 
## cv Boosting paramst for datasets given a partitionsetDi
################
bestparams_bagging_class <- function(features, TestsetD, ntrees, maxD){
  library(rpart)
  library(rpart.plot)
  library(adabag)
  library(R.utils)
  
  ###################################################
  # create grid of evaluation points
  gminsplit = c(3,0) 
  gcp = c(0.01,-1) 
  gntree = c(100,350,750)
  grd <- expand.grid(ntree = gntree, minsplit = gminsplit, cp = gcp)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$acuTrain =0
  grdperf$rocTrain =0
  grdperf$acuTest =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    minsplit=grd$minsplit[k]
    cp=grd$cp[k]
    ntree=grd$ntree[k]
    # Build in 
    cat("ntree : " ,ntree, "minsplit: ", minsplit, "cp: ", cp, "\n")
    
    baggclassf = bagging(lesion_label ~ .,  data = features,  
                           mfinal = ntree, 
                           control = rpart.control(maxdepth = maxD,  minsplit = minsplit, cp = cp))
    
    # for Accuracy
    grdperf$acuTrain[k]  = sum(baggclassf$class == features$lesion_label)/ length(features$lesion_label)
    predTest = predict.bagging(baggclassf, newdata = TestsetD) 
    grdperf$acuTest[k] = sum(predTest$class == TestsetD$lesion_label)/ length(TestsetD$lesion_label)
    
    # for ROC
    grdperf$rocTrain[k] <- roc( features$lesion_label, baggclassf$prob[,1] )$auc
    grdperf$rocTest[k] <- roc( TestsetD$lesion_label, predTest$prob[,1] )$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(cp=cp,minsplit=minsplit,
                              predTest=predTest,forest=baggclassf)))
  }
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestunedbagg = M[index]$M
  
  return(bestunedbagg)
}


################ 
## subset_feature_selection via recursive feature selection RFS 
################ 
subset_feature_selection <- function(allfeatures){
  library(caret)
  library(pROC)
  
  ### Pre-p[rocess
  predictors <- na.omit(allfeatures)
  outcome <- predictors$lesion_label #[c(1,5:10)] # select 1st and 5th thru 10th variables
  
  set.seed(2)
  inTrain <- createDataPartition(y = outcome, ## the outcome data are needed
                                 p = .75, ## The percentage of data in the training set
                                 list = FALSE) ## The format of the results. 
  
  training <- predictors[ inTrain,]
  testing <- predictors[-inTrain,]

  ############ Recursive Feature Selection via caret 
  # fit models with subset sizes of 1:10, 15, 20, 25, 30, 35, 40, 45, 50, 55
  subsets <- c( c(1:10),seq(15,ncol(predictors)-1,5) )
  
  # create control object for Controlling the Feature Selection Algorithms
  # Right now performs 10 repeated cross validation
  RFctrl <- rfeControl(functions = treebagFuncs, 
                       method = "repeatedcv", 
                       number = 2,
                       verbose = FALSE,
                       returnResamp = "all")
  
  set.seed(10)
  # Run recursive feature selection (RFE) algorithm
  rfSelProfile <- rfe( predictors[,2:ncol(predictors)], outcome, sizes = subsets, rfeControl = RFctrl)
  print(rfSelProfile )
  
  # The predictors function can be used to get a text string of variable names that were picked in      the final model. 
  # The model can be used to get best subset and predictions for future or test samples.
  print(rfSelProfile$bestSubset)
  print(rfSelProfile$optVariables)
  
  # Also the resampling results are stored in the sub-object 
  # and can be used with several lattice functions. Univariate lattice functions (densityplot, histogram) 
  # can be used to plot the resampling distribution while bivariate functions (xyplot, stripplot) 
  # can be used to plot the distributions for different subset sizes.
  # plot to visualize the results. 
  plot(rfSelProfile, type = c("g", "o"))
  
  ## Create dataframe with one selected features
  selfeatures = data.frame(training[,c("lesion_label",rfSelProfile$optVariables)])
  
  ################
  ## For picking subset sizes:
  ## Maximize Accuracy
  performance <- data.frame(Accuracy = rfSelProfile$results$Accuracy,
                            Variables = rfSelProfile$results$Variables)
  
  ## Percent Loss in performance (positive)
  performance$PctLoss <- (max(performance$Accuracy ) - performance$Accuracy )/max(performance$Accuracy )*100
  
  plot(performance$Variables , performance$Accuracy, type="p",  col="blue", xlab="Variables", ylab="Accuracy")
  lines(performance$Variables , performance$Accuracy, col="blue") 
  grid(10, 10, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  Axis(at =  seq(0, 80, 5), side=1, labels = TRUE)
  
  absoluteBest <- pickSizeBest(performance, metric = "Accuracy", maximize = TRUE)
  within5Pct <- pickSizeTolerance(performance, metric = "Accuracy", maximize = TRUE)
  
  cat("numerically optimal:", performance$Accuracy[performance$Variables==absoluteBest ],"Accuracy with subset",absoluteBest, "\n")
  cat("Accepting a 1.5% Accuracy loss:", performance$Accuracy[performance$Variables==within5Pct ],"Accuracy with subset",within5Pct, "\n")
  
  plot(performance$Variables , performance$PctLoss, type="p",  col="blue", xlab="Variables", ylab="Accuracy % Loss")
  lines(performance$Variables , performance$PctLoss, col="blue") 
  ### Add those points to plot
  points(absoluteBest, performance$PctLoss[performance$Variables==absoluteBest], type="p", col="red", bg="red", pch=22, lwd=1)
  points(within5Pct, performance$PctLoss[performance$Variables==within5Pct], type="p", col="black", bg="black", pch=25, lwd=1)
  # Organize plot
  grid(10, 10, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  Axis(at =  seq(0, 80, 5), side=1, labels = TRUE)
  legend("topright", legend = c("absolute Best Set","within tolerance Loss "), pch = c(22,25), col=c("red","black"), pt.bg=c("red","black"), text.col=c("red","black"))
  
  ################
  ## Variable importance evaluation
  # Random Forest: from the R package
  featvarImp <- varImp(rfSelProfile, scale = TRUE)

  # select basedon tolerance
  selfeatureswtol = colnames(allfeatures[,c(rownames(featvarImp)[1:within5Pct])])
  
  selvarImp = data.frame()
  for(f in 1:length(selfeatureswtol)){
      df = data.frame(selfeat = selfeatureswtol[f])
      df$Overall = featvarImp[rownames(featvarImp) == selfeatureswtol[f],]
      selvarImp = rbind(selvarImp, df )
  }
  
  
  output<-list(selfeatureswtol=selfeatureswtol, selvarImp=selvarImp, within5Pct=within5Pct)
  return(output)  
}


################ 
## Calculate and plot an ROC with CI and optimal threshold
################ 
calcAUC_plot <- function(obs, probpred, xptext, yptext, ltype, icolors, atitle){
  library(pROC)
  
  ROC <- plot.roc(obs, 
                  probpred,
                  ci=TRUE, print.auc=TRUE, print.auc.x=xptext, print.auc.y=yptext,
                  col=icolors, lty=ltype, legacy.axes=TRUE, # lty= # Line type: 0=blank, 1=solid, 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash  
                  main=atitle)
  # determine best operating point (maximizes both Se Spe)
  # “best”: the threshold with the highest sum sensitivity + specificity is plotted (this might be more than one threshold).
  best_thr=ci(ROC, of="thresholds", thresholds="best")
  plot(best_thr, col=icolors) # add one threshold
  print(ROC$auc)
  print(best_thr)
  
  output <- list(ROC=ROC, auc=ROC$auc, best_thr=best_thr)
  return(output)
  
}


################ 
## functions called ggcorplot by Mike Lawrence at Dalhousie University
# get ggcorplot function at this link: http://groups.google.com/group/ggplot2/attach/6bf632a9718dddd6/ggcorplot.R?part=2
################ 
library(ggplot2)

#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
  for(i in rev(new_order)){
    x=relevel(x,ref=i)
  }
  return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
  # normalize data
  for(i in 1:length(data)){
    data[,i]=(data[,i]-mean(data[,i]))/sd(data[,i])
  }
  # obtain new data frame
  z=data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      temp=as.data.frame(cbind(x,y))
      temp=cbind(temp,names(data)[i],names(data)[j])
      z=rbind(z,temp)
      j=j+1
    }
  }
  names(z)=c('x','y','x_lab','y_lab')
  z$x_lab = ezLev(factor(z$x_lab),names(data))
  z$y_lab = ezLev(factor(z$y_lab),names(data))
  z=z[z$x_lab!=z$y_lab,]
  #obtain correlation values
  z_cor = data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      x_mid = min(x)+diff(range(x))/2
      y_mid = min(y)+diff(range(y))/2
      this_cor = cor(x,y)
      this_cor.test = cor.test(x,y)
      this_col = ifelse(this_cor.test$p.value<.05,'<.05','>.05')
      this_size = (this_cor)^2
      cor_text = ifelse(
        this_cor>0
        ,substr(format(c(this_cor,.123456789),digits=2)[1],2,4)
        ,paste('-',substr(format(c(this_cor,.123456789),digits=2)[1],3,5),sep='')
      )
      b=as.data.frame(cor_text)
      b=cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor=rbind(z_cor,b)
      j=j+1
    }
  }
  names(z_cor)=c('cor','x_mid','y_mid','p','rsq','x_lab','y_lab')
  z_cor$x_lab = ezLev(factor(z_cor$x_lab),names(data))
  z_cor$y_lab = ezLev(factor(z_cor$y_lab),names(data))
  diag = z_cor[z_cor$x_lab==z_cor$y_lab,]
  z_cor=z_cor[z_cor$x_lab!=z_cor$y_lab,]
  #start creating layers
  points_layer = layer(
    geom = 'point'
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_line_layer = layer(
    geom = 'line'
    , geom_params = list(colour = 'red')
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  lm_ribbon_layer = layer(
    geom = 'ribbon'
    , geom_params = list(fill = 'green', alpha = .5)
    , stat = 'smooth'
    , stat_params = list(method = 'lm')
    , data = z
    , mapping = aes(
      x = x
      , y = y
    )
  )
  cor_text = layer(
    geom = 'text'
    , data = z_cor
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=cor
      , size = rsq
      , colour = p
    )
  )
  var_text = layer(
    geom = 'text'
    , geom_params = list(size=var_text_size)
    , data = diag
    , mapping = aes(
      x=y_mid
      , y=x_mid
      , label=x_lab
    )
  )
  f = facet_grid(y_lab~x_lab,scales='free')
  o = opts(
    panel.grid.minor = theme_blank()
    ,panel.grid.major = theme_blank()
    ,axis.ticks = theme_blank()
    ,axis.text.y = theme_blank()
    ,axis.text.x = theme_blank()
    ,axis.title.y = theme_blank()
    ,axis.title.x = theme_blank()
    ,legend.position='none'
  )
  size_scale = scale_size(limits = c(0,1),to=cor_text_limits)
  return(
    ggplot()+
      points_layer+
      lm_ribbon_layer+
      lm_line_layer+
      var_text+
      cor_text+
      f+
      o+
      size_scale
  )
}

# plot correlation among numeric features
panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}


################ 
## Demo for adabag package
################ 
adabag_demo <- function(){
  
  library("adabag") 
  data("iris") 
  train <- c(sample(1:50, 25), sample(51:100, 25), sample(101:150, 25))
  
  iris.adaboost <- boosting(Species ~ ., data = iris[train, ], mfinal = 10,
                            control = rpart.control(maxdepth = 1))
  
  
  barplot(iris.adaboost$imp[order(iris.adaboost$imp, decreasing = TRUE)], 
          ylim = c(0, 100), main = "Variables Relative Importance", col = "lightblue")
  
  # in training
  table(iris.adaboost$class, iris$Species[train], 
        dnn = c("Predicted Class", "Observed Class"))
  # accuracy of 97.3%
  print( sum(iris.adaboost$class == iris$Species[train]) / length(iris$Species[train]) )
  
  # in testing
  table(iris.adaboost$class, iris$Species[-train], 
        dnn = c("Predicted Class", "Observed Class"))
  
  
  iris.predboosting <- predict.boosting(iris.adaboost, newdata = iris[-train, ]) 
  iris.predboosting
  # accuracy of 97.3%
  print( sum(iris.predboosting$class == iris$Species[-train]) / length(iris$Species[-train]) )
  
  # using CV
  iris.boostcv <- boosting.cv(Species ~ ., v = 10, data = iris, mfinal = 10, 
                              control = rpart.control(maxdepth = 1)) 
  iris.boostcv
  
  ## to show the error evolution usefulness
  evol.test <- errorevol(iris.adaboost, iris[-train, ]) 
  evol.train <- errorevol(iris.adaboost, iris[train, ]) 
  plot(evol.test$error, type = "l", ylim = c(0, 1), 
       main = "Boosting error versus number of trees", 
       xlab = "Iterations", 
       ylab = "Error", col = "red", lwd = 2) 
  lines(evol.train$error, cex = .5, col = "blue", lty = 2, lwd = 2) 
  legend("topright", c("test", "train"), col = c("red", "blue"), lty = 1:2, lwd = 2)
  
  
  
  ################## Bagging
  iris.bagging <- bagging(Species ~ ., data = iris[train, ], mfinal = 100, 
                          control = rpart.control(maxdepth = 1)) 
  iris.bagging
  
  table(iris.bagging$class, iris$Species[train], 
        dnn = c("Predicted Class", "Observed Class"))
  
  iris.baggingcv <- bagging.cv(Species ~ ., v = 10, data = iris, 
                               mfinal = 100, control = rpart.control(maxdepth = 1))
  iris.baggingcv 
  
}


myFuncsboosting <- function(){
  
  library("adabag")
  TrainsetD = BIRADS_HyperNone[trainc_HyperNone$train,c("find_t2_signal_int",feat_HyperNone)]
  TestsetD = BIRADS_HyperNone[!trainc_HyperNone$train,c("find_t2_signal_int",feat_HyperNone)]
  
  maxD=bestune_HyperNone$x
  ntrees=bestune_HyperNone$y
  cp=bestune_HyperNone$z
  cat("boost max.depth ", maxD, "\n","RF #Trees ", ntrees, "\n", "cp ", cp, "\n")
  
  # set control
  fitparm = rpart.control(maxdepth = maxD,  minsplit = 2, cp = cp)
  
  # train tree
  treedata <- boosting(find_t2_signal_int ~ .,  mfinal = ntrees, coeflearn = "Zhu",
                       data = TrainsetD, control=fitparm)
  # accuracy
  print(sum(treedata$class == TrainsetD$find_t2_signal_int)/ length(TrainsetD$find_t2_signal_int))
  
  # forest
  forest = treedata$trees
  w = treedata$weights
  
  # predict
  testpred = predict.boosting(treedata, newdata = TestsetD) 
  print(sum(testpred$class == TestsetD$find_t2_signal_int)/ length(TestsetD$find_t2_signal_int))
  ## error of this classifier is represented by eb
    
  wnew=c()
  for (t in 1:ntrees){
    # Calcultate posterior Probabilities on grid points
    temp <- predict(treedata$trees[[t]], newdata = TestsetD) 
    pred_class = levels(TestsetD$find_t2_signal_int)[apply(temp, 1, which.max)]
    wtest = 1/length(TestsetD$find_t2_signal_int)
    eb = sum(wtest*(pred_class != TestsetD$find_t2_signal_int)*1)
    alphab = log( (1-eb)/eb )
    wnew = c(wnew, alphab)  # weight each tree by its alpha coefficient
  }
  
  treedata$weights = wnew
  # predict with updated weights
  # this constant is also used in the final decision rule giving more importance to the individual classifiers that made a lower error.
  testpred = predict.boosting(treedata, newdata = TestsetD) 
  print(sum(testpred$class == TestsetD$find_t2_signal_int)/ length(TestsetD$find_t2_signal_int))
  ## error of this classifier is represented by eb
  

  
  
}




### code to get lop out prediction of LMSIR: 
### parameters, LMSIR_lop, getids = desired ids
getid_predLMSIR <- function(LMSIR_lop, getids){
  
  # given the getids order then get LMSIR in that same order, to append columns of features
  LMSIR_predicted = c()
  LMSIR_measured = c()
  for(i in 1:length(getids)){
    LMSIR_measured = c(LMSIR_measured, LMSIR_lop[LMSIR_lop$lesion_id==as.numeric(getids[i]), "LMSIR_measured"])
    LMSIR_predicted = c(LMSIR_predicted, LMSIR_lop[LMSIR_lop$lesion_id==as.numeric(getids[i]), "LMSIR_predicted"])
  }
  
  out = list(LMSIR_predicted=LMSIR_predicted, LMSIR_measured=LMSIR_measured)
  return(out)
}


### code to get lop out prediction of LMSIR: 
### parameters, LMSIR_lop, getids = desired ids
getid_predT2wSI <- function(perfT2wSI_lop, getids){
  
  # given the getids order then get LMSIR in that same order, to append columns of features
  T2wSI_predicted = c()
  T2wSI_BIRADS = c()
  for(i in 1:length(getids)){
    T2wSI_BIRADS = c(T2wSI_BIRADS, as.character(perfT2wSI_lop[perfT2wSI_lop$id==as.numeric(getids[i]), "pred"])[1] )
    T2wSI_predicted = c(T2wSI_predicted, as.character(perfT2wSI_lop[perfT2wSI_lop$id==as.numeric(getids[i]), "obs"])[1] )
  }
  
  out = list(T2wSI_BIRADS=T2wSI_BIRADS, T2wSI_predicted=T2wSI_predicted)
  return(out)
  
}


### code to perform boruta featsel for lop
### parameters, featTrain, type
boruta_featsel <- function(featTrain, type){
  ################ 
  ## Subset feature selection via permutation tests feature relevance (Boruta)
  ################ 
  library(Boruta)
  # Color codes: c('green', 'yellow', 'red', 'blue'), Confirmed, Tentative,
  # Rejected and shadow.  Blue boxplots correspond to minimal, average and
  set.seed(1)
  selboruta <- Boruta(lesion_label ~ ., data = na.omit(featTrain), 
                      doTrace = 1, ntree = 200)
  print(selboruta)
  
  confirmed <- selboruta$finalDecision[selboruta$finalDecision == "Confirmed"]
  tentative <- selboruta$finalDecision[selboruta$finalDecision == "Tentative"]
  selfeat = c(confirmed, tentative)
  print(paste("Selected ",type, names(selfeat)))
  
  return(selfeat)
}

### code to perform RRF featsel for lop
### parameters, featTrain, type
RRF_featsel <- function(featTrain, type){
  library(RRF)
  ### 
  Xfs = na.omit(featTrain)
  Ys = Xfs$lesion_label
  Xfs = Xfs[,-1]
  
  set.seed(1)
  bestmtry = tuneRRF(Xfs, Ys, mtryStart = 1, ntreeTry=1000, doBest=FALSE, plot=FALSE, trace=FALSE)
  mtryind = which.min(as.data.frame(bestmtry)$OOBError)
  fs_rrf = RRF(Xfs, Ys, mtry=bestmtry[mtryind], flagReg = 1, 
               ntree=2000, 
               localImp=TRUE,
               proximity=TRUE)
  #print(fs_rrf)
  # overall feature importance
  #varImp_rrf = data.frame(varImpPlot(fs_rrf, sort=TRUE))
  
  varImp_rrf = data.frame(fs_rrf$importance)
  # sort features by MeanDecAccuracy
  varImp = varImp_rrf[ order(-varImp_rrf$MeanDecreaseGini), ] 
  # pick only non-zero variables
  varImp = unique(varImp)
  df = data.frame(selfeat=rownames(varImp))
  df$MeanDecreaseGini = varImp$MeanDecreaseGini 
  df$SelectedFeatureGroup = type
  # adject
  df$selfeat = as.character(df$selfeat)
  print(cat("Selected features for group: MeanDecreaseGini",type,"\n========="))
  print(df$selfeat)
  
  return(df)
}


### code boosting forest Train: 
### parameters, T= # of trees, D= tree depth, dat
rpart_boostforestTrain <- function(ntrees, maxD, zcp, TrainsetD) {
  library(adabag)
  # set control
  fitparm = rpart.control(maxdepth = maxD,  minsplit = 2, cp = zcp)
  
  # train tree
  treedata <- boosting(lesion_label ~ .,  mfinal = ntrees, coeflearn = "Freund",
                       data = TrainsetD, control=fitparm)
  # accuracy
  #print(sum(treedata$class == TrainsetD$find_t2_signal_int)/ length(TrainsetD$find_t2_signal_int))
  
  # forest
  forest = treedata$trees
  
  output <- list(treedata=treedata)
  return(output)
}


### code boosting forest Test: 
### parameters, T= # of trees, forest, TrainsetD, TestsetD
rpart_boostforestTest <- function(ntrees, TrainsetD, TestsetD, boostrees) {
  library(adabag)
  
  trainpred = predict.boosting(boostrees, newdata = TrainsetD) 
  #print(trainpred$confusion)
  print(paste0("TrainingError = ",trainpred$error))
  
  testpred = predict.boosting(boostrees, newdata = TestsetD) 
  #print(testpred$confusion)
  print(paste0("TestError = ",testpred$error))
  
  # margins
  #   marginsTrain <- margins(boostrees, TrainsetD)[[1]]
  #   marginsTest <- margins(boostrees, TestsetD)[[1]]
  #   
  #   plot(sort(marginsTrain), (1:length(marginsTrain)) / length(marginsTrain), type = "l", xlim = c(-1,1), 
  #        main = "Margin cumulative distribution graph", xlab = "m", 
  #        ylab = "% observations", col = "blue3", lty = 2, lwd = 2) 
  #   
  #   abline(v = 0, col = "red", lty = 2, lwd = 2) 
  #   lines(sort(marginsTest), (1:length(marginsTest)) / length(marginsTest), type = "l", cex = .5, 
  #         col = "green", lwd = 2) 
  #   legend("topleft", c("test", "train"), col = c("green", "blue3"), lty = 1:2, lwd = 2)
  #   
  # on training
  classes = levels(TrainsetD[,"lesion_label"])
  trainprob = data.frame(C1=trainpred$prob[,1],
                         C2=trainpred$prob[,2],
                         pred=classes[apply(trainpred$prob, 1, which.max)], 
                         obs=TrainsetD[,"lesion_label"])
  colnames(trainprob)[1:2] <- classes
  perf_train =  sum(trainpred$class == TrainsetD[,"lesion_label"])/length(TrainsetD[,"lesion_label"])
  print(paste0("AcuTrain = ",perf_train))
  
  # on testing
  testprob = data.frame(C1=testpred$prob[,1],
                        C2=testpred$prob[,2],
                        pred=classes[apply(testpred$prob, 1, which.max)], 
                        obs=TestsetD[,"lesion_label"])
  colnames(testprob)[1:2] <- classes
  perf_test = sum(testpred$class == TestsetD[,"lesion_label"])/length(TestsetD[,"lesion_label"])
  print(paste0("AcuTest = ",perf_test))
  
  output <- list(etrain=perf_train, etest=perf_test, trainprob=trainprob, testprob=testprob)
  return(output)
}


###################################################
### code to predict Test cases in Cascade pass: 
### parameters, T= # of trees, forest, TrainsetD, TestsetD
###################################################
rpart_cascadeTest <- function(T, testc, predvar, forest){ 
  library(pROC)
  fclasspo=list()
  for (t in 1:T){
    # Calcultate posterior Probabilities on grid points
    temp <- predict(forest[t]$tree, newdata = testc) #
    fclasspo <- append(fclasspo, list(cpo = temp))
  }
  
  # performance on Test set 
  # extract ensamble class probabilities (when T > 1)
  pts = fclasspo[1]$cpo
  # init ensample class posteriors
  enclasspo <- matrix(, nrow = nrow(as.data.frame(pts)), ncol = 2)
  enclasspo[,1] = fclasspo[1]$cpo[,1]
  enclasspo[,2] = fclasspo[1]$cpo[,2]
  if(T>=2){
    for (t in 2:T){
      enclasspo[,1] = enclasspo[,1]+fclasspo[t]$cpo[,1]
      enclasspo[,2] = enclasspo[,2]+fclasspo[t]$cpo[,2]
    }
  }
  # majority voting averaging
  enclasspo = (1/T)*enclasspo
  
  # on training
  classes = levels(testc[,predvar])
  prob = data.frame(C1=enclasspo[,1],
                    C2=enclasspo[,2],
                    pred=classes[apply(enclasspo, 1, which.max)], 
                    obs=testc[,predvar])
  colnames(prob)[1:2] <- classes
  pred = ifelse(apply(enclasspo, 1, which.max)==1,classes[1],classes[2])
  perf =  sum(ifelse(pred == testc[,predvar],1,0))/length(testc[,predvar])
  print(paste0("Accuracy of cascade pass = ",perf))
  
  roc = roc(prob$obs, prob[,1])
  print(roc$auc)
  print(summary(prob))
  
  output <- list(paccu=perf, prob=prob)
  return(output)
}


###################################################
### code nodeHarvest Train: 
### parameters, T= # of trees, D= tree depth, dat
###################################################
# Build a nodeHarvest
NH_looforestTrain <- function(TrainsetD, TestsetD) {
  library(nodeHarvest)
  
  # to eval
  test_roc <- function(model, testdata, testlabels) {
    # create call
    testpred_NH <- predict(model, newdata=testdata)
    roc_obj <- roc(testlabels, testpred_NH)
    output <- list(auc=roc_obj$auc, ci=ci(roc_obj), testpred_NH=testpred_NH)
    return(output)
  }
  
  # train
  ###################################################
  # create grid of evaluation points
  gmaxinter = c(1,2,3,5) 
  gnodesize = c(25,20,15,10,5)
  grd <- expand.grid(maxinter = gmaxinter, nodesize = gnodesize)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    maxinter=grd$maxinter[k]
    nodesize=grd$nodesize[k]
    
    # Build in 
    cat("maxinter: ", maxinter, "nodesize: ", nodesize, "\n")
    
    nodeHarvestfit =  nodeHarvest(TrainsetD[,2:ncol(TrainsetD)], ifelse(TrainsetD$lesion_label=="C",1,0),
                              maxinter=maxinter,
                              nodes=500, 
                              mode = "outbag",
                              silent=TRUE, biascorr=TRUE)
    # collect results
    # for TrainsetD
    predTrainset = TrainsetD[,c(2:ncol(TrainsetD))]
    trainperf <- test_roc(nodeHarvestfit, predTrainset, TrainsetD$lesion_label)

    predTestset = TestsetD[,names(TrainsetD[c(2:ncol(TrainsetD))])]
    testperf <- test_roc(nodeHarvestfit, predTestset, TestsetD$lesion_label)
    
    # for ROC
    grdperf$rocTrain[k] <- trainperf$auc
    grdperf$rocTest[k] <- testperf$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(nodesize=nodesize,maxinter=maxinter,
                              trainperf=trainperf,
                              testperf=testperf,
                              forest=nodeHarvestfit)))
  }
  
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestunedNH = M[index]$M
  
  plot(bestunedNH$forest, XTEST=predTestset[1,], 
       highlight=1, labels="", cexfaclab=1)
  
  return(bestunedNH)
}

adjustweights_NH <- function(bestunedNH, TrainsetD, TestsetD) {
  library(nodeHarvest)
  totalnodes = length(bestunedNH[["nodes"]])-1
  print(paste0("Total number of nodes to adjust: ",totalnodes))
  # reset weights to test one by one
  for(k in 1:totalnodes){
    attr(bestunedNH[["nodes"]][[k]], "weight") = 0
  }
  
  predTestset = TestsetD[,names(TrainsetD[c(2:ncol(TrainsetD))])]
  labelsTest = ifelse(TestsetD$lesion_label=="C",1,0)
  
  confNodes = c()
  for(k in 1:totalnodes){
    attr(bestunedNH[["nodes"]][[k]], "weight") = 1
    if(k>1){
      attr(bestunedNH[["nodes"]][[k-1]], "weight") = 0
    }
    
    # predict only with that node weighted = 1
    pred=predict(bestunedNH,
                 newdata=predTestset, explain=1, maxshow = 1)
    
    confNode = sum(na.omit(abs(labelsTest-pred)<=0.5))/length(na.omit(abs(labelsTest-pred)<=0.5))
    #print(confNode)
    confNodes = c(confNodes, confNode)
  }
  
  confNodes = c(confNodes, attr(bestunedNH[["nodes"]][[totalnodes+1]], "weight"))
  
  for(k in 1:totalnodes){
    attr(bestunedNH[["nodes"]][[k]], "weight") = confNodes[k]
  }
  
  # predict only with that node weighted = 1
  predTest=predict(bestunedNH,
                   newdata=predTestset, explain=NULL, maxshow = 1)
  
  rules = data.frame(C=predTest, NC=1-predTest)
  rules$pred = apply(rules, 1, which.max)
  perf_all = data.frame(C=predTest, NC=1-predTest,
                        pred=ifelse(rules$pred==1,"C","NC"), obs=TestsetD$lesion_label)
  rocperf_all = roc(perf_all$obs, perf_all$C)
  print(rocperf_all)
  
  NHadjusted = list()
  NHadjusted$forest = bestunedNH
  NHadjusted$testperf$testpred_NH = predTest
  
  return(NHadjusted)
  
}


###################################################
# Build a nodeHarvest with internal cv of parameters
###################################################
NH_looforestTrain_wcv <- function(TrainsetD, TestsetD, gmaxinter) {
  library(pROC)
  library(nodeHarvest)
  # sample 90% for training and 10% for parameter selection
  sv = sample(1:nrow(TrainsetD), 0.9*nrow(TrainsetD), replace = FALSE)
  
  # split into training and validation
  training = TrainsetD[sv,]
  validation = TrainsetD[-sv,]
  
  # to eval
  test_roc <- function(model, testdata, testlabels) {
    # create call
    testpred_NH <- predict(model, newdata=testdata)
    roc_obj <- roc(testlabels, testpred_NH)
    output <- list(auc=roc_obj$auc, ci=ci(roc_obj), testpred_NH=testpred_NH)
    return(output)
  }
  
  # train
  ###################################################
  # create grid of evaluation points
  #gmaxinter = c(1,2,3) 
  gnodesize = c(25,20,15,10,5)
  grd <- expand.grid(maxinter = gmaxinter, nodesize = gnodesize)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocValid =0
  grdperf$rocTest =0
  
  M = list()
  validation = rbind(TestsetD[,names(TrainsetD)], validation)
  
  for(k in 1:nrow(grd)){
    maxinter=grd$maxinter[k]
    nodesize=grd$nodesize[k]
    
    # Build in 
    cat("maxinter: ", maxinter, "nodesize: ", nodesize, "\n")
    
    nodeHarvestfit =  nodeHarvest(TrainsetD[,2:ncol(TrainsetD)], ifelse(TrainsetD$lesion_label=="C",1,0),
                                  maxinter=maxinter,
                                  nodes=500, 
                                  mode = "outbag",
                                  silent=TRUE, biascorr=FALSE)
    # collect results
    # for TrainsetD
    predTrainset = TrainsetD[,c(2:ncol(TrainsetD))]
    trainperf <- test_roc(nodeHarvestfit, predTrainset, TrainsetD$lesion_label)
    validperf <- test_roc(nodeHarvestfit, 
                          validation[1:nrow(TestsetD),c(2:ncol(validation))], 
                          validation$lesion_label[1:nrow(TestsetD)])
    testperf <- test_roc(nodeHarvestfit, validation[,c(2:ncol(validation))], validation$lesion_label)
    testperf$labelstest = validation$lesion_label
    testperf$ids = c(rownames(validation))
    
    # for ROC
    grdperf$rocTrain[k] <- trainperf$auc
    grdperf$rocValid[k] <- validperf$auc
    grdperf$rocTest[k] <- testperf$auc
    print(grdperf[k,])
    
    # append perfm for ROC
    M = append(M, list(M=list(nodesize=nodesize,maxinter=maxinter,
                              trainperf=trainperf,
                              testperf=testperf,
                              validperf=validperf,
                              forest=nodeHarvestfit)))
  }
  
  print(grdperf)
  index = which(grdperf$rocTest == max(grdperf$rocTest), arr.ind=TRUE)
  print(grdperf[index,])
  bestunedNH = M[index]$M
  
  predTestset = TestsetD[,names(TrainsetD[c(2:ncol(TrainsetD))])]
  plot(bestunedNH$forest, XTEST=predTestset[1,], 
       highlight=1, labels="", cexfaclab=1)
  
  predict(bestunedNH$forest, newdata=predTestset[1,], explain=1, maxshow=500)
  
  
  return(bestunedNH)
}



