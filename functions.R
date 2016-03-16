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
  gntrees = c(50,100,250,350)
  gminsplit = c(3,1) 
  gcp = c(0.01,-1) 
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
calcAUC_plot <- function(obs, probpred, xptext, yptext, icolors, atitle){
  library(pROC)
  ROC <- plot.roc(obs, 
                  probpred,
                  ci=TRUE, print.auc=TRUE, print.auc.x=xptext, print.auc.y=yptext,
                  col=icolors, lty=1, legacy.axes=TRUE,
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
### code forest Train: 
### parameters, T= # of trees, D= tree depth, dat
###################################################

# Build a nodeHarvest
NH_looforestTrain <- function(TrainsetD, TestsetD) {
  library(klaR)
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
  gnodes = c(25,50,100,150,200,250,400)
  gnodesize = c(5,3,1)
  grd <- expand.grid(maxinter = gmaxinter, nodes = gnodes, nodesize=gnodesize)
  
  # for oneshot
  grdperf = data.frame(grd)
  grdperf$rocTrain =0
  grdperf$rocTest =0
  
  M = list()
  for(k in 1:nrow(grd)){
    Nnodes=grd$nodes[k]
    maxinter=grd$maxinter[k]
    nodesize=grd$nodesize[k]
    
    # Build in 
    cat("nodes : " ,Nnodes, "maxinter: ", maxinter, "nodesize: ", nodesize, "\n")
    
    nodeHarvestfit =  nodeHarvest(TrainsetD[,2:ncol(TrainsetD)], ifelse(TrainsetD[,1]=="C",1,0),
                              maxinter=maxinter,
                              nodes=Nnodes, 
                              mode = "outbag",
                              silent=TRUE, biascorr=FALSE)
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
    M = append(M, list(M=list(nodes=Nnodes,maxinter=maxinter,
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

