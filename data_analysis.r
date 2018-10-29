####################################################################################################################
# Modeling of ADME relevant endpoints from AstraZeneca data set: A Case study
# Author: Swapnil Chavan
# Date: 27/10/2018

####################################################################################################################
# Load required libraries

all_required_packages <- c("data.table","plyr","RSQLite","caret","DAAG","randomForest","e1071","FNN",
                           "miscTools","ggplot2","e1071","ggpubr","hydroGOF")
for (package in all_required_packages){
  if(!require(package, character.only=T, quietly=T)){
    install.packages(package)
  }
}

####################################################################################################################
# Import from .txt

#my_data <- fread('C:/Swapnil/AZ_case_study/bioactivity-18_12_03_49.txt')

#OR

# Import from SQLITE
# In case one doesn't want to download chembl data, they can directly load 'chembl_data.Rdata' file provided in 
# repository

cur_dir <- getwd() 
url_of_file <- c('http://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_24_1/chembl_24_1_sqlite.tar.gz')
download.file(url_of_file, 'chembl_db_downloaded')
untar("chembl_db_downloaded")

sqlite    <- dbDriver("SQLite")
file_full_name <- paste(cur_dir,'/chembl_24/chembl_24_sqlite/chembl_24.db',sep = '')
exampledb <- dbConnect(sqlite,file_full_name)

# -------------------------- Extract data using MYSQL command -----------------------------

memory.limit(size=56000)

# 1.logD
res <- dbSendQuery(exampledb, "SELECT * FROM ACTIVITIES AS a 
                   LEFT OUTER JOIN ASSAYS AS b ON a.assay_id=b.assay_id 
                   LEFT OUTER JOIN ASSAY_PARAMETERS AS c ON b.assay_id=c.assay_id
                   LEFT OUTER JOIN COMPOUND_PROPERTIES AS d ON a.molregno=d.molregno
                   LEFT OUTER JOIN COMPOUND_STRUCTURES AS e ON a.molregno=e.molregno
                   LEFT OUTER JOIN MOLECULE_DICTIONARY AS f ON a.molregno=f.molregno
                   WHERE b.chembl_id=='CHEMBL3301363'")
chembl_logd <- dbFetch(res)
dbClearResult(res)
chembl_logd <- chembl_logd[,colSums(is.na(chembl_logd))<nrow(chembl_logd)]

# 2.clint
res <- dbSendQuery(exampledb, "SELECT * FROM ACTIVITIES AS a 
                   LEFT OUTER JOIN ASSAYS AS b ON a.assay_id=b.assay_id 
                   LEFT OUTER JOIN ASSAY_PARAMETERS AS c ON b.assay_id=c.assay_id
                   LEFT OUTER JOIN COMPOUND_PROPERTIES AS d ON a.molregno=d.molregno
                   LEFT OUTER JOIN COMPOUND_STRUCTURES AS e ON a.molregno=e.molregno
                   LEFT OUTER JOIN MOLECULE_DICTIONARY AS f ON a.molregno=f.molregno
                   WHERE b.chembl_id=='CHEMBL3301370'")
chembl_clint <- dbFetch(res)
dbClearResult(res)
chembl_clint <- chembl_clint[,colSums(is.na(chembl_clint))<nrow(chembl_clint)]

# 3.hppb
res <- dbSendQuery(exampledb, "SELECT * FROM ACTIVITIES AS a 
                   LEFT OUTER JOIN ASSAYS AS b ON a.assay_id=b.assay_id 
                   LEFT OUTER JOIN ASSAY_PARAMETERS AS c ON b.assay_id=c.assay_id
                   LEFT OUTER JOIN COMPOUND_PROPERTIES AS d ON a.molregno=d.molregno
                   LEFT OUTER JOIN COMPOUND_STRUCTURES AS e ON a.molregno=e.molregno
                   LEFT OUTER JOIN MOLECULE_DICTIONARY AS f ON a.molregno=f.molregno
                   WHERE b.chembl_id=='CHEMBL3301365'")
chembl_hppb <- dbFetch(res)
dbClearResult(res)
chembl_hppb <- chembl_hppb[,colSums(is.na(chembl_hppb))<nrow(chembl_hppb)]

# 4.solu
res <- dbSendQuery(exampledb, "SELECT * FROM ACTIVITIES AS a 
                   LEFT OUTER JOIN ASSAYS AS b ON a.assay_id=b.assay_id 
                   LEFT OUTER JOIN ASSAY_PARAMETERS AS c ON b.assay_id=c.assay_id
                   LEFT OUTER JOIN COMPOUND_PROPERTIES AS d ON a.molregno=d.molregno
                   LEFT OUTER JOIN COMPOUND_STRUCTURES AS e ON a.molregno=e.molregno
                   LEFT OUTER JOIN MOLECULE_DICTIONARY AS f ON a.molregno=f.molregno
                   WHERE b.chembl_id=='CHEMBL3301364'")
chembl_solu <- dbFetch(res)
dbClearResult(res)
chembl_solu <- chembl_solu[,colSums(is.na(chembl_solu))<nrow(chembl_solu)]

rm(list=setdiff(ls(),ls()[grepl('chembl',ls())]))

#save.image(file = 'chembl_data.Rdata')

####################################################################################################################
# In case one doesn't want to download chembl data, they can directly load 'chembl_data.Rdata' file provided in 
# repository

#load("C:/Swapnil/AZ_case_study/chembl_db/chembl_data.Rdata")


####################################################################################################################
# Build model for logD

#Find duplicates
temp_tab <- chembl_logd[c(which(duplicated(tolower(chembl_logd$canonical_smiles))),which(duplicated(tolower(chembl_logd$canonical_smiles),fromLast = T))),]
chembl_logd_sorted <- chembl_logd[!tolower(chembl_logd$canonical_smiles) %in% tolower(temp_tab$canonical_smiles),]
rm(temp_tab)

#salts and mixtures
chembl_logd_sorted <- chembl_logd_sorted[-grep('\\.',chembl_logd_sorted$canonical_smiles),]
chembl_logd_sorted <- chembl_logd_sorted[order(chembl_logd_sorted$molregno),]

# Data analysis
# Normality test and plots
hist(chembl_logd_sorted$published_value)
ggdensity(chembl_logd_sorted$published_value,main = "plot",xlab = "logD")
shapiro.test(chembl_logd_sorted$published_value)

write.table(chembl_logd_sorted$canonical_smiles,'smiles_4187_logd.csv',row.names = F,col.names = F,quote = F)
write.table(chembl_logd_sorted$published_value,'logd_set_4187_endpoints.csv',row.names = F,col.names = F)
#saveRDS(chembl_logd_sorted$published_value,'logd_set_4187_endpoints.rds')
#saveRDS(chembl_logd_sorted,'logd_set_4187_chem_table.rds')

# File 'smiles_4187_logd.csv' will be used as input for PaDEL software in order to calculate descriptors.

####################################################################################################################
#  Calculate descriptors using PaDEL software outside of Rstudio and then import those descriptors in Rstudio.
#  File containing descriptors 'descr_4187_chem_logd_set.csv' has been provided in github repository  
####################################################################################################################

# Load descriptors that were calculated from PaDEL software

set.seed(3456)

x_data <- fread('C:/Swapnil/AZ_case_study/desc_calc/descr_4187_chem_logd_set.csv',stringsAsFactors = F)
row.names(x_data) <- c(x_data$Name)
x_data <- x_data[,-'Name']
y_data <- fread('C:/Swapnil/AZ_case_study/logd_set_4187_endpoints.csv')

# Remove descriptors with missing values
x_data <- x_data[,colSums(is.na(x_data))<nrow(x_data),with=F]    # No missing values

# Remove descriptors with near zero variance
nzv <- nearZeroVar(x_data)
x_data <- x_data[, -nzv,with=F]             # 346 descriptors are omitted
rm(nzv)

# Data scaling and centering
x_data <- scale(x_data,center = T, scale = T)

# Remove highly correlated descriptors
correlationMatrix <- cor(x_data)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9)
x_data <- x_data[,-highlyCorrelated]                                  #185 descriptors removed
x_data <- as.data.frame(x_data)
rm(correlationMatrix,highlyCorrelated)


# Data splitting in training and test set
trainIndex <- createDataPartition(y_data$V1, p = .8,list = FALSE,times = 1)

x_train <- as.data.frame(x_data[trainIndex,])
x_test <- as.data.frame(x_data[-trainIndex,])
y_train <- y_data$V1[trainIndex]
y_test <- y_data$V1[-trainIndex]

####################################################################################################################
# Multiple linear regression method
# 1. lm

rm(list = setdiff(ls(),c('x_train','x_test','y_train','y_test')))
x_train_backup <- x_train
y_train_backup <- y_train

# Outlier detection
temp_mod <- lm(y_train ~ ., data=x_train)
cooksd <- cooks.distance(temp_mod)
plot(cooksd, pch="*", cex=2, main="Influential chemicals by Cooks distance")
abline(h = 20*mean(cooksd, na.rm=T), col="red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>50*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  
outlier_index <- as.vector(which(cooksd>20*mean(cooksd, na.rm=T)))
x_train <- x_train[-outlier_index,]
y_train <- y_train[-outlier_index]


# Perform 5 fold cross validation
ncrssv <- 5
row.names(x_train) <- c(1:nrow(x_train))
cross_val_ref_table <- as.data.table(cbind(group_id = sample(unique(row.names(x_train))), group_no = 1:ncrssv))
all_models <- list()
all_models_res <- list()

ptm <- proc.time()
for(cv in 1:ncrssv){
  ytrn <- y_train[which(!row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id])]
  ytst <- y_train[which(row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id])]
  xtrn <- x_train[which(!row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id]),]
  xtst <- x_train[which(row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id]),]
  lmFit <- lm(ytrn~., data=xtrn)
  Pred_int_test <- predict(lmFit, xtst)
  actuals_preds_test <- data.frame(cbind(obs=ytst, pred=Pred_int_test))
  temp_res <- defaultSummary(actuals_preds_test)
  training_name     <- paste('model_cv',cv,sep = '_')
  all_models[[training_name]] <- lmFit
  all_models_res[[training_name]] <- temp_res
  print(training_name)
}
proc.time() - ptm

# collect model output in a table
res_all <- data.frame(cbind(
  Sr=1:length(all_models_res),
  Model_name=unlist(names(all_models_res)),
  Group_no=gsub('.*cv_','\\1',names(all_models_res)),
  RMSE=unlist(lapply(1:length(all_models_res),function(x) all_models_res[[x]][1])),
  Rsquared= unlist(lapply(1:length(all_models_res),function(x) all_models_res[[x]][2]))
))

res_all$Rsquared <- as.numeric(as.character(res_all$Rsquared))
res_all$RMSE <- as.numeric(as.character(res_all$RMSE))
res_all$Model_name <- as.character(res_all$Model_name)

# calculate average of Rsquared and RMSE over all 5 folds of CV
res_all[1,'avg_r2'] <- c(mean(res_all$Rsquared))
res_all[1,'avg_RMSE'] <- c(mean(res_all$RMSE))


# best model selection based on highest Rsquared
signi_model_name <- res_all$Model_name[which(res_all$Rsquared==max(res_all$Rsquared[!is.na(res_all$Rsquared)]))]
best_model <- all_models[[eval(signi_model_name)]]

# External test set prediction
test_pred <- predict(best_model, newdata = x_test)
actuals_preds_test <- data.frame(cbind(obs=y_test, pred=test_pred))
defaultSummary(actuals_preds_test)

# Plot observed vs predicted
#xyplot(actuals_preds_test$obs ~ actuals_preds_test$pred,type = c("p", "g"),xlab = "Predicted", ylab = "Observed")
p <- ggplot(aes(x=actual, y=pred),data=data.frame(actual=actuals_preds_test$obs, pred=actuals_preds_test$pred))
p + geom_point() + geom_abline(color="red") + ggtitle(paste("Linear regression: r^2=", round(defaultSummary(actuals_preds_test)[2],2), sep=""))

####################################################################################################################
# Instance based method
# 2. kNN regression

rm(list = setdiff(ls(),c('x_test','y_test','x_train_backup','y_train_backup')))
x_train <- x_train_backup
y_train <- y_train_backup

# Perform 5 fold cross validation for k value ranging from 1 to 10
ncrssv <- 5
row.names(x_train) <- c(1:nrow(x_train))
cross_val_ref_table <- as.data.table(cbind(group_id = sample(unique(row.names(x_train))), group_no = 1:ncrssv))
all_models <- list()
all_models_res <- list()

ptm <- proc.time()
for(k_in in 1:10){
  for(cv in 1:ncrssv){
    ytrn <- y_train[which(!row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id])]
    ytst <- y_train[which(row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id])]
    xtrn <- x_train[which(!row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id]),]
    xtst <- x_train[which(row.names(x_train) %in% cross_val_ref_table[group_no==cv,group_id]),]
    knn_fit <- knnreg(xtrn, ytrn, k = k_in)
    test_pred <- predict(knn_fit, newdata = xtst)
    actuals_preds_test <- data.frame(cbind(obs=ytst, pred=test_pred))
    temp_res <- defaultSummary(actuals_preds_test)
    training_name     <- paste('model_cv',cv,'k',k_in,sep = '_')
    all_models[[training_name]] <- knn_fit
    all_models_res[[training_name]] <- temp_res
    print(training_name)
  }
}
proc.time() - ptm

# collect model output in a table
res_all <- data.frame(cbind(
  Sr=1:length(all_models_res),
  Model_name=unlist(names(all_models_res)),
  Group_no=gsub('model_cv_([0-9]{1})_.*','\\1',names(all_models_res)),
  RMSE=unlist(lapply(1:length(all_models_res),function(x) all_models_res[[x]][1])),
  Rsquared= unlist(lapply(1:length(all_models_res),function(x) all_models_res[[x]][2])),
  k_value=gsub('.*_k_','\\1',names(all_models_res))
  ))

res_all$Rsquared <- as.numeric(as.character(res_all$Rsquared))
res_all$RMSE <- as.numeric(as.character(res_all$RMSE))
res_all$k_value <- as.numeric(as.character(res_all$k_value))
res_all$Model_name <- as.character(res_all$Model_name)

all_groups <- seq(1,50,5)             # calculate average of Rsquared and RMSE over all 5 folds of CV
for(j in all_groups){
  res_all[j:j+4,'avg_r2'] <- c(mean(res_all$Rsquared[j:j+4]))
  res_all[j:j+4,'avg_RMSE'] <- c(mean(res_all$RMSE[j:j+4]))
}

# Plot RMSE vs k and Rsquared vs K, in order to figure out optimal value of k
xyplot(Rsquared ~ k_value , group=Group_no, data = res_all,xlim=c(seq(0,11,1)),auto.key=list(title="CV", x = 0.99, y =0.8, cex=1.0))
xyplot(RMSE ~ k_value , group=Group_no, data = res_all,xlim=c(seq(0,11,1)),auto.key=list(title="CV", x = 0.99, y =0.8, cex=1.0))

# best model selection based on highest Rsquared
signi_model_name <- res_all$Model_name[which(res_all$avg_r2==max(res_all$avg_r2[!is.na(res_all$avg_r2)]))]
best_model <- all_models[[eval(signi_model_name)]]

# External test set prediction
test_pred <- predict(best_model, newdata = x_test)
actuals_preds_test <- data.frame(cbind(obs=y_test, pred=test_pred))
defaultSummary(actuals_preds_test)

####################################################################################################################
# Random subspace method based
# 3. RF

rm(list = setdiff(ls(),c('x_test','y_test','x_train_backup','y_train_backup')))
x_train <- x_train_backup
y_train <- y_train_backup

# Find optimal ntree
temp_model <- randomForest(y = y_train, x = x_train, ytest = y_test, xtest = x_test, do.trace = TRUE,importance = TRUE,
                                      keep.inbag = TRUE,keep.forest = TRUE,replace = TRUE)
plot(temp_model)
# 100 trees are enough

# Lets find optimal mtry
t <- tuneRF(x_train,y_train, stepFactor = 1,plot = T, ntreeTry = 100, trace = T, improve = 0.1)
# mtry = 67

# Final model
rf_fit <-  randomForest(y = y_train, x = x_train, do.trace = TRUE,importance = TRUE,ntree = 100,mtry = 67, 
                        keep.inbag = TRUE,keep.forest = TRUE,replace = TRUE)
# training set predictions
r2_train <- rSquared(y_train, y_train - rf_fit$predicted)
mse_train <- mean((y_train - rf_fit$predicted)^2)
rmse_train <- sqrt(mse_train)


r2_test <- rSquared(y_test, y_test - predict(rf_fit, x_test))
mse_test <- mean((y_test - predict(rf_fit, x_test))^2)
rmse_test <- sqrt(mse_test)

p <- ggplot(aes(x=actual, y=pred),data=data.frame(actual=y_test, pred=predict(rf_fit, x_test)))
p + geom_point() + geom_abline(color="red") + ggtitle(paste("RF results: r^2=", round(r2_test,2), sep=""))

####################################################################################################################
# Kernel based method
# 4. Support vector regression

rm(list = setdiff(ls(),c('x_test','y_test','x_train_backup','y_train_backup')))
x_train <- x_train_backup
y_train <- y_train_backup

#Regression with SVM
input_table <- cbind(x_train,y_train)
modelsvm <- svm(y_train~.,data=x_train)

postResample(predict(modelsvm, x_train),y_train)

#Predict using SVM regression
pred_test <- predict(modelsvm, x_test)

postResample(pred_test,y_test)

# End