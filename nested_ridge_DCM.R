library(mlbench)
library(readxl)
library(MASS)
library(reshape)
library(ggplot2)
library(dplyr)
library(gtools)
library(stringr)
library(glmnet)

# clinical data import
df <- read_excel("/Users/.../Data/clinical_data.xlsx")

# connectivity data import and labeling
DCMB <- read_excel("/Users/.../Data/DCM_B_model1.xlsx")
DCMA <- read_excel("/Users/.../Data/DCM_A_model1.xlsx")

names(DCMB) <- gsub("_aff", "_iles", names(DCMB))
names(DCMA) <- gsub("_aff", "_iles", names(DCMA))

names(DCMA) <- gsub("_unaff", "_cles", names(DCMA))
names(DCMB) <- gsub("_unaff", "_cles", names(DCMB))

names(DCMB) <- gsub("DCMB_", "DCMB ", names(DCMB))
names(DCMA) <- gsub("DCMA_", "DCMA ", names(DCMA))

# select connectivity of subsamples (healthy and patient group)
DCM_A_and_B = cbind(DCMB,DCMA)
df_pat = df[14:26,]
DCMB_pat = DCMB[14:26,]
DCMB_hc = DCMB[1:13,]
DCMA_pat = DCMA[14:26,]
DCMA_hc = DCMA[1:13,]
DCM_AnB_pat = DCM_A_and_B[14:26,]

# define function testing for significantly positive or negative connections, corrected for multiple comparisons
mult_comp_connectome <- function(data = DCM_mat){
  pvals = NULL
  test_df = data
  for (col in names(test_df)){
    test <- t.test(test_df[[col]])
    pvals = rbind(pvals, data.frame(test$estimate, test$p.value))
  }
  row.names(pvals) <- names(test_df)
  BH = p.adjust(pvals$test.p.value, "BH")
  BH = data.frame(BH)
  BH_mat = cbind(pvals, BH)
  return(BH_mat)
}

# define dependend variable
y = df_pat$aips_Vel_RD
#y = df_pat$dpmc_Vel_RD
#y = df_pat$m_Vel_RD
y = as.matrix(y)

n <- dim(DCMA_pat)[1]

# Prepare nested Cross Validation using the LASSO
predictions2 <- c()
sig_cons <- c()
connection_counter <- c()

# for data handling purposes, find the maximum of potentially predictive connections (all significant connections in patient subgroups or healthy control group)

for (i in 1:n) {
  training_DCMA_pat = as.data.frame(DCMA_pat[-i,])
  training_DCMB_pat = as.data.frame(DCMB_pat[-i,])
  BH_DCMA_pat_i <- mult_comp_connectome(data = training_DCMA_pat)
  BH_DCMB_pat_i <- mult_comp_connectome(data = training_DCMB_pat)
  BH_DCM_AnB_pat_i <- rbind(BH_DCMA_pat_i,BH_DCMB_pat_i)
  sig_cons <- cbind(sig_cons,BH_DCM_AnB_pat_i[,"BH"])
  
}
sig_cons = as.data.frame(sig_cons)
rownames(sig_cons) <- rownames(BH_DCM_AnB_pat_i)
any_sig_cons_pat <- as.matrix(rowSums(sig_cons < 0.05))
colnames(any_sig_cons_pat)[1] <- "count"
any_sig_cons_pat = as.data.frame(any_sig_cons_pat[which(rowSums(any_sig_cons_pat) > 0),])

# additionaly, check significant connections in healthy controls
BH_DCMB_hc <- mult_comp_connectome(data = DCMB_hc)
BH_DCMA_hc <- mult_comp_connectome(data = DCMA_hc)
BH_DCM_AnB_hc <- rbind(BH_DCMA_hc,BH_DCMB_hc)
BH_hc_sig = row.names(BH_DCM_AnB_hc[BH_DCM_AnB_hc[, "BH"] < 0.05,])
BH_pat_sig_any = rownames(any_sig_cons_pat)
col.num_hc_any <- which(colnames(DCM_AnB_pat) %in% BH_hc_sig)
col.num_pat_any <- which(colnames(DCM_AnB_pat) %in% BH_pat_sig_any)

# summarize connections significant in HC OR Pat
col.num_both_any = c(col.num_hc_any,col.num_pat_any)
col.num_both_df_any = data.frame(col.num_both_any)
col.num_both_df_any = col.num_both_df_any[!duplicated(col.num_both_df_any), ]
DCM_AnB_pat_in_hc_or_patsig_any <- as.matrix(DCM_AnB_pat[,sort(c(col.num_both_df_any))])

pred_connections <- as.data.frame(matrix(0, nrow = ncol(DCM_AnB_pat_in_hc_or_patsig_any), ncol = 1))
rownames(pred_connections) <- colnames(DCM_AnB_pat_in_hc_or_patsig_any)

# begin nested cross-validation loops
for (i in 1:n) {
  training_DCMA_pat = as.data.frame(DCMA_pat[-i,])
  training_DCMB_pat = as.data.frame(DCMB_pat[-i,])
  training_DCM_AnB_pat = as.data.frame(DCM_AnB_pat[-i,])

  training_y = as.matrix(y[-i])
  
  test_DCMA_pat = as.data.frame(DCMA_pat[i,])
  test_DCMB_pat = as.data.frame(DCMB_pat[i,])
  test_DCM_AnB_pat = as.data.frame(DCM_AnB_pat[i,])

  BH_DCMB_pat_i = mult_comp_connectome(data = training_DCMB_pat)
  BH_DCMA_pat_i = mult_comp_connectome(data = training_DCMA_pat)
  BH_DCM_AnB_pat_i = rbind(BH_DCMA_pat_i,BH_DCMB_pat_i)
  
  # find connections significant in HC OR Pat for training- and test-subset
  BH_pat_sig_i = row.names(BH_DCM_AnB_pat_i[BH_DCM_AnB_pat_i[, "BH"] < 0.05,])
  BH_hc_sig = row.names(BH_DCM_AnB_hc[BH_DCM_AnB_hc[, "BH"] < 0.05,])
  col.num_hc = which(colnames(training_DCM_AnB_pat) %in% BH_hc_sig)
  col.num_pat_i = which(colnames(training_DCM_AnB_pat) %in% BH_pat_sig_i)
  col.num_both = c(col.num_hc,col.num_pat_i)
  col.num_both_df = data.frame(col.num_both)
  col.num_both_df = col.num_both_df[!duplicated(col.num_both_df), ]
  train_DCM_AnB_pat_in_hc_or_patsig_i = as.matrix(training_DCM_AnB_pat[,sort(c(col.num_both_df))])
  test_DCM_AnB_pat_in_hc_or_patsig_i = as.matrix(test_DCM_AnB_pat[,sort(c(col.num_both_df))])

  #cross-validation for training subset with alpha set to 1
  lassoResults_i<-cv.glmnet(train_DCM_AnB_pat_in_hc_or_patsig_i,training_y,alpha=1,nfolds=12)
  fit_i <- lassoResults_i$glmnet.fit
  
  #find lambda with least mean cross validated error and compute prediction of TMS effects
  bestlambda_i<-lassoResults_i$lambda.min
  y_predicted_i <- predict(fit_i, s = bestlambda_i, newx = test_DCM_AnB_pat_in_hc_or_patsig_i)

  #add to predictions from other subsamples of leave-one-out cross validation
  predictions2 <- c(predictions2, y_predicted_i)
  results_i<-predict(lassoResults_i,s=bestlambda_i,type="coefficients")
  resmat_i = as.data.frame(as.matrix(results_i)[-1,])
  colnames(resmat_i)[1] <- paste(i)

  choicePred_i<-rownames(results_i)[which(results_i !=0)]
  connection_counter <- append(connection_counter,choicePred_i)
}

colnames(pred_connections)[1] <- "count"
for (j in 1:nrow(pred_connections)) {
  pred_connections[j,"count"] = length(grep(rownames(pred_connections)[j], connection_counter))
}

cv_picked_connections = as.matrix(pred_connections)
cv_picked_connections = as.data.frame(cv_picked_connections)
colnames(cv_picked_connections)[1] <- "selection_frequency"
cv_picked_connections$selection_frequency <- cv_picked_connections$selection_frequency / (n)

cv_picked_connections$type <- substr(rownames(cv_picked_connections), 0, 4)

pearson_correlation <- cor.test(predictions2,y, method = "pearson", use = "complete.obs")

cv_picked_connections_reduced<-cv_picked_connections[!(cv_picked_connections$selection_frequency==0),]

#plot selected connections
tiff('cv_selections.tiff', units="in", width=6.5, height=3, res=300)
ggplot(cv_picked_connections_reduced, aes(x = reorder(row.names(cv_picked_connections_reduced), selection_frequency), y = selection_frequency, fill = factor(type))) +
  scale_fill_manual(values = c("#a1d99b", "#31a354")) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  ggtitle("Connections consistently selected by the LASSO")
dev.off()

#plot model predictions
new_df <- data.frame(y, predictions2)
colnames(new_df)[1] <- "Observed_TMS_Effect"
names(new_df)[2] <- "Predicted_Effect"
tiff('cv_predicted_observed.tiff', units="in", width=6.5, height=4, res=300)
ggplot(new_df, aes(x=Predicted_Effect, y=Observed_TMS_Effect)) +
  geom_point(color = I("#666666"), size=3) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color = I("#999999"),linetype="dotted",size=1) +
  theme_minimal() +
  ggtitle(paste("Connectivity-dependent TMS effect: Pearson r: ", round(pearson_correlation$estimate,2), " p = ", round(pearson_correlation$p.value,4)))
dev.off()
