if (!require("glmnet")) install.packages("glmnet")
library('glmnet');
if (!require("pROC")) install.packages("pROC");
library('pROC');

SEED <- 1               # choose rand n seed for this run
set.seed(SEED)

filename<-"HUMAN_ES.txt"
dir_name<-"HIGH_AUC_PPV_FILTERED_FUNCTIONS/"
#dir_name<- "SUBSET_FUNCTIONS"

data<-read.delim(filename);
data[is.na(data)]<-0; #set missing data to 0. R's regression models don't always work well with missing data.

rownames(data)<-data$idr_name

gene2row<-list();
gene2idr<-list();
for (i in 1:length(data$idr_name)) {
  gene<-strsplit(toString(data$idr_name[[i]]),'_')[[1]][1];
  idr<-as.numeric(strsplit(toString(data$idr_name[[i]]),'_')[[1]][3])
  gene2row[[gene]]<-c(gene2row[[gene]],i);
  gene2idr[[gene]]<-c(gene2idr[[gene]],idr);
  
}

X<-as.matrix(data[,3:length(data[1,])]);
#X<-scale(X) ##standarize feature matrix
initw<-NULL; #will store initial weights
multi<-NULL;
for (g in 1:length(gene2row)) {
  initw[gene2row[[g]][1]]<-1;
  if (length(gene2row[[g]])>1) { #if there is more than one idr for this gene ...
    multi<-c(multi,g);
    initw[gene2row[[g]]] <- 1/length(gene2row[[g]]); #downweight them all equally
  }
}


iva_GO_files<-list.files(dir_name)

##single category
#targets<- read.delim(paste(dir_name,iva_GO_files[34],sep="/"))
#rownames(targets)<-targets$idr_name
#targets<-targets[rownames(data),]

##all categories
targets<-data.frame(matrix(NA,nr=dim(X)[1],nc=length(iva_GO_files)+1))
targets[1]<-read.delim(paste(dir_name,iva_GO_files[1],sep="/"))$idr_name ##get the names from the first file
for (i in 1:length(iva_GO_files)) { ##get all the annotations
  targets[i+1]<-read.delim(paste(dir_name,iva_GO_files[i],sep="/"))[3] 
}
colnames(targets)<-c("idr_name",iva_GO_files) ##set the names of the columns
rownames(targets)<-targets$idr_name ##need this to get the order to match!
targets<-targets[rownames(data),]

fit_FAIDR_to_split<-function(train_ind,initw,lam) {
  
  w<-initw;
  z<-NULL;
  #infer IDR posterior probabilities in logistic regression model.
  niter<-20;
  pr<-NULL;
  mod<-glmnet(X[train_ind,],y[train_ind],family="gaussian",weights=w[train_ind]); #intial model
  for (i in 1:niter) {
    for (j in 1:5) { #M-step: fit the weighted logistic regression by iterative least squares
      pr<-1/(1+exp(-predict(mod,X,s=lam)));
      w1<-pr*(1-pr);
      w1[w1==0]<-exp(-20); ### avoid 0s
      z<-predict(mod,X,s=lam) + (y - pr)/w1;
      mod<-glmnet(X[train_ind,],z[train_ind],family="gaussian",weights=(w1*w)[train_ind])
    }
    #E-step: recalculate weights given latest model
    for (g in multi) {
      rows<-gene2row[[g]];
      w[rows] <- y[rows]*pr[rows] + (1-y[rows])*(1-pr[rows])
      w[rows] <- w[rows]/sum(w[rows]);
    }
  }
  return(mod)
}


lam <- 0.2;
N_splits <- 6;
Balance <- 3;
# Specificity <- 0.99; # No longer needed for thresholding, though kept if referenced elsewhere (commented out to avoid confusion)
Tstats<-matrix(nr=(length(X[1,])+1), nc=length(targets));
PostIDR<-matrix(nr=length(X[,1]), nc=length(targets));
Pred<-matrix(nr=length(gene2row), nc=length(targets));
geney<-matrix(nr=length(gene2row), nc=length(targets));

IDR_level_pr<-data.frame(matrix(nr=length(X[,1]), nc=length(targets)));
classi_info<-data.frame(matrix(nr=length(targets)-1,nc=7));
classi_info[1]<-iva_GO_files
#for (t in 2:5) {
for (t in 2:length(targets)) {
  y<-targets[,t];
  
  #split_assign <- NULL;
  gene_split_assign <-cut(sample(1:length(gene2row)), breaks = N_splits, labels = FALSE);
  for (g in 1:length(gene2row)) { 
    geney[g,t]<-y[gene2row[[g]][1]];
    #split_assign[gene2row[[g]]] <- gene_split_assign[g]; 
  }
  
  unbiased_pr<-NULL
  held_out_y<-NULL
  for (this_split in 1:N_splits) { 
    test_pos<-intersect(which(y==1),unlist(gene2row[which(gene_split_assign==this_split)])) ##use all the positives
    #test_neg<-intersect(which(y==0),which(split_assign==this_split)) ###use all the negatives
    #test_neg<-sample(intersect(which(y==0),which(split_assign==this_split)),Balance*length(test_pos)) ##sample IDRs
    neg_genes<-sample( intersect(which(geney[,t]==0),which(gene_split_assign==this_split)) , sum(geney[which(gene_split_assign==this_split),t])*Balance ) ##sample negative genes
    test_neg<-unlist(gene2row[ neg_genes ]) ##get their IDRs
    
    train_pos<-intersect(which(y==1),unlist(gene2row[which(gene_split_assign!=this_split)]))##use all the positives
    
    neg_genes<-sample( intersect(which(geney[,t]==0),which(gene_split_assign!=this_split)) , sum(geney[which(gene_split_assign!=this_split),t])*Balance ) ##sample negative genes
    train_neg<- unlist(gene2row[ neg_genes ])##get their IDRs
    #train_neg<-intersect(which(y==0),which(split_assign!=this_split)) ###use all the negatives
    train<-c(train_pos,train_neg)
    
    split_mod<-fit_FAIDR_to_split(train,initw,lam)
    split_pr <- 1/(1+exp(-predict(split_mod,X,s=lam)));
    unbiased_pr[test_pos]<-split_pr[test_pos]
    #held_out_y[test_pos]<-1
    unbiased_pr[test_neg]<-split_pr[test_neg]
    #held_out_y[test_neg]<-0
  }  
  ##get the unbiased FPR  
  roc_obj<-roc(as.numeric(y[which(!is.na(unbiased_pr))]),as.numeric(unbiased_pr[which(!is.na(unbiased_pr))]))
  plot(roc_obj)
  #print("unbiased AUC is")
  #auc(roc_obj)
  
  # --- CHANGED SECTION START ---
  # Calculate G-means (Geometric Mean of Sensitivity and Specificity)
  # Note: roc_obj$sensitivities and roc_obj$specificities are safer than indices [[3]] or [[2]]
  g_means <- sqrt(roc_obj$sensitivities * roc_obj$specificities)
  
  # Find the index where G-mean is maximized
  g_means_max_index <- which.max(g_means)
  
  # Set the threshold corresponding to the max G-mean
  thresh <- roc_obj$thresholds[g_means_max_index]
  # --- CHANGED SECTION END ---
  
  classi_info[t-1,2] <- thresh
  
  classi_info[t-1,5] <- as.numeric(auc(roc_obj))
  classi_info[t-1,7] <- sum(y)
  
  #get an unbiased prediction at the protein level
  unbiased_w <- initw;
  for (g in multi) {
    rows<-gene2row[[g]];
    unbiased_w[rows] <- y[rows]*unbiased_pr[rows] + (1-y[rows])*(1-unbiased_pr[rows])
    unbiased_w[rows] <- unbiased_w[rows]/sum(unbiased_w[rows]);
  }
  unbiased_gene_pr<-NULL
  for (g in 1:length(gene2row)) {
    
    #predictions for each gene are the sum of the idrs in the gene
    unbiased_gene_pr[g]<-sum(unbiased_pr[gene2row[[g]]]*unbiased_w[gene2row[[g]]]);
  }
  
  classi_info[t-1,6] <- as.numeric(auc(roc(
    as.numeric(geney[which(!is.na(unbiased_gene_pr)),t]),
    as.numeric(unbiased_gene_pr[which(!is.na(unbiased_gene_pr))])
  )))
  
  ##fit the model to all the training data, but get predictions for everything.
  
  train_pos <- which(y==1) #use all the positives
  #train_neg<-sample(which(y==0),Balance*length(train_pos)) 
  train_neg<- unlist(gene2row[ sample(which(geney[,t]==0) , sum(geney[,t])*Balance) ]) ##sample negatives and get their IDRs
  
  #train_neg <- which(y==0) ##use all negatives
  train<-c(train_pos,train_neg)
  
  
  w<-initw;
  z<-NULL;
  #infer IDR posterior probabilities in logistic regression model.
  niter<-20;
  pr<-NULL;
  mod<-glmnet(X[train,],y[train],family="gaussian",weights=w[train]); #intial model
  for (i in 1:niter) {
    for (j in 1:5) { #M-step: fit the weighted logistic regression by iterative least squares
      pr<-1/(1+exp(-predict(mod,X,s=lam)));
      w1<-pr*(1-pr);
      w1[w1==0]<-exp(-20); ### avoid 0s
      z<-predict(mod,X,s=lam) + (y - pr)/w1;
      mod<-glmnet(X[train,],z[train],family="gaussian",weights=(w1*w)[train])
    }
    #E-step: recalculate weights given latest model
    for (g in multi) {
      rows<-gene2row[[g]];
      w[rows] <- y[rows]*pr[rows] + (1-y[rows])*(1-pr[rows])
      w[rows] <- w[rows]/sum(w[rows]);
    }
  }
  #final weights using a formula that doesn't depend on y
  for (g in multi) {
    rows<-gene2row[[g]];
    w[rows] <- (pr[rows]+0.05)/sum(pr[rows]+0.05) #post(x) given that this is a positive example (y=1), smoothed so small pr are ignored
  }
  PostIDR[,t]<-w;
  #store the final predictions
  for (g in 1:length(gene2row)) {
    geney[g,t]<-y[gene2row[[g]][1]];
    #predictions for each gene are the sum of the idrs in the gene
    Pred[g,t]<-sum(pr[gene2row[[g]]]*w[gene2row[[g]]]);
  }
  
  #find the non-zero coefficients
  cols<-abs(as.numeric(coef(mod,s=lam))[2:(length(X[1,])+1)])>0;
  
  if (sum(cols)>0) {
    subX<-X[,cols];
    Tstats[c(1,1+which(cols)),t]<-summary(glm(z~subX, weights=w*w1))$coefficients[,3];
  }
  #get the number of positive IDR-level predictions
  classi_info[t-1,3] <- length(intersect(which(y==1),which(pr>thresh)))
  classi_info[t-1,4] <- length(intersect(which(y==0),which(pr>thresh)))
  IDR_level_pr[t]<-pr;
}

IDR_level_pr[1]<-targets[1]
colnames(IDR_level_pr)<-colnames(targets)

out_probs <- sprintf("IDR_probability_eq1_bal_control_seed%03d.txt", SEED) # tag the outputs with the seed number ID

write.table(IDR_level_pr, "IDR_probability_eq1_bal_control.txt", sep="\t",row.names = FALSE)
colnames(classi_info)<-c("category","thresh","above_thresh_in_annotated","above_thresh_not_in_annotated","IDR_AUC_CV","Prot_AUC_CV","N_annotated_IDRs")
write.table(classi_info, "IDR_classification_bal_control.txt", sep="\t",row.names = FALSE)