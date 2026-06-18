
if (!require("glmnet")) install.packages("glmnet")
library('glmnet');
if (!require("pROC")) install.packages("pROC");
library('pROC');

SEED <- 1
set.seed(SEED)

filename<-"HUMAN_ES.txt"
target_file<-"HIGH_AUC_PPV_FILTERED_FUNCTIONS/GO0005730_nucleolus.tsv"

data<-read.delim(filename);
data[is.na(data)]<-0; #set missing data to 0. R's regression models don't always work well with missing data.

rownames(data)<-data$idr_name

gene2row<-list();
for (i in 1:length(data$idr_name)) {
  gene<-strsplit(toString(data$idr_name[[i]]),'_')[[1]][1];
  gene2row[[gene]]<-c(gene2row[[gene]],i);
}

X<-as.matrix(data[,3:length(data[1,])]);

initw<-NULL; #will store initial weights
multi<-NULL;
for (g in 1:length(gene2row)) {
  initw[gene2row[[g]][1]]<-1;
  if (length(gene2row[[g]])>1) { #if there is more than one idr for this gene ...
    multi<-c(multi,g);
    initw[gene2row[[g]]] <- 1/length(gene2row[[g]]); #downweight them all equally
  }
}

targets<-read.delim(target_file)
rownames(targets)<-targets$idr_name
targets<-targets[rownames(data),]
y<-targets[,2]
category_name<-sub("\\.(txt|tsv)$", "", basename(target_file))

fit_FAIDR<-function(train_ind,initw,lam) {
  w<-initw;
  z<-NULL;
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
  list(mod=mod, pr=pr, w=w, w1=w1, z=z)
}


lam <- 0.2;
N_splits <- 6;
N_reps<-as.integer(Sys.getenv("N_REPS", unset = "100"))
Balance <- 3;
Specificity <- 0.99;
Tstats<-matrix(nr=(length(X[1,])+1), nc=1);
PostIDR<-numeric(length(X[,1]));
Pred<-numeric(length(gene2row));
geney<-numeric(length(gene2row));

IDR_level_pr<-data.frame(idr_name=targets$idr_name, pr=0)
IDR_level_pred_counts<-data.frame(idr_name=targets$idr_name, pred=0)
classi_info<-data.frame(matrix(0,nr=1,nc=7));
classi_info[1,1]<-category_name
for (g in 1:length(gene2row)) { geney[g]<-y[gene2row[[g]][1]]; }

auc_cv_reps<-numeric(N_reps)
roc_cv_list<-vector("list", N_reps)
for (rep in 1:N_reps) {
  gene_split_assign <-cut(sample(1:length(gene2row)), breaks = N_splits, labels = FALSE);
  unbiased_pr<-NULL
  for (this_split in 1:N_splits) { 
    test_pos<-intersect(which(y==1),unlist(gene2row[which(gene_split_assign==this_split)])) ##use all the positives
    neg_genes<-sample( intersect(which(geney==0),which(gene_split_assign==this_split)) , sum(geney[which(gene_split_assign==this_split)])*Balance ) ##sample negative genes
    test_neg<-unlist(gene2row[ neg_genes ]) ##get their IDRs
    
    train_pos<-intersect(which(y==1),unlist(gene2row[which(gene_split_assign!=this_split)]))##use all the positives
  
    neg_genes<-sample( intersect(which(geney==0),which(gene_split_assign!=this_split)) , sum(geney[which(gene_split_assign!=this_split)])*Balance ) ##sample negative genes
    train_neg<- unlist(gene2row[ neg_genes ])##get their IDRs
    train<-c(train_pos,train_neg)
    
    split_fit<-fit_FAIDR(train,initw,lam)
    split_pr <- 1/(1+exp(-predict(split_fit$mod,X,s=lam)));
    unbiased_pr[test_pos]<-split_pr[test_pos]
    unbiased_pr[test_neg]<-split_pr[test_neg]
  }  
  roc_cv_list[[rep]]<-roc(as.numeric(y[which(!is.na(unbiased_pr))]),as.numeric(unbiased_pr[which(!is.na(unbiased_pr))]))
  auc_cv_reps[rep]<-as.numeric(auc(roc_cv_list[[rep]]))
}
classi_info[1,5] <- mean(auc_cv_reps)
if (Sys.getenv("SKIP_PLOT", unset = "") != "1") {
plot(roc_cv_list[[1]], col=adjustcolor("steelblue", 0.2), lwd=0.8, legacy.axes=TRUE,
     main=sprintf("CV ROC (%d reps), mean AUC = %.3f", N_reps, mean(auc_cv_reps)))
for (rep in 2:N_reps) {
  lines(roc_cv_list[[rep]], col=adjustcolor("steelblue", 0.2), lwd=0.8)
}
}

for (rep in 1:N_reps) {
  ##fit the model to all the training data, but get predictions for everything.
  
  train_pos <- which(y==1) #use all the positives 
  gene_neg <- sample(which(geney==0) , sum(geney)*Balance) #sample negatives
  train_neg<- unlist(gene2row[gene_neg]) ## get their IDRs
  train<-c(train_pos,train_neg)
  fit<-fit_FAIDR(train,initw,lam)
  mod<-fit$mod; pr<-fit$pr; w<-fit$w; w1<-fit$w1; z<-fit$z
  #final weights using a formula that doesn't depend on y
  for (g in multi) {
    rows<-gene2row[[g]];
    w[rows] <- (pr[rows]+0.05)/sum(pr[rows]+0.05) #post(x) given that this is a positive example (y=1), smoothed so small pr are ignored
  }
  PostIDR<-w;
  for (g in 1:length(gene2row)) {
    Pred[g]<-sum(pr[gene2row[[g]]]*w[gene2row[[g]]]);
  }
  
  #find the non-zero coefficients
  cols<-abs(as.numeric(coef(mod,s=lam))[2:(length(X[1,])+1)])>0;
  
  if (sum(cols)>0) {
    subX<-X[,cols];
    Tstats[c(1,1+which(cols)),1]<-summary(glm(z~subX, weights=w*w1))$coefficients[,3];
  }
  
  roc_train<-roc(as.numeric(y[train]),as.numeric(pr[train]))
  FPR_index <- min(which(roc_train[[3]]>Specificity)) ##get the index when the FPR < 0.01 (ROC is in Specificity=1-FPR)
  thresh<-roc_train[[4]][FPR_index]
  #get the number of positive IDR-level predictions
  classi_info[1,3] <- classi_info[1,3] + length(intersect(which(y==1),which(pr>thresh)))/N_reps;
  classi_info[1,4] <- classi_info[1,4] + length(intersect(which(y==0),which(pr>thresh)))/N_reps;
  classi_info[1,2] <- classi_info[1,2] + thresh/N_reps;
  classi_info[1,6] <- classi_info[1,6] + as.numeric(auc(roc(
    as.numeric(geney[c(which(geney==1),gene_neg)]),
    as.numeric(Pred[c(which(geney==1),gene_neg)])
  )))/N_reps;
  
  IDR_level_pr$pr<- IDR_level_pr$pr + pr/N_reps;
  IDR_level_pred_counts$pred<- IDR_level_pred_counts$pred + (pr>thresh)/N_reps;
}
classi_info[1,7] <- sum(y)
colnames(IDR_level_pr)<-c("idr_name", category_name)
colnames(IDR_level_pred_counts)<-c("idr_name", category_name)
write.table(IDR_level_pred_counts, "IDR_above_cut_bal_control_no_cv.txt", sep="\t",row.names = FALSE)
colnames(classi_info)<-c("category","thresh","above_thresh_in_annotated","above_thresh_not_in_annotated","IDR_AUC_CV","Prot_AUC_CV","N_annotated_IDRs")
write.table(classi_info, "IDR_classification_bal_control_no_cv.txt", sep="\t",row.names = FALSE)

