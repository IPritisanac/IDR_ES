

if (!require("glmnet")) install.packages("glmnet")
library('glmnet');
if (!require("pROC")) install.packages("pROC");
library('pROC');
args<-commandArgs(trailingOnly=TRUE)
filename<-"HUMAN_ES.txt"
go_file<-"GO0005730_nucleolus.txt"
dir_name<-"HIGH_AUC_PPV_FILTERED_FUNCTIONS"
if (length(args)>=1) go_file<-args[1]
if (length(args)>=2) filename<-args[2]
if (length(args)>=3) dir_name<-args[3]
save_roc_plots<-toupper(Sys.getenv("FAIDR_SAVE_ROC_PLOTS","FALSE"))=="TRUE"
out_dir<-file.path("FAIDR_output", sub("\\.txt$", "", go_file))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

data<-read.delim(filename);
data[is.na(data)]<-0; #set missing data to 0. R's regression models don't always work well with missing data.
data<-data[!is.na(data$idr_name) & nzchar(as.character(data$idr_name)),];

rownames(data)<-data$idr_name

X<-as.matrix(data[,3:length(data[1,])]);
#X<-scale(X) ##standarize feature matrix

iva_GO_files<-go_file

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

## group IDRs by protein (NAME from label file, not parsed idr_name)
go_df<-read.delim(paste(dir_name,iva_GO_files[1],sep="/"))
protein_name<-go_df$NAME[match(data$idr_name, go_df$idr_name)];
gene2row<-list();
for (i in 1:nrow(data)) {
  gene<-as.character(protein_name[i]);
  gene2row[[gene]]<-c(gene2row[[gene]],i);
}
initw<-NULL; #will store initial weights
multi<-NULL;
for (g in 1:length(gene2row)) {
  initw[gene2row[[g]][1]]<-1;
  if (length(gene2row[[g]])>1) { #if there is more than one idr for this gene ...
    multi<-c(multi,g);
    initw[gene2row[[g]]] <- 1/length(gene2row[[g]]); #downweight them all equally
  }
}

lam <- 0.2;
#N_splits <- 6;
N_reps<-100
Balance <- 3;
Specificity <- 0.99;
feat_names<-colnames(X)
tstats_dir<-file.path(out_dir, "Tstats")
dir.create(tstats_dir, recursive=TRUE, showWarnings=FALSE)
Tstats<-matrix(nr=(length(X[1,])+1), nc=length(targets));
PostIDR<-matrix(nr=length(X[,1]), nc=length(targets));
Pred<-matrix(nr=length(gene2row), nc=length(targets));
geney<-matrix(nr=length(gene2row), nc=length(targets));

IDR_level_pr<-data.frame(matrix(0,nr=length(X[,1]), nc=length(targets)));
IDR_level_pred_counts<-data.frame(matrix(0,nr=length(X[,1]), nc=length(targets)));
classi_info<-data.frame(matrix(0,nr=length(targets)-1,nc=7));
classi_info[1]<-iva_GO_files
for (rep in 1:N_reps) {
for (t in 2:length(targets)) {
  y<-targets[,t];
  for (g in 1:length(gene2row)) {  geney[g,t]<-y[gene2row[[g]][1]];  }

  ##fit the model to all the training data, but get predictions for everything.
  
  train_pos <- which(y==1) #use all the positives
  #train_neg<-sample(which(y==0),Balance*length(train_pos)) 
  gene_neg <- sample(which(geney[,t]==0) , sum(geney[,t])*Balance) #sample negatives
  train_neg<- unlist(gene2row[gene_neg]) ## get their IDRs
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
  
  t_vec<-rep(NA, length(feat_names)+1)
  names(t_vec)<-c("(Intercept)", feat_names)
  if (sum(cols)>0) {
    subX<-X[,cols];
    coef_glm<-summary(glm(z~subX, weights=w*w1))$coefficients
    idx<-c(1, 1+which(cols))
    t_vec[idx]<-coef_glm[,3]
    Tstats[idx,t]<-coef_glm[,3]
  }
  write.table(data.frame(feature=names(t_vec), t_stat=t_vec),
    file.path(tstats_dir, sprintf("Tstats_rep_%03d.txt", rep)), sep="\t", row.names=FALSE)
  
  roc_obj<-roc(as.numeric(y[train]),as.numeric(pr[train]))
  if (save_roc_plots) plot(roc_obj)
  FPR_index <- min(which(roc_obj[[3]]>Specificity)) ##get the index when the FPR < 0.01 (ROC is in Specificity=1-FPR)
  thresh<-roc_obj[[4]][FPR_index]
  #get the number of positive IDR-level predictions
  classi_info[t-1,3] <- classi_info[t-1,3] + length(intersect(which(y==1),which(pr>thresh)))/N_reps;
  classi_info[t-1,4] <- classi_info[t-1,4] + length(intersect(which(y==0),which(pr>thresh)))/N_reps;
  classi_info[t-1,2] <- classi_info[t-1,2] + thresh/N_reps;
  classi_info[t-1,5] <- classi_info[t-1,5] + as.numeric(auc(roc_obj))/N_reps;
  classi_info[t-1,7] <- classi_info[t-1,7] + sum(y)/N_reps;
  classi_info[t-1,6] <- classi_info[t-1,6] + as.numeric(auc(roc(
    as.numeric(geney[c(which(geney[,t]==1),gene_neg),t]),
    as.numeric(Pred[c(which(geney[,t]==1),gene_neg),t])
  )))/N_reps;
  
  IDR_level_pr[t]<- IDR_level_pr[t] + pr/N_reps;
  IDR_level_pred_counts[t]<- IDR_level_pred_counts[t] + (pr>thresh)/N_reps;
}
}
IDR_level_pr[1]<-targets[1]
colnames(IDR_level_pr)<-colnames(targets)
IDR_level_pred_counts[1]<-targets[1]
colnames(IDR_level_pred_counts)<-colnames(targets)
write.table(IDR_level_pred_counts, file.path(out_dir, "IDR_above_cut_bal_control_no_cv.txt"), sep="\t",row.names = FALSE)
colnames(classi_info)<-c("category","thresh","above_thresh_in_annotated","above_thresh_not_in_annotated","IDR_AUC_CV","Prot_AUC_CV","N_annotated_IDRs")
write.table(classi_info, file.path(out_dir, "IDR_classification_bal_control_no_cv.txt"), sep="\t",row.names = FALSE)

