"""
__author__: Iva Pritisanac
__date__: Oct 11, 2021
"""

#
# DESCRIPTION #
#
#    When run as main -- for a set of IDRs (with UNIPROT IDs) of interest
#    Needs input files:
#        i)   query cluster file (.cdt) of interest (second column must contains Uniprot ID)
#        ii)  background proteome file (.cdt) (second column must contains Uniprot ID)
#        iii) GO annotation (e.g. goa_human.gaf)
#        iv)  GO information details (go.obo)
#
#    On the cluster of interest does:
#    1. GO enrichment analysis
#    2. Feature enrichment analysis
#        For 1 and 2 outputs:
#       > O/E fold enrichment and two-sided p-value (by default with bonferroni correction)
#    To run in bulk on a directory (folder) that contains one or many clusters (.cdt files) of interest:
#    python3 PATH/TO/PROTEOME/CDT/FILE PATH/TO/QUERY/DIRECTORY/WITH/CLUSTER/FILES
#
import os,sys
from scipy import stats
from statsmodels.stats import multitest
import numpy as np

class StatsTests():

    def fisher_test(self,pos,neg):
        oddsratio,pvalue=stats.fisher_exact([[pos[0],pos[1]],[neg[0],neg[1]]])

        return oddsratio,pvalue

    def compute_fisher(self,pos_n1,neg_n1,pos_n2,neg_n2):
        pos_tuple=(pos_n1,pos_n2)
        neg_tuple=(neg_n1,neg_n2)
        ratio,p_value=self.fisher_test(pos_tuple,neg_tuple)

        return p_value

# class that runs GO enrichment analysis on a set of genes (uniprot_ids) of interest
# the enrichment is computed relative to the IDR-ome (the set of all uniprot IDs that are included in the IDR-ome analysis)
class GoAnalysis(StatsTests):
    def __init__(self,go_obo_file,go_annot_file,background_file):

        self.get_go_labels(go_obo_file)
        # get all the background (proteome) data when initiating the class
        self.background_list=self.read_query(background_file,1)
        #self.get_go(go_annot_file,self.background_list) #old GO annotations from GO website
        self.parse_go_doc(go_annot_file) # to get GOs from PANTHER file

    def get_go_labels(self,go_obo_path):
        fin=open(go_obo_path,'r')
        collect=False
        self.go_labels={}
        for line in fin:
            stripped=line.strip()
            splitted=stripped.split('\t')
            if splitted[0].startswith("id: GO:"):
                collect=True
                go_number=splitted[0].split(":")[-1]
                self.go_labels.setdefault(go_number,{})
            if collect:
                if splitted[0].startswith("name:"):
                    name=splitted[0].split(":")[-1]
                    self.go_labels[go_number]["name"]=name
                if splitted[0].startswith("namespace:"):
                    gocat=splitted[0].split(":")[-1]
                    self.go_labels[go_number]["namespace"]=gocat
            if splitted[0]=="[Term]":
                collect=False

        if len(self.go_labels.keys())==0:
            print("WARNING: NO GO TERMS READ IN FROM\nCHECK FILE FORMAT %s"%(go_obo_path))

    # indx gives the column of the text file with uniprot IDs, by default the second (1st indx)
    def read_query(self,infile,indx=1):
        uniprot_list=[]
        fin=open(infile,"r")
        for line in fin:
            splitted=line.split()
            if len(splitted)==0:
                continue
            #uniprotid=splitted[indx]
            try:
                uniprotid=splitted[indx].split("_")[0]
                #uniprotid=splitted[indx].split(":")[0]
                if uniprotid=="IDRID": # skip header lines
                    continue
                if uniprotid=="0.895822":
                    #print("??",uniprotid)
                    continue
                uniprot_list.append(uniprotid)
            except IndexError:
                continue

        uniprot_list=list(set(uniprot_list))
        #print(uniprot_list)
        #sys.exit(0)

        return uniprot_list

    # qindx = 1 for uniprotids, qindx = 2 for geneids
    def get_go(self,infile,uniprot_list,qindx=1):
        go_terms_uniprot={}
        uniprot_go_terms={}
        fin=open(infile,"r")
        cnt=0
        for line in fin:
            splitted=line.split()
            if len(splitted)==0 or splitted[0][0]=="!":
                continue
            if splitted[qindx] in set(uniprot_list):
                #if (splitted[3]=="colocalizes_with" or splitted[3]=="NOT" or splitted[3]=="contributes_to" or splitted[3]=="NOT|colocalizes_with" or splitted[3]=="NOT|contributes_to") and ("GO:" in splitted[4]):
                if "GO:" in splitted[3]:
                    go_term=splitted[3]
                elif ("GO:" not in splitted[3]) and ("GO:" in splitted[4]):
                    go_term=splitted[3]+"_"+splitted[4]
                else:
                    print("UNRECOGNIZED ENTRY IN goa_human.gaf\n%s"%(line))
                go_terms_uniprot.setdefault(go_term,[]).append(splitted[qindx])
                uniprot_go_terms.setdefault(splitted[qindx],[]).append(go_term)
        ## REMOVE DUPLICATES or multi entries
        self.go_terms_uniprot=self.remove_duplicates_dict(go_terms_uniprot)
        self.uniprot_go_terms=self.remove_duplicates_dict(uniprot_go_terms)


    def get_go_query(self,query_list):
        go_terms_query={}
        for uniprot in query_list:
            if uniprot in self.uniprot_go_terms.keys():
                go_terms=self.uniprot_go_terms[uniprot]
                for go_term in go_terms:
                    go_terms_query.setdefault(go_term,[]).append(uniprot)
            else:
                print("SKIPPING ENTRY %s IN QUERY UNIPROT -- NO GO TERM DATA AVAILABLE"%(uniprot))

        go_terms_query=self.remove_duplicates_dict(go_terms_query)

        #print(go_terms_query)

        return go_terms_query

    def remove_duplicates_dict(self,data_dict):
        for key,value in data_dict.items():
            data_dict[key]=list(set(value))

        return data_dict

    def parse_go_doc(self,panther_go_file):
        self.uniprot_go_terms={}
        fin=open(panther_go_file,"r")
        for line in fin:
            stripped=line.strip()
            if len(stripped)==0:
                continue
            splitted=stripped.split("\t")
            uniprot=splitted[1]
            for i in range(len(splitted)):
                if "GO" in splitted[i]:
                    go_terms=splitted[i].split(";")
                    self.uniprot_go_terms.setdefault(uniprot,go_terms)

            #print(len(splitted))
            #print(splitted)
            #if stripped.startswith("id:"):
            #    key=stripped.split(": ")[1]
            #    go_doc.setdefault(key,{})
            #else:
            #    try:
            #        nkey=stripped.split(": ")[0]
            #        nval=stripped.split(": ")[1]
            #        go_doc[key].setdefault(nkey,[]).append(nval)
            #    except IndexError:
            #        print("ERROR IN LINE %s"%(line))
        self.go_terms_uniprot={}
        for key,value in self.uniprot_go_terms.items():
            for val in value:
                go_split=val.split(";")
                for entry in go_split:
                    self.go_terms_uniprot.setdefault(entry,[]).append(key)

        #print(self.uniprot_go_terms)
        #print("\n\n")
        #print(self.go_terms_uniprot)
        #sys.exit()


    def compute_enrichment(self,query_go_dict,total_query,background_go_dict,total_background,outputf,correct='bonferroni'):
        fout=open(outputf,'w')
        gos=[]
        fisher_pvalues=[]
        enrichments=[]
        nbs=[]
        nqs=[]
        for go,queries in query_go_dict.items():

            nq=len(set(queries))
            try:
                nb=len(set(background_go_dict[go]))
            except KeyError:
                print("GO TERM %s NOT FOUND IN THE BACKGROUND\nSKIPPING..."%(go))
                continue
            expected=float(nb)/float(total_background)
            observed=float(nq)/float(total_query)
            fold_enrichment=observed/expected
            fpval=self.compute_fisher(nq,total_query-nq,nb,total_background-nb)
            gos.append(go)
            fisher_pvalues.append(fpval)
            enrichments.append(fold_enrichment)
            nbs.append(nb)
            nqs.append(nq)

        ## DO A P-VAL CORRECTION FOR MULTIPLE STATISTICAL TESTS
        try:
            #print(fisher_pvalues)
            acceptance,p_values_corr,alpha_sidak,alpha_bon=multitest.multipletests(fisher_pvalues, alpha=0.05, method=correct, is_sorted=False, returnsorted=False)

            fout.write("GO_TERM\tFOLD_ENRICHMENT\tp-val(%s corrected)"%(correct))
            for i in range(len(acceptance)):
                if acceptance[i] and (float(enrichments[i])>=2.5): # export only the GO terms with enrichment over 2.5x relative to expectation
                    go=gos[i].split(":")[-1]
                    try:
                        go_name=self.go_labels[go]["name"]
                        go_cat=self.go_labels[go]["namespace"]
                        fout.write("\n%s\t%s\t%s\t%.3f\t%.3e\t%s\t%s"%(gos[i],go_name,go_cat,enrichments[i],p_values_corr[i],str(nbs[i])+"/"+str(total_background),str(nqs[i])+"/"+str(total_query)))

                    except KeyError:
                        fout.write("\n%s\t%s\t%s\t%.3f\t%.3e\t%s\t%s"%(gos[i],"NA","NA",enrichments[i],p_values_corr[i],str(nbs[i])+"/"+str(total_background),str(nqs[i])+"/"+str(total_query)))

                    #fout.write("\n%s\t%.3f\t%.3e"%(gos[i],enrichments[i],p_values_corr[i]))
            fout.close()
        except (KeyError, ZeroDivisionError) as error:
            print(error)
            print("ERROR FOUND IN p-val ADJUSTMENT FOR %s"%(outputf))


class FeatAnalysis(StatsTests):

    def __init__(self):
        pass

    def read_cluster_file(self,data_file,indx_rows=1,indx_cols=2): # indx_rows and cols - indx at which to start taking values
        darray=np.genfromtxt(data_file, delimiter='\t') #, skip_header, skip_footer, converters, missing_values, filling_values, usecols, names, excludelist, deletechars, replace_space, autostrip, case_sensitive, defaultfmt, unpack, usemask, loose, invalid_raise, max_rows, encoding)
        #features=np.genfromtxt(data_file, delimiter='\t',usecols=[0],dtype=str)
        cnt=0
        fin=open(data_file,'r')
        for line in fin:
            #print(line)
            cnt+=1
            stripped=line.strip()
            myarray=stripped.split("\t") # KEEP \t delimiter HERE OR ELSE FEATURES GET MESSED UP
            if cnt==1:
                break
        features=np.array(myarray[indx_cols:],dtype='str')
        #print(features)
        data=darray[indx_rows:,indx_cols:]
        data=np.nan_to_num(data,nan=0)
        #sys.exit()
        #print(data)


        return data,features

    def get_feat_counts(self,data,features,threshold=5):
        indices_pos=np.argwhere(data>=float(threshold)) # get indices of rows (IDRs) and columns (features) that have Zscores above a threshold
        indices_neg=np.argwhere(data<=float(threshold)*-1)
        total_rows=data.shape[0]

        #for indx in indices_pos:
        #    print(indx[0],indx[1],features[indx[1]],data[indx[0],indx[1]])
            #print()

        #for i in range(len(features)):
        #    print(i,features[i])
        #sys.exit()

        feat_pos={}
        feat_neg={}
        for i in range(len(features)): # make an empty list for every feature
            feat_pos.setdefault(features[i],[])
            feat_neg.setdefault(features[i],[])

        for indx in indices_pos[:,1]: # these are indices of columns with pos Zscores above t
            feat=features[indx]
            feat_pos[feat].append(indx) # append everytime the column (feat) above threshold

        for indx in indices_neg[:,1]: # these are indices of columns with neg Z scores above t
            feat=features[indx]
            feat_neg[feat].append(indx) # append everytime the column (feat) above threshold

        feat_count={}
        for feat,occ in feat_pos.items():
            feat_count.setdefault(feat,(len(occ),len(feat_neg[feat])))

        return feat_count,total_rows

    # given the cluster file
    # get a numpy array -- find which features (columns) >= threshold or <= threshold
    def find_enriched_features(self,df,ef,feats,outputf,threshold=10,correct='bonferroni'): # ef - expectation dataframe, df - query data frame
        qcounts,qtotal=self.get_feat_counts(df,feats,threshold) # query counts and total
        ecounts,etotal=self.get_feat_counts(ef,feats,threshold) # expectation counts and total
        fout=open(outputf,'w')
        fisher_pvalues_pos=[]
        fisher_pvalues_neg=[]
        enrichments_pos=[]
        enrichments_neg=[]
        #enriched_pos=[]
        #enriched_neg=[]
        feats=[]
        for key,value in qcounts.items():

            qpos=float(value[0])
            qneg=float(value[1])
            try:
                epos=ecounts[key][0]
                eneg=ecounts[key][1]
            except KeyError:
                print("ERROR:Missing feat %s for proteome background"%(key))
                continue
            feats.append(key)

            expected_pos=float(epos)/etotal
            expected_neg=float(eneg)/etotal
            observed_pos=float(qpos)/qtotal
            observed_neg=float(qneg)/qtotal
            try:
                fold_enrich_pos=observed_pos/expected_pos
                #print("pos",key,fold_enrich_pos)

                enrichments_pos.append(fold_enrich_pos)
                fpval_pos=self.compute_fisher(qpos,qtotal-qpos,epos,etotal-epos)
                fisher_pvalues_pos.append(fpval_pos)
            except ZeroDivisionError: # skip any features that do not have any Zscore>thresh in reference
                #continue
                #print(key)
                enrichments_pos.append(0)
                fisher_pvalues_pos.append(100)

            try:
                fold_enrich_neg=observed_neg/expected_neg
                #print("neg",key,fold_enrich_neg)
                enrichments_neg.append(fold_enrich_neg)
                fpval_neg=self.compute_fisher(qneg,qtotal-qneg,eneg,etotal-eneg)
                fisher_pvalues_neg.append(fpval_neg)
            except ZeroDivisionError: # skip any features
                #continue
                #print(key)
                enrichments_neg.append(0)
                fisher_pvalues_neg.append(100)

        ## DO A P-VAL CORRECTION FOR MULTIPLE STATISTICAL TESTS
        pos_acceptance,pos_p_values_corr,pos_alpha_sidak,pos_alpha_bon=multitest.multipletests(fisher_pvalues_pos, alpha=0.05, method=correct, is_sorted=False, returnsorted=False)

        fout.write("POS_ZSCORES\tFOLD_ENRICHMENT\tp-val(%s corrected)"%(correct))
        for i in range(len(pos_acceptance)):
            if pos_acceptance[i] and (enrichments_pos[i] >= threshold):
                fout.write("\n%s\t%.3e\t%.3e"%(feats[i],enrichments_pos[i],pos_p_values_corr[i]))
        #fout.close()
        ## DO A P-VAL CORRECTION FOR MULTIPLE STATISTICAL TESTS
        neg_acceptance,neg_p_values_corr,neg_alpha_sidak,neg_alpha_bon=multitest.multipletests(fisher_pvalues_neg, alpha=0.05, method=correct, is_sorted=False, returnsorted=False)

        fout.write("\n\nNEG_ZSCORES\tFOLD_ENRICHMENT\tp-val(%s corrected)"%(correct))
        for i in range(len(neg_acceptance)):
            if neg_acceptance[i] and (enrichments_neg[i] >= threshold):
                fout.write("\n%s\t%.3f\t%.3e"%(feats[i],enrichments_neg[i], neg_p_values_corr[i]))
        fout.close()

if __name__=="__main__":
    backfile=sys.argv[1] #PATH TO FILE WITH FULL IDR-ome CLUSTERING IN THE CDT FILE FORMAT
    qdir=sys.argv[2] #PATH TO DIR/FOLDER WITH SELECTED CLUSTER .CDT FILES
    goannotfile="PANTHER_GO_ANNOTATIONS.txt" #"goa_human.gaf"
    GOA=GoAnalysis("go.obo",goannotfile,backfile)
    cnt_f=0
    for inputf in os.listdir(qdir):
        if inputf.endswith('.cdt') and (inputf[:-4]+"_GO_ANALYSIS.txt" not in os.listdir(qdir)):
            fpath=qdir+os.sep+inputf
            print("RUNNING GO ANALYSIS ON %s"%(fpath))
            qlist=GOA.read_query(fpath)
            if len(set(qlist))<5:
                print("SKIPPING DUE TO SMALL SAMPLE SIZE")
                continue

            goquery=GOA.get_go_query(qlist)
            cnt_f+=1
            GOA.compute_enrichment(goquery,len(qlist),GOA.go_terms_uniprot,len(GOA.background_list),fpath[:-4]+"_GO_ANALYSIS.txt")
            FE=FeatAnalysis()
            bf,bfeat=FE.read_cluster_file(backfile,2,4) # read in background (IDRome) - the whole cdt file
            cf,cfeat=FE.read_cluster_file(fpath,2,4) # read in specific cluster with a subset of IDRs of interest
            FE.find_enriched_features(cf,bf,bfeat,fpath[:-4]+"_FEAT_ENRICH_ANALYSIS.out",5)
    print("TOTAL ELIGIBLE CLUSTERS %s"%(cnt_f))
