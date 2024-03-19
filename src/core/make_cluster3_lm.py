import os,sys
import networkx as nx

"""
Write out a linkage matrix from CLUSTER3.0 .gtr tree file
store Node: parent:children structure
Test dendogram method with this linkage matrix
"""


# needed imports
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram,leaves_list,linkage,fcluster
import numpy as np


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

def plot_from_file(linkage_matrix,cthresh,leaves_labels,outfile):
    # Generate and Plot the corresponding dendrogram
    ddict=dendrogram(linkage_matrix,orientation='left',color_threshold=float(cthresh),above_threshold_color='#bcbddc',get_leaves=True,show_leaf_counts=True)
    plt.savefig(outfile[:-4]+"_DENDROGRAM.pdf",format='pdf')
    #plt.show()
    fout=open(outfile,'w')
    leaves_order=list(ddict['ivl']) # goes from bottom to top
    for indx in reversed(leaves_order): # now goes from top to bottom
        ix=int(indx)
        fout.write("%s\n"%(leaves_labels[ix]))
    fout.close()

def read_idr_labels(idrs_file):
    clustered_idrs_order=[]
    fin=open(idrs_file,"r")
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        clustered_idrs_order.append(splitted[0])

    return clustered_idrs_order

def label_clusters(cdt_file,gene_cluster_dict,distthresh):
    fin=open(cdt_file,"r")
    fout=open(cdt_file[:-4]+"_CLUSTERS_DIST_"+str(distthresh)+".cdt.txt","w")
    cluster_ids=[value for value in gene_cluster_dict.values()]
    cluster_ids=list(set(cluster_ids))
    data_cluster_table=[]
    for i in range(len(cluster_ids)):
        data_cluster_table.append([])

    print(len(data_cluster_table))
    gene_uniprotid={}
    header=[]
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split("\t")
        if splitted[0].startswith("GENE"):
            fout.write("%s\t%s"%(splitted[0],str(splitted[1])+"_CLUST_"+str(gene_cluster_dict[splitted[0]])))
            cluster_id=int(gene_cluster_dict[splitted[0]])
            gene_uniprotid.setdefault(splitted[0],splitted[1])
            #print(cluster_id)
            data_cluster_table[cluster_id-1].append(line)
            for i in range(len(splitted[2:])):
                fout.write("\t%s"%(splitted[i+2]))
            fout.write("\n")
        else:
            header.append(line)
            fout.write(line)
    ## WRITE OUT GENE:UNIPROT IDR CODE AND EXport
    for i in range(len(data_cluster_table)):

        fout=open("CLUSTER_"+str(i)+"_"+str(distthresh)+".out.cdt","w")
        for head_line in header:
            fout.write(head_line)
        for entry in data_cluster_table[i]:
            fout.write(entry)
    return gene_uniprotid

#USE GENE UNIPROTID CODE TO EXPORT FOR GO ANALYSIS
def write_out_clusters(cluster_gene_dict,distthresh,gene_uniprotid):
    for cluster,genes in cluster_gene_dict.items():
        name_cluster=str(cluster)
        fout=open("CLUSTER_"+name_cluster+"_DIST_THRESH_"+str(distthresh)+".txt","w")
        for gene in genes:
            fout.write("%s\t%s\n"%(gene_uniprotid[gene],gene))



gtrfile=sys.argv[1] # CLUSTER3.0 output .gtr file - tree file
max_d = float(sys.argv[2]) # distance for thresholding
cdtfile=sys.argv[3] # CLUSTER3.0 output .cdt file

fin=open(gtrfile,"r")
#nodes_dict={} # stores parents' children node structure
T=nx.DiGraph() # directed networkx graph
all_genes=[]
all_nodes=[]
for line in fin:
    stripped=line.strip()
    splitted=stripped.split()
    all_nodes.append(splitted[0])
    if not "NODE" in splitted[1]:
        all_genes.append(splitted[1])
    if not "NODE" in splitted[2]:
        all_genes.append(splitted[2])

all_genes=sorted(all_genes,key=lambda x: int(x[4:-1]))
all_nodes=sorted(all_nodes,key=lambda x: int(x[4:-1]))

genes_nodes=all_genes+all_nodes
indx_node_dict={i:genes_nodes[i] for i in range(0,len(genes_nodes))}
#print(indx_node_dict)
fin=open(gtrfile,"r")
for line in fin:
    stripped=line.strip()
    splitted=stripped.split()
    parent=genes_nodes.index(splitted[0])
    child1=genes_nodes.index(splitted[1])
    child2=genes_nodes.index(splitted[2])
    T.add_node(parent)
    T.add_edge(parent, child1)
    T.add_edge(parent, child2)
    #nodes_dict.setdefault(key,value)
fin.close()

print(T.number_of_nodes())
print(nx.is_tree(T))

leaf_nodes = [node for node in T.nodes() if T.in_degree(node)!=0 and T.out_degree(node)==0]
print(len(leaf_nodes))
#print(leaf_nodes)

fin=open(gtrfile,"r")
Z=[]
for line in fin:
    stripped=line.strip()
    splitted=stripped.split()
    parent=genes_nodes.index(splitted[0])
    child1=genes_nodes.index(splitted[1])
    child2=genes_nodes.index(splitted[2])
    dnodes=nx.descendants(T,parent) # get all nodes descending from this node

    dleaves=[node for node in dnodes if T.in_degree(node)!=0 and T.out_degree(node)==0]
    #print(len(dleaves))
    dist=1-float(splitted[3])

    t=[child1,child2,dist,len(dleaves)]
    Z.append(t)
Zn=np.array(Z)
#print(Zn)
dendict=dendrogram(Zn,no_plot=True)
leaves=leaves_list(Zn)

cluster_annotations=fcluster(Zn,max_d,criterion='distance')

gene_cluster={}
cluster_gene={}
for i in range(len(leaves)):
    #print(leaves[i],genes_nodes[leaves[i]],cluster_annotations[leaves[i]])
    gene_cluster.setdefault(genes_nodes[leaves[i]],cluster_annotations[leaves[i]])
    cluster_gene.setdefault(cluster_annotations[leaves[i]],[]).append(genes_nodes[leaves[i]])

geneuniprot_dict=label_clusters(cdtfile,gene_cluster,max_d)

#visualize dendrogram
#dendrogram(Zn,
#    truncate_mode='lastp',  # show only the last p merged clusters#
#    p=100,  # show only the last p merged clusters
#    leaf_rotation=90.,
#    leaf_font_size=12.,
#    show_contracted=True,  # to get a distribution impression in truncated branches
#)
#plt.show()

#fancy_dendrogram(
#    Zn,
#    truncate_mode='lastp',
#    p=1000,
#    leaf_rotation=90.,
#    leaf_font_size=10.,
#    show_contracted=True,
#    annotate_above=10,  # useful in small plots so annotations don't overlap
#    max_d=max_d,
#)
#plt.show()
