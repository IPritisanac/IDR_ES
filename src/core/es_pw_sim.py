import io
import os,sys
import math
import random
import re
import time
import numpy as np
import scipy.stats
from src.utils.data import ParseInput
from src.core.sequence_features import SequenceFeatures

class Alignment():

	def __init__(self):
		self.names = [] #list of names of sequences in the alignment
		self.seq = {} #the actual sequence indexed by their names
		self.info = "" #any other information

	#read in mfa from a fasta file
	def read_mfa(self,file):
		f = open(file, 'r') #open the file
		for line in f:
			line=line.rstrip("\n\r")
			if ('>' in line):
				self.names.append(line) # add this header line to the list
			else:
				seqname=self.names[-1] # this is the most recently found header line
				if (not seqname in self.seq):
					self.seq[seqname] = line
				else:
					self.seq[seqname] += line

		f.close()

	def print_mfa(self): #print fasta format to the screen
		for s in self.names:
			if ('>' in s):
				 print(s)
			else:
				print ('>' + s)
			print(self.seq[s])

	def print_maf(self): #print maf format to the screen
						#assumes the score line has been stored as the info for the alignment
		print(self.info),;
		for s in self.names:
			print ("s "),;print(s),;print(self.seq[s])
		print ('\n')

# Class that computes Z-scores for features defined in sequence_features
class EVSign():

	def __init__(self,inputf):
		#self.dirpath=os.path.dirname(os.path.abspath(inputf)) # extract directory name out of the input filename
		if not os.path.exists("output"): # if the output results directory does not exist on this path
			os.makedirs("output")   # create the directory
		self.outdir="output"   # set the output directory
		  # create instance of the parser class to parse the input file
		self.parse_param_files(inputf)
		self.parse_variables(inputf)
		self.get_features()
		self.print_var=False # do not print output variances unless this is set to True

	## run parser and return only file paths of various parameter files
	def parse_param_files(self,input_file):
		PP=ParseInput() # create an instance of a Parser class
		try:
			PP.parse(input_file)
		except Exception as e:
			print("ERROR: Could not parse input file:\n%s"%e)
			sys.exit(1)
		try:
			PP.check_variables()
		except Exception as e:
			print("ERROR in input variables!\n")
			print(e)
			sys.exit(1)

		PP.read_in_files()
		# Parse parameters

		self.repeats=PP.repeats
		self.motifs=PP.motifs
		self.motifs_expect=PP.motifs_expect
		self.aln_dir=PP.aln_dir
		self.aa_freq=PP.aa_freq
		self.max_indel_size=PP.max_indel_size
		self.indel_prob=PP.indel_prob

	## run parser and return only variables -- different modes, integer parameter values etc.
	def parse_variables(self,input_file):
		PV=ParseInput()  # create a instance of the parser class to parse the input file
		try:
			PV.parse(input_file)
		except Exception as e:
			print("ERROR: Could not parse input file:\n%s"%e)
			sys.exit(1)
		try:
			PV.check_variables()
		except Exception as e:
			print("ERROR in input variables!\n")
			print(e)
			sys.exit(1)

		self.use_indels=PV.use_indels
		self.n_sim=PV.n_sim
		self.ref_n=PV.ref_n
		self.ref_name=PV.ref_name

		self.min_sd=PV.min_sd
		self.len_min=PV.len_min
		self.len_factor=PV.len_factor
		self.dist_ratio=PV.dist_ratio
		self.dist_tot=PV.dist_tot

	def get_features(self,):
		SF=SequenceFeatures()
		self.aafeats=SF.aafeats
		self.aafeats_names=SF.aafeats_names
		self.seqfeats=SF.seqfeats
		self.seqfeats_names=SF.seqfeats_names


	def P_JC69_aa(self,A,B,d): # probability of changing from B to A in evolutionary distnace d
		if (A is B):
			return 0.05*(1 + float(19)*math.exp(-float(20)*d/float(19)))
		else:
			return 0.05*(1 - math.exp(-float(20)*d/float(19)))

	# computes distance between two aligned sequences
	def JC69_aa_dis(self,s1,s2): #distance between two aligned sequences
		nid=0
		ndiff=0
		ratio = float(19)/float(20)
		for pos in range(0,len(s1)):
			if ('-' in s1[pos]+s2[pos]):
				pass
			elif (s1[pos] == s2[pos]):
				nid+=1
			else:
				ndiff+=1
		if (nid == 0):
			return float(10)
		elif (ndiff == 0):
			return float(0)
		else:
			frac = float(ndiff)/float(nid+ndiff)

			if (frac>=ratio):
				return float(10)
			else:
				return -ratio*math.log(float(1) - frac/ratio)

	def P_F81_aa(self,A,B,d,freq):
		ss=float(0)
		for x in freq:
			#print(x)
			ss+=freq[x]*freq[x]
		beta=float(1)/(float(1)-ss)
		if (A == B):
			return math.exp(-beta*d) + freq[A]*(1 - math.exp(-beta*d))
		else:
			return freq[A]*(1 - math.exp(-beta*d))

	def fastF81_aa(self,A,B,d,fA,beta):
		if (A == B):
			return math.exp(-beta*d) + fA*(1 - math.exp(-beta*d))
		else:
			return fA*(1 - math.exp(-beta*d))


	#computes distance between two aligned sequences, assuming s1 evolved from s2
	## uses JC69 model as an initial guess
	# does brute force 1D ML estimation
	def F81_aa_dis(self,s1,s2,freq):
		initdis=self.JC69_aa_dis(s1,s2) ## use JC69 as an initial guess
		#brute force 1D ML estimation
		N=150 #how many values to try
		ss=float(0) ###for faster f81 calculations
		for x in freq.values():
			ss += x*x
		beta=float(1)/(float(1)-ss)
		alnpairs = zip(list(s1),list(s2))

		ilogL=float(0)
		for aa1, aa2 in alnpairs:
			if ((aa1 == '-') or (aa2 == '-')):
				continue
			ilogL += math.log(self.fastF81_aa(aa1,aa2,initdis,freq[aa1],beta))

		bestd=initdis
		bestL=ilogL
		for i in range(0,N):
			d=initdis + 0.01*initdis*(float(N)/float(2) - float(i))
			logL=0
			for aa1, aa2 in alnpairs:
				if ((aa1 == '-') or (aa2 == '-')):
					continue
				logL += math.log(self.fastF81_aa(aa1,aa2,d,freq[aa1],beta))

			if(logL>bestL):
				bestL=logL
				bestd=d

		return bestd

	# defines a model same as JC, but with beta=3/4
	def fastF81_aa_dis(self,s1,s2,freq):
		ss=float(0)
		for x in freq.values(): ss += x*x
		beta=float(1)/(float(1)-ss)
		nid=0
		ndiff=0
		alnpairs = zip(list(s1),list(s2))
		for aa1, aa2 in alnpairs:
			if ((aa1 == '-') or (aa2 == '-')):
				pass
			elif (aa1 == aa2):
				nid+=1
			else:
				ndiff+=1
		if (nid == 0):
			return float(10)
		elif (ndiff == 0):
			return float(0)
		frac = float(ndiff)/float(nid+ndiff)
		if (frac*beta>float(1)):
			return float(10)
		else:
			frac = float(ndiff)/float(nid+ndiff)
			return -math.log(float(1) - frac*beta)/beta

	# defines a heuristic used to get a set of sequences with good distance distribution from REFERENCE
	# avoids use of 'redundant' sequences that have same/similar distance to reference
	def seq_choosing_heuristic(self,aln,ref,idrfreq,D_RAT,D_TOT,L_FACT): ##try to get dis_tot, but don't use 'redundant' sequences
		sp_dis = {}
		refs=aln.seq[ref].replace("-","") ### reference has to be checked before this is run...
		for s in aln.names:
			if (s == ref):
				continue
			ugseq=aln.seq[s].replace("-","")
			if ('X' in ugseq) or ('B' in ugseq) or ('Z' in ugseq):
				continue #no bad data
			if ( float(L_FACT*len(ugseq)) < float(len(refs)) ):
				continue #this sequence is too short
			if ( float(len(ugseq)) > float(L_FACT*len(refs)) ):
				continue #this sequence is too long
			dis = self.fastF81_aa_dis(aln.seq[s],aln.seq[ref],idrfreq)
			if ((dis>=float(10)) or (dis<=float(0))):
				continue #distance looks bad
			sp_dis[dis]=s  ### include this sequence

		totd=float(0)
		seqlist = []
		for d in sorted(list(sp_dis.keys())):
			if (len(seqlist)==0):
				seqlist.append(sp_dis[d])
				totd=totd+d
			elif (totd<D_TOT):
				mindis=D_TOT
				for s in seqlist:
					###distance of this sequence to the previous ons added
					dis = self.fastF81_aa_dis(aln.seq[s],aln.seq[sp_dis[d]],idrfreq)
					if (dis<mindis):
						mindis=dis

				## include this sequence if it is D_RATIO further away from the reference than it is from any other sequence
				if (D_RAT*mindis > d):
					seqlist.append(sp_dis[d])
					totd=totd+d

		return seqlist

	def sim_multi_prot_given_anc(self,seq,d,freq,n): ##simulate n sequences
		SM ={}
		aalist = list("ACDEFGHIKLMNPQRSTVWY")
		#for faster F81 calculations
		ss=float(0)
		for x in freq:
			ss+=freq[x]*freq[x]
		beta=float(1)/(float(1)-ss)
		for a in aalist: # precompute the substitution matrix
			SM[a] = [self.fastF81_aa(b,a,d,freq[b],beta) for b in aalist]
		newseq= []
		seqlist=list(seq)

		for i in range(0,n):
			try:
				sseq="".join([aalist[self.random_int(SM[b])] for b in seqlist])
				newseq.append(sseq)
			except TypeError:

				print("ERROR IN CREATING SEQ N %s"%(i))
				sys.exit()
				#continue

		return(newseq)

	# create an indel of particular size, taking into the account AA probabilities
	def indel_seq(self,freq_diso,indsize):
		aa_keys=[key for key in freq_diso.keys()]
		int_key={i:aa_keys[i] for i in range(len(aa_keys))}
		aafreq=[freq_diso[key] for key in aa_keys]
		indel=''
		for i in range(0,indsize):
			aaint=np.random.choice(np.arange(0,len(aa_keys)),p=aafreq)
			indel+=int_key[aaint]

		return indel

	# get probabilities for a range of indel sizes
	#z=1.5 # exponent for the power law distribution
	#size_range=range(1,len(seq)) or range(1,1001) # indel size probability distribution
	def indel_prob(self,size_range,z):
		from mpmath import nsum, exp, inf
		indel_probs=[] # probability of an indel given its size
		for size in size_range:
			p_indel=size**-z/nsum(lambda x: x**-z, [1, inf])
			indel_probs.append(p_indel)

		probsum=sum(indel_probs)
		norm_probs=[float(prob)/float(probsum) for prob in indel_probs] # normalize probabilities so they sum up to 1

		return norm_probs

	# read in precomputed and normalized probabilities of indels of different size
	# saves computing time
	def read_indel_prob(self,indel_file):
		fin=open(indel_file,'r')
		freq=[]
		for line in fin:
			stripped=line.strip()
			splitted=stripped.split()
			size=int(splitted[0])
			freq.append(float(splitted[1]))
			#if size==int(1/3.*seq_len): # when the size approaches sequence length -- stop
			#	break
		sumfreq=sum(freq)
		normfreq=[float(f)/float(sumfreq) for f in freq]
		fin.close()

		return size+1,normfreq

	def introduce_indels(self,seq,pos_array,size_array,type_array):
		cnt1=-1
		for i in range(len(pos_array)):
			if pos_array[i]>0:
				for j in range(0,pos_array[i]):
					cnt1+=1
					size_indx=cnt1
					isize=size_array[size_indx]
					itype=type_array[size_indx]
					if itype=='ins':
						iseq=self.indel_seq(self.aa_freq,isize) # sequence of an indel to insert
						seq=seq[:i]+iseq+seq[i:] # insert the indel at random position in the sequence
					else:
						seq=seq[:i]+seq[i+isize:] # delete an indel starting at random position in the sequence
		return seq

	# d - distance between ref seq and a homologous sequence
	def simple_indel(self,start_seq,max_indel,size_prob,d):
		indels_pos=self.indel_poisson(d,start_seq) # positions of indels determined using poisson process
		samples=int(np.sum(indels_pos))
		indels_type=['ins','del']
		isizes=np.random.choice(np.arange(1,max_indel),samples,p=size_prob)
		itypes=np.random.choice(indels_type,samples) #choose deletion or insertion with equal probability

		mod_seq=self.introduce_indels(start_seq,indels_pos,isizes,itypes) # modified sequence with indels introduced

		if len(mod_seq)==0: # if the entire sequence ablated
			return start_seq # return a sequence with no indels
		else:
			return mod_seq

	def indel_poisson(self,d,seq):
		from scipy.stats import poisson
		indel_rate=1/20.*d #indel rate
		r = poisson.rvs(indel_rate, size=len(seq))

		return r

	def count_aa(self,seq):
		return {aa:seq.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}



	##//KNOWN BUG// -- if AA probabilities given to limited precision
	##		-- does not return a number 0-19 but 'NoneType'
	def random_int(self,pr):
		#print(pr)
		cumpr=float(0)
		rn = random.random()
		i=int(0);
		for p in pr:
			cumpr += p
			i+=1
			if (cumpr>=rn):
				#print(i-1)
				return (i-1)

	def expected_repeat_residues(self,r,freq,L): ##expected repeat residue content in sequence of length L
		aalist=list("ACDEFGHIKLMNPQRSTVWY")
		matchpr=float(0) ## probability that a residue matches this repeat.
		for a in r:
			if (a in aalist):
				matchpr += freq[a]
		startpr=matchpr*matchpr ## probability that a repeat starts (i.e., two residues need to match
		repeatpr = [float(0)] #repeat can't start at the first poisition
		for pos in range(1,L-1):
			repeatpr.append( startpr + repeatpr[pos-1]*matchpr )    ##looks like some kind of geometric series...
		repeatpr.append( repeatpr[L-2]*matchpr ) #...or at the last position

		return sum(repeatpr)

	def sim_prot_given_anc(self,seq,d,freq):
		newseq=[]
		aalist = list("ACDEFGHIKLMNPQRSTVWY")
		#for faster F81 calculations
		ss=float(0)
		for x in freq:
			#print(x)
			ss+=freq[x]*freq[x]
		beta=float(1)/(float(1)-ss)

		for b in list(seq):
			thispr=float(0)
			rn=random.random()
			for a in aalist:
				thispr += fastF81_aa(a,b,d,freq[a],beta)
				if (thispr>=rn):
					newseq.append(a)
					break;
		return("".join(newseq))

	def repeats_norm(self,r,seq):
		repregex=r.split("m=")[-1].replace('"','')
		er= self.expected_repeat_residues(r, self.aa_freq, len(seq))
		repeatlength=0
		pat = re.compile(repregex)
		for m in pat.finditer(seq):
			repeatlength += len(m.group(0))
		rep_value=float(repeatlength) - er

		return rep_value

	def motifs_norm(self,m,seq):
		en= self.motifs_expect[m]
		if (not self.motifs[m][-1] == "$"):
			en*= float(len(seq)-len(self.motifs[m])) ##C-terminal motifs don't get this
			#eprint (self.motifs[m],"en=",en)
		pat = re.compile("".join(self.motifs[m]))
		mot_value=float(len(pat.findall(seq)))-en

		return mot_value

	#@param sim_sd:standard deviation over simulated obs_values
	#@param sim_values: array of feature values for simulated sequences
	#@param obs_values: array of feature values for real (observed/extant) sequences
	#@param obs_boolean: array of feature presence True/False for real (observed) sequences
	def test_limit_cases(self,feat,refseq,sim_values,obs_values,obs_boolean):
		obs_values=np.array(obs_values)
		sim_values=np.array(sim_values)
		obs_boolean=np.array(obs_boolean)
		#print(obs_values)
		#print(sim_values)
		#print(obs_boolean)
		if (np.mean(sim_values)<np.mean(obs_values)):
			# if feature observed in >= 90% extant sequences
			if (np.sum(obs_boolean>0)>=int(0.9*len(obs_values))):
				return True,15,-15
			# if feature observed in few extant sequences
			elif (np.sum(obs_boolean>0)>int(0.1*len(obs_values))) and (np.sum(obs_boolean>0)<int(0.5*len(obs_boolean))):
				return True,5,5
			elif (np.sum(obs_boolean>0)>int(0.5*len(obs_values))) and (np.sum(obs_boolean>0)<int(0.9*len(obs_boolean))):
				return True,10,10
			else:
				return False,0,0
		else:
			return False,0,0


	def eprint(self,*args):
		sys.stderr.write(str(args)+"\n")

	def mean(self,x):
		tot =float(0)
		n=0
		for val in x:
			if (not val == ""):
				tot += float(val)
				n+=1
		if (n>0):
			return tot/float(n)
		else:
			return ""

	def stdev(self,x):
		tot =float(0)
		n=0
		for val in x:
			if (not val == ""):
				tot += float(val)*float(val)
				n+=1
		var=0
		if (n>0):
			var = tot/float(n) - self.mean(x)**2
		if (var>=0):
			return math.sqrt( var )
		else:
			return ""

	def log(self,x): ##vectorized log function
		if (type(x) is list):
			return [log(val) for val in x]
		elif (x>0):
			return math.log(float(x))
		else:
			return ""

	# main loop, computes evolutionary Z-scores for alignment files in a directory specified in input_file.txt
	def compute_es_dir(self):
		start_time=time.time()
		obs = {}
		obs_b={}
		Z = {}
		featnamesorted = (self.aafeats_names + self.seqfeats_names +  list(self.repeats.keys()) + list(self.motifs.keys()))
		output_file_path = self.outdir+os.sep+"ES"+"_"+os.path.basename(self.aln_dir.rstrip("/"))+".out.txt"

		if os.path.exists(output_file_path):  # Check if file exists
			os.remove(output_file_path)       # Delete the file if it exists
			print(f"Deleted existing file: {output_file_path}")

		# open output file in append mode
		outres=open(output_file_path,"a")
		outres.write("IDR_ID\t")
		outres.write("\t".join([str(f)+"_meanZ" for f in featnamesorted])) # write header in output file
		outres.write("\n")

		for mfa in os.listdir(self.aln_dir):
			ALN = Alignment()
			try:
				ALN.read_mfa(self.aln_dir+os.sep+mfa)
			except:
				print("WARNING: Skipping file:%s"%(mfa))
				continue # skip a file if problems with reading

		    ###find the reference sequence if possible
			if (not self.ref_name==""):
				REF_NUM=0 ###set back to default in case this alignment doesn't have the reference
				for i in range(0,len(ALN.names)):
					if (self.ref_name in ALN.names[i]): # if you can find reference name or part of the name
						REF_NUM=i # set ref_num to the reference
						break

			refs = ALN.seq[ALN.names[REF_NUM]].replace("-","")

			if ((len(refs)<30) or ('X' in refs)): # include only IDRs with 30aa or more.
				continue

			print(str(len(ALN.names))+' sequences read in '+mfa+' reference is '+ALN.names[REF_NUM]+" length="+str(len(refs)))

			##normalization to help avoid zeros in the variances
			ex_rep_res = {r:self.expected_repeat_residues(self.repeats[r][0],self.aa_freq, len(refs)) for r in list(self.repeats.keys())}

			ugseqs = {}
			sequences_to_use = self.seq_choosing_heuristic(ALN,ALN.names[REF_NUM],self.aa_freq,self.dist_ratio,self.dist_tot,self.len_factor) ### this should do all the quality control

			print(str(len(sequences_to_use))+' sequences passed filtering')

			for s in sequences_to_use:
				ugseqs[s]=ALN.seq[s].replace("-","")

			if (len(list(ugseqs.keys()))<10): ##only include alignments with at least 10 orthologous sequences, length>10
				continue

			Z[mfa]={}
			SIM = {}
			for i in range(0,self.n_sim):
				SIM[i] = Alignment()

			for s in ALN.names:
				if (s == ALN.names[REF_NUM]):
					continue

				if (not s in ugseqs): #quality control
					continue
				dis=self.F81_aa_dis(ALN.seq[s],ALN.seq[ALN.names[REF_NUM]],self.aa_freq)
				if ((dis>=float(10)) or (dis<=float(0))): #probably diverged or identical
					print('df81='+str(dis)+' skipping '+s)
					continue

				i=0
				for rs in self.sim_multi_prot_given_anc(refs,dis,self.aa_freq,self.n_sim):
					#-- add indels to each of these sequences
					if self.use_indels:
						rs=self.simple_indel(rs,self.max_indel_size,self.indel_prob,dis)
					SIM[i].seq[s]=rs
					SIM[i].names.append(s)
					i+=1

			print(str(self.n_sim)+ " simulations completed for " + mfa)

			sm ={}
			sv ={}
			obs[mfa] = {}
			obs_b[mfa] = {}

			for f in self.aafeats_names + list(self.repeats.keys()) + list(self.motifs.keys()) + self.seqfeats_names:

				sm[f] = []
				sv[f] = []
				obs[mfa][f] = [] #keep the observations
				obs_b[mfa][f+"_boolean"]=[]

			### get the observed values
			for s in ALN.names:
				if (not s in ugseqs): #quality control
					continue
				aas = self.count_aa(ugseqs[s])
				aas = {a:float(aas[a]) for a in aas.keys()}
				for j in range(0,len(self.aafeats)):
					val =self.aafeats[j](aas)
					obs[mfa][self.aafeats_names[j]].append(val)
				for r in list(self.repeats.keys()):
					repeatlength=0
					er= self.expected_repeat_residues(self.repeats[r][0], self.aa_freq, len(ugseqs[s]))
					repregex=r.split("m=")[-1].replace('"','')
					pat = re.compile(repregex)
					for m in pat.finditer(ugseqs[s]):
						repeatlength += len(m.group(0))
					obs[mfa][r].append(float(repeatlength) - er)
				for m in list(self.motifs.keys()):
					en= self.motifs_expect[m]
					if (not self.motifs[m][-1] == "$"):
						en*= float(len(ugseqs[s])-len(self.motifs[m])) ##C-terminal motifs don't get this
					pat = re.compile("".join(self.motifs[m]))
					obs[mfa][m].append(float(len(pat.findall(ugseqs[s])))-en)
					if len(pat.findall(ugseqs[s]))>0:
						obs_b[mfa][m+"_boolean"].append(True)
				for j in range(0,len(self.seqfeats)):
					val =self.seqfeats[j](ugseqs[s])
					obs[mfa][self.seqfeats_names[j]].append(val)
		    ### get the mean and variance for the simulations
			for i in list(SIM.keys()):
				if (len(SIM[i].names)<5):
					continue ###do some quality control on the simulations
				for f in self.aafeats_names + list(self.repeats.keys()) + list(self.motifs.keys())+ self.seqfeats_names:
					sm[f].append(float(0)) #do not keep the values from the simulations
					sv[f].append(float(0))
				for s in SIM[i].names: #OK to have different numbers of sequences in each simulation
					aas = self.count_aa(SIM[i].seq[s])
					aas = {a:float(aas[a]) for a in aas.keys()}
					if not self.use_indels:
						length_ratio = float(len(ugseqs[s]))/float(len(refs))
						# aas take into account lenght difference **if indels not used**
						aas = {a:float(aas[a])*length_ratio for a in aas.keys()}
					for j in range(0,len(self.aafeats)):
						val =self.aafeats[j](aas)
						sm[self.aafeats_names[j]][-1]+= val
						sv[self.aafeats_names[j]][-1]+= val*val
					for j in range(0,len(self.seqfeats)):
						val=self.seqfeats[j](SIM[i].seq[s])
						sm[self.seqfeats_names[j]][-1]+= val
						sv[self.seqfeats_names[j]][-1]+= val*val
					for r in list(self.repeats.keys()):
						repeatlength=0
						repregex=r.split("m=")[-1].replace('"','')
						pat = re.compile(repregex)
						for m in pat.finditer(SIM[i].seq[s]):
							repeatlength += len(m.group(0))
						val=(float(repeatlength) - ex_rep_res[r])
						if not self.use_indels:
							val=(float(repeatlength) - ex_rep_res[r])*length_ratio
						sm[r][-1]+= val
						sv[r][-1]+= val*val
					for m in list(self.motifs.keys()):
						pat = re.compile("".join(self.motifs[m]))
						en=self.motifs_expect[m]
						if (not self.motifs[m][-1] == "$"):
							#en *= float(len(refs)-len(self.motifs[m]))
							en*=float(len(SIM[i].seq[s])-len(self.motifs[m]))
						val=(float(len(pat.findall(SIM[i].seq[s]))) - en)
						if not self.use_indels:
							val=(float(len(pat.findall(SIM[i].seq[s]))) - en)*length_ratio

						sm[m][-1]+= val
						sv[m][-1]+= val*val
				for f in self.aafeats_names + list(self.repeats.keys()) + list(self.motifs.keys()) + self.seqfeats_names:
					sm[f][-1] /= float(len(SIM[i].names))
					sv[f][-1] /= float(len(SIM[i].names))
					sv[f][-1] -= sm[f][-1]*sm[f][-1]
					if (sv[f][-1]>float(0)):
						sv[f][-1] = math.log(sv[f][-1])
					else:
						sv[f][-1] = np.nan

			for f in self.aafeats_names + list(self.repeats.keys()) + list(self.motifs.keys()) + self.seqfeats_names:
				# FOR MEAN Z-scores
				simsd=np.nanstd(sm[f])
				if ((simsd=="") or (not simsd > self.min_sd)):
					simsd=self.min_sd
					mz = (np.nanmean(obs[mfa][f]) - np.nanmean(sm[f]))/simsd
					# test limit cases
					if (f in self.motifs.keys()) and (mz>100 or mz<-100):
						limit_test,Zscore_mean,Zscore_var=self.test_limit_cases(f,refs,sm[f],obs[mfa][f],obs_b[mfa][f+"_boolean"])
						if limit_test: # giving an arbitrarily large mean and var score
							Z[mfa][f+" meanZ"]=Zscore_mean
							Z[mfa][f+" varZ"]=Zscore_var
							continue
						else:
							pass
					else:
						Z[mfa][f+" meanZ"]=mz
				else:
					mz = (np.nanmean(obs[mfa][f]) - np.nanmean(sm[f]))/simsd
					Z[mfa][f+" meanZ"]=mz

				# FOR VARIANCE Z-scores
				#obssd=np.nanstd(obs[mfa][f])
				#simsd=np.nanstd(sv[f])
				#elogv=np.nanmean(sv[f])

				#if ((simsd=="") or (not simsd > self.min_sd)):
				#	simsd=self.min_sd
				#if (elogv<float(2)*math.log(self.min_sd)):
					#adding a small number here
				#	elogv=float(2)*math.log(self.min_sd)
				#if ((obssd=="") or (not obssd > self.min_sd)):
				#	obssd=self.min_sd

				#vz = (float(2)*math.log(obssd) - elogv)/simsd
				#Z[mfa][f+" varZ"]=vz

			outres.write(mfa +"\t")
			outres.write("\t".join([str(Z[mfa][f+" meanZ"]) for f in featnamesorted]))
			if self.print_var: # print variance Z-scores only if set to True in constructor
				outres.write("\t")
				outres.write("\t".join([str(Z[mfa][f+" varZ"]) for f in featnamesorted]))
			outres.write("\n")
		print("Completed computation of ES for IDRs in dir %s"%(self.aln_dir))	#print("TOOK %s SECONDS"%(time.time()-start_time))
