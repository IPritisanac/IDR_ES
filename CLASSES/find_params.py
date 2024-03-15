import os,sys


indir=sys.argv[1]
mean={}
var_sim={}
var_obs={}

for f in os.listdir(indir):
    if f.startswith("RES_PARAM_"):
        fpath=indir+os.sep+f
        fin=open(fpath,'r')
        for line in fin:
            splitted=line.split("\t")
            if splitted[0]=="VAR_SIM:":
                var_sim.setdefault(splitted[1],[]).append([float(splitted[2]),float(splitted[3]),abs(float(splitted[4]))])
            elif splitted[0]=="VAR_OBS:":
                var_obs.setdefault(splitted[1],[]).append([float(splitted[2]),float(splitted[3]),abs(float(splitted[4]))])
            elif splitted[0]=="MEAN:":
                mean.setdefault(splitted[1],[]).append([float(splitted[2]),float(splitted[3]),abs(float(splitted[4]))])
            else:
                print(line)
import numpy as np
#print(var_sim["net_charge"])
import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 3, sharey=True, tight_layout=True)
# We can set the number of bins with the *bins* keyword argument.
physchem=['net_charge', 'WF_complexity', 'KL_hydropathy', 'isoelectric_point', 'FCR', 'ED_ratio', 'RK_ratio', 'SCD', 'my_kappa', 'my_omega']
repeats=['id=REP_Q2', 'id=REP_N2', 'id=REP_S2', 'id=REP_G2', 'id=REP_E2', 'id=REP_D2', 'id=REP_K2', 'id=REP_R2', 'id=REP_P2', 'id=REP_QN2', 'id=REP_RG2', 'id=REP_PR', 'id=REP_FG2', 'id=REP_SG2', 'id=REP_SR2', 'id=REP_KAP2', 'id=REP_PTS2']
slims_prob=['LIG_WRPW_1', 'LIG_PAM2_1','LIG_PDZ_Class_1', 'LIG_EF_ALG2_ABM_1','LIG_HOMEOBOX','LIG_PDZ_Wminus1_1','TRG_ER_FFAT_1','LIG_EH_1']
slims=['CLV_C14_Caspase3-7', 'DEG_APCC_KENBOX_2', 'DEG_Kelch_Keap1_1', 'DEG_SCF_TRCP1_1', 'DOC_ANK_TNKS_1', 'DOC_CYCLIN_RxL_1', 'DOC_MAPK_gen_1', 'DOC_MAPK_JIP1_4', 'DOC_MAPK_MEF2A_6', 'DOC_MAPK_NFAT4_5', 'DOC_PP1_RVXF_1', 'DOC_PP2A_B56_1', 'DOC_PP4_FxxP_1', 'DOC_WW_Pin1_4', 'LIG_14-3-3_CanoR_1', 'LIG_CaM_IQ_9', 'LIG_CtBP_PxDLS_1', 'LIG_EF_ALG2_ABM_1', 'LIG_EH_1', 'LIG_HCF-1_HBM_1', 'LIG_HOMEOBOX', 'LIG_KEPE_2', 'LIG_KLC1_WD_1', 'LIG_LIR_Gen_1', 'LIG_PCNA_PIPBox_1', 'LIG_PDZ_Class_1', 'LIG_PDZ_Wminus1_1', 'LIG_PTAP_UEV_1', 'LIG_PTB_Apo_2', 'LIG_PTB_Phospho_1', 'LIG_Rb_LxCxE_1', 'LIG_SH2_CRK', 'LIG_SH2_GRB2like', 'LIG_SH2_NCK_1', 'LIG_SH2_SRC', 'LIG_SH2_STAP1', 'LIG_SH2_STAT5', 'LIG_SH3_2', 'LIG_SUMO_SIM_anti_2', 'LIG_SUMO_SIM_par_1', 'LIG_WRPW_1', 'MOD_CDK_SPK_2', 'MOD_CDK_SPxK_1', 'MOD_CDK_SPxxK_3', 'MOD_CK1_1', 'MOD_CK2_1', 'MOD_DYRK1A_RPxSP_1', 'MOD_GSK3_1', 'MOD_N-GLC_1', 'MOD_NMyristoyl', 'MOD_PIKK_1', 'MOD_PKA_1', 'MOD_PKA_2', 'MOD_PKB_1', 'MOD_Plk_1', 'MOD_ProDKin_1', 'MOD_SUMO_for_1', 'MOD_SUMO_rev_2', 'TRG_ER_diArg_1', 'TRG_ER_FFAT_1', 'TRG_LysEnd_APsAcLL_1', 'TRG_NES_CRM1_1', 'TRG_NLS_MonoExtN_4']

homorep1=['A_homorep', 'C_homorep', 'D_homorep', 'E_homorep', 'F_homorep', 'G_homorep', 'H_homorep', 'I_homorep', 'K_homorep']
homorep2=['L_homorep', 'M_homorep', 'N_homorep', 'P_homorep', 'Q_homorep', 'R_homorep', 'S_homorep', 'T_homorep', 'V_homorep']
homorep3=['W_homorep', 'Y_homorep','ELASTIN_LIKE', 'FGDF']
phasesep=['FG_rich', 'PY', 'FRG', 'SGFYSG', 'PG_rich', 'ELASTIN_LIKE', 'FGDF', 'R_plus_Y', 'REP_RGG']
aacontent=['AA_S', 'AA_P', 'AA_T', 'AA_A', 'AA_H', 'AA_Q', 'AA_N', 'AA_G', 'AA_R', 'AA_C', 'AA_D', 'AA_E', 'AA_F', 'AA_I', 'AA_K', 'AA_L', 'AA_M', 'AA_V', 'AA_W', 'AA_Y']
comp=['acidic', 'basic', 'aliphatic', 'polar_fraction', 'chain_expanding', 'aromatic', 'disorder_promoting']
i=-1
n_bins=100
for key,value in var_obs.items():
    #if not key.split("m=")[0] in repeats:

    if not key in slims_prob:
        continue

    i+=1
    value=np.array(value)
    average=np.mean(value,axis=0)
    stdev=np.std(value,axis=0)
    if i<3:
        x=0
        y=i
    elif i>=3 and i<6:
        x=1
        y=i-3
    else:
        x=2
        y=i-6
    print(x,y)
    axs[x,y].hist(value[:,0], bins=n_bins)
    axs[x,y].set_xlabel(key)
    axs[x,y].set_title(r'$\mu=%.3E$, $\sigma=%.3E$'%(average[0],stdev[0]))

    if i==14:
        break
plt.show()
