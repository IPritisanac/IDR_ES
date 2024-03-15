from es_pw_sim import EVSign
#from es_pw_sim_testing import EVSign
import time
import os,sys

if (len(sys.argv)==0):
    print ("usage: run_es.py input_file.txt")

input_file=sys.argv[1]
ES=EVSign(input_file)
ES.compute_es_dir()
