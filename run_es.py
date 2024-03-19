from src.core.es_pw_sim import EVSign
import time
import os,sys

def main(input_file):
    # Instantiate EVSign

    ES=EVSign(input_file)
    ES.compute_es_dir()


if __name__ == "__main__":
    #input_file=sys.argv[1]
    if (len(sys.argv)==0):
        print ("usage: run_es.py input_file.txt")
    main(sys.argv[1])
