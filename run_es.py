from src.core.es_pw_sim import EVSign
import time
import os, sys


def main(input_file, protocol):
    protocol = protocol.upper()
    if protocol not in ("ES", "FS"):
        print("ERROR: Unknown protocol '%s'. Expected 'ES' or 'FS'." % protocol)
        sys.exit(1)

    # Instantiate EVSign (also parses input file and validates the configured paths)
    ES = EVSign(input_file)

    if protocol == "ES":
        # compute evolutionary signatures given MSAs from orthologous IDR sequences
        if ES.aln_dir == "":
            print("ERROR: protocol 'ES' requires 'align_dir' to be set in %s" % input_file)
            sys.exit(1)
        ES.compute_es_dir()
    else:  # FS
        # compute feature signatures given IDRs in a single FASTA file
        if ES.fasta_dir == "":
            print("ERROR: protocol 'FS' requires 'fasta_dir' to be set in %s" % input_file)
            sys.exit(1)
        ES.compute_fs_seq()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: run_es.py input_file.txt <ES|FS>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
