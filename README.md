# IDR_ES
A python package for calculation of evolutionary signatures from homologues sequence of intrinsically disordered regions (IDRs) of proteins. Evolutionary signatures of an IDR are defined as an array/vector of evolutionary Z-scores for molecular features that can be easily computed from the protein sequence of the IDR and a set of its respective homologous sequences. For each molecular feature, a mean Z-score is computed over the set of homologs. The mean Z-score evaluates how much the mean value of the distribution of a feature over the homologous sequences varies relative to a null hypothesis.

**Main features**

o importing of FASTA file(s) with full-length protein sequences for the query species --  (e.g. 'HUMAN'); supplied by the user from the Uniprot database

o querying and filtering ENSEMBL database for suitable orthologous sequences -- based on Uniprot to ENSEMBL GENE ID mapping and Uniprot fasta files (supplied   by the user; retrivable from the Uniprot database)

o alignment of selected sequences using query as a reference -- call to the external program (e.g. MAFFT)

o extraction of IDR segments from the alignments -- based on a three column file with uniprot ID and disordered regions coordinates (supplied by the user)

o calculation of a proxy evolutionary distance (F81) between a query IDR sequence and each of the orthologous IDR sequences -- this is **not** based on a phylogenetic tree

o calculation of a null hypothesis based on the proxy evolutionary distances -- using a set of user-modifyable parameters

o calculation of evolutionary Z-scores for each of the knowledge-based features included in the current version -- can be extended by the user to include additional features

**Installation/dependencies**

The code can be run on Linux/Mac or Windows that have installed appropriate Python and site-packages versions (SciPy, NumPy). Installation steps are detailed below:

1. Install Python 3: If Python 3 is not already installed on your system, you can download and install it from the official Python website: https://www.python.org/downloads/

2. Install pip: Most Python distributions come with pip pre-installed, but you can check by running

         pip --version

in your terminal/command prompt. If it's not installed, you can follow the instructions for installation here: https://pip.pypa.io/en/stable/installation/

5. Install NumPy, SciPy: You can install  NumPy & SciPy using pip by running the following commands in your terminal/command prompt:

          pip install numpy scipy

6. Run run_es.py: Once you have Python installed along with the necessary site-packages, you can run run_es.py from the command line. **NOTE** When running the script, provide the full path to input_file.txt i.e., src/utils/input_file.txt. Type in your terminal/command prompt:

          python run_es.py src/utils/input_file.txt

Alternatively, you can use Anaconda, as detailed below.

1. Install Anaconda: Download the Anaconda distribution suitable for your operating system from the official Anaconda website: https://www.anaconda.com/products/distribution

2. Create a new environment (optional): Anaconda allows you to create isolated Python environments, which can be helpful for managing dependencies. You can create a new environment using the following command:

        conda create --name myenv python=3

3. Replace myenv with the desired name for your environment.

4. Activate the environment (if created): If you created a new environment, activate it using the following command:

  On Windows:

       conda activate myenv

  On macOS and Linux:

      source activate myenv

5. Install packages: Anaconda comes with many pre-installed packages, including NumPy and SciPy. However, you can install additional packages using conda or pip. To install NumPy and SciPy, you can run:

      conda install numpy scipy

6. Run run_es.py: Once you've installed conda and the necessary packages, go to "CLASSES" directory (i.e. IDR_ES-main/CLASSES/). Then type in your terminal/command prompt:
   
      python run_es.py src/utils/input_file.txt

Note: If you've created and activated a new environment, make sure it's activated before running the script!

**Output**
The ES calculation will create an output directory with an output file in a directory from which run_es.py is run. The output file lists IDR files from the directory specified under "align_dir" in the input file and the respective Z-scores, which are listed on the same line as the identifier of the IDR alignment file.

**Authors**
Iva Pritisanac (iva.pritisanac[at]gmail.com), Medical University of Graz 
Alan Moses (alan.moses[at]utoronto.ca), University of Toronto 
Julie Forman-Kay (forman[at]sickkids.ca), The Hospital for Sick Children

**Known bugs**

method random_int in src/core/es_pw_sim.py --> if AA probabilities given to limited precision (.2/.4 floating points) -- does not return an integer but 'NoneType' Future treatment: Exception upstream Current treatment: Provide AA probabilities with long precision
