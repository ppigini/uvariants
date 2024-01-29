

CLINVAR ANALYSIS BY UVARIANTS


Uvariants is a short python script that can be used for:

• analyzing “benign” and “pathogenic” mutations in proximity of donor splice sites from ClinVar database;

• extracting “pathogenic” variants that could potentially be treated by U1 therapy.

 
###########################################################################


1.Requirements

The following languages and packages must be installed prior to run the pipeline: 

• Python (executable from macOS or Linux environments)

• Python packages: “os”, “itertools”

The pipeline has been tested using an 8-core CPU and 32 GB of RAM.


###########################################################################

 
2.Input files

The pipeline requires the complete ClinVar database as input, which can be downloaded from the following link: 

• https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

Decompress the file and extract the features by running the following command: 

  % grep -E '<ClinVarSet ID=|  <Title>|    <ClinVarAccession Acc=|      <Description>|        <SequenceLocation Assembly="GRCh38"' path_to_input_file > path_to_uvariants/input/ClinVar_features.txt

Replace "path_to_input_file" with the directory of the ClinVar database. Replace "path_to_uvariants" with the directory of the uvariants package.
 
The pipeline also requires the human genome annotation and assembly files as inputs. You can download the human annotation and assembly files (GRCh38) using the following links: 

• annotation: https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz

• assembly: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

Store these files in the respective directories "path_to_uvariants/input/annotation.gtf" and "path_to_uvariants/input/assembly.fna". Replace "path_to_uvariants" with the directory of the uvariants package. In case different sources are used, ensure that the annotation and assembly files are in the same format as those found in the above links.


###########################################################################

 
3.Executing the pipeline

The pipeline can be executed with the following commands:

  % cd /path_to_uvariants ; chmod a+x uvariants.py ; python uvariants.py

Replace "path_to_uvariants" with the directory of the uvariants package. The analysis results can be found in the “path_to_uvariants/results/” folder as:

• 1) a tab format file containing the distribution analysis (“path_to_uvariants/results/distribution.txt”);

• 2) a tab format file containing the list of “pathogenic” variants (and relative exon coordinates and sequences) that could be potentially treated by U1 therapy (“path_to_uvariants/results/exons.txt”).

Note that each execution will first delete all the files in the “path_to_uvariants/results/” folder and will then generate new files.
