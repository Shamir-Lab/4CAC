# 4CAC

4CAC is a four-class classifier designed to classify contigs in mixed metagenome assemblies as phages, plasmids, prokaryotes (bacteria and archaea), microeukaryote (fungi and protozoa) or uncertain. 4CAC generates an initial four-way classification using XGBoost algorithms and further improves the classification using the assembly graph.


# Installation

4CAC requires python 3.7 or above to run. You will need the following python dependencies to run 4CAC and related support scripts.

* scikit-learn
* xgboost

To run 4CAC, please download the respository as follows.
```sh
git clone https://github.com/Shamir-Lab/4CAC.git
cd 4CAC
```

# Usage

## 1. Input 

4CAC requires the following input files:

**(1) Contig file in "fasta" format:** a set of contigs to be classified. 

**(2) Assembly grah file in "gfa" format:** the assembly graph generated by metaSPAdes or metaFlye when assembling reads to generate the input contigs. 

**(3) A path file has path information for each contig,** such as `scaffolds.path` in metaSPAdes assembly and `assembly_info.txt` in metaFlye assembly.

For contigs assembled from short reads by metaSPAdes, files `scaffolds.fasta`, `assembly_graph_with_scaffolds.gfa`, and `scaffolds.path` can be used as input.
For contigs assembled from long reads by metaFlye, files `assembly.fasta`, `assembly_graph.gfa`, `assembly_info.txt` can be used as input.

## 2. Running XGBoost classifier

4CAC generates an initial four-way classification using XGBoost algorithms.

 ```sh
 python classify_xgb.py -f <contig file> 
 ``` 
 
The command line options for this script are:

`-f/--fasta`: The contig file to be classified

`-o/--outfile`: The name of the output file. If not specified, \<input filename>.probs_xgb_4class.out

`-p/--num_processes`: The number of processes to use. Default=8

## 3. Running 4CAC classifier

4CAC further improves the XGBoost classification using the assembly graph.

 ```sh
 python classify_4CAC.py --assembler <metaSPAdes/metaFlye> --asmdir <output directory of assembly graph and path file> 
  ``` 
  
The command line options for this script are:

`--assembler`: The assembler used to generate the contig file. metaSPAdes or metaFlye

`--asmdir`: The directory of the assembly graph and the path file. 

`--classdir`: The directory of the XGBoost four-class classification file. If not specified, default = `--asmdir`

`--outdir`: The output directory. If not specified, default = `--asmdir`

## 4. An example
A small test dataset could be found under the `test_data` folder. Note that the `ground_truth_class.fasta` file under test_data is not necessary for running 4CAC. It can be used to compare the XGBoost classification and 4CAC classification with the ground truth if users are interested.


# Contacts

In case of any questions or suggestions please feel free to contact lianrong.pu@gmail.com
  
 
 
