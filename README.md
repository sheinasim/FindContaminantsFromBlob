# FindContaminantsFromBlob 
Print list of contigs identified as "Not Arthropod" from using the Blobtools taxrule: bestsumorder

This script is designed to be run on the output of the blobtools pipeline following adding hits from blast and or diamond blast. 

Dependencies:

* Python3 
* Python pandas package
* Python json package
* Python os package
* Python re package
* Python numpy package
* Python scipy package
* Python Biopython package

> export PATH=$PATH:[PATH TO FindContaminantsFromBlob]  

### Usage
  
> python findContaminants.py pathToBlobdir pathToAssembly (.fasta format)  pathToBuscoTable 

If no arguments are provided, the script will return help message.

## Outputs

* ./{blobdir prefix}\_contaminants.tsv 
* ./{blobdir prefix}\_keepers.tsv 
* ./{blobdir prefix}\_all.tsv

### Citation

If this script is useful to you, please cite the following in your publication:

```
@software{FindContaminantsFromBlob,
  author = {Sim, Sheina B.},
  title = {FindContaminantsFromBlob},
  url = {https://github.com/sheinasim/FindContaminantsFromBlob}
}
```

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
