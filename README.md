# FindContaminantsFromBlob 
Print list of contigs identified as "Not Arthropod" from using the Blobtools taxrule: bestsumorder

Dependencies:

* Python3 
* Python pandas package
* Python json package

> export PATH=$PATH:[PATH TO FindContaminantsFromBlob]  

### Usage
  
> python findContaminants.py pathToBlobdir  

If no arguments are provided, the script will return help message.

## Outputs

* {blobdir prefix}\_contaminants.txt 
* {blobdir prefix}\_keepers.txt 

### Citation

If this script is useful to you, please cite the following in your publication:

```
@software{FindContaminantsFromBlob,
  author = {Sim, Sheina B.},
  title = {FindContaminantsFromBlob},
  url = {https://github.com/sheinasim/HiFiAdapterFilt}
}
```

Sheina B. Sim  
USDA-ARS  
US Pacific Basin Agricultural Research Service  
Hilo, Hawaii, 96720 USA  
sheina.sim@usda.gov  

This script is in the public domain in the United States per 17 U.S.C. § 105
