# Scripts used for analysis of A.cardenasii introgression in Line8

## Requirements
* minimap2 

## Steps
1. Align Arahy.09 and Arahy.09_alt with Tifrunner Arahy.09 and 
   A.cardenasii Chr09 using minimap2

## Input files
* Chromosome 9 fastas for Line8, Tifrunner, and A.cardenasii
  * File names used below
    * Line8 Arahy.09: Line8_09.fa
    * Line8 Arahy.09_alt: Line8_09alt.fa
    * Tifrunner Arahy.09: Tifrunner_09.fa
    * A.cardenasii Chr09: Acardenasii_09.fa

## Alignments using minimap2
```
minimap2 -cx asm5 --eqx -t 20 Tifrunner_09.fa Line8_09.fa > L809_v_Tif09.paf

minimap2 -cx asm5 --eqx -t 20 Tifrunner_09.fa Line8_09alt.fa > \
L809alt_v_Tif09.paf

minimap2 -cx asm5 --eqx -t 20 Acardenasii_09.fa Line8_09.fa > \
L809_v_Acard09.paf

minimap2 -cx asm5 --eqx -t 20 Acardenasii_09.fa Line8_09alt.fa > \
L809alt_v_Acard09.paf
```


