# countFFRs
a perl script to simulate species delimitation using haplowebs

This is perl script was used in Dellicour S, Flot J-F (2015) Delimiting species-poor datasets using single molecular markers: a study of barcode gaps, haplowebs and GMYC. Systematic Biology 64:900-908 (https://doi.org/10.1093/sysbio/syu130).

Usage: perl countffrs.pl -i <input file in sequential FASTA format> -n <number of individuals to be sampled randomly from each simulated species> -s <number of simulated species; the program will divide the aligment into s equal blocks of sequences and consider each of them as coming from a different species, meaning that the haplotypes of heterozygotes will always be drawn from the same block> -r <number of replicate samplings> -v <verbiose> -f <save test dataset as FASTA> -g <save test dataset as gzipped FASTA> -R <generate Roehl output file for Network>

To generate .list files (to be used as inputs for the R scripts to calculate %oversplitting, %overlumping and %success):
Simcoal:
for i in `ls *fasta` ; do countffrs.pl -i $i -n 100 -s 1 -r 10 -g > $i.ffrs ; done
for i in `ls *ffrs` ; do csplit -f $i $i 1 {10} ; done ; rm *00; rm *11 ; sed -i 's/;/\n/g' *ffrs[0-9]*
for k in `ls *ffrs[0-9]*` ; do awk '{for (i=1; i<=NF; i++) {print $i"\t"NR; print $i"\t"NR}}' $k | sort -n > $k.list ; done

DendroPy 3 species:
for i in `ls *.fasta` ; do countffrs.pl -i $i -n 50 -s 3 -r 10 -g > $i.ffrs ; done
for i in `ls *ffrs` ; do csplit -f $i $i 1 {10} ; done ; rm *00; rm *11 ; sed -i 's/;/\n/g' *ffrs[0-9]*
for k in `ls *ffrs[0-9]*` ; do awk '{for (i=1; i<=NF; i++) {print $i"\t"NR; print $i"\t"NR}}' $k | sort -n > $k.list ; done

Dendropy 6 species:
for i in `ls *.fasta` ; do countffrs.pl -i $i -n 30 -s 6 -r 10 -g > $i.ffrs ; done
for i in `ls *ffrs` ; do csplit -f $i $i 1 {10} ; done ; rm *00; rm *11 ; sed -i 's/;/\n/g' *ffrs[0-9]*
for k in `ls *ffrs[0-9]*` ; do awk '{for (i=1; i<=NF; i++) {print $i"\t"NR; print $i"\t"NR}}' $k | sort -n > $k.list ; done

