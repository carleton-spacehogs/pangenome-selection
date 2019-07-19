###BELOW IS TO GET THE PN/PS CALCULATIONS
for i in SAAV_Profile/*;
do file=$(echo $i | cut -d"/" -f 2);
contig=$(echo $file | cut -d"_" -f 1);
scv=$(echo $file | cut -d"_" -f 1-3);
anvi-script-calculate-pn-ps-ratio -a $i -b SCV_Profile/${scv}_SCV.txt -c ${contig}_contigs.db -m 10 -o pn_ps_03_25/${scv}