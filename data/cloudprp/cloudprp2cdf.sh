filename=$1

gawk 'FNR==6{nreff=$1};/Nleg/{if($4>maxmom){maxmom=$4}};END{print nreff,maxmom}' $filename > tmpdim.dat

gawk 'FNR==1{print "netCDF created by cloudprp2cdf.sh"; print "source: "FILENAME; print "header information:"}; FNR<=6 {print}; FNR==6 {exit}' $filename > tmphed.dat

gawk '/Nleg/{printf("%s, ", $1)>"tmpref.dat";printf("%s, ", $2)>"tmpext.dat";printf("%s, ", $3)>"tmpssa.dat";printf("%s, ", $4)>"tmpnmo.dat";}' $filename 

gawk '/Nleg/{nmom=$4+1;counter=0; while (counter<nmom) {getline; for (i=1;i<=NF; i++) printf("%s ",$i); counter+=NF}; printf("\n")}' $filename > tmpmom.dat

gawk 'BEGIN{maxmom='`gawk '{print $2}' tmpdim.dat`'}{for (i=1; i<=maxmom+1; i++) printf("%.5f, ",$i+0); printf("\n");}' tmpmom.dat > tmpmom.cdl

gawk -f cloudprp2cdf.awk > $filename.cdl

ncgen -o$filename.cdf $filename.cdl

rm -f tmphed.dat tmpdim.dat tmpref.dat tmpext.dat tmpssa.dat tmpnmo.dat tmpmom.dat tmpmom.cdl $filename.cdl
