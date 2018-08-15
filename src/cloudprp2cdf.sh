#! /bin/sh

# cloudprp2cdf.sh converts an ASCII Mie table 
# (output from Frank Evans' cloudprp) to a netCDF file that 
# can be read by read_mie_table() and hence by uvspec.
# The name of the file is set in the variable filename below.

allwater=$(ls ./wc.*.mie)

for filename in $allwater
do

  gawk 'FNR==6{nreff=$1};/Nleg/{if($4>maxmom){maxmom=$4}};/Nphase/{nphase=$5};END{print nreff,maxmom+1,nphase}' $filename > tmpdim.dat

  gawk 'FNR==1{print "netCDF created by cloudprp2cdf.sh"; print "source: "FILENAME; print "header information:"}; FNR<=6 {print}; FNR==6 {exit}' $filename > tmphed.dat

  gawk '/Nleg/{printf("%s, ", $1)>"tmpref.dat";printf("%s, ", $2)>"tmpext.dat";printf("%s, ", $3)>"tmpssa.dat";printf("%s, ", $4)>"tmpnmo.dat";}' $filename 

  gawk '/Nleg/{nmom=$4; nphase=$5; counter=0; while (counter<nmom*nphase) {getline; for (i=1;i<=NF; i++) printf("%s ",$i); counter+=NF}; printf("\n")}' $filename > tmpmom.dat
   
  gawk 'BEGIN{maxmom='`gawk '{print $2}' tmpdim.dat`'; nphase='`gawk '{print $3}' tmpdim.dat`'}{for (i=1; i<=(maxmom*nphase); i++) printf("%.5f, ",$i+0); printf("\n");}' tmpmom.dat > tmpmom.cdl

  gawk -f ./cloudprp2cdf.awk > $filename.cdl

  ncgen -o $filename.cdf $filename.cdl

  rm -f tmphed.dat tmpdim.dat tmpref.dat tmpext.dat tmpssa.dat tmpnmo.dat tmpmom.dat tmpmom.cdl $filename.cdl

done


allice=$(ls ./ic.*.mie)
    
for filename in $allice
do
   gawk 'FNR==6{nreff=$1};/Nleg/{if($4>maxmom){maxmom=$4}};/Nphase/{nphase=$5};END{print nreff,maxmom+1,nphase}' $filename > tmpdim.dat

  gawk 'FNR==1{print "netCDF created by cloudprp2cdf.sh"; print "source: "FILENAME; print "header information:"}; FNR<=6 {print}; FNR==6 {exit}' $filename > tmphed.dat

  gawk '/Nleg/{printf("%s, ", $1)>"tmpref.dat";printf("%s, ", $2)>"tmpext.dat";printf("%s, ", $3)>"tmpssa.dat";printf("%s, ", $4)>"tmpnmo.dat";}' $filename 

  gawk '/Nleg/{nmom=$4; nphase=$5; counter=0; while (counter<nmom*nphase) {getline; for (i=1;i<=NF; i++) printf("%s ",$i); counter+=NF}; printf("\n")}' $filename > tmpmom.dat
    
  gawk 'BEGIN{maxmom='`gawk '{print $2}' tmpdim.dat`'; nphase='`gawk '{print $3}' tmpdim.dat`'}{for (i=1; i<=(maxmom*nphase); i++) printf("%.5f, ",$i+0); printf("\n");}' tmpmom.dat > tmpmom.cdl

  gawk -f  /export/home/emde_cl/libRadtran/src/cloudprp2cdf.awk > $filename.cdl

  ncgen -o $filename.cdf $filename.cdl

  rm -f tmphed.dat tmpdim.dat tmpref.dat tmpext.dat tmpssa.dat tmpnmo.dat tmpmom.dat tmpmom.cdl $filename.cdl

done
