# companion of cloudprp2cdf.sh

BEGIN {
  getline < "tmpdim.dat";
  nreff = $1;
  maxmom = $2;
  nphase = $3;

  print "netcdf mie.cdf {";
  print "dimensions:";
  print "    reff = "nreff";";
  print "    maxmom = "maxmom";";
  print "    nphase = "nphase";";
  print "";
  print "variables:";
  print "    double reff(reff), ext(reff), ssa(reff);";
  print "    int nmom(reff);";
  print "    float pmom(reff,nphase,maxmom);";
  print "    double wavelen, nre, nim, alpha, reffmin, reffmax;";
  print "    char distflag;";
  print "";
  print "    reff:long_name = \"effective radius\";";
  print "    reff:units = \"micron\";";
  print "    ext:long_name = \"extinction\";";
  print "    ext:units = \"1/km\";";
  print "    ssa:long_name = \"single scattering albedo\";";
  print "    ssa:units = \"\";";
  print "    nmom:long_name = \"number of Legendre moments\";";
  print "    nmom:units = \"\";";
  print "    pmom:long_name = \"moments of the phase matrix\";";
  print "    pmom:units = \"\";";
  print "    pmom:_FillValue = 0.f;";
  print "";
  while ((getline line < "tmphed.dat") > 0)
      print "    :line_"counter++" = \""line"\";";
  close ("tmphed.dat");
  print "";
  print "data:";
  getline < "tmphed.dat";
  getline < "tmphed.dat";
  getline < "tmphed.dat";
  getline < "tmphed.dat";
  getline < "tmphed.dat";
  print "    wavelen = ",$1";";
  getline < "tmphed.dat";
  print "    nre = ",$1";";
  print "    nim = ",$2";";
  getline < "tmphed.dat";
  print "    distflag = \""$1"\";";
  getline < "tmphed.dat";
  print "    alpha = ",$1";";
  getline < "tmphed.dat";
  print "    reffmin = ",$2";";
  print "    reffmax = ",$3";";
  close ("tmphed.dat");

  getline < "tmpref.dat";
  reff = substr($0,1,length($0)-2);
  print "    reff = "reff";";
  print "";
  getline < "tmpext.dat";
  ext = substr($0,1,length($0)-2);
  print "    ext = "ext";";
  print "";
  getline < "tmpssa.dat";
  ssa = substr($0,1,length($0)-2);
  print "    ssa = "ssa";";
  print "";
  getline < "tmpnmo.dat";
  nmom = substr($0,1,length($0)-2);
  print "    nmom = "nmom";";
  printf ("pmom = ");
  for (i=1; i<=nreff; i++) {
    getline < "tmpmom.cdl";
    if (i<nreff)
      pmom=$0;
    else
      pmom=substr($0,1,length($0)-2)";";
    print pmom;
  }
  print "}";
}
