c     
c     aerosol profile data, nicked from modtran (see relevant 
c     references from AFGL) and cleaned up a little.
c
c   zaerd              The altitude at which the densities are given
c  
c 0-2km
c   
c   hzviz(lc,iviz)     Haze visibility profiles 0-2km,
c                      iviz = 1 visibility 50km
c                      iviz = 2 visibility 23km
c                      iviz = 3 visibility 10km
c                      iviz = 4 visibility 5km
c                      iviz = 5 visibility 2km
c >2-10km
c                      fawi50=fall/winter   50km vis
c                      fawi23=fall/winter    23km vis
c                      spsu50=spring/summer  50km vis
c                      spsu23=spring/summer  23km vis
c >10-30km
c                      bastfw=background stratospheric   fall/winter
c                      vumofw=moderate volcanic          fall/winter
c                      hivufw=high volcanic              fall/winter
c                      exvufw=extreme volcanic           fall/winter
c                      bastss,vumoss,hivuss,exvuss=      spring/summer
c >30-100km
c                      upnatm=normal upper atmospheric
c                      vutono=transition from volcanic to normal
c                      vutoex=transition from volcanic to extreme
c                      exupat=extreme upper atmospheric
c
c
      REAL zaerd(34),
     $     hzviz(35,5), fawi50(34), fawi23(34),
     $     spsu50(34), spsu23(34), bastfw(34), vumofw(34),
     $     hivufw(34),exvufw(34),bastss(34), vumoss(34),
     $     hivuss(34),exvuss(34),upnatm(34),vutono(34),
     $     vutoex(34),exupat(34)
c
      data zaerd/
     $    0.,    1.,    2.,    3.,    4.,
     $     5.,    6.,    7.,    8.,
     $    9.,   10.,   11.,   12.,   13.,
     $     14.,   15.,   16.,   17.,
     $   18.,   19.,   20.,   21.,   22.,
     $     23.,   24.,   25.,   30.,
     $   35.,   40.,   45.,   50.,   70.,
     $     100.,99999./
c
       data hzviz(1,1),hzviz(1,2),hzviz(1,3),hzviz(1,4),hzviz(1,5)/
     $      6.62e-02,  1.58e-01,  3.79e-01,  7.70e-01,  1.94e+00/
       data hzviz(2,1),hzviz(2,2),hzviz(2,3),hzviz(2,4),hzviz(2,5)/
     $      4.15e-02,  9.91e-02,  3.79e-01,  7.70e-01,  1.94e+00/
       data hzviz(3,1),hzviz(3,2),hzviz(3,3),hzviz(3,4),hzviz(3,5)/
     $      2.60e-02,  6.21e-02,  6.21e-02,  6.21e-02,  6.21e-02/
c
      data fawi50  /3*0.,
     $ 1.14e-02, 6.43e-03, 4.85e-03, 3.54e-03, 2.31e-03, 1.41e-03,
     $ 9.80e-04,7.87e-04,23*0./
c
      data fawi23              /3*0.,
     1 2.72e-02, 1.20e-02, 4.85e-03, 3.54e-03, 2.31e-03, 1.41e-03,
     2 9.80e-04,7.87e-04,23*0./
c
      data  spsu50              / 3*0.,
     1 1.46e-02, 1.02e-02, 9.31e-03, 7.71e-03, 6.23e-03, 3.37e-03,
     2 1.82e-03  ,1.14e-03,23*0./
c
      data  spsu23              / 3*0.,
     1 3.46e-02, 1.85e-02, 9.31e-03, 7.71e-03, 6.23e-03, 3.37e-03,
     2 1.82e-03  ,1.14e-03,23*0./
c
      data bastfw       /11*0.,
     1           7.14e-04, 6.64e-04, 6.23e-04, 6.45e-04, 6.43e-04,
     2 6.41e-04, 6.00e-04, 5.62e-04, 4.91e-04, 4.23e-04, 3.52e-04,
     3 2.95e-04, 2.42e-04, 1.90e-04, 1.50e-04, 3.32e-05 ,7*0./
c
      data    vumofw       /11*0.,
     1           1.79e-03, 2.21e-03, 2.75e-03, 2.89e-03, 2.92e-03,
     2 2.73e-03, 2.46e-03, 2.10e-03, 1.71e-03, 1.35e-03, 1.09e-03,
     3 8.60e-04, 6.60e-04, 5.15e-04, 4.09e-04, 7.60e-05 ,7*0./
c
      data    hivufw       /11*0.,
     1           2.31e-03, 3.25e-03, 4.52e-03, 6.40e-03, 7.81e-03,
     2 9.42e-03, 1.07e-02, 1.10e-02, 8.60e-03, 5.10e-03, 2.70e-03,
     3 1.46e-03, 8.90e-04, 5.80e-04, 4.09e-04, 7.60e-05 ,7*0./
c
      data    exvufw       /11*0.,
     1           2.31e-03, 3.25e-03, 4.52e-03, 6.40e-03, 1.01e-02,
     2 2.35e-02, 6.10e-02, 1.00e-01, 4.00e-02, 9.15e-03, 3.13e-03,
     3 1.46e-03, 8.90e-04, 5.80e-04, 4.09e-04, 7.60e-05 ,7*0./
c
      data    bastss       /11*0.,
     1           7.99e-04, 6.41e-04, 5.17e-04, 4.42e-04, 3.95e-04,
     2 3.82e-04, 4.25e-04, 5.20e-04, 5.81e-04, 5.89e-04, 5.02e-04,
     3 4.20e-04, 3.00e-04, 1.98e-04, 1.31e-04, 3.32e-05 ,7*0./
c
      data    vumoss       /11*0.,
     1           2.12e-03, 2.45e-03, 2.80e-03, 2.89e-03, 2.92e-03,
     2 2.73e-03, 2.46e-03, 2.10e-03, 1.71e-03, 1.35e-03, 1.09e-03,
     3 8.60e-04, 6.60e-04, 5.15e-04, 4.09e-04, 7.60e-05 ,7*0./
c
      data    hivuss       /11*0.,
     1           2.12e-03, 2.45e-03, 2.80e-03, 3.60e-03, 5.23e-03,
     2 8.11e-03, 1.20e-02, 1.52e-02, 1.53e-02, 1.17e-02, 7.09e-03,
     3 4.50e-03, 2.40e-03, 1.28e-03, 7.76e-04, 7.60e-05 ,7*0./
c
      data    exvuss       /11*0.,
     1           2.12e-03, 2.45e-03, 2.80e-03, 3.60e-03, 5.23e-03,
     2 8.11e-03, 1.27e-02, 2.32e-02, 4.85e-02, 1.00e-01, 5.50e-02,
     3 6.10e-03, 2.40e-03, 1.28e-03, 7.76e-04, 7.60e-05 ,7*0./
c
      data upnatm       /26*0.,
     1 3.32e-05, 1.64e-05, 7.99e-06, 4.01e-06, 2.10e-06, 1.60e-07,
     2 9.31e-10, 0.      /
c
      data vutono       /26*0.,
     1 7.60e-05, 2.45e-05, 7.99e-06, 4.01e-06, 2.10e-06, 1.60e-07,
     2 9.31e-10, 0.      /
c
      data vutoex       /26*0.,
     1 7.60e-05, 7.20e-05, 6.95e-05, 6.60e-05, 5.04e-05, 1.03e-05,
     2 4.50e-07, 0.      /
c
      data exupat       /26*0.,
     1 3.32e-05, 4.25e-05, 5.59e-05, 6.60e-05, 5.04e-05, 1.03e-05,
     2 4.50e-07, 0.      /
c
