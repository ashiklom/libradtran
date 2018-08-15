C *** MODULE radparams ***
      integer nvx
      parameter (nvx = 105)
      integer nv1x
      parameter (nv1x = nvx + 1)
      integer ndfsx, mdfsx
      parameter (ndfsx = nvx, mdfsx = nvx + 1)
      integer ndfs4x, mbx, mbsx
      parameter (ndfs4x = 4*ndfsx, mbx = 18, mbsx = 6)
      integer mbirx, ncx, mby
      parameter (mbirx = 12, ncx = 8, mby = 10)
      integer naer, nrh, mxat, mxac
      parameter (naer = 18, nrh = 8, mxat = 7, mxac = 3)
 
!#INPUTS
      logical fourssl, foursir
      integer nv, nv1, mb, mbs, mbir, nc, ndfs, mdfs, ndfs4, ndfs2
 
      integer iaform, n_atau, ivd, ifg, nac
      integer itps(3)
      real a_wlis(7,3), a_taus(7,3)
      real aprofs(100,3)
 
      logical edding, quadra, hemisp, mquadr
      common /radparams/a_wlis, a_taus, aprofs, nv, nv1, mb, mbs, mbir, 
     1   nc, ndfs, mdfs, ndfs4, ndfs2, iaform, n_atau, ivd, ifg, nac, 
     2   itps, fourssl, foursir, edding, quadra, hemisp, mquadr
C *** END MODULE radparams ***
