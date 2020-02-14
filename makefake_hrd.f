*************************************************************
* MAKEFAKE_HRD - generates Gaussian deviates given an input *
* mean and standard deviation. Employs Numerical Recipes    *
* code "gasdev". The program MAKEFAKE is customized here to *
* output a pair of Gaussian deviates (presumably logT,logL  *
* for HRDs).                                                *
*                                                           *
* > f77 -o makefake_hrd makefake_hrd.f                      *
* > ./makefake_hrd                                          *
*                                                           *
* A. Martin  5/23/2006                                      *
*************************************************************

      program makefake_hrd
      implicit double precision (a-h,l-m,o-z)
      character*32 filen1
      double precision x,ux,y,uy,xpt,ypt
      integer idum,i,nmax

      write(*,*)' '
      write(*,*)' Lets make synthetic HRD data! '
      write(*,*)' '
      write(*,*)'NUMBER of synthetic data pairs? Name of OUTPUT file?'
      read(*,*)nmax,filen1
      write(*,*)'Input NEGATIVE SEED INTEGER'
      read(*,*)idum
      open(10,file=filen1,status='new')
      write(*,*)'MEAN, STANDARD DEVIATION for log(Teff):'
      read(*,*)x,ux
      write(*,*)'MEAN, STANDARD DEVIATION for log(L/Lo):'
      read(*,*)y,uy
      do i = 1,nmax,1
         xpt = x + ux*gasdev(idum)
         ypt = y + uy*gasdev(idum)
         write(10,*)xpt,ypt
      end do
      end

******************************************************************
* RANN2 - Numerical recipes random number generator.             *
* W.H. Press says about ran2: "We think that, within the limits  *
* of its floating-point precision, ran2 provides perfect random  *
* numbers; a practical definition of "perfect" is that we will   *
* pay $1000 to the first reader that convinces us otherwise.     *
******************************************************************
      function rann2(idum)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real rann2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *  ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     *  ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum .le. 0) then
         idum = max(-idum,1)
         idum2 = idum
         do j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1) - k*ir1
            if (idum .lt. 0) idum = idum + im1
            if (j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1) - k*ir1
      if (idum .lt. 0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2) - k*ir2
      if (idum2 .lt. 0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy + imm1
      rann2 = min(am*iy,rnmx)
      return
      end

*************************************************
* GASDEV - from Numerical Recipes FORTRAN p.280 *
*************************************************
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES rann2
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*rann2(idum)-1.
        v2=2.*rann2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
