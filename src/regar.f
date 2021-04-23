C=======================================================================
      subroutine s_bdesfe(xk,n,m,n0,bcov,xy,uhat,st,tauef,xtx,c,uu, 
     +     ww,aux,ipiv) 
*----------------------------------------------------------------------- 
*     This subroutine computes the covariance matrix of the regression 
*     coefficients. 
* 
*     Input: 
*          xk  : constant of the function rho1: rho1(u)=rhof(u/xk) 
*          n   : number of observations 
*          m   : number of independent variables 
*          n0  : length of the vector of initial values 
*          xy  : matrix containing the values of the filtered x's  
*          uhat: series of filtered residuals  
*          st  : series of scales of the uhat's 
* 
*     Output: 
*          bcov : covariance matrix of the regression coefficients 
*          tauef: inverse of the efficiency of the tau-estimate 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z)  
      dimension bcov(m,m),xy(n,m+1),uhat(n),st(n) 
      dimension xtx(m,m),c(m,m),uu(n),ww(n),aux(2*n) 
      integer ipiv(m)
      data zero/0.d0/
*----------------------------------------------------------------------- 
* 
*     We compute standardized residuals. 
* 
      do i=n0+1,n 
         uu(i)=uhat(i)/st(i) 
      enddo 
      call s_calsfe(uu,n,n0,sigmm,aux(1),aux(n+1)) 
      sum1=zero 
      sum2=zero 
      sum3=zero 
* 
*     The function rho2(u) of the tau estimate is rhof(u) and the function 
*     rho1(u) is rhof(u/xk). 
*     In sum1 we store the estimate (n-n0) E(rho2), in sum2 we store  
*     the estimate of (n-n0) E(psi2) and in sum3 we store the estimate  
*     of (n-n0) E(psi1). 
* 
      do i=n0+1,n 
         auxx=uu(i)/sigmm 
         sum1=sum1+s_rhoffe(auxx) 
         sum2=sum2+s_psiffe(auxx)*auxx 
         sum3=sum3+s_psiffe(auxx/xk)*(auxx/xk) 
      enddo 
* 
*     w0 is the weight of rho1 in the tau-estimate. 
* 
      w0=(2.0d0*sum1-sum2)/sum3 
      sum1=zero 
      sum2=zero 
      wwac=zero 
*  
*     In sum1 we store the estimate of (n-n0)^2 E((w0*psi1+psi2)^2). 
*     In sum2 we store the estimate of (n-n0) E(w0psi1'+psi2'). 
*     ww are the implicit weighs of the tau-estimate if it is thought as 
*     an iterative reweighted least squares estimate. 
*     In wwac we store an estimate of E(ww). 
*     a is the efficiency factor of the tau-estimate. 
* 
      do i=n0+1,n 
         auxx=uu(i)/sigmm 
         wwv=w0*s_psiffe(auxx/xk)/xk+s_psiffe(auxx) 
         sum1=sum1+wwv**2 
         ww(i-n0)=wwv/auxx   
         wwac=wwac+ww(i-n0) 
         sum2= sum2+(w0*s_dpsife(auxx/xk)/(xk**2)+s_dpsife(auxx)) 
      enddo 
      wwac=wwac/dble(n-n0) 
      tauef=dble(n-n0)*sum1/(sum2*sum2) 
      if (m.gt.0) then 
* 
*     In matrix xy we store the filtered x's and xtx is a robust 
*     estimate of xy' xy. 
*    
         do i=1,m 
            do j=1,m 
               sum=zero 
               do k=1,n-n0 
                  sum=sum+xy(k,i)*xy(k,j)*ww(k)/(st(k+n0)*st(k+n0)) 
               enddo 
               xtx(i,j)=sum/wwac 
            enddo 
         enddo    
* 
*     We obtain the inverse of xtx. 
*        
         call s_rinvfe(xtx,bcov,m,m,c,ipiv) 
         do i=1,m 
            do ii=1,m 
               bcov(i,ii)=(sigmm**2)*tauef*bcov(i,ii) 
            enddo 
         enddo 
      endif 
      return 
      end
C=======================================================================
      subroutine s_calsfe(u,n,nfirst,sout,v,aux) 
*-----------------------------------------------------------------------
*     This subroutine computes a robust M-estimate of the scale of the  
*     filtered residuals using the bisquared rho function. 
* 
*     Input:   
*           u: series of residuals 
*           n: lenght of the series 
*           nfirst: the scale estimator is computed using u(i) 
*                   from i=nfirst+1 to i=n  
* 
*     Output: 
*           sout: computed scale estimate 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      dimension u(n),v(n),aux(n) 
      data xk/.405d0/eeps/.000000001d0/maxit/10000/
*-----------------------------------------------------------------------
      b=1.625d0 
      n1=n-nfirst 
      do i=1,n 
         v(i)=dabs(u(i)) 
      enddo 
*  
*     A standardized MAD of |u| is computed.  
* 
      call s_mednfe(v(nfirst+1),n1,xmed,aux)  
      xmed=xmed/.6745d0 
      if(xmed.lt.1.d-20) xmed=1.d-20         
      sant=1.d0
*  
*     The M-estimator of scale is computed, and standardized so 
*     that MAD(|u|)=1 if u is N(0,1). The old value of the M-scale found  
*     by the algorithm is stored in sant and the new one in snu.  
*     The first value of sant=1 
* 
      do iter=1,maxit 
         sum=0.0d0 
         do i=1+nfirst,n 
            sum=sum+s_rhoffe(v(i)/(sant*xmed*xk)) 
         enddo 
         sant2=sant*sant 
         snu2=sant2*sum/(dble(n-nfirst)*b)
         snu=dsqrt(snu2) 
         dif=(snu-sant)/sant 
*  
*     If dif<ep we stop the iterations and return the original scale. 
* 
         if (dabs(dif).lt.eeps) then 
            sout=snu*xmed            
            return 
         endif 
* 
*     We make the old value equal to the new one. 
* 
         sant=snu 
      enddo 
      sout=snu*xmed 
      return 
      end 
C=======================================================================
      subroutine s_corsfe(x,bopt,n,m,idif,isp,nsd,zcor,yhat,depshat,aux, 
     +     epshat,work3)
*----------------------------------------------------------------------- 
*     This subroutine computes a series zcor whose correlogram is a  
*     robust correlogram of the differenced residuals of a REGARIMA model. 
* 
*     Input: 
*           x    : matrix of independent variables 
*           bopt : vector of regression coefficients 
*           n    : number of observations 
*           m    : number of independent variables 
*           idif : number of ordinary differences 
*           isp  : seasonal period 
*           nsd  : number of seasonal differences         
*           yhat : vector containing the cleaned input series 
* 
*     Output: 
*           zcor : series whose correlogram is a robust correlogram of the 
*                  differenced residuals 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      dimension x(n,m),bopt(m),zcor(n),yhat(n) 
      dimension depshat(n),aux(n),epshat(n),work3(2*n)
      data zero/0.d0/
*----------------------------------------------------------------------- 
* 
*     We compute the filtered errors. 
* 
      do i=1,n 
         sum=zero
         do j=1,m 
            sum=sum+x(i,j)*bopt(j) 
         enddo 
         epshat(i)= yhat(i)-sum 
      enddo 
* 
*     We compute the regular differences of the filtered errors. 
* 
      if (idif.eq.1) then        
         do i=2,n 
            depshat(i-1)=epshat(i)-epshat(i-1) 
         enddo 
      elseif(idif.eq.2) then        
         do i=3,n 
            depshat(i-2)=epshat(i)-2*epshat(i-1)+epshat(i-2) 
         enddo 
      elseif(idif.eq.0) then 
         do i=1,n 
            depshat(i)=epshat(i) 
         enddo 
      endif  
* 
*     We compute the seasonal differences of the filtered errors. 
* 
      if (nsd.eq.1) then  
         do i=isp+1,n 
            depshat(i-isp)=depshat(i)-depshat(i-isp) 
         enddo 
      elseif (nsd.eq.2) then 
         do i=2*isp+2,n 
            depshat(i-2*isp-1)=depshat(i-1)-2*depshat(i-isp-1)+ 
     +           depshat(i-2*isp-1) 
         enddo 
      endif  
      n1=n-idif-isp*nsd 
      do j=1,n1 
         aux(j)=1 
      enddo 
* 
*     We compute the series zcor whose correlogram is a robust correlogram 
*     of depshat. 
* 
      call s_rcorfe(depshat,aux,n1,0,zcor,work3) 
      return 
      end 
C=======================================================================
      subroutine s_dlpafe(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1, 
     +     wa2,dwarf)
*----------------------------------------------------------------------- 
*     Given an m by n matrix a, an n by n nonsingular diagonal 
*     matrix d, an m-vector b, and a positive number delta, 
*     the problem is to determine a value for the parameter 
*     par such that if x solves the system 
* 
*           a*x = b ,     sqrt(par)*d*x = 0 , 
* 
*     in the least squares sense, and dxnorm is the euclidean 
*     norm of d*x, then either par is zero and 
* 
*           (dxnorm-delta) .le. 0.1*delta , 
* 
*     or par is positive and 
* 
*           abs(dxnorm-delta) .le. 0.1*delta . 
* 
*     This subroutine completes the solution of the problem 
*     if it is provided with the necessary information from the 
*     qr factorization, with column pivoting, of a. That is, if 
*     a*p = q*r, where p is a permutation matrix, q has orthogonal 
*     columns, and r is an upper triangular matrix with diagonal 
*     elements of nonincreasing magnitude, then s_dlpafe expects 
*     the full upper triangle of r, the permutation matrix p, 
*     and the first n components of (q transpose)*b. On output 
*     s_dlpafe also provides an upper triangular matrix s such that 
* 
*            t   t                   t 
*           p *(a *a + par*d*d)*p = s *s . 
* 
*     s is employed within s_dlpafe and may be of separate interest. 
* 
*     Only a few iterations are generally needed for convergence 
*     of the algorithm. If, however, the limit of 10 iterations 
*     is reached, then the output par will contain the best 
*     value obtained so far. 
* 
*     Arguments
* 
*       n is a positive integer input variable set to the order of r. 
* 
*       r is an n by n array. On input the full upper triangle 
*         must contain the full upper triangle of the matrix r. 
*         On output the full upper triangle is unaltered, and the 
*         strict lower triangle contains the strict upper triangle 
*         (transposed) of the upper triangular matrix s. 
* 
*       ipvt is an integer input array of length n which defines the 
*            permutation matrix p such that a*p = q*r. Column j of p 
*            is column ipvt(j) of the identity matrix. 
* 
*       diag is an input array of length n which must contain the 
*            diagonal elements of the matrix d. 
* 
*       qtb is an input array of length n which must contain the first 
*           n elements of the vector (q transpose)*b. 
* 
*       delta is a positive input variable which specifies an upper 
*             bound on the euclidean norm of d*x. 
* 
*       par is a nonnegative variable. On input par contains an 
*           initial estimate of the Levenberg-Marquardt parameter. 
*           On output par contains the final estimate. 
* 
*       x is an output array of length n which contains the least 
*         squares solution of the system a*x = b, sqrt(par)*d*x = 0, 
*         for the output par. 
* 
*       sdiag is an output array of length n which contains the 
*             diagonal elements of the upper triangular matrix s. 
* 
*       wa1 and wa2 are work arrays of length n. 
* 
*       dwarf is the smallest positive magnitude. 
*-----------------------------------------------------------------------   
      implicit double precision (a-h,o-z)
      integer ipvt(n) 
      dimension r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),wa2(n)  
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/ 
*-----------------------------------------------------------------------   
*
*     Compute and store in x the Gauss-Newton direction. If the 
*     jacobian is rank-deficient, obtain a least squares solution. 
*     
      nsing = n 
      do 10 j = 1, n 
         wa1(j) = qtb(j) 
         if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1 
         if (nsing .lt. n) wa1(j) = zero 
 10   continue 
      if (nsing .lt. 1) go to 50 
      do 40 k = 1, nsing 
         j = nsing - k + 1 
         wa1(j) = wa1(j)/r(j,j) 
         temp = wa1(j) 
         jm1 = j - 1 
         if (jm1 .lt. 1) go to 30 
         do 20 i = 1, jm1 
            wa1(i) = wa1(i) - r(i,j)*temp 
 20      continue 
 30      continue 
 40   continue 
 50   continue 
      do 60 j = 1, n 
         l = ipvt(j) 
         x(l) = wa1(j) 
 60   continue 
* 
*     Initialize the iteration counter. 
*     Evaluate the function at the origin, and test 
*     for acceptance of the Gauss-Newton direction. 
* 
      iter = 0 
      do 70 j = 1, n 
         wa2(j) = diag(j)*x(j) 
 70   continue 
      dxnorm = s_dnrmfe(n,wa2) 
      fp = dxnorm - delta 
      if (fp .le. p1*delta) go to 220 
* 
*     If the jacobian is not rank deficient, the Newton 
*     step provides a lower bound, parl, for the zero of 
*     the function. Otherwise set this bound to zero. 
* 
      parl = zero 
      if (nsing .lt. n) go to 120 
      do 80 j = 1, n 
         l = ipvt(j) 
         wa1(j) = diag(l)*(wa2(l)/dxnorm) 
 80   continue 
      do 110 j = 1, n 
         sum = zero 
         jm1 = j - 1 
         if (jm1 .lt. 1) go to 100 
         do 90 i = 1, jm1 
            sum = sum + r(i,j)*wa1(i) 
 90      continue 
 100     continue 
         wa1(j) = (wa1(j) - sum)/r(j,j) 
 110  continue 
      temp = s_dnrmfe(n,wa1) 
      parl = ((fp/delta)/temp)/temp 
 120  continue 
* 
*     Calculate an upper bound, paru, for the zero of the function. 
* 
      do 140 j = 1, n 
         sum = zero 
         do 130 i = 1, j 
            sum = sum + r(i,j)*qtb(i) 
 130     continue 
         l = ipvt(j) 
         wa1(j) = sum/diag(l) 
 140  continue 
      gnorm = s_dnrmfe(n,wa1) 
      paru = gnorm/delta 
      if (paru .eq. zero) paru = dwarf/dmin1(delta,p1) 
* 
*     If the input par lies outside of the interval (parl,paru), 
*     set par to the closer endpoint. 
* 
      par = dmax1(par,parl) 
      par = dmin1(par,paru) 
      if (par .eq. zero) par = gnorm/dxnorm 
* 
*     Beginning of an iteration. 
* 
 150  continue 
      iter = iter + 1 
* 
*     Evaluate the function at the current value of par. 
* 
      if (par .eq. zero) par = dmax1(dwarf,p001*paru) 
      temp = dsqrt(par) 
      do 160 j = 1, n 
         wa1(j) = temp*diag(j) 
 160  continue 
      call s_dqrsfe(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2) 
      do 170 j = 1, n 
         wa2(j) = diag(j)*x(j) 
 170  continue 
      dxnorm = s_dnrmfe(n,wa2) 
      temp = fp 
      fp = dxnorm - delta 
* 
*     If the function is small enough, accept the current value 
*     of par. Also test for the exceptional cases where parl 
*     is zero or the number of iterations has reached 10. 
* 
      if (dabs(fp) .le. p1*delta 
     +     .or. parl .eq. zero .and. fp .le. temp 
     +     .and. temp .lt. zero .or. iter .eq. 10) go to 220 
* 
*     Compute the Newton correction. 
* 
      do 180 j = 1, n 
         l = ipvt(j) 
         wa1(j) = diag(l)*(wa2(l)/dxnorm) 
 180  continue 
      do 210 j = 1, n 
         wa1(j) = wa1(j)/sdiag(j) 
         temp = wa1(j) 
         jp1 = j + 1 
         if (n .lt. jp1) go to 200 
         do 190 i = jp1, n 
            wa1(i) = wa1(i) - r(i,j)*temp 
 190     continue 
 200     continue 
 210  continue 
      temp = s_dnrmfe(n,wa1) 
      parc = ((fp/delta)/temp)/temp 
* 
*     Depending on the sign of the function, update parl or paru. 
* 
      if (fp .gt. zero) parl = dmax1(parl,par) 
      if (fp .lt. zero) paru = dmin1(paru,par) 
* 
*     Compute an improved estimate for par. 
* 
      par = dmax1(parl,par+parc) 
* 
*     End of an iteration. 
* 
      go to 150 
 220  continue 
* 
*     Termination. 
* 
      if (iter .eq. 0) par = zero 
      return 
      end 
C=======================================================================
      function s_dnrmfe(n,x)
*----------------------------------------------------------------------- 
*     Given an n-vector x, this function calculates the  
*     euclidean norm of x.  
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(n)
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
*-----------------------------------------------------------------------
      s1 = zero  
      s2 = zero  
      s3 = zero  
      x1max = zero  
      x3max = zero  
      floatn = dble(n)
      agiant = rgiant/floatn  
      do 90 i = 1, n  
         xabs = dabs(x(i))  
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70  
         if (xabs .le. rdwarf) go to 30  
*  
*     Sum for large components.  
*  
         if (xabs .le. x1max) go to 10  
         s1 = one + s1*(x1max/xabs)**2  
         x1max = xabs  
         go to 20  
 10      continue  
         s1 = s1 + (xabs/x1max)**2  
 20      continue  
         go to 60  
 30      continue  
*  
*     Sum for small components.  
*     
         if (xabs .le. x3max) go to 40  
         s3 = one + s3*(x3max/xabs)**2  
         x3max = xabs  
         go to 50  
 40      continue  
         if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2  
 50      continue  
 60      continue  
         go to 80  
 70      continue  
*  
*     Sum for intermediate components.  
*     
         s2 = s2 + xabs**2  
 80      continue  
 90   continue  
*  
*     Calculation of norm.  
*  
      if (s1 .eq. zero) go to 100  
      s_dnrmfe = x1max*dsqrt(s1+(s2/x1max)/x1max)  
      go to 130  
 100  continue  
      if (s2 .eq. zero) go to 110  
      if (s2 .ge. x3max)  
     +     s_dnrmfe = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))  
      if (s2 .lt. x3max)  
     +     s_dnrmfe = dsqrt(x3max*((s2/x3max)+(x3max*s3)))  
      go to 120  
 110  continue  
      s_dnrmfe = x3max*dsqrt(s3)  
 120  continue  
 130  continue  
      return    
      end  
C=======================================================================
      function s_dpsife(x) 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
*-----------------------------------------------------------------------
      g1= -1.944d0 
      g2=  1.728d0 
      g3=  -.312d0 
      g4=   .016d0 
      ax=dabs(x) 
      if (ax.gt.3.0d0) then 
         s_dpsife=0.0d0 
      else if(ax.le.2.d0) then 
         s_dpsife=1.d0 
c     else if(ax.gt.2.d0.and.ax.le.3.d0) then 
      else
         s_dpsife=7.d0*g4*(x**6)+5.d0*g3*(x**4)+3.d0*g2*(x**2)+g1 
      endif 
      return  
      end     
C=======================================================================
      subroutine s_dqrffe(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa,
     +     epsmch)
*-----------------------------------------------------------------------
*     This subroutine uses Householder transformations with column 
*     pivoting (optional) to compute a qr factorization of the 
*     m by n matrix a. It determines an orthogonal matrix q,
*     a permutation matrix p, and an upper trapezoidal 
*     matrix r with diagonal elements of nonincreasing magnitude, 
*     such that a*p = q*r. The Householder transformation for 
*     column k, k = 1,2,...,min(m,n), is of the form 
*                            t 
*           i - (1/u(k))*u*u 
* 
*     where u has zeros in the first k-1 positions. The form of 
*     this transformation and the method of pivoting first 
*     appeared in the corresponding Linpack subroutine. 
* 
*     Arguments
* 
*       a is an m by n array. On input a contains the matrix for 
*         which the qr factorization is to be computed. On output 
*         the strict upper trapezoidal part of a contains the strict 
*         upper trapezoidal part of r, and the lower trapezoidal 
*         part of a contains a factored form of q (the non-trivial 
*         elements of the u vectors described above). 
* 
*       pivot is a logical input variable. If pivot is set true, 
*             then column pivoting is enforced. If pivot is set false, 
*             then no column pivoting is done. 
* 
*       ipvt is an integer output array of length lipvt. ipvt 
*            defines the permutation matrix p such that a*p = q*r. 
*            column j of p is column ipvt(j) of the identity matrix. 
*            If pivot is false, ipvt is not referenced. 
* 
*       lipvt is a positive integer input variable. If pivot is false, 
*             then lipvt may be as small as 1. If pivot is true, then 
*             lipvt must be at least n. 
* 
*       rdiag is an output array of length n which contains the 
*             diagonal elements of r. 
* 
*       acnorm is an output array of length n which contains the 
*               norms of the corresponding columns of the input matrix a. 
*               If this information is not needed, then acnorm can coincide 
*               with rdiag. 
* 
*       wa is a work array of length n. If pivot is false, then wa 
*          can coincide with rdiag. 
*
*       epsmch is the machine precision. 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical pivot  
      integer ipvt(lipvt) 
      dimension a(lda,n),rdiag(n),acnorm(n),wa(n)  
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/ 
*-----------------------------------------------------------------------  
* 
*     Compute the initial column norms and initialize several arrays. 
* 
      do 10 j = 1, n 
         acnorm(j) = s_dnrmfe(m,a(1,j)) 
         rdiag(j) = acnorm(j) 
         wa(j) = rdiag(j) 
         if (pivot) ipvt(j) = j 
 10   continue 
* 
*     Reduce a to r with Householder transformations. 
* 
      minmn = min0(m,n) 
      do 110 j = 1, minmn 
         if (.not.pivot) go to 40 
* 
*     Bring the column of largest norm into the pivot position. 
* 
         kmax = j 
         do 20 k = j, n 
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k 
 20      continue 
         if (kmax .eq. j) go to 40 
         do 30 i = 1, m 
            temp = a(i,j) 
            a(i,j) = a(i,kmax) 
            a(i,kmax) = temp 
 30      continue 
         rdiag(kmax) = rdiag(j) 
         wa(kmax) = wa(j) 
         k = ipvt(j) 
         ipvt(j) = ipvt(kmax) 
         ipvt(kmax) = k 
 40      continue 
* 
*     Compute the Householder transformation to reduce the 
*     j-th column of a to a multiple of the j-th unit vector. 
*     
         ajnorm = s_dnrmfe(m-j+1,a(j,j)) 
         if (ajnorm .eq. zero) go to 100 
         if (a(j,j) .lt. zero) ajnorm = -ajnorm 
         do 50 i = j, m 
            a(i,j) = a(i,j)/ajnorm 
 50      continue 
         a(j,j) = a(j,j) + one 
* 
*     Apply the transformation to the remaining columns 
*     and update the norms. 
* 
         jp1 = j + 1 
         if (n .lt. jp1) go to 100 
         do 90 k = jp1, n 
            sum = zero 
            do 60 i = j, m 
               sum = sum + a(i,j)*a(i,k) 
 60         continue 
            temp = sum/a(j,j) 
            do 70 i = j, m 
               a(i,k) = a(i,k) - temp*a(i,j) 
 70         continue 
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80 
            temp = a(j,k)/rdiag(k) 
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2)) 
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80 
            rdiag(k) = s_dnrmfe(m-j,a(jp1,k)) 
            wa(k) = rdiag(k) 
 80         continue 
 90      continue 
 100     continue 
         rdiag(j) = -ajnorm 
 110  continue 
      return 
      end 
C=======================================================================
      subroutine s_dqrsfe(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa) 
*-----------------------------------------------------------------------   
*     Given an m by n matrix a, an n by n diagonal matrix d, 
*     and an m-vector b, the problem is to determine an x which 
*     solves the system 
*           a*x = b ,     d*x = 0 , 
*     in the least squares sense. 
* 
*     This subroutine completes the solution of the problem 
*     if it is provided with the necessary information from the 
*     qr factorization, with column pivoting, of a. That is, if 
*     a*p = q*r, where p is a permutation matrix, q has orthogonal 
*     columns, and r is an upper triangular matrix with diagonal 
*     elements of nonincreasing magnitude, then dqrsolv expects 
*     the full upper triangle of r, the permutation matrix p, 
*     and the first n components of (q transpose)*b. The system 
*     a*x = b, d*x = 0, is then equivalent to 
* 
*                  t       t 
*           r*z = q *b ,  p *d*p*z = 0 , 
* 
*     where x = p*z. If this system does not have full rank, 
*     then a least squares solution is obtained. On output dqrsolv 
*     also provides an upper triangular matrix s such that 
* 
*            t   t               t 
*           p *(a *a + d*d)*p = s *s . 
* 
*     s is computed within dqrsolv and may be of separate interest. 
* 
*     Arguments
*
*       n is a positive integer input variable set to the order of r. 
* 
*       r is an n by n array. On input the full upper triangle 
*         must contain the full upper triangle of the matrix r. 
*         On output the full upper triangle is unaltered, and the 
*         strict lower triangle contains the strict upper triangle 
*         (transposed) of the upper triangular matrix s. 
* 
*       ldr is a positive integer input variable not less than n 
*           which specifies the leading dimension of the array r. 
* 
*       ipvt is an integer input array of length n which defines the 
*            permutation matrix p such that a*p = q*r. Column j of p 
*            is column ipvt(j) of the identity matrix. 
* 
*       diag is an input array of length n which must contain the 
*            diagonal elements of the matrix d. 
* 
*       qtb is an input array of length n which must contain the first 
*           n elements of the vector (q transpose)*b. 
* 
*       x is an output array of length n which contains the least 
*         squares solution of the system a*x = b, d*x = 0. 
* 
*       sdiag is an output array of length n which contains the 
*             diagonal elements of the upper triangular matrix s. 
* 
*       wa is a work array of length n. 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z)
      integer ipvt(n) 
      dimension r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)  
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/ 
*-----------------------------------------------------------------------  
* 
*     Copy r and (q transpose)*b to preserve input and initialize s. 
*     In particular, save the diagonal elements of r in x. 
* 
      do 20 j = 1, n 
         do 10 i = j, n 
            r(i,j) = r(j,i) 
 10      continue 
         x(j) = r(j,j) 
         wa(j) = qtb(j) 
 20   continue 
* 
*     Eliminate the diagonal matrix d using a givens rotation. 
* 
      do 100 j = 1, n 
* 
*     Prepare the row of d to be eliminated, locating the 
*     diagonal element using p from the qr factorization. 
* 
         l = ipvt(j) 
         if (diag(l) .eq. zero) go to 90 
         do 30 k = j, n 
            sdiag(k) = zero 
 30      continue 
         sdiag(j) = diag(l) 
*     
*     The transformations to eliminate the row of d 
*     modify only a single element of (q transpose)*b 
*     beyond the first n, which is initially zero. 
* 
         qtbpj = zero 
         do 80 k = j, n 
* 
*     Determine a givens rotation which eliminates the 
*     appropriate element in the current row of d. 
* 
            if (sdiag(k) .eq. zero) go to 70 
            if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40 
            cotan = r(k,k)/sdiag(k) 
            sin = p5/dsqrt(p25+p25*cotan**2) 
            cos = sin*cotan 
            go to 50 
 40         continue 
            tan = sdiag(k)/r(k,k) 
            cos = p5/dsqrt(p25+p25*tan**2) 
            sin = cos*tan 
 50         continue 
* 
*     Compute the modified diagonal element of r and 
*     the modified element of ((q transpose)*b,0). 
*     
            r(k,k) = cos*r(k,k) + sin*sdiag(k) 
            temp = cos*wa(k) + sin*qtbpj 
            qtbpj = -sin*wa(k) + cos*qtbpj 
            wa(k) = temp 
*     
*     Accumulate the tranformation in the row of s. 
*     
            kp1 = k + 1 
            if (n .lt. kp1) go to 70 
            do 60 i = kp1, n 
               temp = cos*r(i,k) + sin*sdiag(i) 
               sdiag(i) = -sin*r(i,k) + cos*sdiag(i) 
               r(i,k) = temp 
 60         continue 
 70         continue 
 80      continue 
 90      continue 
* 
*     Store the diagonal element of s and restore 
*     the corresponding diagonal element of r. 
*     
         sdiag(j) = r(j,j) 
         r(j,j) = x(j) 
 100  continue 
* 
*     Solve the triangular system for z. If the system is 
*     singular, then obtain a least squares solution. 
* 
      nsing = n 
      do 110 j = 1, n 
         if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1 
         if (nsing .lt. n) wa(j) = zero 
 110  continue 
      if (nsing .lt. 1) go to 150 
      do 140 k = 1, nsing 
         j = nsing - k + 1 
         sum = zero 
         jp1 = j + 1 
         if (nsing .lt. jp1) go to 130 
         do 120 i = jp1, nsing 
            sum = sum + r(i,j)*wa(i) 
 120     continue 
 130     continue 
         wa(j) = (wa(j) - sum)/sdiag(j) 
 140  continue 
 150  continue 
* 
*     Permute the components of z back to components of x. 
* 
      do 160 j = 1, n 
         l = ipvt(j) 
         x(l) = wa(j) 
 160  continue 
      return 
      end 
C=======================================================================
      subroutine s_durbfe(rho,lp,partial,ier,phi,ndim2) 
*-----------------------------------------------------------------------
*     This subroutine computes lp partial autocorrelations given  
*     lp autocorrelations using the Durbin algorithm.  
* 
*     Input: 
*          rho  : vector containing the autocorrelations 
*          lp   : length of rho 
*          ndim2: max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                 dimension the auxiliary array phi 
*           
*     Output: 
*          partial: vector containing the lp partial autocorrelations 
*          ier    : error code. If ier=0 all the partial autocorrelations  
*                   are not larger than 1. If ier=1 at least one partial  
*                   autocorrelation is greater or equal than 1. 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      dimension rho(lp),partial(lp),phi(ndim2,ndim2) 
*-----------------------------------------------------------------------
      ier=0 
      phi(1,1)=rho(1) 
      do i=2,lp 
         a=rho(i) 
         do m=1,i-1 
            a=a-rho(m)*phi(i-1,i-m) 
         enddo 
         b=1.0d0 
         do m=1,i-1 
            b=b-rho(m)*phi(i-1,m) 
         enddo 
         phi(i,i)=a/b 
         do j=1,i-1 
            phi(i,j)=phi(i-1,j)-phi(i-1,i-j)*phi(i,i) 
         enddo 
      enddo 
      do i=1,lp 
         partial(i)=phi(i,i)  
      enddo 
* 
*     If |partial(i)|>1, ier=1. 
*               
      do i=1,lp 
         if(dabs(partial(i)).gt.1.0d0) ier=1 
      enddo 
      return 
      end 
C=======================================================================
      subroutine s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,np,nq,n0, 
     +     indth,npo,sigman,sigmau,npred,x,y,xy,yhat,cck,uhat,epshat,
     +     st,epspred,w,auxm,poldif,ndim1,ndim2,work,nw,work4,nw4,
     +     iwork4,niw4,work5,nw5,iwork5,niw5) 
*-----------------------------------------------------------------------     
*     This subroutine distributes the work vector work of lengths nw 
*     between different arrays and vectors to be used by s_fnc1fe, and  
*     then calls it. 
*-----------------------------------------------------------------------    
      implicit double precision (a-h,o-z) 
      dimension par(ndim1),f(n),x(n,m),y(n),xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred) 
      dimension epspred(n+npred),w(n+npred),poldif(isp*nsd+idif+1) 
      dimension auxm(n+npred,ndim2),work(nw),work4(nw4),work5(nw5) 
      integer npo(np+nq),iwork4(niw4),iwork5(niw5)       
*-----------------------------------------------------------------------    
      n1=1 
      n2=ndim2 
      n3=n2+nq 
      n4=n3+np+1 
      n5=n4+ndim2+1 
      n6=n5+m 
      n7=n6+n 
      n8=n7+ndim2+1 
      n9=n8+ndim2+1 
      n10=n9+ndim2 
      n11=n10+ndim2 
      n12=n11+ndim2 
      n13=n12+n 
      call s_fnc1fe(n,npar,par,f,iflag,idif,isp,nsd,m,np,nq,n0,indth, 
     +     npo,sigman,sigmau,npred,x,y,xy,yhat,uhat,epshat,st, 
     +     epspred,w,auxm,poldif,ndim1,ndim2,work(n1),work(n2+1), 
     +     work(n3+1),work(n4+1),work(n5+1),work(n6+1), 
     +     work(n7+1),work(n8+1),work(n9+1),work(n10+1), 
     +     work(n11+1),work(n12+1),work(n13+1),cck,work4,nw4, 
     +     iwork4,niw4,work5,nw5,iwork5,niw5) 
      return 
      end 
C=======================================================================
      subroutine s_fc12fe(phi,theta,thetas,n,beta,cck,idif,isp,nsd,m,
     +     np,nq,n0,indth,x,y,sigman,sigmau,vtau,sigini,tau,xy,yhat,
     +     uhat,epshat,st,epspred,w,auxm,npred,poldif,ndim2,work,nw,
     +     work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
*-----------------------------------------------------------------------    
*     This subroutine distributes the work vector work of lengths nw 
*     between different arrays and vectors to be used by s_fnc2fe, and  
*     then calls it. 
*-----------------------------------------------------------------------     
      implicit double precision (a-h,o-z)    
      dimension phi(ndim2),theta(nq),beta(m),x(n,m),y(n),tau(0:ndim2) 
      dimension xy(n,m+1),yhat(n),uhat(n+npred),epshat(n+npred), 
     +     st(n+npred),epspred(n+npred),w(n+npred), 
     +     auxm(n+npred,ndim2),poldif(isp*nsd+idif+1) 
      dimension work(nw),work4(nw4),work5(nw5) 
      integer iwork4(niw4),iwork5(niw5) 
*-----------------------------------------------------------------------    
      n1=1 
      n2=n 
      n3=n2+ndim2+1 
      n4=n3+ndim2+1 
      n5=n4+n 
      n6=n5+ndim2+1 
      n7=n6+4*n 
      call s_fnc2fe(phi,theta,thetas,n,beta,cck,idif,isp,nsd,m,np,nq,
     +     n0,indth,x,y,sigman,sigmau,vtau,sigini,tau,xy,yhat,uhat, 
     +     epshat,st,epspred,w,auxm,npred,poldif,ndim2,work(n1), 
     +     work(n2+1),work(n3+1),work(n4+1),work(n5+1), 
     +     work(n6+1),work(n7+1),work4,nw4,iwork4,niw4,work5, 
     +     nw5,iwork5,niw5) 
      return 
      end 
C=======================================================================
      subroutine s_flt1fe(x,y,n,m,idif,isp,nsd,phi,beta,theta, 
     +     thetas,k,iq,sigmau,indth,n0,tau,sigmadif, 
     +     indfil,rho,cck,npred,ypure,xy,yhat,uhat, 
     +     epshat,st,epspred,w,auxm,ndim2,work,nw, 
     +     work5,nw5,iwork5,niw5) 
*-----------------------------------------------------------------------    
*     This subroutine distributes the work vectors work of length nw 
*     between severak arrays and vectors to be used by filter, and then  
*     calls it. 
*      
*     The first 32 arguments of this subroutine are described in filter. 
*     work5 is a real auxiliary vector of length nw5 and iwork5 is an 
*     integer auxiliary vector of length niw5. 
*-----------------------------------------------------------------------    
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phi(ndim2+1),beta(m),theta(iq) 
      dimension rho(0:ndim2),tau(0:ndim2),ypure(n),xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred),w(n+npred),auxm(n+npred,ndim2) 
      dimension work(nw),work5(nw5) 
      dimension iwork5(niw5)
*-----------------------------------------------------------------------     
      n1=1 
      n2=ndim2*ndim2 
      n3=n2+ndim2*ndim2 
      n4=n3+ndim2 
      n5=n4+ndim2 
      n6=n5+ndim2 
      n7=n6+ndim2 
      n8=n7+ndim2*ndim2 
      n9=n8+ndim2 
      n10=n9+ndim2*ndim2 
      n11=n10+ndim2 
      n12=n11+ndim2*ndim2 
      n13=n12+iq+1 
      call s_fltrfe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,k,iq, 
     +     sigmau,indth,n0,tau,sigmadif,indfil,rho,cck,npred, 
     +     ypure,xy,yhat,uhat,epshat,st,epspred,w,auxm,ndim2, 
     +     work(n1),work(n2+1),work(n3+1),work(n4+1),work(n5+1), 
     +     work(n6+1),work(n7+1),work(n8+1),work(n9+1), 
     +     work(n10+1),work(n11+1),work(n12+1),work(n13+1), 
     +     work5,nw5,iwork5,niw5) 
      return 
      end 
C======================================================================= 
      subroutine s_flt2fe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,
     +     k,iq,sigmau,indth,n0,tau,sigmadif,indfil,rho,cck,iout,
     +     yhat,uhat,epshat,st,lscan,w,auxm,xy,m1,thetanew,xm,xp,
     +     alfapred,r,v,av,alfafi,rv,alfaux,xaux,v1,v2,ypure,
     +     idim,work3,idimw3,iwork,idimiw)
*----------------------------------------------------------------------- 
*     This subroutine produces a filtered version of the input series 
*     y(t). It also produces a series of filtered residuals and a series 
*     of standard deviations of the prediction errors.
*
*     input  
*             x        : matrix of independent variables.
*             y        : input series 
*             n        : number of observations
*             m        : number of independent variables
*             idif     : number of ordinary differences
*             isp      : seasonal period
*             nsd      : number of seasonal differences
*             phi      : vector of coefficients of the AR models
*             beta     : vector of coeff. of the independent variables        
*             theta    : vector of ordinary MA coefficients
*             thetas   : seasonal moving average coefficient
*             k        : order of the non stationary AR model to be used 
*             q        : order of the ordinary MA model 
*             sigmau   : initial scale
*             indth    : 0 - no seasonal moving average component
*                        1 - seasonal moving average term included
*             n0       : length of the vector of initial values
*             tau      : vector containing e(u(t-j+1) z(t))/e(z(t)^2)    
*             sigmadif : scale of the differenced series
*             indfil   : filtering selection
*                        1 - filtering of y and x with fixed w and st
*                        0 - filtering of epshat
*             rho      : autocorrelation function of the differenced series
*             iout     : the filtering process is conducted for the 
*                        observations t=1,...,iout-1, then we put
*                        w(t)=1 for t=iout,iout+1,...,iout+h and
*                        finally we continue the filtering process
*                        from then to the last observation.  if iout=0 
*                        or iout=n+1, we filter all
*                        the observations.
*              
*     output
*             yhat   : filtered series 
*             uhat   : series of filtered residuals 
*             epshat : series of filtered errors
*             st     : series of standard deviations of the  prediction
*                      errors
*             lscan  : vector  containing  the  position of  possible 
*                      level shifts 
*             w      : series of weights
*             auxm   : matrix containing the first columns of matrix xm 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      dimension x(n,m),y(n),xy(n,m1+1),beta(m),phi(idim+1),theta(iq)
      dimension rho(0:idim),tau(0:iq+isp*indth),w(n),auxm(n,idim)
      dimension yhat(n),uhat(n),epshat(n),st(n),thetanew(iq+isp*indth)
      dimension xm(idim,idim),xp(idim,idim),ypure(n),work3(idimw3)
      dimension alfapred(idim),r(idim),v(idim),av(idim,idim)
      dimension alfafi(idim),rv(idim,idim),alfaux(idim) 
      dimension xaux(idim,idim),v1(iq+isp*indth+1),v2(iq+isp*indth+1)
      integer h,lscan(n),iwork(idimiw)
      data control,zero/2.d5,0.d0/
*-----------------------------------------------------------------------

c initialize variables defined within if statements
      res = 0.d0      

      h=2
      if (iout.eq.0) iout=n+1
      do i=1,n
        lscan(i)=0
      enddo
      ilscan=0
      do i=1,idim
         r(i)=zero
         do j=1,idim
            rv(i,j)=zero
         enddo
      enddo         
      do i=1,iq+isp*indth+1
         v1(i)=zero
         v2(i)=zero
      enddo
      v1(1)=1.d0
      do i=1,iq
         v1(i+1)=-theta(i)
      enddo
      v2(1)=1.d0
      if (indth.eq.1) v2(isp+1)=-thetas
      nv2=indth*isp 
      call s_polyfe(v1,iq,v2,nv2,r,nrdim)  
*     
*     we construct vector thetanew to be used by s_rinife.
*
      do i=1,nrdim
         thetanew(i)=-r(i+1)
      enddo
      lfin=max0(k,nrdim+1)
      do i=k+1, lfin
         phi(i)=zero
      enddo
      do i=1,lfin
         do j=1,lfin
            rv(i,j)=r(i)*r(j)
         enddo
      enddo
*
*     call to subroutine ini which returns initial values of vector 
*     alfafi and matrix xp.
*     previously, we must construct the series ypure= y-x'*beta and the 
*     vector of autocorrelations rho. 
*
*     if indfil=0 it begins the loop to filter x and y.
*
      nloop=1
      if (indfil.eq.1) nloop=m+1
      do mm=1,nloop
         do i=1,n
            if (indfil.eq.0) then      
               ypure(i)=y(i)
               do j=1,m
                  ypure(i)=ypure(i)-x(i,j)*beta(j)
               enddo
            else
               if (mm.eq.nloop) then    
                  ypure(i)=y(i)
               else
                  ypure(i)=x(i,mm)
               endif
            endif
         enddo
*     
*     We compute the position in the work vector work3 for the 
*     arrays used in subroutine ini
*
         iw3covu=1
         iw3uhat=iw3covu+idim*idim
         iw3uuha=iw3uhat+idim
         iw3w=iw3uuha+idim
         iw3rhom=iw3w+idim
         iw3b=iw3rhom+idim*idim
         iw3alfa=iw3b+idim*idim
         iw3rhoi=iw3alfa+idim
         call s_rinife(ypure,n,lfin,phi,thetanew,isp,rho,tau,alfafi,
     +        idif,nsd,xp,nrdim,n0,sigmau,sigmadif,idim,
     +        work3(iw3alfa),work3(iw3rhoi),work3(iw3covu),
     +        work3(iw3uhat),work3(iw3uuha),work3(iw3w),
     +        work3(iw3rhom),work3(iw3b),iwork)        
*                          
*     numout will count the number of succesive outliers.
*     jout will be the index of the first outlier in the group.
*
         jout=0           
         numout=0 
*     
*     inicialization of yhat(i), w(i) and st(i) for i=1,n0.
*
         if (indfil.eq.0) then
            do it=1,n0
               yhat(it)=y(it)
               epshat(it)=ypure(it)
               uhat(it)=zero
            enddo
            w(n0)=1.0d0
         else
            if (mm.eq.nloop) then
               do it=1,n0
                  yhat(it)=y(it)
                  epshat(it)=y(it)
               enddo
            else 
               do it=1,n0
                  yhat(it)=x(it,mm)
                  epshat(it)=x(it,mm)
               enddo
            endif
         endif
*
*     cycle of the observations.
*
         it=n0+1
         do while (it.le.n)
*                    
*     calculus of alfapred and uhat(it). 
*
            if (indth.eq.1.or.iq.gt.0) then
               do ir=1,lfin-1
                  alfapred(ir)=phi(ir)*alfafi(1)+alfafi(ir+1)
               enddo
               alfapred(lfin)=phi(lfin)*alfafi(1)
               uhat(it)=ypure(it)-alfapred(1)
            else
               do ir=1,lfin-1
                  alfapred(ir+1)=alfafi(ir)
               enddo
               sum2=zero
               do ir=1,lfin
                  sum2=sum2+phi(ir)*alfafi(ir)
               enddo
               alfapred(1)=sum2
               uhat(it)=ypure(it)-alfapred(1)
            endif
*
*     updating of matrix xm (if indfil=0).
*     
            if (indfil.eq.0) then 
               if (indth.eq.0.and.iq.eq.0) then
                  do i=2,k           
                     do j=2,k           
                        xm(i,j)=xp(i-1,j-1)
                     enddo
                  enddo
                  do i=2,k
                     sum=zero
                     do j=1,k   
                        sum=sum+phi(j)*xp(i-1,j)
                     enddo
                     xm(1,i) = sum  
                     xm(i,1)=xm(1,i) 
                  enddo
                  sum=sigmau**2
                  do i=1,k
                     do j=1,k
                        sum=sum+phi(i)*phi(j)*xp(i,j)
                     enddo
                  enddo
                  xm(1,1)=sum
               else
                  do i=1,lfin
                     do j=1,lfin
                        xm(i,j)=phi(i)*phi(j)*xp(1,1)
                        if(i.lt.lfin) xm(i,j)=xm(i,j)+phi(j)*xp(i+1,1)
                        if(j.lt.lfin) xm(i,j)=xm(i,j)+phi(i)*xp(1,j+1)
                        if(i.lt.lfin.and.j.lt.lfin) 
     1                       xm(i,j)=xm(i,j)+xp(i+1,j+1)
                     enddo
                  enddo
                  do i=1,lfin
                     do j=1,lfin
                        xm(i,j)=xm(i,j)+rv(i,j)*sigmau**2
                     enddo
                  enddo 
               endif   
*     
*     we build vector v (first column of matrix xm) and we store the first
*     column of xm in the corresponding column of auxm.
*
               do i=1,lfin
                  auxm(it,i)=xm(i,1)
                  v(i)=xm(i,1)
               enddo
*     
*     we compute st.
*
               if (v(1).le.zero) then
                  v(1)=sigmau
               endif
               st(it)=dsqrt(v(1))
*     
*     we compute res and w(it).
*
               hf=zero
               res=uhat(it)/st(it) 
               if ((it.lt.iout).or.(it.gt.iout+h)) then
                  hf=cck*s_psiffe(res/cck)
                  w(it)=1.0d0
                  if(res.ne.zero) w(it)=hf/res
               else
                  w(it)=1.0
               endif
*     
*     we recover from auxm the corresponding first column of xm 
*     and store it in v (if indfil=1).
*    
            else
               do i=1,lfin
                  v(i)=auxm(it,i)
               enddo
            endif 
*     
*     we compute alfafi.
*
            do i=1,lfin
               alfafi(i)=alfapred(i)+w(it)*v(i)*uhat(it)/(st(it)**2)
            enddo
*
*     we compute yhat(it) and epshat(it) and we update matrix xm if indfil=0. 
* 
            if (indfil.eq.0) then
               yhat(it)=alfafi(1)+y(it)-ypure(it) 
               epshat(it)=alfafi(1)
               do i=1,lfin
                  do j=1,lfin
                     av(i,j)=v(i)*v(j)
                  enddo
               enddo
               do i=1,lfin
                  do j=1,lfin
                     xp(i,j)=xm(i,j)-av(i,j)*w(it)/st(it)**2
                  enddo
               enddo
*
*     control of sequence of outliers and possible level shifts.
*      
               if(dabs(res).gt.control) then
                  numout=numout+1                           
                  if (numout.eq.1) then
                     jout=it 
                     do i1=1,lfin
                        do i2=1,lfin
                           xaux(i1,i2)=xm(i1,i2)
                        end do
                        alfaux(i1)=alfapred(i1)
                     end do
                  endif
                  if (numout.gt.2) then
                     numout=0
                     ilscan=ilscan+1
                     lscan(ilscan)=jout
                     yhat(jout)=y(jout)
                     w(jout)=1.0d0
                     it=jout
                     do i1=1,lfin
                        do i2=1,lfin
                           xm(i1,i2)=xaux(i1,i2)
                        enddo
                     enddo 
                     do i=1,lfin
                        v(i)=xm(1,i)
                        alfafi(i)=alfaux(i)+v(i)*uhat(it)/st(it)**2
                     enddo
                     if (v(1).lt.zero) then
                        v(1)=sigmau
                     endif
                     st(it)=dsqrt(v(1))
                     epshat(it)=alfafi(1)
                     do i=1,lfin
                        do j=1,lfin
                           av(i,j)=v(i)*v(j)
                        enddo
                     enddo
                     do i=1,lfin
                        do j=1,lfin
                           xp(i,j)=xm(i,j)-av(i,j)/st(it)**2
                        enddo
                     enddo
                  endif   
               else  
                  numout=0
                  jout=0
               end if
            endif
*
*     we return uhat(it) in column mm of matrix xy 
*     

            if (indfil.eq.1) xy(it-n0,mm)=uhat(it)
*
*     end of the cycle of the observations.      
*
            it=it+1
         enddo
*     
*     end of the cycle over y and x, in case indfil=1.
*
      enddo
      return                                    
      end     
C=======================================================================
      subroutine s_fltrfe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,k, 
     +     iq,sigmau,indth,n0,tau,sigmadif,indfil,rho,cck, 
     +     npred,ypure,xy,yhat,uhat,epshat,st,epspred,w, 
     +     auxm,ndim2,xm,xp,thetanew,alfapred,r,v,av, 
     +     alfafi,rv,alfaux,xaux,v1,v2,work,nw,iwork,niw) 
*-----------------------------------------------------------------------
*     This subroutine produces a filtered version of the input series  
*     y(t). It also produces a series of filtered residuals and a 
*     series of standard deviations of the prediction errors. 
* 
*     Input:   
*             x        : matrix of independent variables. 
*             y        : input series  
*             n        : number of observations 
*             m        : number of independent variables 
*             idif     : number of ordinary differences 
*             isp      : seasonal period 
*             nsd      : number of seasonal differences 
*             phi      : vector of coefficients of the AR models 
*             beta     : vector of coeff. of the independent variables         
*             theta    : vector of ordinary MA coefficients 
*             thetas   : seasonal moving average coefficient 
*             k        : order of the non stationary AR model to be used  
*             iq       : order of the ordinary MA model  
*             sigmau   : initial scale 
*             indth    : 0 - no seasonal moving average component 
*                        1 - seasonal moving average term included 
*             n0       : length of the vector of initial values 
*             tau      : vector containing E(u(t-j+1) z(t))/E(z(t)^2)     
*             sigmadif : scale of the differenced series 
*             indfil   : filtering selection 
*                        1 - filtering of y and x with fixed w and st 
*                        0 - filtering of epshat 
*             rho      : autocorrelation function of the differenced series 
*             cck      : bandwidth of the robust filter 
*             npred    : number of predicted regression errors 
*             ndim2    : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), 
*                        required to dimension several arrays. 
*              
*     Output: 
*             ypure  : vector containing the series y-x'*beta 
*             xy     : if indfil=1, this matrix contains the values of 
*                      the vector uhat corresponding to the filtering of 
*                      the x's and y. 
*             yhat   : filtered series  
*             uhat   : series of filtered residuals  
*             epshat : series of filtered errors 
*             st     : series of standard deviations of the  prediction 
*                      errors 
*             epspred: series of predicted filtered errors 
*             w      : series of weights 
*             auxm   : matrix containing the first columns of matrix xm  
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),xy(n,m+1),phi(ndim2+1),beta(m),theta(iq) 
      dimension rho(0:ndim2),tau(0:ndim2),ypure(n),v1(iq+1),v2(isp+1)  
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred),w(n+npred),auxm(n+npred,ndim2) 
      dimension xm(ndim2,ndim2),xp(ndim2,ndim2),thetanew(ndim2), 
     +     alfapred(ndim2),r(ndim2),v(ndim2),av(ndim2,ndim2),
     +     alfafi(ndim2),rv(ndim2,ndim2),alfaux(ndim2),xaux(ndim2,ndim2)
      dimension work(nw) 
      integer iwork(niw)
      data zero, one/0.d0, 1.d0/
*-----------------------------------------------------------------------

c initialize variables defined within if statements
      res = 0.d0

      control=2.5d0 
*     
*     We initialize r, rv, v1, v2. 
* 
      do i=1,ndim2 
         r(i)=zero 
         do j=1,ndim2 
            rv(i,j)=zero 
         enddo 
      enddo
      do i=1,iq+1 
         v1(i)=zero
      enddo 
      do i=1,isp+1 
         v2(i)=zero 
      enddo 
* 
*    In v1 we allocate the regular MA polynomial,in v2 the seasonal MA 
*    polynomial and in r the product MA polynomial. 
* 
      v1(1)=one
      do i=1,iq 
         v1(i+1)=-theta(i) 
      enddo 
      v2(1)=one 
      if (indth.eq.1) v2(isp+1)=-thetas 
      nv2=indth*isp  
      call s_polyfe(v1,iq,v2,nv2,r,nrdim)   
* 
*     We construct vector thetanew of moving average coefficients 
*     corresponding to r, to be used by s_rinife. 
* 
      do i=1,nrdim 
         thetanew(i)=-r(i+1) 
      enddo 
* 
*     We redefine the order of the regressive part to lfin. This is 
*     required for the space state representation, where the 
*     order of the autoregressive part is at least the order of the 
*     moving average part +1. 
*       
      lfin=max0(k,nrdim+1)
      do i=k+1,lfin 
         phi(i)=zero
      enddo 
* 
*     We compute the matrix r*r'. This appears also in the state space  
*     representation. 
* 
      do i=1,lfin 
         do ii=1,lfin 
            rv(i,ii)=r(i)*r(ii) 
         enddo 
      enddo    
*     
*     We call subroutine s_rinife which returns initial values of vector  
*     alfafi and matrix xp. Previously, we must construct the series  
*     ypure= y-x'*beta and the vector of autocorrelations rho.  
*  
*     If indfil=1 it begins the loop to filter all the x' and y. 
*     If indfil=0 it filters y-beta'x. 
* 
      nloop=1 
      if (indfil.eq.1) nloop=m+1 
      do mm=1,nloop 
         do i=1,n 
            if (indfil.eq.0) then       
               ypure(i)=y(i) 
               do j=1,m 
                  ypure(i)=ypure(i)-x(i,j)*beta(j) 
               enddo 
            else 
               if (mm.eq.nloop) then     
                  ypure(i)=y(i) 
               else 
                  ypure(i)=x(i,mm) 
               endif 
            endif 
         enddo 
*     
*     s_rinife computes initial values for the state space vector alfafi 
*     and the prediction covariace matrix xp, and n0: the number of first 
*     observations which can not be filtered. 
* 
         n2=ndim2 
         n3=n2+ndim2*ndim2 
         n4=n3+ndim2*ndim2 
         n5=n4+ndim2 
         n6=n5+ndim2 
         n7=n6+ndim2 
         n8=n7+ndim2*ndim2 
         call s_rinife(ypure,n,lfin,phi,thetanew,isp,rho,tau,alfafi,
     +        idif,nsd,xp,nrdim,n0,sigmau,sigmadif,ndim2,work(1), 
     +        work(n2+1),work(n3+1),work(n4+1),work(n5+1), 
     +        work(n6+1),work(n7+1),work(n8+1),iwork) 
*                           
*     numout  will count the number of outliers in the current batch. 
*     jout will be the index of the first outlier in the present batch. 
* 
         jout=0            
         numout=0  
* 
*     Initialization of yhat(i), w(i) and st(i) for i=1,n0. 
* 
         if (indfil.eq.0) then 
            do it=1,n0 
               yhat(it)=y(it) 
               epshat(it)=ypure(it) 
               uhat(it)=zero 
            enddo 
            w(n0)=one 
         else 
            if (mm.eq.nloop) then 
               do it=1,n0 
                  yhat(it)=y(it) 
                  epshat(it)=y(it) 
               enddo 
            else  
               do it=1,n0 
                  yhat(it)=x(it,mm) 
                  epshat(it)=x(it,mm) 
               enddo 
            endif 
         endif 
* 
*     Cycle of the observations. 
* 
         it=n0+1 
         do while (it.le.n+npred) 
*                     
*     Calculus of alfapred when the model is not pure autoregressive. 
* 
            if (indth.eq.1.or.iq.gt.0) then 
               do ir=1,lfin-1 
                  alfapred(ir)=phi(ir)*alfafi(1)+alfafi(ir+1) 
               enddo 
               alfapred(lfin)=phi(lfin)*alfafi(1) 
            else 
* 
*     Calculus of alfapred when the model is pure autoregressive.  
* 
               do ir=1,lfin-1 
                  alfapred(ir+1)=alfafi(ir) 
               enddo 
               sum2=zero 
               do ir=1,lfin 
                  sum2=sum2+phi(ir)*alfafi(ir) 
               enddo 
               alfapred(1)=sum2 
            endif 
* 
*     Calculus of the predicted value of the ARIMA series. 
* 
            epspred(it)=alfapred(1) 
* 
*     Calculus of uhat, the filtered innovation. 
* 
            if(it.le.n) then 
               uhat(it)=ypure(it)-alfapred(1) 
            else 
               uhat(it)=zero 
            endif 
* 
*     Updating of matrix xm when the model is pure autoregressive 
*     (if indfil=0). 
*       
            if (indfil.eq.0) then  
               if (indth.eq.0.and.iq.eq.0) then 
                  do i=2,k            
                     do j=2,k            
                        xm(i,j)=xp(i-1,j-1) 
                     enddo 
                  enddo 
                  do i=2,k 
                     sum=zero 
                     do j=1,k    
                        sum=sum+phi(j)*xp(i-1,j) 
                     enddo 
                     xm(1,i)=sum   
                     xm(i,1)=xm(1,i)  
                  enddo 
                  sum=sigmau**2 
                  do i=1,k 
                     do j=1,k 
                        sum=sum+phi(i)*phi(j)*xp(i,j) 
                     enddo 
                  enddo 
                  xm(1,1)=sum 
               else 
* 
*     Updating of matrix xm  when the model is not pure autoregressive 
*     (if indfil=0). 
* 
                  do i=1,lfin 
                     do j=1,lfin 
                        xm(i,j)=phi(i)*phi(j)*xp(1,1) 
                        if(i.lt.lfin) xm(i,j)=xm(i,j)+phi(j)*xp(i+1,1) 
                        if(j.lt.lfin) xm(i,j)=xm(i,j)+phi(i)*xp(1,j+1) 
                        if(i.lt.lfin.and.j.lt.lfin)  
     +                       xm(i,j)=xm(i,j)+xp(i+1,j+1) 
                     enddo 
                  enddo 
                  do i=1,lfin 
                     do j=1,lfin 
                        xm(i,j)=xm(i,j)+rv(i,j)*sigmau**2 
                     enddo 
                  enddo  
               endif    
*                
*     We build vector v (first column of matrix xm) and we store the first 
*     column of xm in the corresponding column of auxm. 
* 
               do i=1,lfin 
                  auxm(it,i)=xm(i,1) 
                  v(i)=xm(i,1) 
               enddo 
* 
*     We compute st: the prediction error scale of observation i. 
* 
               if (v(1).le.zero) then 
                  v(1)=sigmau 
               endif 
               st(it)=dsqrt(v(1)) 
* 
*     We compute w(it): the weight that wil be given at the  
*     observed process, res: the standardized innovation and hf: the  
*     standardized pseudo innovation. 
* 
               hf=zero 
               res=uhat(it)/st(it)  
               hf=cck*s_psiffe(res/cck) 
               w(it)=one 
               if (dabs(res).ge.0.0000000001d0) w(it)=hf/res 
* 
*     For predicted observations w=0. 
* 
               if (it.gt.n) w(it)=zero 
*  
*     We recover from auxm the corresponding first column of xm  
*     and store it in v (if indfil=1). 
*     
            else 
               do i=1,lfin 
                  v(i)=auxm(it,i) 
               enddo 
            endif  
* 
*     We compute alfafi. 
* 
            do i=1,lfin 
               alfafi(i)=alfapred(i)+w(it)*v(i)*uhat(it)/(st(it)**2) 
            enddo 
* 
*     We compute the cleaned series yhat(it), the filtered 
*     regression error epshat(it) and update matrix xm (if indfil=0).  
*  
            if (indfil.eq.0) then 
               if(it.le.n) yhat(it)=alfafi(1)+y(it)-ypure(it)  
               epshat(it)=alfafi(1) 
               do i=1,lfin 
                  do ii=1,lfin 
                     av(i,ii)=v(i)*v(ii) 
                  enddo 
               enddo    
               do i=1,lfin 
                  do j=1,lfin 
                     xp(i,j)=xm(i,j)-av(i,j)*w(it)/st(it)**2 
                  enddo 
               enddo 
* 
*     Control of sequence of outliers and possible level shifts. 
* 
*     We increase the number of outliers of the current batch.  
* 
               if (dabs(res).gt.control) then 
                  numout=numout+1                         
* 
*     If we start a new batch of outliers, we store the matrix xm and  
*     alfapred. 
*    
* 
                  if (numout.eq.1) then 
                     jout=it  
                     do i1=1,lfin 
                        do i2=1,lfin 
                           xaux(i1,i2)=xm(i1,i2) 
                        enddo 
                        alfaux(i1)=alfapred(i1) 
                     enddo 
                  endif 
* 
*     If the number of outliers in the current batch is larger than 2 we 
*     go back to the first outlier of the batch and redefine it as no outlier.
* 
                  if (numout.gt.2) then 
                     numout=0 
                     yhat(jout)=y(jout) 
                     w(jout)=one 
                     it=jout 
                     do i1=1,lfin 
                        do i2=1,lfin 
                           xm(i1,i2)=xaux(i1,i2) 
                        enddo 
                     enddo  
*     
*     We redefine alfafi using alfapred stored in alfaux, now using w(it)=1 
* 
                     do i=1,lfin 
                        v(i)=xm(1,i) 
                        alfafi(i)=alfaux(i)+v(i)*uhat(it)/st(it)**2 
                     enddo 
                     if (v(1).lt.zero) v(1)=sigmau 
                     st(it)=dsqrt(v(1)) 
                     epshat(it)=alfafi(1) 
                     do i=1,lfin 
                        do ii=1,lfin 
                           av(i,ii)=v(i)*v(ii) 
                        enddo 
                     enddo 
                     do i=1,lfin 
                        do j=1,lfin 
                           xp(i,j)=xm(i,j)-av(i,j)/st(it)**2 
                        enddo 
                     enddo 
                  endif    
               else   
                  numout=0 
                  jout=0 
               end if 
            endif 
* 
*     If indfil=1 we return uhat(it) in column mm of matrix xy.  
* 
            if (indfil.eq.1) xy(it-n0,mm)=uhat(it) 
* 
*     End of the cycle of the observations.       
* 
            it=it+1 
         enddo 
* 
*     End of the cycle over y and x, in case indfil=1. 
* 
      enddo 
      return                                     
      end 
C=======================================================================
      subroutine s_fnc1fe(n,npar,par,f,iflag,idif,isp,nsd,m,np,nq,n0, 
     +     indth,npo,sigman,sigmau,npred,x,y,xy,yhat,uhat,epshat,st,
     +     epspred,w,auxm,poldif,ndim1,ndim2,phi,theta,phiaux,phiaux2,
     +     beta,uaux,rho,tau,para,para1,theprod,ypure,aux,cck,work4, 
     +     nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
*-----------------------------------------------------------------------
*     This subroutine computes a vector f(1),..., f(n) 
*     for a given set of transformed parameters par(1),...,par(npar). 
*     The minimization of f(1)^2+...+f(n)^2 is equivalent to 
*     the minimization of the tau-pseudo-likelihood. 
* 
*     The order of the parameters within the vector par is: 
* 
*        1.  par(1),...,par(np): transformed AR parameters 
*        2.  par(np+1),...,par(np+nq): transformed regular MA parameters 
*        3.  par(np+nq+1): transformed seasonal MA parameter 
*        4.  par(np+nq+indth+1),...,par(np+nq+indth+m): regression coeff
*        5   par(np+nq+indth+m+1) the filter bandwidth 
* 
*     All the regular AR and MA parameters as well as the seasonal MA 
*     parameter are transformed so that they lie in the interval (-1,1). 
* 
*     Input: 
*           n     : number of observations 
*           npar  : number of transformed parameters 
*           par   : transformed parameters 
*           idif  : number of regular differences 
*           isp   : seasonal period 
*           nsd   : number of seasonal differences 
*           m     : number of regression independent variables 
*           np    : order of the autoregressive polinomial 
*           nq    : order of the moving average polinomial 
*           n0    : number of initial non filtered observations 
*           indth : 1 indicates that the model has a seasonal MA parameter 
*           npo   : vector containig the orders of the AR and MA polynomials. 
*                   In the present version npo=(1,...,np,1,...,nq). 
*           x     : matrix of independent regression variables 
*           y     : response variable 
*           ndim1 : ip+iqfin+m+indth+1, required to dimension several arrays. 
*           ndim2 : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to  
*                   dimension several arrays 
*           poldif: vector containing the coefficients of the differences  
*                   polynomial 
* 
*     Output          
*           f      : vector whose sum of squares gives the value of the 
*                    tau-pseudo-likelihood 
*           sigman : innovations scale 
*           sigmau : corrected innovations scale 
*           yhat   : cleaned y 
*           uhat   : innovation residuals 
*           epshat : regression error residuals 
*           st     : scales of epspred 
*           epspred: predicted regression errors 
*           w      : series of weights, required by filter 
*           auxm   : matrix required by filter 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      dimension par(ndim1),f(n),x(n,m),y(n),xy(n,m+1),ypure(n),aux(4*n)
      dimension yhat(n),uhat(n+npred),epshat(n+npred),w(n+npred), 
     +     st(n+npred),epspred(n+npred),auxm(n+npred,ndim2)
      dimension poldif(isp*nsd+idif+1),theprod(ndim2),uaux(n)
      dimension phi(ndim2),theta(nq),phiaux(np+1),phiaux2(ndim2+1), 
     +     beta(m),rho(0:ndim2),tau(0:ndim2),para(ndim2),para1(ndim2) 
      dimension work4(nw4),work5(nw5) 
      integer iwork4(niw4),iwork5(niw5),npo(np+nq) 
      data zero,one/0.d0,1.d0/
*-----------------------------------------------------------------------
*     
*     We antitransform the npar parameters from par, using subroutine 
*     s_tranfe. 
* 
      iflag = iflag
      npar = npar
      ndif=isp*nsd+idif 
      do i=1,ndim2 
         theprod(i)=zero 
      enddo 
      do i=1,nq 
         theta(i)=zero
      enddo 
      do i=1,ndim2 
         phi(i)=zero
      enddo  
      call s_tranfe(par,ndim1,ndim2,np,nq,indth,m,para,para1,work4, 
     +     phi,theta,thetas,beta)       
* 
*     We compute the  scale  (sigini) of the stationary ARMA component 
*     of the regression model, that is, of the regression errors after  
*     differencing.  
* 
      sigini=s_xmadfe(x,y,beta,m,n,aux(1),aux(n+1),aux(2*n+1), 
     +     poldif,ndif) 
* 
*     npaux and nqaux are the maximum of the powers with non zero 
*     AR and MA coefficients respectively. In this version npaux=np and  
*     nqaux=nq. 
* 
      if (np.eq.0) then  
         npaux=0 
      else 
         npaux=npo(np) 
      endif 
      if (nq.eq.0) then  
         nqaux=0 
      else 
         nqaux=npo(np+nq) 
      endif 
* 
*     lfin gives the smallest order of the stationary autoregressive model 
*     compatible with the state space representation of the ARIMA model. 
* 
      lfin=max0(npaux+idif+isp*nsd,nqaux+isp*indth+1) 
      lfin=lfin-idif-isp*nsd 
* 
*     We call subroutine s_sys2fe which returns the autocorrelations, rho,  
*     and the correlation between the lagged innovations and the ARMA 
*     series, tau.       
* 
      n2=ndim2 
      n3=n2+(ndim2+1)*(ndim2+1) 
      n4=n3+(ndim2+1) 
      call s_sys2fe(phi,theta,thetas,lfin,nqaux,isp,indth,rho,tau, 
     +     work4(1),work4(n2+1),work4(n3+1),work4(n4+1), 
     +     iwork4,ndim2) 
* 
*     theprod is the vector containing the MA parameters corresponding to  
*     the product of the two MA polynomial operators. 
* 
      nqq=nqaux+indth*isp 
      do i=1,nqaux 
         theprod(i)=theta(i) 
      enddo 
      if (indth.ne.0) then 
         theprod(isp)=thetas 
         do i=1,nqaux 
            theprod(isp+i)=-theta(i)*thetas 
         enddo 
      endif 
* 
*     aa is the ratio between the innovation scale and the scale of the 
*     stationary ARMA process. 
* 
      aa=zero 
      do i=1,lfin 
         do j=1,lfin 
            aa=aa+phi(i)*phi(j)*rho(abs(i-j)) 
         enddo 
      enddo 
      bb=zero
      do i=1,lfin 
         do j=i,nqq 
            bb=bb+phi(i)*theprod(j)*tau(abs(j-i)) 
         enddo 
      enddo 
      cc=one
      do i=1,nqq 
         cc=cc+theprod(i)**2 
      enddo 
      dd=(one-aa+bb)/cc 
* 
*     We compute the innovation scale sigmau. 
* 
      sigmau=sigini*dsqrt(dd) 
      do i=1,np+1 
         phiaux(i)=zero 
      enddo  
      phiaux(1)=one 
      do i=1,np 
         phiaux(npo(i)+1)=-phi(npo(i)) 
      enddo  
* 
*     We compute the non stationary autoregressive polynomial. 
* 
      call s_polyfe(phiaux,npaux,poldif,ndif,phiaux2,k)       
      do i=1,k 
         phiaux2(i)=-phiaux2(i+1) 
      enddo 
      do i=k+1,ndim2 
         phiaux2(i)=zero 
      enddo 
* 
*     We filter the original series y using the given parameters and obtain 
*     the filtered innovations uhat and their scales st. 
*     
      call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux2,beta,theta,thetas,k, 
     +     nqaux,sigmau,indth,n0,tau,sigini,0,rho,cck,0,ypure, 
     +     xy,yhat,uhat,epshat,st,epspred,w,auxm,ndim2,work4, 
     +     nw4,work5,nw5,iwork5,niw5) 
* 
*     We start computing  pseudo residuals f such that the sum of the 
*     squares is the tau-pseudo likelihood  of the filtered residuals, 
*     except by a monotone transformation. We begin by initializing f. 
* 
      do i=1,n0 
         f(i)=zero 
      enddo 
* 
*     We compute the standardized filtered innovations. 
* 
      do i=n0+1,n 
         uaux(i)=uhat(i)/st(i) 
      enddo         
* 
*     We compute the initial M-scale of the filtered innovations 
*     used by the tau scale. 
* 
      call s_calsfe(uaux,n,n0,sout,aux(1),aux(n+1)) 
* 
*     Sigman is the corrected innovation scale. 
* 
      sigman=sigmau*sout    
      prod=zero
      do i=n0+1,n 
         prod= prod + 2.0d0 * dlog(st(i)) 
      enddo 
      prod=prod/dble(n-n0) 
      prod=dexp(prod) 
      prod=prod*sout*sout 
      do i=n0+1,n 
         f(i)=uaux(i)/sout 
         f(i)=s_rhoffe(f(i)) 
         f(i)=dsqrt(prod*f(i)) 
      enddo 
      return 
      end 
C=======================================================================
      subroutine s_fnc2fe(phi,theta,thetas,n,beta,cck,idif,isp,nsd,m,np, 
     +     nq,n0,indth,x,y,sigman,sigmau,vtau,sigini,tau,xy,yhat,uhat,
     +     epshat,st,epspred,w,auxm,npred,poldif,ndim2,f,phiaux,
     +     phiaux2,uaux,rho,aux,ypure,work4,nw4,iwork4,niw4,work5,nw5,
     +     iwork5,niw5) 
*-----------------------------------------------------------------------
*     This subroutine computes the tau pseudo likelihood for a given  
*     set of parameters. The difference with s_fc11fe is that in the  
*     latter the AR and MA parameters are transformed, so that 
*     they are restricted to lie in the interval (-1,1). 
* 
*     Input: 
*          phi   : vector of ordinary AR coefficients 
*          theta : vector of ordinary MA coefficients 
*          thetas: seasonal MA coefficient 
*          n     : length of the series 
*          beta  : vector of regression coefficient 
*          cck   : bandwidth of the robust filter 
*          idif  : number of ordinary differences 
*          isp   : seasonal period 
*          nsd   : number of seasonal differences 
*          m     : length of beta 
*          np    : length of phi 
*          nq    : length of theta 
*          n0    : length of the vector of initial values 
*          indth : 0 - no seasonal moving average component 
*                  1 - seasonal moving average term included 
*          x     : matrix of independent variables 
*          y     : observed series 
*          sigmau: scale of the filtered residuals 
*          ndim2 : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                  dimension several arrays 
* 
*     Output: 
*          sigman: corrected innovation scale 
*          vtau  : value of the tau-pseudo-likelihood 
*          sigini: scale of the differenced regression errors 
*          tau   : vector containing E(u(t-j+1) z(t))/E(z(t)^2)     
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      dimension phi(ndim2),theta(nq),beta(m),x(n,m),y(n),tau(0:ndim2) 
      dimension xy(n,m+1),yhat(n),uhat(n+npred),epshat(n+npred), 
     +     st(n+npred),epspred(n+npred),w(n+npred),f(n),uaux(n), 
     +     auxm(n+npred,ndim2),poldif(isp*nsd+idif+1),aux(4*n)
      dimension phiaux(ndim2+1),phiaux2(ndim2+1),rho(0:ndim2),ypure(n) 
      dimension work4(nw4),work5(nw5) 
      integer iwork4(niw4),iwork5(niw5)
      data zero/0.d0/
*-----------------------------------------------------------------------
* 
*     We clean the last elements of vectors theta and phi. 
* 
      ndif=isp*nsd+idif
*     do i=nq+1 
      do i=1,nq 
         theta(i)=zero
      enddo 
      do i=np+1,ndim2 
         phi(i)=zero
      enddo  
* 
*     We compute the scale (sigini) of the stationary ARMA component 
*     of the regression model, that is, of the regression errors after  
*     differencing.  
* 
      sigini=s_xmadfe(x,y,beta,m,n,aux(1),aux(n+1),aux(2*n+1),
     +     poldif,ndif) 
* 
*     npaux and nqaux are the maximum of the powers with non zero AR and MA  
*     coefficients respectively. In this version npaux=np and nqaux=nq. 
* 
      npaux=np 
      nqaux=nq 
* 
*     lfin gives the smallest order of the stationary autoregressive model 
*     compatible with the state space representation of the ARIMA model. 
* 
      lfin=max0(npaux+idif+isp*nsd,nqaux+isp*indth+1) 
      lfin=lfin-idif-isp*nsd 
* 
*     We call subroutine s_sys2fe which returns the autocorrelations, rho,  
*     and the correlation between the lagged innovations and the ARMA 
*     series, tau.       
* 
      n1=ndim2 
      n2=n1+(ndim2+1)*(ndim2+1) 
      n3=n2+(ndim2+1) 
      call s_sys2fe(phi,theta,thetas,lfin,nqaux,isp,indth,rho,tau, 
     +     work4(1),work4(n1+1),work4(n2+1),work4(n3+1),iwork4,ndim2) 
      do i=1,ndim2 
         phiaux(i)=zero 
      enddo  
      phiaux(1)=1.0d0 
      do i=1,np 
         phiaux(i+1)=-phi(i) 
      enddo  
* 
*     We compute the non-stationary AR polynomial to be used by filter. 
* 
      call s_polyfe(phiaux,npaux,poldif,ndif,phiaux2,k)
      do i=1,k 
         phiaux2(i)=-phiaux2(i+1) 
      enddo 
      do i=k+1,ndim2 
         phiaux2(i)=zero 
      enddo 
* 
*     We filter the original series y using the given parameters and obtain 
*     the filtered innovations uhat and their scales st. 
* 
      call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux2,beta,theta,thetas,
     +     k,nqaux,sigmau,indth,n0,tau,sigini,0,rho,cck,0,ypure, 
     +     xy,yhat,uhat,epshat,st,epspred,w,auxm,ndim2,work4, 
     +     nw4,work5,nw5,iwork5,niw5) 
      do i=1,n0 
         f(i)=zero 
      enddo 
* 
*     We compute the scale of the standarized filtered residuals. 
* 
      do i=n0+1,n 
         uaux(i)=uhat(i)/st(i) 
      enddo
      call s_calsfe(uaux,n,n0,sout,aux(1),aux(n+1)) 
* 
*     Sigman is the corrected innovation scale. 
* 
      sigman=sigmau*sout    
* 
*     We compute the tau pseudo likelihood criterion vtau. 
* 
      prod=zero
      do i=n0+1,n 
         prod= prod + 2.0d0 * dlog(st(i)) 
      enddo 
      do i=n0+1,n 
         f(i)=uaux(i)/sout 
         f(i)=s_rhoffe(f(i)) 
      enddo 
      vtau=zero 
      do i=n0+1,n 
         vtau=vtau+f(i) 
      enddo 
      vtau=vtau/dble(n-n0) 
      vtau=(sout*sout)*vtau/(.488d0) 
      vtau=dble(n-n0)*dlog(vtau)+prod
      return 
      end
C=======================================================================
      subroutine s_gd11fe(x,y,interc,n,m,idif,isp,nsd,iqfin,phi,beta, 
     +     phidif,sigini,sigmau,akai,cck,parold,ssqold, 
     +     xy,yhat,uhat,epshat,st,epspred,npred,poldif, 
     +     w,auxm,ndim1,ndim2,work2,nw2,iwork2,niw2, 
     +     work3,nw3,iwork3,niw3,work4,nw4,iwork4,niw4, 
     +     work5,nw5,iwork5,niw5,utol,maxfev,epsmch,dwarf) 
*-----------------------------------------------------------------------   
*     This subroutine distributes the work vectors work2 and iwork2, 
*     of lengths nw2 and niw2 respectively, between different arrays  
*     and vectors to be used by s_grd1fe, and then calls it. 
*-----------------------------------------------------------------------    
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phi(ndim2,ndim2),beta(m),parold(ndim1)
      dimension xy(n,m+1),yhat(n),uhat(n+npred),epshat(n+npred)
      dimension st(n+npred),epspred(n+npred),w(n+npred), 
     +     auxm(n+npred,ndim2),poldif(isp*nsd+idif+1)  
      dimension work2(nw2),work3(nw3),work4(nw4),work5(nw5) 
      dimension iwork2(niw2),iwork3(niw3),iwork4(niw4),iwork5(niw5) 
*-----------------------------------------------------------------------   
      n1=1 
      n2=m 
      n3=n2+m 
      n4=n3+idif+isp*nsd+1 
      n5=n4+ndim2 
      n6=n5+ndim2 
      n7=n6+ndim2+1 
      n8=n7+ndim2+1 
      n9=n8+ndim2+1 
      n10=n9+n 
      n11=n10+iqfin 
      n12=n11+n*ndim1 
      n13=n12+ndim1 
      n14=n13+n 
      n15=n14+ndim1 
      n16=n15+ndim1 
      n17=n16+ndim1 
      n18=n17+n 
      n19=n18+ndim1 
      n20=n19+ndim1 
      n21=n20+4*n 
      call s_grd1fe(x,y,interc,n,m,idif,isp,nsd,iqfin,phi,beta,phidif, 
     +     sigini,sigmau,akai,cck,parold,ssqold,xy,yhat,uhat, 
     +     epshat,st,epspred,w,auxm,npred,poldif,ndim1,ndim2, 
     +     work2(n1),work2(n2+1),work2(n3+1),work2(n4+1), 
     +     work2(n5+1),work2(n6+1),work2(n7+1),work2(n8+1), 
     +     work2(n9+1),work2(n10+1),work2(n11+1),work2(n12+1), 
     +     work2(n13+1),work2(n14+1),work2(n15+1),work2(n16+1), 
     +     work2(n17+1),work2(n18+1),work2(n19+1),work2(n20+1), 
     +     work2(n21+1),iwork2(1),iwork2(ndim1+1),work3,nw3, 
     +     iwork3,niw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5, 
     +     niw5,utol,maxfev,epsmch,dwarf) 
      return 
      end 
C=======================================================================
      subroutine s_gdk1fe(x,y,n,m,ip,idif,isp,nsd,iqfin,phi,beta, 
     +     thetas,phidif,k,sigini,sigmau,akai,cck,parold,ssqold, 
     +     xy,yhat,uhat,epshat,st,epspred,w,auxm,npred,poldif,ndim1,
     +     ndim2,work2,nw2,iwork2,niw2,work3,nw3,iwork3,niw3,work4,nw4,
     +     iwork4,niw4,work5,nw5,iwork5,niw5,utol,maxfev,epsmch,dwarf) 
*-----------------------------------------------------------------------   
*     This subroutine distributes the work vectors work2 and iwork2, 
*     of lengths nw2 and niw2 respectively, between different arrays  
*     and vectors to be used by s_grdkfe, and then calls it. 
*-----------------------------------------------------------------------   
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phi(ndim2,ndim2),beta(m),phidif(ip+1,ip+1), 
     +     parold(ndim1),xy(n,m+1),yhat(n),uhat(n+npred), 
     +     epshat(n+npred),st(n+npred),epspred(n+npred),w(n+npred), 
     +     auxm(n+npred,ndim2),poldif(idif+isp*nsd+1) 
      dimension work2(nw2),work3(nw3),work4(nw4),work5(nw5) 
      integer iwork2(niw2),iwork3(niw3),iwork4(niw4),iwork5(niw5) 
*-----------------------------------------------------------------------   
      n1=1 
      n2=iqfin 
      n3=n2+ndim2 
      n4=n3+n 
      n5=n4+n 
      n6=n5+n 
      n7=n6+ndim2+1 
      n8=n7+ip 
      n9=n8+ip 
      n10=n9+ip 
      n11=n10+ndim1 
      n12=n11+ndim2+1 
      n13=n12+ndim1*n 
      n14=n13+ndim1 
      n15=n14+n 
      n16=n15+ndim1 
      n17=n16+ndim1 
      n18=n17+ndim1 
      n19=n18+n 
      n20=n19+ndim1 
      n21=n20+ip+1 
      n22=n21+ndim2+1
      call s_grdkfe(x,y,n,m,ip,idif,isp,nsd,iqfin,phi,beta,thetas,
     +     phidif,k,sigini,sigmau,akai,cck,parold,ssqold,xy,yhat,uhat, 
     +     epshat,st,epspred,w,auxm,npred,poldif,ndim1,ndim2, 
     +     work2(n1),work2(n2+1),work2(n3+1),work2(n4+1), 
     +     work2(n5+1),work2(n6+1),work2(n7+1),work2(n8+1), 
     +     work2(n9+1),work2(n10+1),work2(n11+1),work2(n12+1), 
     +     work2(n13+1),work2(n14+1),work2(n15+1),work2(n16+1), 
     +     work2(n17+1),work2(n18+1),work2(n19+1),work2(n20+1), 
     +     work2(n21+1),work2(n22+1),iwork2(1),iwork2(ip+iqfin+1),
     +     work3,nw3,iwork3,niw3,work4,nw4,iwork4,niw4,work5,
     +     nw5,iwork5,niw5,utol,maxfev,epsmch,dwarf) 
      return 
      end 
C=======================================================================
      subroutine s_gdt1fe(x,y,n,m,ip,idif,isp,nsd,ipfin,iqfin,beta, 
     +     thetas,phidif1,kopt,sigini,sopt,tau,n0,cckopt,xy,yhat,uhat,
     +     epshat,st,epspred,w,auxm,npred,poldif,thetaopt,ndim2,work2,
     +     nw2,work3,nw3,iwork3,niw3,work4,nw4,iwork4,niw4, 
     +     work5,nw5,iwork5,niw5) 
*-----------------------------------------------------------------------    
*     This subroutine distributes the work vector work2 of lengths nw2  
*     between different arrays and vectors to be used by s_gdthfe, and then  
*     calls it. 
*-----------------------------------------------------------------------     
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),beta(m),phidif1(ip),tau(0:ndim2),xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred),w(n+npred),auxm(n+npred,ndim2), 
     +     poldif(idif+isp*nsd+1),thetaopt(iqfin+1)  
      dimension work2(nw2),work3(nw3),work4(nw4),work5(nw5) 
      integer iwork3(niw3),iwork4(niw4),iwork5(niw5)
*-----------------------------------------------------------------------     
      n1=1 
      n2=ndim2 
      n3=n2+ip+1 
      n4=n3+ip 
      n5=n4+ndim2 
      n6=n5+ndim2 
      n7=n6+ndim2+1 
      call s_gdthfe(x,y,n,m,ip,idif,isp,nsd,ipfin,iqfin,beta,thetas, 
     +     phidif1,kopt,sigini,sopt,tau,n0,cckopt,xy,yhat,uhat, 
     +     epshat,st,epspred,w,auxm,npred,poldif,thetaopt,ndim2, 
     +     work2(n1),work2(n2+1),work2(n3+1),work2(n4+1), 
     +     work2(n5+1),work2(n6+1),work2(n7+1),work3,nw3,iwork3, 
     +     niw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
      return 
      end 
C=======================================================================
      subroutine s_gdthfe(x,y,n,m,ip,idif,isp,nsd,ipfin,iqfin,beta,
     +     thetas,phidif1,k,sigini,sigmau,tau,n0,cck,xy,yhat, 
     +     uhat,epshat,st,epspred,w,auxm,npred,poldif, 
     +     theta,ndim2,thetabes,phiau1,phidbes,phithe, 
     +     phiaux,rho1,rho,work3,nw3,iwork3,niw3, 
     +     work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5)    
*-----------------------------------------------------------------------   
*     This subroutine computes the initial estimates for the ARMA  
*     coefficients  when there is a seasonal MA parameter 
* 
*     Input: 
*           x      : regression matrix 
*           y      : original series 
*           n      : length of the series 
*           m      : number of regression variables 
*           idif   : number of ordinary differences 
*           isp    : seasonal period 
*           nsd    : number of seasonal differences 
*           iqfin  : order of the MA polynomial 
*           beta   : on input contains the estimated regression coefficients  
*                    of the precedent model 
*           thetas : on input contains the previous estimated seasonal MA  
*                    coefficient 
*           phidif1: on input contains the estimated coefficients of the  
*                    precedent stationary model 
*           k      : order of the precedent stationary AR model 
*           sigini : scale of the differenced regression residuals, that is,  
*                    of the regression errors after differencing 
*           sigmau : estimated scale of the precedent estimated AR model 
*           tau    : vector containing e(u(t-j+1) z(t))/e(z(t)^2)     
*           n0     : number of initial non filtered observasions 
*           cck:   : bandwidth of the robust filter 
*           ndim2  : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                    dimension several arrays 
*           poldif : vector containing the coefficients of the differences 
*                    polynomial            
* 
*     Output: 
*           phi    : on output contains the estimated AR coefficients of the  
*                    present model 
*           beta   : on output contains the estimated regression coefficients  
*                    of the present model 
*           thetas : on output contains the new estimated seasonal MA  
*                    coefficient 
*           phidif1: on output contains the estimated coefficients of the  
*                    present stationary model 
*           sigmau : on output contains the estimated scale of the new 
*                    model
*-----------------------------------------------------------------------   
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),beta(m),phidif1(ip),tau(0:ndim2),xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred),w(n+npred),auxm(n+npred,ndim2), 
     +     poldif(idif+isp*nsd+1),thetabes(ndim2),xgrid(21),delta(3),
     +     phiau1(ip+1),phidbes(ip),phithe(ndim2),the(2),phiaux(ndim2), 
     +     theta(iqfin),xxxgrid(21),dddelta(3),rho1(0:ndim2), 
     +     rho(0:ndim2),work3(nw3),work4(nw4),work5(nw5) 
      integer iwork3(niw3),iwork4(niw4),iwork5(niw5),iiwork(300) 
      data xxxgrid / -0.6d0,0.4d0,0.2d0,0.8d0,0.0d0,0.1d0,-.3d0,0.3d0, 
     +     -0.5d0,0.9d0,-0.2d0,-0.7d0,0.6d0,0.98d0,-0.4d0, 
     +     -0.98d0,0.5d0,0.7d0,-0.1d0,-0.9d0,-0.8d0/ 
      data dddelta /.05d0,.025d0,.0125d0 /
      data zero/0.d0/

c initialize variables defined within if statements
      vtaumin = 0.d0
      themej = 0.d0
      sigmin = 0.d0

*-----------------------------------------------------------------------   
* 
*     We  will minimize the tau-filtered pseudo likelihood. 
*     This is done in two stages. In the first stage we use a grid 
*     in the interval [-1,1]. In the second stage we use three steps  
*     of bisection. The program has two blocks corresponding to these 
*     two stages. 
* 
      do ih=1,21 
         xgrid(ih)=xxxgrid(ih) 
      enddo 
      do ih=1,3 
         delta(ih)=dddelta(ih) 
      enddo 
      do nn=1,ndim2 
         phiaux(nn)=zero
      enddo 
      do nn=1,k 
         phiaux(nn)=phidif1(nn)  
      enddo 
* 
*     Using the selected AR model, s_sys2fe computes the autocorrelations rho  
*     and the covariances tau between the stationary  process and the  
*     innovations. 
* 
      n1=ndim2 
      n2=n1+(ndim2+1)*(ndim2+1) 
      n3=n2+(ndim2+1) 
      call s_sys2fe(phiaux,theta,thetas,k,0,0,0,rho1,tau,work3(1), 
     +     work3(n1+1),work3(n2+1),work3(n3+1),iwork3,ndim2) 
* 
*     Beginning of block 1. 
* 
*     We start the optimization over a 21 points grid: xgrid. 
* 
      indth=1 
      do l=1,21 
         thetas=xgrid(l) 
* 
*     Variance of the innovation for the current value of thetas. 
* 
         thesig=sigmau 
*           
*       Given the rho's computed with s_sys2fe, s_sys1fe computes a 
*       model with the same rho's, which includes an AR operator  
*       of the same order than the selected one, and a seasonal MA  
*       parameter equal to the current value of thetas. The 
*       autoregressive parameters are in the array phithe. 
* 
         n1=ndim2 
         n2=n1+ip*ip 
         call s_sys1fe(phidif1,ip,thetas,rho1,k,isp,phithe,tau,ier, 
     +        work3(1),work3(n1+1),work3(n2+1),iwork3,ndim2) 
         if (iqfin.gt.0) then 
            n1=ndim2 
            n2=n1+(ndim2+1)*(ndim2+1) 
            n3=n2+(ndim2+1) 
* 
*     s_sys2fe computes the autocorrelations rho and the covariances  
*     tau between the stationary  process and the innovations of the 
*     model computed by s_sys1fe. 
* 
            call s_sys2fe(phithe,theta,thetas,k,0,isp,indth,rho,tau, 
     +           work3(1),work3(n1+1),work3(n2+1),work3(n3+1), 
     +           iwork3,ndim2) 
*           
*     Given the rho's and tau's computed by s_sys2fe, s_sys3fe computes  
*     the coefficients of an ARMA model with the same rho's, and tau's 
*     and with a seasonal MA parameter equal to the current value of  
*     thetas.  
* 
            call s_sys3fe(phithe,k,phiaux,theta,thetas,ipfin,iqfin,isp, 
     +           indth,rho,tau,ndim2,work3(1), 
     +           work3(ndim2*ndim2+1),iwork3)   
         else 
            do i=1,ndim2 
               phiaux(i)=zero 
            enddo 
            do i=1,k 
               phiaux(i)=phithe(i) 
            enddo 
         endif 
* 
*     We check if the AR operator computed by s_sys3fe is stationary. 
*     Otherwise is modified. 
*   
         if (ipfin.gt.0) then  
            call s_yulefe(phiaux,rho,ipfin,work3(1),iiwork,ndim2) 
            call s_durbfe(rho,ipfin,phiau1,ier,work3(1),ndim2) 
            if (ier.eq.1) then 
               do i=1,ipfin 
                  if (phiau1(i).gt. 1.d0) phiau1(i)= .98d0 
                  if (phiau1(i).lt.-1.d0) phiau1(i)=-.98d0 
               enddo 
               call s_invdfe(phiau1,ipfin,phiaux,work3(1),ndim2)   
            endif 
         endif 
* 
*     We check if the MA operator computed by s_sys3fe is invertible. 
*     Otherwise is modified. 
*       
         if (iqfin.gt.0) then 
            call s_yulefe(theta,rho,iqfin,work3(1),iiwork,ndim2) 
            call s_durbfe(rho,iqfin,phiau1,ier,work3(1),ndim2) 
            if (ier.eq.1) then 
               do i=1,iqfin 
                  if (phiau1(i).gt. 1.d0) phiau1(i)= .98d0 
                  if (phiau1(i).lt.-1.d0) phiau1(i)=-.98d0 
               enddo 
               call s_invdfe(phiau1,iqfin,theta, work3(1),ndim2)   
            endif 
         endif 
* 
*     We call s_fc12fe to compute the goal function vtau. 
* 
         kk=k 
         if (iqfin.gt.0) kk=ipfin
         call s_fc12fe(phiaux,theta, thetas,n,beta,cck,idif,isp,nsd, 
     +        m,kk,iqfin,n0,indth,x,y,sigman,thesig,vtau, 
     +        sigini,tau,xy,yhat,uhat,epshat,st,epspred,w, 
     +        auxm,npred,poldif,ndim2,work3,nw3,work4,nw4, 
     +        iwork4,niw4,work5,nw5,iwork5,niw5) 
* 
*     We determine if the value of thetas is a new optimal. 
* 
         if (vtau.lt.vtaumin.or.l.eq.1) then 
            sigmin=sigman 
            vtaumin=vtau 
            do i=1,kk 
               phidbes(i)=phiaux(i) 
            enddo 
            do i=1,iqfin 
               thetabes(i)=theta(i) 
            enddo 
            themej=thetas 
         endif 
      enddo 
* 
*     End of the block 1. 
* 
*     Begining of block 2.  
* 
*     We look for the minimum whith 3 steps of bisection. 
*     Block 2  has identical stucture than block 1. 
* 
      do 1800 ivez=1,3 
         the(1)=themej-delta(ivez) 
         the(2)=themej+delta(ivez) 
         do l=1,2 
            thetas=the(l) 
            thesig=sigmau 
            if (thetas.le.0.99d0.and.thetas.ge.-0.99d0) then 
               n1=ndim2 
               n2=n1+ip*ip
               call s_sys1fe(phidif1,ip,thetas,rho1,k,isp,phithe,tau,
     +              ier,work3(1),work3(n1+1),work3(n2+1),iwork3,ndim2) 
               if (ier.eq.1) go to 1800       
               if (iqfin.gt.0) then 
                  n1=ndim2 
                  n2=n1+(ndim2+1)*(ndim2+1) 
                  n3=n2+(ndim2+1)
                  call s_sys2fe(phithe,theta,thetas,k,0,isp,indth,rho, 
     +                 tau,work3(1),work3(n1+1),work3(n2+1), 
     +                 work3(n3+1),iwork3,ndim2) 
                  call s_sys3fe(phithe,k,phiaux,theta,thetas,ipfin, 
     +                 iqfin,isp,indth,rho,tau,ndim2,work3(1), 
     +                 work3(ndim2*ndim2+1),iwork3)          
               else 
                  do i=1,ndim2 
                     phiaux(i)=zero
                  enddo 
                  do i=1,k 
                     phiaux(i)=phithe(i) 
                  enddo 
               endif 
               if (ipfin.gt.0) then 
                  call s_yulefe(phiaux,rho,ipfin,work3(1),iiwork,ndim2) 
                  call s_durbfe(rho,ipfin,phiau1,ier,work3(1),ndim2) 
                  if (ier.eq.1) then 
                     do i=1,ipfin 
                        if (phiau1(i).gt. 1.d0) phiau1(i)= .98d0 
                        if (phiau1(i).lt.-1.d0) phiau1(i)=-.98d0 
                     enddo
                     call s_invdfe(phiau1,ipfin,phiaux,work3(1),ndim2)   
                  endif 
               endif 
               if (iqfin.gt.0) then 
                  call s_yulefe(theta,rho,iqfin,work3(1),iiwork,ndim2) 
                  call s_durbfe(rho,iqfin,phiau1,ier,work3(1),ndim2) 
                  if (ier.eq.1) then 
                     do i=1,iqfin 
                        if (phiau1(i).gt. 1.d0) phiau1(i)= .98d0 
                        if (phiau1(i).lt.-1.d0) phiau1(i)=-.98d0 
                     enddo 
                     call s_invdfe(phiau1,iqfin,theta, work3(1),ndim2)   
                  endif 
               endif 
               kk=k 
               if (iqfin.gt.0) kk=ipfin 
               call s_fc12fe(phiaux,theta, thetas,n,beta,cck,idif,isp, 
     +              nsd,m,kk,iqfin,n0,indth,x,y,sigman,thesig, 
     +              vtau,sigini,tau,xy,yhat,uhat,epshat,st, 
     +              epspred,w,auxm,npred,poldif,ndim2,work3, 
     +              nw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5)  
               if (vtau.lt.vtaumin) then 
                  sigmin=sigman 
                  vtaumin=vtau 
                  do i=1,kk 
                     phidbes(i)=phiaux(i) 
                  enddo 
                  do i=1,iqfin 
                     thetabes(i)=theta(i) 
                  enddo 
                  themej=thetas 
               endif 
            endif
         enddo 
 1800 continue 
      do i=1,ipfin 
         phidif1(i)= phidbes(i) 
      enddo 
      do i=1,iqfin 
         theta(i)=thetabes(i) 
      enddo 
      thetas=themej 
      sigmau=sigmin 
      return 
      end 
C=======================================================================
      subroutine s_gesvfe(n, nrhs, a, lda, ipiv, b, ldb, info) 
*----------------------------------------------------------------------- 
*     This subroutine computes the solution to a real system of linear 
*     equations 
*                  a * x = b, 
*     where a is an n-by-n matrix and x and b are n-by-nrhs matrices. 
*     
*     The lu decomposition with partial pivoting and row interchanges is 
*     used to factor a as 
*               a = p * l * u, 
*     where p is a permutation matrix, l is unit lower triangular, and 
*     u is upper triangular. The factored form of a is then used to 
*     solve the system of equations a * x = b. 
* 
*     Arguments 
*
*     a       (input/output) double precision array, dimension (lda,n) 
*             On entry, the n-by-n coefficient matrix a. 
*             On exit, the factors l and u from the factorization 
*             a = p*l*u; the unit diagonal elements of l are not stored. 
* 
*     ipiv    (output) integer array, dimension (n) 
*             the pivot indices that define the permutation matrix p; 
*             row i of the matrix was interchanged with row ipiv(i). 
* 
*     b       (input/output) double precision array, dimension (ldb,nrhs) 
*             On entry, the n-by-nrhs matrix of right hand side matrix b. 
*             On exit, if info = 0, the n-by-nrhs solution matrix x. 
* 
*     info    (output) integer 
*              = 0:  successful exit 
*              < 0:  if info = -i, the i-th argument had an illegal value 
*              > 0:  if info = i, u(i,i) is exactly zero. The 
*                    factorization has been completed, but the factor u 
*                    is exactly singular, so the solution could not be 
*                    computed. 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      integer ipiv( * ) 
      dimension a(lda, *), b(ldb, *)
*----------------------------------------------------------------------- 
*
*     Test the input parameters. 
*     
      info = 0 
      if( n.lt.0 ) then 
          info = -1 
      else if( nrhs.lt.0 ) then 
          info = -2 
      else if( lda.lt.max( 1, n ) ) then 
          info = -4 
      else if( ldb.lt.max( 1, n ) ) then 
          info = -7 
      end if 
      if( info.ne.0 ) return 
* 
*     Compute the lu factorization of a. 
* 
      call dgetf2(n, n, a, lda, ipiv, info) 
      if( info.eq.0 ) then 
* 
*     Solve the system a*x = b, overwriting b with x. 
*     
          call dgetrs( 'no transpose', n, nrhs, a, lda, ipiv, b, ldb, 
     +        info ) 
      end if 
      return 
      end 
C=======================================================================
      subroutine s_grd1fe(x,y,interc,n,m,idif,isp,nsd,iqfin,phi,beta, 
     +     phidif1,sigini,sigmau,akai,cck,parold,ssqold, 
     +     xy,yhat,uhat,epshat,st,epspred,w,auxm,npred, 
     +     poldif,ndim1,ndim2,beta0,beta1,prod,phibes, 
     +     phiaux,tau,rho,phiau2,yy,theta,fjac,qtf,f, 
     +     wa1,wa2,wa3,wa4,diag,par,aux2,ypure,ipvt,npo, 
     +     work3,nw3,iwork3,niw3,work4,nw4,iwork4,niw4, 
     +     work5,nw5,iwork5,niw5,utol,maxfev,epsmch,dwarf) 
*----------------------------------------------------------------------- 
*     This subroutine calculates the AR and the regression parameters when 
*     it is asummed that the errors of the regression follow an AR(1) or 
*     ARI(1,d) model. 
* 
*     Input:  
*               x    : matrix of independent variables 
*               y    : response or dependent variable (series which 
*                      will be fitted)  
*               n    : length of the series 
*               m    : number of independent variables in the linear model 
*               idif : number of ordinary differences 
*               isp  : seasonal period 
*               nsd  : number of seasonal differences 
*               iqfin: order of the MA polynomial, required to dimension 
*                      auxiliary arrays 
*               ndim1: ip+iqfin+m+indth+1, required to dimension several 
*                      arrays 
*               ndim2: max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required 
*                      to dimension several arrays 
*               utol : We make ftol=xtol=gtol=utol in the optimizer
*                      soubroutine s_lmdffe. 
*               maxfev: Maximum number of calls to the function
*                       which calculates the pseudo likelihood in s_lmdffe
* 
*     Output: 
*               phi    : Row 1 contains the phi estimator for the 
*                        non-stationary model  
*               beta   : regression coefficients estimators 
*               phidif1: phi estimator for the stationary model of order 1 
*               sigini : scale of the differenced regression errors 
*               sigmau : scale estimator of the filtered residuals 
*               akai   : value of the robust Akaike criterion 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phi(ndim2,ndim2),beta(m),parold(ndim1) 
      dimension xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred),w(n+npred),auxm(n+npred,ndim2) 
      dimension poldif(isp*nsd+idif+1) 
      dimension beta0(m),beta1(m),fi0(21),fi1(7),delt1(3),fi00(3) 
      dimension prod(ndim2+1),phibes(ndim2),phiaux(ndim2), 
     +     polphi(2),tau(0:ndim2),rho(0:ndim2),phiau1(2), 
     +     phiau2(ndim2+1),yy(n),theta(iqfin),fjac(n,ndim1), 
     +     qtf(ndim1),f(n),wa1(ndim1),wa2(ndim1),wa3(ndim1), 
     +     wa4(n),diag(ndim1),par(ndim1),aux2(4*n),ypure(n) 
      dimension ipvt(ndim1),npo(1+iqfin) 
      dimension work3(nw3),work4(nw4),work5(nw5) 
      dimension iwork3(niw3),iwork4(niw4),iwork5(niw5) 
      dimension para(ndim2), para1(ndim2), rh(ndim2+1)
      external s_fc11fe
      data zero,one/0.d0,1.d0/
      data fi0 /0.5d0, -0.1d0,  0.9d0,  0.6d0,  0.8d0,  0.4d0,  0.99d0,
     +          0.3d0,  0.2d0,  0.1d0,  0.0d0,  0.7d0, -0.2d0, -0.30d0,
     +         -0.4d0, -0.5d0, -0.6d0, -0.7d0, -0.8d0, -0.9d0, -0.99d0/  
      data fi1 /0.0d0,  0.1d0, -0.1d0,  0.2d0, -0.2d0,  0.3d0, -0.30d0/ 
      data delt1/.05d0,.025d0,.0125d0/    

c initialize variables defined within if
      vtaumin = 0.d0
      vtaumin1 = 0.d0
      phi1 = 0.d0

*-----------------------------------------------------------------------  
* 
*     s_lmdffe parameters. 
* 
      ftol=zero
      xtol=utol
      gtol=zero
      mode=1 
      nprint=0 
      factor=100.d0 
      ldfjac=n 
*  
*     s_grd1fe is divided in four blocks. The first block minimizes the goal 
*     function over a gross grid. In the second block we do a similar 
*     minimization using a finer grid. In the third block we minimize the 
*     goal function by three steps of minimization. In the fourth block 
*     we minimize the goal function using the Marquard subroutine s_lmdffe. 
* 
*     Here it begins the first block. 
* 
      ndiff=isp*nsd+idif 
      iq=0 
      rho(0)=one
* 
*     We look for the minimum of the filtered residuals' scale along a 
*     grid of phi values which are in fi0. 
* 
      do i=1,21 
         phidif1=fi0(i) 
*     
*     We calculate the product of two polynomials:  
*     1-fi0(i)*B and polds = (1-B)**d * (1-B**isp)**nsd. 
*     ktrue is the degree of the product. Then we obtain the nonstationary 
*     polynomial. 
* 
         if ((idif.eq.0).and.(nsd.eq.0)) then  
            ktrue=1 
            phi(1,1)=fi0(i) 
         else 
            polphi(1)=one 
            polphi(2)=-fi0(i) 
            call s_polyfe(poldif,ndiff,polphi,1,prod,ktrue) 
            do ii=1,ktrue 
               phi(ktrue,ii)=-prod(ii+1)  
            enddo 
         endif 
* 
*     We estimate the regression vector beta. 
*     We create matrix xy filtering the x's and y with the nonstationary  
*     autoregresive polynomial with coefficients phi(ktrue,.). 
* 
         if (m.gt.0) then 
            do it=1,n-ktrue 
               do j=1,m 
                  sum=zero 
                  do ii=1,ktrue 
                     sum=sum+phi(ktrue,ii)*x(it+ktrue-ii,j)            
                  enddo 
                  xy(it,j) = x(it+ktrue,j) - sum  
               enddo 
               sum=zero 
               do ii=1,ktrue 
                  sum=sum+phi(ktrue,ii)*y(it+ktrue-ii)            
               enddo 
               xy(it,m+1)= y(it+ktrue) - sum 
               yy(it)=xy(it,m+1) 
            enddo 
            eps=0.001d0 
            itmax=20 
* 
*     We call the robust regression procedure. 
* 
         call s_rqr1fe(n-ktrue,m,xy,yy,eps,itmax,beta,sumre,
     +           interc,n,work3,nw3,iwork3,niw3) 
         endif 
* 
*     We calculate an initial estimate of the residuals scale. 
* 
         sigini=s_xmadfe(x,y,beta,m,n,aux2(1),aux2(n+1),aux2(2*n+1), 
     +        poldif,ndiff) 
* 
*     We compute the innovations scale. 
* 
         sigmau=sigini*dsqrt(one-phidif1*phidif1) 
         rho(1)=phidif1 
         do iii=1,ndim2 
            phiaux(iii)=zero 
         enddo  
         phiaux(1)=phidif1 
* 
*     We compute the goal function for the current solution. 
* 
         call s_fc12fe(phiaux,theta,thetas,n,beta,one,idif,isp,nsd,m,
     +        1,0,n0,0,x,y,sigman,sigmau,vtau,sigini,tau,xy,yhat, 
     +        uhat,epshat,st,epspred,w,auxm,npred,poldif,ndim2, 
     +        work3,nw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
* 
*     We look for the minimum, comparing the current value with the minimum. 
* 
         if (vtau.lt.vtaumin.or.i.eq.1) then 
            vtaumin=vtau 
            sigma0=sigman 
            phi0=fi0(i)  
            do j=1,m 
               beta0(j)=beta(j) 
            enddo 
            sigini0=sigini 
         endif  
      enddo 
* 
*     The initial grid ends here. We have just computed the initial estimates 
*     phi0, beta0 and sigma0 (an estimate of the filtered residuals' scale). 
* 
*     We proceed in a similar way user a finer grid now given in f1,  
*     but now we filter the y(t) too. 
* 
      do i=1,7 
         phidif1=phi0+fi1(i) 
         if (phidif1 .eq.  one) phidif1= .99d0 
         if (phidif1 .eq. -one) phidif1=-.99d0 
         if (phidif1 .le.  .99d0 .and. phidif1 .ge. -.99d0) then 
*     
*     In the following block we calculate the product of two polynomials: 
*     1-fi0(i)*B and polds = (1-B)**d * (1-B**isp)**nsd. 
*     ktrue is the degree of the product. Then we obtain the nonstationary 
*     polynomial. 
* 
            if ((idif.eq.0).and.(nsd.eq.0)) then  
               ktrue=1 
               phi(1,1)=phidif1 
            else 
               polphi(1)=one
               polphi(2)=-phidif1 
               call s_polyfe(poldif,ndiff,polphi,1,prod,ktrue) 
               do ii=1,ktrue 
                  phi(ktrue,ii)=-prod(ii+1)  
               enddo 
            endif 
*     
*     We filter the observed y(t), setting beta=beta0. 
* 
            rho(1)=phidif1 
            do iii=1,ndim2 
               phiaux(iii)=zero 
            enddo  
            do iii=1,ktrue 
               phiaux(iii)=phi(ktrue,iii) 
            enddo 
* 
*     We filter twice, the first time with nfil=0 and the second with  
*     nfil=1. The first time we filter robustly the residuals  
*     getting the values uhat. The weights used to get uhat are stored on w. 
*     The second time we use a linearized version of the filter using  
*     these weights w to filter the x's and y, obtaining a matrix xy. 
* 
            call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux,beta0,theta, 
     +           thetas,ktrue,iq,sigma0,0,n0,tau,sigini0,0,rho, 
     +           one,0,ypure,xy,yhat,uhat,epshat,st,epspred, 
     +           w,auxm,ndim2,work3,nw3,work5,nw5,iwork5,niw5) 
            if (m.gt.0) then 
               call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux,beta0,theta, 
     +              thetas,ktrue,iq,sigma0,0,n0,tau,sigini0, 
     +              1,rho,one,0,ypure,xy,yhat,uhat,epshat,st, 
     +              epspred,w,auxm,ndim2,work3,nw3,work5,nw5,
     +              iwork5,niw5) 
* 
*     We compute beta. 
* 
* 
               do it=1,n-ktrue 
                  yy(it)=xy(it,m+1) 
               enddo 
* 
*     We estimate the regression vector beta for a linear regression 
*     of y(t)-phi*yhat(t-1) over x(t)-phi*x(t-1). These are in the 
*     matrix xy computed by filter. 
*     We apply the robust regression to the filtered observations. 
* 
               call s_rqr1fe(n-n0,m,xy,yy,eps,itmax,beta,sumre, 
     +              interc,n,work3,nw3,iwork3,niw3) 
* 
*     We compute a robust scale of the differenced residuals. 
* 
               sigini=s_xmadfe(x,y,beta,m,n,aux2(1),aux2(n+1),
     +              aux2(2*n+1),poldif,ndiff) 
            endif 
* 
*     We compute the innovations scale. 
*     
            sigmau=sigini*dsqrt(one-phidif1*phidif1) 
            rho(1)=phidif1 
* 
*     We assign the optimal value of the nonstationary AR polynomial 
*     to phiaux and compute the optimal value of the goal function. 
* 
            do iii=1,ndim2 
               phiaux(iii)=zero 
            enddo  
            do iii=1,ktrue 
               phiaux(iii)=phi(ktrue,iii) 
            enddo 
            phiaux(1)=phidif1 
* 
*     We compute the goal function for the current solution. 
* 
            call s_fc12fe(phiaux,theta, thetas,n,beta,one,idif,isp,
     +           nsd,m,1,0,n0,0,x,y,sigman,sigmau,vtau,sigini,tau, 
     +           xy,yhat,uhat,epshat,st,epspred,w,auxm,npred, 
     +           poldif,ndim2,work3,nw3,work4,nw4,iwork4,niw4, 
     +           work5,nw5,iwork5,niw5) 
* 
*     We  determine if the current solution is better than the old optimal. 
* 
            if (vtau.lt.vtaumin1.or.i.eq.1) then 
               vtaumin1=vtau 
               sigma1=sigman 
               phi1=phidif1 
               do j=1,m 
                  beta1(j)=beta(j) 
               enddo 
               do ii=1,ktrue 
                  phibes(ii)=phi(ktrue,ii)  
               enddo 
               sigini1=sigini 
            endif  
         endif 
      enddo 
* 
*     End of the second block.    
* 
*     Beginning of the third block. 
* 
*     We apply three steps of the bisection algorithm to minimize the  
*     goal function. 
* 
      do ivez=1,3 
* 
*     We set the values of the autoregressive parameter. 
* 
         fi00(1)=phi1-delt1(ivez) 
         fi00(2)=phi1+delt1(ivez) 
         do i=1,2 
            phidif1=fi00(i) 
            if (phidif1 .le. .99d0 .and. phidif1 .ge. -.99d0) then 
* 
*     In the following block we calculate the product of two polynomials: 
*     1-fi00(i)*B and polds = (1-B)**d * (1-B**isp)**nsd. 
*     ktrue is the degree of the product. Then we obtain the nonstationary 
*     polynomial. 
* 
               if ((idif.eq.0).and.(nsd.eq.0)) then  
                  ktrue=1 
                  phi(1,1)=fi00(i) 
               else 
                  polphi(1)=one 
                  polphi(2)=-fi00(i) 
                  call s_polyfe(poldif,ndiff,polphi,1,prod,ktrue) 
                  do ii=1,ktrue 
                     phi(ktrue,ii)=-prod(ii+1)  
                  enddo 
               endif 
               rho(1)=phidif1 
               do iii=1,ndim2 
                  phiaux(iii)=zero 
               enddo  
               do iii=1,ktrue 
                  phiaux(iii)=phi(ktrue,iii) 
               enddo 
* 
*     We filter twice, the first time with nfil=0 and the second time with 
*     with nfil=1. The first time we filter robustly the residuals 
*     getting the values uhat. The weights used to get uhat are stored on w. 
*     The second time we use a linearized version of the filter using  
*     these weights w to filter the x's and the y obtaining a matrix xy. 
* 
               call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux,beta1,theta, 
     +              thetas,ktrue,iq,sigma1,0,n0,tau,sigini1, 
     +              0,rho,one,0,ypure,xy,yhat,uhat,epshat,st, 
     +              epspred,w,auxm,ndim2,work3,nw3,work5,nw5, 
     +              iwork5,niw5)  
               rho(1)=phidif1 
               if (m.gt.0) then 
                  call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux,beta1,
     +                 theta,thetas,ktrue,iq,sigma1,0,n0,tau,sigini1, 
     +                 1,rho,one,0,ypure,xy,yhat,uhat,epshat, 
     +                 st,epspred,w,auxm,ndim2,work3,nw3,work5, 
     +                 nw5,iwork5,niw5)
* 
*     We compute beta. 
* 
* 
                  do it=1,n-ktrue 
                     yy(it)=xy(it,m+1) 
                  enddo 
* 
*     We estimate the regression vector beta for a linear regression 
*     of y(t)-phi*yhat(t-1) over x(t)-phi*x(t-1). These are in the 
*     matrix xy computed by filter. 
*     We apply the robust regression to the filtered observations. 
* 
                  call s_rqr1fe(n-n0,m,xy,yy,eps,itmax,beta,sumre, 
     +                 interc,n,work3,nw3,iwork3,niw3) 
* 
*     We compute a robust scale of the differenced residuals. 
* 
                  sigini=s_xmadfe(x,y,beta,m,n,aux2(1),aux2(n+1), 
     +                 aux2(2*n+1),poldif,ndiff) 
               endif 
* 
*     We compute the innovations scale.     
* 
               sigmau=sigini*dsqrt(one-phidif1*phidif1) 
               rho(1)=phidif1 
* 
*     We assign the optimal value of the nonstationary AR polynomial 
*     to phiaux and compute the optimal value of the goal function. 
* 
               do iii=1,ndim2 
                  phiaux(iii)=zero
               enddo  
               phiaux(1)=phidif1 
* 
*     We compute the goal function for the current solution. 
* 
               call s_fc12fe(phiaux,theta,thetas,n,beta,one,idif,
     +              isp,nsd, m,1,0,n0,0,x,y,sigman,sigmau,vtau,
     +              sigini,tau,xy,yhat,uhat,epshat,st,epspred,w,auxm,
     +              npred,poldif,ndim2,work3,nw3,work4,nw4,iwork4,
     +              niw4,work5,nw5,iwork5,niw5) 
* 
*     We determine if the current solution is better than the old optimal. 
* 
               if (vtau.lt.vtaumin1) then 
                  vtaumin1=vtau 
                  sigma1=sigman 
                  phi1=fi00(i) 
                  do j=1,m 
                     beta1(j)=beta(j) 
                  enddo 
                  do ii=1,ktrue 
                     phibes(ii)=phi(ktrue,ii)  
                  enddo 
                  sigini1=sigini 
               endif  
            endif 
         enddo 
      enddo 
* 
*     We have just computed  phi1, beta1, sigma1 and we allocate them 
*     in their right position. 
* 
      sigmau=sigma1 
      do j=1,m 
         beta(j)=beta1(j) 
      enddo 
      phidif1=phi1 
      sigini=sigini1 
* 
*     End of the third block. 
* 
*     Beginninig of the forth block. 
* 
      do ii=1,ktrue 
         phi(ktrue,ii)=phibes(ii)  
      enddo 
      phiaux(1)=phidif1 
*  
*     We set npar: the number of parameters of the optimization subroutine 
*     s_lmdffe and the filter bandwidth, par(npar), to 1.0. 
* 
      npar=1+m
      ccknew=one
* 
*     We call subroutine s_tranfe which computes partial autocorrelations 
*     and transforms them in variables that may take any real value. 
* 
      call s_trasfe(phiaux,theta,thetas,beta,ndim2,1,0,0,m,para,par,
     +     ndim1,rh,work3,nw3,iwork3,niw3,0) 
* 
*     The vector npo contains the powers of B with non zero coefficients  
*     in the AR and MA operators. In this version 0's are not allowed. 
* 
      do i=1,1+iqfin 
         npo(i)=0 
      enddo 
      npo(1)=1 
      itte=0       
* 
*     The minimization process start. We make a loop with three different 
*     starting values of the robust filter bandwidth. 
* 
      do while (itte.le.2) 
         itte=itte+1 
* 
*     s_lmdffe is the Marquard routine which minimizes the goal function. 
* 
         call s_lmdffe(s_fc11fe,n,npar,par,f,ftol,xtol,gtol,maxfev, 
     +        diag,mode,factor,nprint,info,nfev,fjac,ldfjac, 
     +        ipvt,qtf,wa1,wa2,wa3,wa4,idif,isp,nsd,m,1,0,n0,0, 
     +        npo,0,x,y,sigman,sigmanew,xy,yhat,uhat,epshat,st, 
     +        epspred,w,auxm,poldif,ccknew,ndim1,ndim2,work3,
     +        nw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5,
     +        epsmch,dwarf) 
         ssq=0 
         do jjj=n0+1,n 
            ssq=ssq+f(jjj)*f(jjj) 
         enddo 
*  
*     We compute in ssq the value of the goal function for the 
*     optimal solution. 
* 
*     We call s_fc11fe to compute two estimates of the innovation scale  
*     sigman y sigmanew. The difference between these two estimates is  
*     that sigmanew is computed only using the variance of the differenced  
*     stationary process and sigman corrects this value using the  
*     filtered residuals uhat. 
* 
         call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,1,0,n0,0,npo, 
     +        sigman,sigmanew,0,x,y,xy,yhat, ccknew,uhat, 
     +        epshat,st, epspred,w,auxm,poldif,ndim1,ndim2, 
     +        work3,nw3,work4,nw4,iwork4,niw4,work5,nw5,
     +        iwork5,niw5)  
* 
*     We compare the solution with new bandwith filter with the 
*     current optimal and decide if we should change it. 
* 
         if (itte.eq.1.or. ssqold.gt.ssq)then 
            cck=ccknew
            ssqold=ssq 
            do i=1,npar 
               parold(i)=par(i) 
            enddo 
         endif 
* 
*     We set the new value of the bandwidth filter. 
* 
         vv=(0.8d0**itte)*(sigmau/sigmanew) 
         ccknew=dmin1(vv,one) 
      enddo 
* 
*     We set the optimal parameters. 
* 
      do i=1,npar 
         par(i)=parold(i) 
      enddo 
      do i=1,ndim2 
         phiaux(i)=zero
      enddo  
      call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,1,0,n0,0,npo, 
     +     sigmau,sigmanew,0,x,y,xy,yhat,cck,uhat,epshat,st, 
     +     epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,
     +     work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
      vtau=zero
* 
*     We compute the pseudo tau-likelihood of the optimal solution.  
* 
      do i=1,n 
         vtau=vtau+f(i)*f(i) 
      enddo 
      znnn=dble(n-n0) 
      vtau=znnn*(dlog(vtau)-dlog(0.41d0)-dlog(znnn)) 
*     
*     We antitransform the autoregressive parameter from par, using  
*     subroutine s_tranfe. 
* 
      call s_tranfe(par,ndim1,ndim2,1,0,0,0,para,para1,work3,para, 
     +     theta,thetas,beta)
      phidif1=para(1)
      phiaux(1)=para(1)
* 
*     In the following block we calculate the product of two polynomials: 
*     1-fi00(i)*B and polds = (1-B)**d * (1-B**isp)**nsd. 
*     ktrue is the degree of the product. Then we obtain the nonstationary 
*     polynomial. 
* 
      if ((idif.eq.0).and.(nsd.eq.0)) then 
         do i=1,ktrue 
            phi(ktrue,i)=phidif1 
         enddo 
      else 
         phiau1(1)=one 
         phiau1(2)=-phidif1 
         call s_polyfe(poldif,ndiff,phiau1,1,phiau2,ktrue) 
         do i=1,ktrue 
            phi(ktrue,i)=-phiau2(i+1) 
         enddo 
      endif 
* 
*     We recover the regression parameters from the optimal solution. 
* 
      do i=1,m 
         beta(i)=par(1+i) 
      enddo         
* 
*     We compute the scale of the stationary AR component of the 
*     regression model, that is, of the regression errors after differencing. 
* 
      sigini=s_xmadfe(x,y,beta,m,n,aux2(1),aux2(n+1),aux2(2*n+1),
     +     poldif,ndiff) 
* 
*     We compute the robust AIC criterium. 
* 
      akai=vtau+2 
      return 
      end 
C=======================================================================
      subroutine s_grdkfe(x,y,n,m,ip,idif,isp,nsd,iqfin,phi,beta,
     +     thetas,phidif,k,sigini,sigmau,akai,cck,parold,ssqold, 
     +     xy,yhat,uhat,epshat,st,epspred,w,auxm,npred, 
     +     poldif,ndim1,ndim2,theta,phiaux,ydifh,ydift, 
     +     utilde,tau,rh,para,para1,par,rho,fjac,qtf,f, 
     +     wa1,wa2,wa3,wa4,diag,phiau1,phiau2,aux2,npo, 
     +     ipvt,work3,nw3,iwork3,niw3,work4,nw4,iwork4, 
     +     niw4,work5,nw5,iwork5,niw5,utol,maxfev,epsmch,dwarf) 
*-----------------------------------------------------------------------  
*     This subroutine estimates the AR and the regression coefficients 
*     when it is assumed that the errors follow an AR(k) model. 
* 
*     Input: 
*           x     : regression matrix 
*           y     : original series 
*           n     : length of the series 
*           m     : number of regression variables 
*           ip    : maximum order of the AR model 
*           idif  : number of ordinary differences 
*           isp   : seasonal period 
*           nsd   : number of seasonal differences 
*           iqfin : order of the MA part of the final model 
*           phi   : on input contains the estimated AR coefficients of the  
*                   precedent model 
*           beta  : on input contains the estimated regression coefficients  
*                   of the precedent model 
*           thetas: on input contains the previous estimated seasonal MA  
*                   coefficient 
*           phidif: on input contains the estimated coefficients of the  
*                   precedent stationary model 
*           k:      order of the present stationary AR model 
*           sigini: scale of the differenced  regression residuals 
*           sigmau: estimated scale of the estimated AR(k-1) model 
*           npred : number of predicted regression residuals 
*           poldif: vector containing the coefficients of the differences 
*                   polynomial 
*           ndim1 : ip+iqfin+m+indth+1, required to dimension several arrays 
*           ndim2 : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                   dimension several arrays
*           utol  : We make ftol=xtol=gtol=utol in the optimizer
*                   soubroutine s_lmdffe. 
*           maxfev: Maximum number of calls to the function
*                   which calculates the pseudo likelihood in s_lmdffe
* 
*     Output: 
*           phi:    on output contains the estimated AR coefficients of the  
*                   present model 
*           beta:   on output contains the estimated regression coefficients  
*                   of the present model 
*           thetas: on output contains the new estimated seasonal MA coeff 
*           phidif: on output contains the estimated coefficients of the  
*                   present stationary model 
*           sigmau: on output contains the estimated scale of the new 
*                   model 
*           akai  : robust AIC criterium 
*           cck   : bandwidth of the robust filter 
*           parold: optimal transformed parameter used in s_lmdffe 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phi(ndim2,ndim2),beta(m),phidif(ip+1,ip+1), 
     +     parold(ndim1),xy(n,m+1),yhat(n),uhat(n+npred),st(n+npred),
     +     epshat(n+npred),epspred(n+npred),w(n+npred),theta(iqfin),
     +     auxm(n+npred,ndim2),poldif(isp*nsd+idif+1),phiaux(ndim2),
     +     ydifh(n),ydift(n),utilde(n),tau(0:ndim2),rh(ip+1),para(ip),
     +     para1(ip),par(ndim1),rho(0:ndim1),fjac(n,ndim1),qtf(ndim1),
     +     f(n),wa1(ndim1),wa2(ndim1),wa3(ndim1),wa4(n),diag(ndim1), 
     +     phiau1(ip+1),phiau2(ndim2+1),aux2(4*n)
      dimension work3(nw3),work4(nw4),work5(nw5) 
      integer npo(ip+iqfin),ipvt(ndim1)
      integer iwork3(niw3),iwork4(niw4),iwork5(niw5) 
      external s_fc11fe 
      data zero/0.d0/
*-----------------------------------------------------------------------   
*     We will minimize a robust M-estimate of the scale of the residuals 
*     by computing this scale for each element of a grid on the interval  
*     (0,1). We will choose the minimum among these scales. 
* 
      ftol=zero
      xtol=utol
      gtol=zero
      mode=1 
      nprint=0 
      factor=100.d0 
      ldfjac=n 
      ndiff=isp*nsd+idif 
      iq=0 
      ktrue=k+ndiff 
      ktrue1=ktrue-1 
      if (k.eq.1) then 
         do j=1,ktrue1 
            phi(ktrue1,j)=-poldif(j+1) 
         enddo 
      endif 
      do iii=1,ndim2 
         phiaux(iii)=zero 
      enddo  
      do iii=1,ktrue1 
         phiaux(iii)=phi(ktrue1,iii) 
      enddo 
* 
*     We filter using the precedent model. 
* 
      call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux,beta,theta,thetas, 
     +     ktrue1, iq,sigmau,0,n0,tau,sigini,0,rho,1.d0,0, 
     +     aux2,xy,yhat,uhat,epshat,st,epspred,w,auxm, 
     +     ndim2,work3,nw3,work5,nw5,iwork5,niw5) 
* 
*     We calculate the forward residuals : utilde.    
*     aux2 contains y-x*beta. 
* 
*     ydift and ydifh are estimates of the stationary process obtained by 
*     differencing the regression errors. ydift(j) is completely cleaned of 
*     the influencing outliers and ydifh(j) is cleaned with exception of y(j). 
* 
      do j=ndiff+1,n 
         ydifh(j)=aux2(j) 
         ydift(j)=epshat(j) 
         do jj=1,ndiff 
            ydifh(j)=ydifh(j)+poldif(jj+1)*epshat(j-jj) 
            ydift(j)=ydift(j)+poldif(jj+1)*epshat(j-jj) 
         enddo 
      enddo 
* 
*     utilde are backward residuals. 
* 
      do j=ndiff+1,n-(k-1) 
         utilde(j)=ydifh(j) 
         do jj=1,k-1 
            utilde(j)=utilde(j)-phidif(k-1,jj)*ydift(j+jj) 
         enddo 
      enddo 
* 
*     We compute the ratio between forward and backward residuals.  
* 
      nf=n-ndiff-k 
      lj=0 
      do j=ndiff+k+1,n 
         if (dabs(utilde(j-k)).ge.0.1d-10) then  
            lj=lj+1 
            wa4(lj)=uhat(j)/utilde(j-k) 
         else 
            nf=nf-1 
         endif    
      enddo   
* 
*     We compute an estimate of the partial autocorrelation of order k. 
* 
      call s_mednfe(wa4,nf,xmed,aux2) 
      if(xmed.gt.0.99d0)   xmed= 0.99d0 
      if(xmed.lt.-0.99d0)  xmed=-0.99d0 
* 
*     We compute the robust Durbin-Levinson coefficients. 
* 
      phidif(k,k)=xmed 
      if (k.gt.1) then 
         do i=1,k-1 
            phidif(k,i)=phidif(k-1,i)-phidif(k,k)*phidif(k-1,k-i) 
         enddo 
      endif  
      do i=1,k 
         phiaux(i)=phidif(k,i) 
      enddo 
* 
*     We set npar: the number of parameters of the optimization subroutine 
*     s_lmdffe and the filter bandwidth, par(npar), to 1.0. 
*
      npar=k+m
      ccknew=1.d0
* 
*     We call subroutine s_tranfe which computes partial autocorrelations 
*     and transforms them in variables that may take any real value. 
* 
      call s_trasfe(phiaux,theta,thetas,beta,ndim2,k,0,0,m,para,par,
     +     ndim1,rh,work3,nw3,iwork3,niw3,0) 
* 
*     The vector npo contains the powers of B with non zero coefficients  
*     in the AR and MA operators. In this version 0's are not allowed. 
* 
      do i=1,ip+iqfin 
         npo(i)=0 
      enddo 
      do i=1,k 
         npo(i)=i 
      enddo 
      itte=0 
* 
*     The minimization process starts. We make a loop with three different 
*     starting values of the robust filter bandwith. 
* 
      do while (itte.le.2) 
         itte=itte+1
* 
*     We compute the initial value of s_fc11fe. Since s_fc11fe gives only 
*     pseudoresiduals in vector f, ssq computes the value of the goal  
*     function. 
* 
         call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,k,0,n0,0,npo, 
     +        sigman,sigmanew,0,x,y,xy,yhat,ccknew,uhat,epshat,
     +        st,epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,
     +        work4, nw4,iwork4,niw4,work5,nw5,iwork5,niw5)  
         if (itte.eq.1)then 
            ssq=0 
            do jjj=n0+1,n 
               ssq=ssq+f(jjj)*f(jjj) 
            enddo 
* 
*     We compute the value of the goal function for the AR(k-1) optimal  
*     solution in order to decide which solution should we use to start the 
*     minimization process. 
* 
            do ij=1,npar-k 
               parold(npar-ij+1)=parold(npar-ij) 
            enddo 
            parold(k)=zero
            call s_fc11fe(n,npar,parold,f,iflag,idif,isp,nsd,m,k,0,n0, 
     +           0,npo,sigman,sigmanew,0,x,y,xy,yhat,cck,uhat,epshat,
     +           st,epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,work4,
     +           nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
            ssq1=0 
            do jjj=n0+1,n 
               ssq1=ssq1+f(jjj)*f(jjj) 
            enddo      
            if (ssq1.lt.ssq)then 
               do ij=1,npar 
                  par(ij)=parold(ij) 
               enddo 
               ccknew=cck
            endif 
         endif 
* 
*     s_lmdffe is the Marquard routine which minimizes the goal function. 
* 
         call s_lmdffe(s_fc11fe,n,npar,par,f,ftol,xtol,gtol,maxfev, 
     +        diag,mode,factor,nprint,info,nfev,fjac,ldfjac, 
     +        ipvt,qtf,wa1,wa2,wa3,wa4,idif,isp,nsd,m,k,0,n0,0, 
     +        npo,0,x,y,sigman,sigmanew,xy,yhat,uhat,epshat,st, 
     +        epspred,w,auxm,poldif,ccknew,ndim1,ndim2,work3,nw3,
     +        work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5,epsmch,
     +        dwarf) 
* 
*     We compute in ssq the value of the goal function for the 
*     optimal solution. 
* 
         ssq=zero
         do jjj=n0+1,n 
            ssq=ssq+f(jjj)*f(jjj) 
         enddo 
* 
*     We call s_fc11fe to compute two estimates of the innovation scale: 
*     sigman y simanew. The difference between these two estimates is  
*     that sigmanew is computed only using the variance of the differenced  
*     stationary process and sigman corrects this value using the  
*     filtered residuals uhat. 
*
         call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,k,0,n0,0,npo, 
     +        sigman,sigmanew,0,x,y,xy,yhat,ccknew,uhat,epshat,
     +        st,epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,
     +        work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
* 
*     We compare the solution with new bandwith filter with the 
*     current  optimal and decide if we should change it. 
* 
         if (itte.eq.1.or. ssqold.gt.ssq) then
            cck=ccknew 
            ssqold=ssq 
            do i=1,npar 
               parold(i)=par(i) 
            enddo 
         endif 
* 
*     We set the new value of the bandwith filter. 
* 
         vv=(0.8d0**itte)*(sigmau/sigmanew) 
         ccknew=dmin1(vv,1.d0) 
      enddo 
*        
*     We set up the autoregressive parameters. 
* 
      do i=1,npar 
         par(i)=parold(i) 
      enddo 
      do i=1,ndim2 
         phiaux(i)=zero
      enddo
      call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,k,0,n0,0,npo, 
     +     sigmau,sigmanew,0,x,y,xy,yhat,cck,uhat,epshat,st, 
     +     epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,work4, 
     +     nw4,iwork4,niw4,work5,nw5,iwork5,niw5) 
      vtau=zero
* 
*     We compute the pseudo tau- likelihood of the optimal solution.  
* 
      do i=1,n 
         vtau=vtau+f(i)*f(i) 
      enddo 
      znnn=dble(n-n0) 
      vtau=znnn*(dlog(vtau)-dlog(0.41d0)-dlog(znnn)) 
* 
*     We antitransform the autoregressive parameters from par, using  
*     subroutine s_tranfe. 
* 
      call s_tranfe(par,ndim1,ndim2,k,0,0,0,para,para1,work3,phiaux, 
     +     theta,thetas,beta) 
      do i=1,k 
         phidif(k,i)=phiaux(i) 
*     phiaux(i)=para(i) 
      enddo 
* 
*     We compute the non-stationary autoregressive polynomial operator. 
* 
      if ((idif.eq.0).and.(nsd.eq.0)) then 
         do i=1,ktrue 
            phi(ktrue,i)=phidif(k,i) 
         enddo 
      else 
         phiau1(1)=1.0d0 
         do i=2,k+1 
            phiau1(i) = -phidif(k,i-1) 
         enddo 
         call s_polyfe(poldif,ndiff,phiau1,k,phiau2,ktrue) 
         do i=1,ktrue 
            phi(ktrue,i)=-phiau2(i+1) 
         enddo 
      endif 
* 
*     We recover the regression parameters from the optimal solution. 
* 
      do i=1,m 
         beta(i)=par(k+i) 
      enddo         
* 
*     We compute the scale of the stationary AR component of the 
*     regression model, that is, of the regression errors after differencing.  
* 
      sigini=s_xmadfe(x,y,beta,m,n,aux2(1),aux2(n+1),aux2(2*n+1),
     +     poldif,ndiff) 
* 
*     We compute the robust AIC criterion. 
* 
      akai=vtau+2*k 
      return 
      end 
C=======================================================================
      subroutine s_impcfe(x,y,m,idif,isp,nsd,phi,beta,theta,thetas,
     +     ktrue,q,sigmau,indth,n0,tau,sigmadif,rho,d,n,resid,psils,
     +     psia,waux,lambaux,cck,sigfil,yhat,uhat,epshat,st,lscan,w,
     +     auxm,xy,m1,idim,work2,idimw2,work3,idimw3,iwork,idimiw,
     +     iw2then,iw2xm,iw2xp,iw2alfp,iw2r,iw2v,iw2av,iw2alff,iw2rv,
     +     iw2alfx,iw2xaux,iw2v1,iw2v2,iw2ypur)
*-----------------------------------------------------------------------
*     This subroutine computes the impact of each outlier candidate
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer t,d,q
      integer lscan(n),iwork(idimiw)
      dimension x(n,m),y(n),phi(idim+1),beta(m),theta(q),resid(n)
      dimension rho(0:idim),tau(0:q+isp*indth),psia(n),psils(n)
      dimension yhat(n),uhat(n),epshat(n),st(n),w(n),waux(3,n)
      dimension xy(n,m1+1),auxm(n,idim),work2(idimw2),work3(idimw3)
      double precision lambaux(3,n)
      data zero,one/0.d0,1.d0/
*----------------------------------------------------------------------- 
*     We call s_flt2fe with iout=d (we use the robust filter for the 
*     observations
*     t=1,...,d-1, and t=d+h+1,...n, and the classical filter from  
*     t=d,...,d+h)
*
      call s_flt2fe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,ktrue,
     +     q,sigfil,indth,n0,tau,sigmadif,0,rho,cck,d,yhat,uhat,epshat,
     +     st,lscan,w,auxm,xy,m1,work2(iw2then),work2(iw2xm),
     +     work2(iw2xp),work2(iw2alfp),work2(iw2r),work2(iw2v),
     +     work2(iw2av),work2(iw2alff),work2(iw2rv),work2(iw2alfx),
     +     work2(iw2xaux),work2(iw2v1),work2(iw2v2),work2(iw2ypur),
     +     idim,work3,idimw3,iwork,idimiw)
      do i=1,n
         resid(i)=uhat(i)
      enddo
      do i=n0+1,n
         uhat(i)=uhat(i)/st(i)
      enddo
      call s_calsfe(uhat,n,n0,sigmanew,work2(1),work3)
*     
*     We define the indicators vector for the additive outliers and the 
*     level shifts.
*     
      do t=1,d-1
         psia(t)=zero
         psils(t)=zero
      enddo
      psia(d)=one
      psils(d)=one
      do t=d+1,n
         psia(t)=zero
         psils(t)=one
      enddo
*
*     We call s_flt2fe (with the indicator indfil=1) to apply to the vectors 
*     psia and psils the same transformation that was applied to the vector 
*     eps(t) 
*   
      indfil=1
      call s_flt2fe(psia,psils,n,1,idif,isp,nsd,phi,beta,theta,thetas,
     +     ktrue,q,sigmau,indth,n0,tau,sigmadif,indfil,rho,cck,d,yhat,
     +     uhat,epshat,st,lscan,w,auxm,xy,m1,work2(iw2then),
     +     work2(iw2xm),work2(iw2xp),work2(iw2alfp),work2(iw2r),
     +     work2(iw2v),work2(iw2av),work2(iw2alff),work2(iw2rv),
     +     work2(iw2alfx),work2(iw2xaux),work2(iw2v1),work2(iw2v2),
     +     work2(iw2ypur),idim,work3,idimw3,iwork,idimiw)
*     
*     The transformed psia and psils are placed by subroutine s_flt2fe in
*     matrix xy. Using this we compute waux(i,d) (i=1, the impact of the
*     innovation outlier, i=2 of the additive outlier, i=3 of the level
*     shift) and lambaux(i,d) (waux(i,d)/standard error)
*     
*     Fist we correct the st(i) so that the stardardized residuals have
*     scale estimate = 1
*
      do i=n0+1,n
         st(i)=st(i)*sigmanew
      enddo
      waux(1,d)=resid(d)
      waux(2,d)=zero
      waux(3,d)=zero
      dena=zero
      denls=zero 
      do t=d,n
         vart=st(t)*st(t)
         waux(2,d)=waux(2,d)+resid(t)*xy(t-n0,1)/vart
         waux(3,d)=waux(3,d)+resid(t)*xy(t-n0,2)/vart
         dena=dena+xy(t-n0,1)*xy(t-n0,1)/vart
         denls=denls+xy(t-n0,2)*xy(t-n0,2)/vart
      enddo  
      waux(2,d)=waux(2,d)/dena
      waux(3,d)=waux(3,d)/denls
      vara=one/dena
      varls=one/denls
      lambaux(1,d)=dabs(waux(1,d))/st(d)
      lambaux(2,d)=dabs(waux(2,d))/dsqrt(vara)
      lambaux(3,d)=dabs(waux(3,d))/dsqrt(varls)
      return
      end
C======================================================================= 
      subroutine s_invdfe(partial,lp,phif,phi,ndim2) 
*-----------------------------------------------------------------------
*     This subroutine computes the coefficients of an AR model of order 
*     lp given lp partial autocorrelations, using the Durbin algorithm. 
* 
*     Input: 
*          partial: vector containing the lp partial autocorrelations 
*          lp     : length of partial 
*          ndim2  : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                   dimension the auxiliary array phi 
*           
*     Output: 
*          phif  : vector of coefficients of the AR model 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension partial(lp),phif(lp),phi(ndim2,ndim2) 
*-----------------------------------------------------------------------
      phi(1,1)=partial(1) 
      do i=2,lp 
         phi(i,i)=partial(i)            
         do j=1,i-1 
            phi(i,j)=phi(i-1,j)-phi(i-1,i-j)*phi(i,i) 
         enddo 
      enddo 
      do i=1,lp 
         phif(i)=phi(lp,i)  
      enddo 
      return 
      end
C=======================================================================
      subroutine s_jac2fe(fcn,m,n,x,fvec,fjac,ldfjac,iflag,wa,  
     +     idif,isp,nsd,mm,np,nq,n0,indth,npo,npred,xx,yy,  
     +     sigman,sigmau,xy,yhat,uhat,epshat,st,epspred,  
     +     w,cck,auxm,poldif,ndim1,ndim2,work3,nw3,work4,  
     +     nw4,iwork4,niw4,work5,nw5,iwork5,niw5)
*-----------------------------------------------------------------------
*     This subroutine computes a forward-difference approximation  
*     to the m by n jacobian matrix associated with a specified  
*     problem of m functions in n variables.  
*  
*     Arguments
*  
*      fcn is the name of the user-supplied subroutine which  
*          calculates the functions. fcn must be declared  
*          in an external statement in the user calling  
*          program, and should be written as follows.  
*  
*  
*        subroutine fcn(m,n,x,fvec,iflag)  
*        integer m,n,iflag  
*        double precision x(n),fvec(m)  
*        ----------  
*        Calculate the functions at x and  
*        Return this vector in fvec.  
*        ----------  
*        return  
*        end  
*  
*        The value of iflag should not be changed by fcn unless  
*        the user wants to terminate execution of s_jac2fe.  
*        In this case set iflag to a negative integer.  
*  
*      m is a positive integer input variable set to the number  
*        of functions.  
*  
*      n is a positive integer input variable set to the number  
*        of variables. n must not exceed m.  
*  
*      x is an input array of length n.  
*  
*      fvec is an input array of length m which must contain the  
*           functions evaluated at x.  
*  
*      fjac is an output m by n array which contains the  
*           approximation to the jacobian matrix evaluated at x.  
*  
*      ldfjac is a positive integer input variable not less than m  
*             which specifies the leading dimension of the array fjac.  
*  
*      iflag is an integer variable which can be used to terminate  
*            the execution of s_jac2fe. See description of fcn.  
*    
*      wa is a work array of length m.  
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      integer npo(np+nq),iwork4(niw4),iwork5(niw5)
      dimension x(n),fvec(m),fjac(ldfjac,n),xx(m,mm),yy(m),xy(m,mm+1) 
      dimension yhat(m),uhat(m+npred),epshat(m+npred),st(m+npred), 
     +     epspred(m+npred),auxm(m+npred,ndim2),poldif(isp*nsd+idif+1)
      dimension w(m+npred),wa(m),work3(nw3),work4(nw4),work5(nw5)
      data beta11/53.1036811792759202d0/, 
     +     beta22/2.36701536161196047d0/, 
     +     beta33/0.66580545056039918d-12/, 
     +     beta44/1.d-7/
      data dqtol/0.32927225399136d-9/
      data zero/0.d0/
      external fcn
*----------------------------------------------------------------------- 
      do 20 j = 1, n  
         temp = x(j)
         c=dabs(fvec(j)) 
         if(c .eq. zero) go to 200 
         id=int(dlog10(c))
         if(id .gt. 0) id=id+1 
         go to 300 
 200     id=0 
 300     continue 
         a=dabs(x(j)) 
         if(a .eq. zero) go to 400 
         ib=int(dlog10(a))
         if(ib .eq. 0) ib=ib+1 
         go to 500 
 400     ib=0 
 500     continue 
         if(id.gt.-1) then 
            ae=(beta11**ib) * (beta22**id) * beta33 
            delta=ae+beta44*a 
         else 
            delta=dqtol*dmax1(0.1d0,a) 
         end if 
         x(j) = temp + delta
         call fcn(m,n,x,wa,iflag,idif,isp,nsd,mm,np,nq,n0,indth,npo,  
     +        sigman,sigmau,npred,xx,yy,xy,yhat,cck,uhat,epshat,  
     +        st,epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,  
     +        work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5)  
         if (iflag .lt. 0) go to 30  
         x(j) = temp  
         do 10 i = 1, m 
            fjac(i,j) = (wa(i) - fvec(i))/delta  
 10      continue  
 20   continue  
 30   continue  
      return    
      end  
C=======================================================================
      subroutine s_lmdffe(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,  
     +     diag,mode,factor,nprint,info,nfev,fjac,ldfjac,  
     +     ipvt,qtf,wa1,wa2,wa3,wa4,idif,isp,nsd,mm,np,nq,  
     +     n0,indth,npo,npred,xx,yy,sigman,sigmau,xy,yhat,  
     +     uhat,epshat,st,epspred,w,auxm,poldif,cck,  
     +     ndim1,ndim2,work3,nw3,work4,nw4,iwork4,niw4,  
     +     work5,nw5,iwork5,niw5,epsmch,dwarf)  
*----------------------------------------------------------------------- 
*     This subroutine minimizes the sum of the squares of  
*     m nonlinear functions in n variables by a modification of  
*     the Levenberg-Marquardt algorithm. The user must provide a  
*     subroutine which calculates the functions. The jacobian is  
*     then calculated by a forward-difference approximation.  
*  
*     The subroutine statement is  
*  
*     subroutine s_lmdffe(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,  
*                      diag,mode,factor,nprint,info,nfev,fjac,  
*                      ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)  
*  
*     where  
*  
*         fcn is the name of the user-supplied subroutine which  
*             calculates the functions. fcn must be declared  
*             in an external statement in the user calling  
*             program, and should be written as follows.  
*  
*     NOTE: in our program fcn is s_fc11fe and has more arguments.  
*           See the description of s_fnc1fe.  
*        
*         subroutine fcn(m,n,x,fvec,iflag)  
*         integer m,n,iflag  
*         double precision x(n),fvec(m)  
*         ----------  
*         calculate the functions at x and  
*         return this vector in fvec.  
*         ----------  
*         return  
*         end  
*  
*         The value of iflag should not be changed by fcn unless  
*         the user wants to terminate execution of s_lmdffe.  
*         in this case set iflag to a negative integer.  
*  
*       m is a positive integer input variable set to the number  
*         of functions.  
*  
*       n is a positive integer input variable set to the number  
*         of variables. n must not exceed m.  
*  
*       x is an array of length n. on input x must contain  
*         an initial estimate of the solution vector. on output x  
*         contains the final estimate of the solution vector.  
*  
*       fvec is an output array of length m which contains  
*            the functions evaluated at the output x.  
*  
*       ftol is a nonnegative input variable. termination  
*            occurs when both the actual and predicted relative  
*            reductions in the sum of squares are at most ftol.  
*            Therefore, ftol measures the relative error desired  
*            in the sum of squares.  
*  
*       xtol is a nonnegative input variable. termination  
*            occurs when the relative error between two consecutive  
*            iterates is at most xtol. Therefore, xtol measures the  
*            relative error desired in the approximate solution.  
*  
*       gtol is a nonnegative input variable. Termination  
*            occurs when the cosine of the angle between fvec and  
*            any column of the jacobian is at most gtol in absolute  
*            value. Therefore, gtol measures the orthogonality  
*            desired between the function vector and the columns  
*            of the jacobian.  
*  
*       maxfev is a positive integer input variable. Termination  
*              occurs when the number of calls to fcn is at least  
*              maxfev by the end of an iteration.  
*  
*  
*       diag is an array of length n. If mode = 1 (see  
*            below), diag is internally set. If mode = 2, diag  
*            must contain positive entries that serve as  
*            multiplicative scale factors for the variables.  
*  
*       mode is an integer input variable. If mode = 1, the  
*            variables will be scaled internally. If mode = 2,  
*            the scaling is specified by the input diag. Other  
*            values of mode are equivalent to mode = 1.  
*  
*       factor is a positive input variable used in determining the  
*              initial step bound. This bound is set to the product of  
*              factor and the euclidean norm of diag*x if nonzero, or else  
*              to factor itself. In most cases factor should lie in the  
*              interval (.1,100.). 100. is a generally recommended value.  
*  
*       nprint is an integer input variable that enables controlled  
*              printing of iterates if it is positive. In this case,  
*              fcn is called with iflag = 0 at the beginning of the first  
*              iteration and every nprint iterations thereafter and  
*              immediately prior to return, with x and fve* available  
*              for printing. If nprint is not positive, no special calls  
*              of fcn with iflag = 0 are made.  
*  
*       info is an integer output variable. If the user has  
*            terminated execution, info is set to the (negative)  
*            value of iflag. See description of fcn. Otherwise,  
*            info is set as follows.  
*  
*            info = 0  improper input parameters.  
*  
*            info = 1  both actual and predicted relative reductions  
*                   in the sum of squares are at most ftol.  
*  
*            info = 2  relative error between two consecutive iterates  
*                   is at most xtol.  
*  
*            info = 3  conditions for info = 1 and info = 2 both hold.  
*  
*            info = 4  the cosine of the angle between fvec and any  
*                      column of the jacobian is at most gtol in  
*                      absolute value.  
*  
*            info = 5  number of calls to fcn has reached or  
*                      exceeded maxfev.  
*  
*            info = 6  ftol is too small. No further reduction in  
*                      the sum of squares is possible.  
*  
*            info = 7  xtol is too small. No further improvement in  
*                      the approximate solution x is possible.  
*  
*            info = 8  gtol is too small. fvec is orthogonal to the  
*                      columns of the jacobian to machine precision.  
*  
*       nfev is an integer output variable set to the number of  
*            calls to fcn.  
*  
*       fjac is an output m by n array. The upper n by n submatrix  
*            of fjac contains an upper triangular matrix r with  
*            diagonal elements of nonincreasing magnitude such that  
*  
*                t     t           t  
*               p *(ja* *jac)*p = r *r,  
*  
*            where p is a permutation matrix and jac is the final  
*            calculated jacobian. Column j of p is column ipvt(j)  
*            (see below) of the identity matrix. The lower trapezoidal  
*            part of fjac contains information generated during  
*            the computation of r.  
*  
*       ldfjac is a positive integer input variable not less than m  
*              which specifies the leading dimension of the array fjac.  
*  
*       ipvt is an integer output array of length n. ipvt  
*            defines a permutation matrix p such that jac*p = q*r,  
*            where jac is the final calculated jacobian, q is  
*            orthogonal (not stored), and r is upper triangular  
*            with diagonal elements of nonincreasing magnitude.  
*            Column j of p is column ipvt(j) of the identity matrix.  
*  
*       qtf is an output array of length n which contains  
*           the first n elements of the vector (q transpose)*fvec.  
*  
*       wa1, wa2, and wa3 are work arrays of length n.  
*  
*       wa4 is a work array of length m.  
*  
*       epsmch is the machine precision.  
*
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      integer ipvt(n),npo(np+nq),iwork4(niw4),iwork5(niw5)
      dimension x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n)  
      dimension xx(m,mm),yy(m),xy(m,mm+1),poldif(isp*nsd+idif+1)
      dimension yhat(m),uhat(m+npred),epshat(m+npred),st(m+npred), 
     +     epspred(m+npred),auxm(m+npred,ndim2) 
      dimension w(m+npred),wa1(n),wa2(n),wa3(n),wa4(m)
      dimension work3(nw3),work4(nw4),work5(nw5)
      data one,p1,p5,p25,p75,p0001,zero  
     +     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/  
      external fcn
*----------------------------------------------------------------------- 
      info = 0  
      iflag = 0  
      nfev = 0  
*  
*     Check the input parameters for errors.  
*  
c      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m  
c     +     .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero  
c     +     .or. maxfev .le. 0 .or. factor .le. zero) go to 300  
      if (mode .ne. 2) go to 20  
      do 10 j = 1, n  
         if (diag(j) .le. zero) go to 300  
 10   continue  
 20   continue  
*  
*     Evaluate the function at the starting point  
*     and calculate its norm.  
*  
      iflag = 1
      call fcn(m,n,x,fvec,iflag,idif,isp,nsd,mm,np,nq,n0,indth,npo,  
     +     sigman,sigmau,npred,xx,yy,xy,yhat,cck,uhat,epshat,st,  
     +     epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,work4,  
     +     nw4,iwork4,niw4,work5,nw5,iwork5,niw5)  
      nfev = 1  
      if (iflag .lt. 0) go to 300  
      fnorm = s_dnrmfe(m,fvec)  
*  
*     Initialize Levenberg-Marquardt parameter and iteration counter.  
*  
      par = zero  
      iter = 1  
*  
*     Beginning of the outer loop.  
*  
 30   continue  
*  
*     Calculate the jacobian matrix.  
*     
      iflag = 2  
      call s_jac2fe(fcn,m,n,x,fvec,fjac,ldfjac,iflag,wa4,  
     +     idif,isp,nsd,mm,np,nq,n0,indth,npo,npred,xx,yy,  
     +     sigman,sigmau,xy,yhat,uhat,epshat,st,epspred,  
     +     w,cck,auxm,poldif,ndim1,ndim2,work3,nw3,work4,nw4,  
     +     iwork4,niw4,work5,nw5,iwork5,niw5)  
      nfev = nfev + n  
      if (iflag .lt. 0) go to 300  
*  
*     If requested, call fcn to enable printing of iterates.  
*  
      if (nprint .le. 0) go to 40  
      iflag = 0  
      if (mod(iter-1,nprint) .eq. 0) call fcn(m,n,x,fvec,iflag,  
     +     idif,isp,nsd,mm,np,nq,n0,indth,npo,sigman,  
     +     sigmau,npred,xx,yy,xy,yhat,cck,uhat,epshat,st,  
     +     epspred,w,auxm,poldif,ndim1,ndim2,work3,nw3,  
     +     work4,nw4,iwork4,niw4,work5,nw5,iwork5,niw5)  
      if (iflag .lt. 0) go to 300  
 40   continue  
*  
*     Compute the qr factorization of the jacobian.  
*     
      call s_dqrffe(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3,epsmch)  
*     
*     On the first iteration and if mode is 1, scale according  
*     to the norms of the columns of the initial jacobian.  
*     
      if (iter .ne. 1) go to 80  
      if (mode .eq. 2) go to 60  
      do 50 j = 1, n  
         diag(j) = wa2(j)  
         if (wa2(j) .eq. zero) diag(j) = one  
 50   continue  
 60   continue  
*     
*     On the first iteration, calculate the norm of the scaled x  
*     and initialize the step bound delta.  
*     
      do 70 j = 1, n  
         wa3(j) = diag(j)*x(j)  
 70   continue  
      xnorm = s_dnrmfe(n,wa3)  
      delta = factor*xnorm  
      if (delta .eq. zero) delta = factor  
 80   continue  
*  
*     Form (q transpose)*fvec and store the first n components in  
*     qtf.  
*  
      do 90 i = 1, m  
         wa4(i) = fvec(i)  
 90   continue  
      do 130 j = 1, n  
         if (fjac(j,j) .eq. zero) go to 120  
         sum = zero  
         do 100 i = j, m  
            sum = sum + fjac(i,j)*wa4(i)  
 100     continue  
         temp = -sum/fjac(j,j)  
         do 110 i = j, m  
            wa4(i) = wa4(i) + fjac(i,j)*temp  
 110     continue  
 120     continue  
         fjac(j,j) = wa1(j)  
         qtf(j) = wa4(j)  
 130  continue  
*  
*     Compute the norm of the scaled gradient.  
*     
      gnorm = zero  
      if (fnorm .eq. zero) go to 170  
      do 160 j = 1, n  
         l = ipvt(j)  
         if (wa2(l) .eq. zero) go to 150  
         sum = zero  
         do 140 i = 1, j  
            sum = sum + fjac(i,j)*(qtf(i)/fnorm)  
 140     continue  
         gnorm = dmax1(gnorm,dabs(sum/wa2(l)))  
 150     continue  
 160  continue  
 170  continue  
*  
*     Test for convergence of the gradient norm.  
*  
      if (gnorm .le. gtol) info = 4  
      if (info .ne. 0) go to 300  
*     
*     Rescale if necessary.  
*  
      if (mode .eq. 2) go to 190  
      do 180 j = 1, n  
         diag(j) = dmax1(diag(j),wa2(j))  
 180  continue  
 190  continue  
*  
*     Beginning of the inner loop.  
*     
 200  continue  
*  
*     Determine the Levenberg-Marquardt parameter.  
*     
      call s_dlpafe(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,  
     +     wa3,wa4,dwarf)  
*  
*     Store the direction p and x + p. calculate the norm of p.  
*     
      do 210 j = 1, n  
         wa1(j) = -wa1(j)  
         wa2(j) = x(j) + wa1(j)  
         wa3(j) = diag(j)*wa1(j)  
 210  continue  
      pnorm = s_dnrmfe(n,wa3)  
*  
*     On the first iteration, adjust the initial step bound.  
*     
      if (iter .eq. 1) delta = dmin1(delta,pnorm)  
*     
*     Evaluate the function at x + p and calculate its norm.  
*     
      iflag = 1  
      call fcn(m,n,wa2,wa4,iflag,idif,isp,nsd,mm,np,nq,n0,  
     +     indth,npo,sigman,sigmau,npred,xx,yy,xy,yhat,cck,  
     +     uhat,epshat,st,epspred,w,auxm,poldif,ndim1,ndim2,  
     +     work3,nw3,work4,nw4,iwork4,niw4,work5,nw5,iwork5,  
     +     niw5)
      nfev = nfev + 1  
      if (iflag .lt. 0) go to 300  
      fnorm1 = s_dnrmfe(m,wa4)  
*     
*     Compute the scaled actual reduction.  
*     
      actred = -one  
      if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2  
*  
*     Compute the scaled predicted reduction and  
*     the scaled directional derivative.  
*     
      do 230 j = 1, n  
         wa3(j) = zero  
         l = ipvt(j)  
         temp = wa1(l)  
         do 220 i = 1, j  
            wa3(i) = wa3(i) + fjac(i,j)*temp  
 220     continue  
 230  continue  
      temp1 = s_dnrmfe(n,wa3)/fnorm  
      temp2 = (dsqrt(par)*pnorm)/fnorm  
      prered = temp1**2 + temp2**2/p5  
      dirder = -(temp1**2 + temp2**2)  
*     
*     Compute the ratio of the actual to the predicted  
*     reduction.  
*     
      ratio = zero  
      if (prered .ne. zero) ratio = actred/prered  
*     
*     Update the step bound.  
*     
      if (ratio .gt. p25) go to 240  
      if (actred .ge. zero) temp = p5  
      if (actred .lt. zero)  
     +     temp = p5*dirder/(dirder + p5*actred)  
      if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1  
      delta = temp*dmin1(delta,pnorm/p1)  
      par = par/temp  
      go to 260  
 240  continue  
      if (par .ne. zero .and. ratio .lt. p75) go to 250  
      delta = pnorm/p5  
      par = p5*par  
 250  continue  
 260  continue  
*     
*     Test for successful iteration.  
*     
      if (ratio .lt. p0001) go to 290  
*     
*     Successful iteration. Update x, fvec, and their norms.  
*     
      do 270 j = 1, n  
         x(j) = wa2(j)  
         wa2(j) = diag(j)*x(j)  
 270  continue  
      do 280 i = 1, m  
         fvec(i) = wa4(i)  
 280  continue  
      xnorm = s_dnrmfe(n,wa2)  
      fnorm = fnorm1  
      iter = iter + 1  
 290  continue  
*     
*     Tests for convergence.  
*     
      if (dabs(actred) .le. ftol .and. prered .le. ftol  
     +     .and. p5*ratio .le. one) info = 1  
      if (delta .le. xtol*xnorm) info = 2  
      if (dabs(actred) .le. ftol .and. prered .le. ftol  
     +     .and. p5*ratio .le. one .and. info .eq. 2) info = 3  
      if (info .ne. 0) go to 300  
*  
*     Tests for termination and stringent tolerances.  
*     
      if (nfev .ge. maxfev) info = 5  
      if (dabs(actred) .le. epsmch .and. prered .le. epsmch  
     +     .and. p5*ratio .le. one) info = 6  
      if (delta .le. epsmch*xnorm) info = 7  
      if (gnorm .le. epsmch) info = 8  
      if (info .ne. 0) go to 300  
*     
*     End of the inner loop. Repeat if iteration unsuccessful.  
*     
      if (ratio .lt. p0001) go to 200  
*     
*     End of the outer loop.  
*     
      go to 30  
 300  continue  
*  
*     Termination, either normal or user imposed.  
*  
      if (iflag .lt. 0) info = iflag  
      iflag = 0  
      if (nprint .gt. 0) call fcn (m,n,x,fvec,iflag,idif,isp,nsd,mm,np,  
     +     nq,n0,indth,npo,sigman,sigmau,npred,xx,yy,xy,  
     +     yhat,cck,uhat,epshat,st,epspred,w,auxm,poldif,  
     +     ndim2,ndim2,work3,nw3,work4,nw4,iwork4,niw4,  
     +     work5,nw5,iwork5,niw5)  
      return 
      end  
C=======================================================================
      subroutine s_mednfe(z,n,xmed,aux)                                    
*-----------------------------------------------------------------------
*     This routine computes xmed=median(z(i))                          
* 
*     Input: 
*           z : vector containing the input data 
*           n : length of z    
* 
*     Output: 
*           xmed: median of (z(1),....,z(n))  
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension z(n),aux(n) 
*-----------------------------------------------------------------------
      do i=1,n                                             
         aux(i)=z(i)                                                    
      enddo 
      call s_sortfe(aux,n,1)
      i=n/2                                                          
      k=i*2                                                  
      xmed=aux(i+1)                            
      if (k.ge.n) xmed=(xmed+aux(i))/2.d0
      return 
      end
C=======================================================================
      subroutine s_out2fe(x,y,n,m,idif,isp,nsd,k,q,indth,beta,phidif,
     +     theta,thetas,sigmadif,indio,cck,sigfil,critv,nout,indtipo,
     +     t0,wout,lambda,sigmau0,sigmau,yhat,uhat,epshat,st,lscan,w,
     +     auxm,phi,thetapro,pols,poldif,phiau1,phiau2,resid,indcan,
     +     yaux,psia,psils,waux,lambaux,rho,tau,xy,m1,auxil,idim,work2,
     +     idimw2,work3,idimw3,iwork,idimiw,ierror,n0) 
*----------------------------------------------------------------------- 
*     This subroutine detects innovation outliers, additive outliers and 
*     level shifts. The method used is only valid in the case idif (number 
*     of ordinary differences) > 0 
* 
*     input   
*          x            : matrix of independent variables. 
*          y            : input series  
*          n            : number of observations 
*          m            : number of independent variables 
*          idif         : number of ordinary differences 
*          isp          : seasonal period 
*          nsd          : number of seasonal differences 
*          k            : order of the AR model to be used  
*          q            : order of the ordinary MA  model  
*          indth        : 0 - no seasonal moving average component 
*                         1 - seasonal moving average term included 
*          beta         : vector of regression coefficientes         
*          phidif       : vector of coefficients of the AR model 
*          theta        : vector of ordinary MA coefficients 
*          thetas       : seasonal moving average coefficient 
*          sigmadif     : scale of the differenced series 
*          indio        : innovation outlier indicator, 
*                         if indio=1, the subroutine looks for innovation 
*                         outliers, additive outliers and level shifts. 
*                         if indio=0, only looks for ao and ls 
*          cck          : the constant for the psi function used in  
*                         subroutine s_flt2fe 
*          sigfil       : an estimator of the scale parameter used by 
*                         subroutine s_flt2fe 
*          critv        : critical value for detecting outliers 
* 
*   output: 
* 
*          nout         : number of outliers found 
*          t0(i)        : number of the observation where the i-th outlier is 
*                         found. 
*          indtipo(i)   : if indtipo=0, no outlier is detected, 
*                         if indtipo=1, an innovation outlier is found, 
*                         if indtipo=2, an additive outlier is detected, 
*                         if indtipo=3, a level shift is detected  
*          wout(i)      : impact of the i-th outlier  
*          lambda(i)    : wout(i) / (standard error) 
*          sigmau0      : initial value of the scale parameter 
*          sigmau       : final value of the scale parameter 
*          ierror       : error indicator 
*                         if ierror=1 the subroutine stops the search 
*                            for outliers, because many were detected 
*                         if ierror=2 the subroutine stops the search 
*                            for outliers, because it found outliers in 
*                            the same position twice running
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      integer q,qtrue,d
      integer t0(n),indtipo(n),indcan(n),lscan(n),iwork(idimiw)
      dimension x(n,m),y(n),beta(m),phidif(idim),theta(q),wout(n)
      dimension yhat(n),uhat(n),epshat(n),st(n),w(n),auxm(n,idim)
      dimension work2(idimw2),work3(idimw3),xy(n,m1+1)
      dimension phi(idim+1),thetapro(q+isp*indth)
      dimension pols(idim+1),poldif(idim+1),phiau1(idim+1),
     +     phiau2(idim+1),resid(n),yaux(n),psia(n),psils(n),
     +     waux(3,n),auxil(3,n),rho(0:idim),tau(0:q+isp*indth) 
      double precision lambda(n),lambaux(3,n)
      data zero,one/0.d0,1.d0/
*----------------------------------------------------------------------- 
*     We compute the vector phi containing the coefficients of the AR 
*     part of the model. Phi is computed by multiplying the polynomial 
*     phidif, containing the coefficients of the stationary part of the 
*     AR model with the polynomials that represents the ordinary and 
*     seasonal differences). 
*----------------------------------------------------------------------- 
      ierror=0
      ndif = idif+isp*nsd 
      ktrue = k+ndif 
      call s_pindfe(idif,nsd,isp,poldif,ndif) 
      if((idif.eq.0).and.(nsd.eq.0)) then 
         do i=1,ktrue 
            phi(i)=phidif(i) 
         enddo 
      else 
         phiau1(1)=one
         do i=2,k+1 
            phiau1(i) = -phidif(i-1) 
         enddo 
         call s_polyfe(poldif,ndif,phiau1,k,phiau2,ktrue) 
         do i=1,ktrue 
            phi(i)=-phiau2(i+1) 
         enddo 
      endif 
*  
*     We construct the vector thetapro containig the coefficients of the 
*     MA part of the model. Thetapro is computed as the product of the  
*     polynomial theta (the ordinary MA model) and the seasonal theta. 
*     
      if(indth.eq.0) then 
         pols(1)=one
         do i=2,isp+1 
            pols(i)=zero 
         enddo 
      else 
         pols(1)=one
         do i=2,isp 
          pols(i)=zero 
       enddo 
       pols(isp+1)=-thetas 
      endif   
      phiau1(1)=one       
      do i=2,q+1 
         phiau1(i) = -theta(i-1) 
      enddo 
      call s_polyfe(pols,isp*indth,phiau1,q,phiau2,qtrue) 
      do i=1,qtrue 
         thetapro(i)=-phiau2(i+1) 
      enddo 
*     
*     We call subroutine s_sys2fe to compute the vector of the 
*     autocorrelations rho and the vector tau, which are used in the  
*     subroutine s_flt2fe 
*     
      lfin=max0(k+idif+isp*nsd,qtrue+1) 
      lfin=lfin-idif-isp*nsd          
      do i=k+1,lfin 
         phidif(i)=zero 
      enddo 
* 
*     We compute the position in the work vector work2 and iwork for the  
*     vector used in subroutine s_sys2fe 
* 
      maxqtru=q+isp*indth 
      iw2thea=1 
      iw2coef=iw2thea+idim 
      iw2a=iw2coef+idim 
      iw2b=iw2a+(idim+1)*(idim+1) 
      call s_sys2fe(phidif,theta,thetas,lfin,q,isp,indth,rho,tau, 
     +     work2(iw2coef),work2(iw2a),work2(iw2b),  
     +     work2(iw2thea),iwork(1),idim) 
*     
*     We filter the residuals 
* 
*     We compute the position in the work vector work2 for the  
*     vector used in subroutine s_flt2fe 
* 
      maxqtru=q+isp*indth 
      idim1=maxqtru 
      iw2then=1 
      iw2xm=iw2then+idim1 
      iw2xp=iw2xm+idim*idim 
      iw2alfp=iw2xp+idim*idim 
      iw2r=iw2alfp+idim 
      iw2v=iw2r+idim 
      iw2av=iw2v+idim 
      iw2alff=iw2av+idim*idim 
      iw2rv=iw2alff+idim 
      iw2alfx=iw2rv+idim*idim 
      iw2xaux=iw2alfx+idim 
      iw2v1=iw2xaux+idim*idim 
      iw2v2=iw2v1+idim1+1 
      iw2ypur=iw2v2+idim1+1  
      indfil=0 
      iout=0 
      call s_flt2fe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,ktrue,q,
     +     sigfil,indth,n0,tau,sigmadif,indfil,rho,cck,iout,yhat,uhat,
     +     epshat,st,lscan,w,auxm,xy,m1,work2(iw2then),work2(iw2xm),
     +     work2(iw2xp),work2(iw2alfp),work2(iw2r),work2(iw2v),
     +     work2(iw2av),work2(iw2alff),work2(iw2rv),work2(iw2alfx),
     +     work2(iw2xaux),work2(iw2v1),work2(iw2v2),work2(iw2ypur), 
     +     idim,work3,idimw3,iwork,idimiw) 
* 
*     We compute the scale parameter sigmau  
*     
      do i=n0+1,n 
         uhat(i)=uhat(i)/st(i) 
      enddo 
      call s_calsfe(uhat,n,n0,sigmanew,work2(1),work3) 
      sigmau = sigmanew*sigfil 
      sigmau0 = sigmau 
* 
*     We begin the search for outliers. 
*     
      indica=0 
      nout=0 
      do while(indica.eq.0)     
* 
*     We look for the observations with large residuals 
* 
         indfil=0 
         iout=0 
*     
*     First we call s_flt2fe with iout=0 (which means that we apply the 
*     robust filter to all the observations) to compute the residuals 
*  
         call s_flt2fe(x,y,n,m,idif,isp,nsd,phi,beta,theta,thetas,
     +        ktrue,q,sigfil,indth,n0,tau,sigmadif,indfil,rho,cck,iout,
     +        yhat,uhat,epshat,st,lscan,w,auxm,xy,m1,work2(iw2then),
     +        work2(iw2xm),work2(iw2xp),work2(iw2alfp),work2(iw2r),
     +        work2(iw2v),work2(iw2av),work2(iw2alff),work2(iw2rv),
     +        work2(iw2alfx),work2(iw2xaux),work2(iw2v1),work2(iw2v2),
     +        work2(iw2ypur),idim,work3,idimw3,iwork,idimiw) 
         do i=n0+1,n 
            uhat(i)=uhat(i)/st(i) 
         enddo 
         call s_calsfe(uhat,n,n0,sigmanew,work2(1),work3) 
*     
*     We correct the standardized residuals multiplying them by a factor 
*     so that they have scales estimate = 1 
* 
         do i=n0+1,n 
            uhat(i)=uhat(i)/sigmanew 
         enddo 
*
*     We compute xmaxres the maximum of the absolute value of the  
*     standardized residuals 
*     
         xmaxres=dabs(uhat(n0+1)) 
         do i=n0+2,n 
            resst=dabs(uhat(i)) 
            if (resst.gt.xmaxres) xmaxres=resst 
         enddo 
*     
*     We compute ncan = number of standardized residuals which absolute 
*     value are greater than 2 and greater than xmaxrex*0.7.          
*     The indices of this residuals are allocated in the vector indcan.  
*     
         ncan=0 
         do i=n0+1,n 
            resst=dabs(uhat(i)) 
            if ((resst.gt.2.d0).and.(resst.gt.xmaxres*0.7d0)) then 
               ncan=ncan+1 
               indcan(ncan)=i 
            endif 
         enddo 
*     
*     For each observation with standardized residual greater than 2 
*     and greater than xmaxres*0.7, we compute its impact and its  
*     lambda value. 
*     
         do ii=1,ncan 
            d=indcan(ii) 
            call s_impcfe(x,y,m,idif,isp,nsd,phi,beta,theta,thetas,
     +           ktrue,q,sigmau,indth,n0,tau,sigmadif,rho,d,n,resid,
     +           psils,psia,waux,lambaux,cck,sigfil,yhat,uhat,epshat,
     +           st,lscan,w,auxm,xy,m1,idim,work2,idimw2,work3,idimw3,
     +           iwork,idimiw,iw2then,iw2xm,iw2xp,iw2alfp,iw2r,iw2v,
     +           iw2av,iw2alff,iw2rv,iw2alfx,iw2xaux,iw2v1,iw2v2,
     +           iw2ypur) 
         enddo  
*     
*     We calculate the maximum of the  lambda's. The indice of the 
*     maximum is the outlier candidate 
*     
         xmaxl=zero
         lambaux(3,n)=zero 
         i1=1 
         if (indio.eq.0) i1=2 
         do ii=1,ncan 
            d=indcan(ii) 
            do iii=i1,3 
               if (lambaux(iii,d).ge.xmaxl) then 
                  xmaxl=lambaux(iii,d) 
                  xmaxw=waux(iii,d) 
                  ind=d 
                  itipo=iii 
               endif 
            enddo 
         enddo 
*     
*     We look if the outlier candidate is in the same observation that 
*     any of the outliers detected before  
* 
         indrep=0 
         if (nout.ne.0) then 
            do i=1,nout 
               if (ind.eq.t0(i)) then 
                  indrep=i 
               endif 
            enddo 
         endif 
*     
*     We look if the outlier candidate has a lambda value greater than the 
*     critical value. If this is the case and if the candidate has not 
*     been detected before, we have found a new outlier! So we remove the  
*     effect of this outlier from the series. 
* 
         if(xmaxl.ge.critv) then 
            if(indrep.eq.0) then 
               call s_remvfe(itipo,ind,xmaxw,n,ktrue,phi,qtrue, 
     +              thetapro,y,yaux,0,auxil,maxqtru,idim) 
*
*     We compute a new sigmau 
*     
               indfil=0 
               iout=0
               call s_flt2fe(x,yaux,n,m,idif,isp,nsd,phi,beta,theta,
     +              thetas,ktrue,q,sigfil,indth,n0,tau,sigmadif, 
     +              indfil,rho,cck,iout,yhat,uhat,epshat,st,lscan,
     +              w,auxm,xy,m1,work2(iw2then),work2(iw2xm),
     +              work2(iw2xp),work2(iw2alfp),work2(iw2r),work2(iw2v),
     +              work2(iw2av),work2(iw2alff),work2(iw2rv),
     +              work2(iw2alfx),work2(iw2xaux), 
     +              work2(iw2v1),work2(iw2v2),work2(iw2ypur), 
     +              idim,work3,idimw3,iwork,idimiw)   
               do i=n0+1,n 
                  uhat(i)=uhat(i)/st(i) 
               enddo 
               call s_calsfe(uhat,n,n0,sigmanew,work2(1),work3) 
               sigmau = sigmanew*sigfil  
* 
*     As we have found a new outlier, we save the information  
*     about this outlier in the output vectors t0, indtipo, wout, 
*     and lambda and then we continue the search for outliers 
*     
               do i=1,n 
                  y(i)=yaux(i) 
               enddo 
               nout=nout+1 
               if (nout.ge.n/2) then 
                  ierror=1 
                  return 
               endif 
               indtipo(nout)=itipo 
               t0(nout)=ind 
               wout(nout)=xmaxw 
               lambda(nout)=xmaxl  
            else 
* 
*     This is the case in which the outlier candidate is in the same position 
*     that an outlier detected before,  so we rest the correction made before 
*   
               if (indrep.eq.nout) then 
                  ierror=2 
                  return 
               endif 
               call s_remvfe(indtipo(indrep),t0(indrep),wout(indrep), 
     +              n,ktrue,phi,qtrue,thetapro,y,yaux,1,auxil, 
     +              maxqtru,idim) 
*
*     we compute a new sigmau 
*     
               indfil=0 
               iout=0 
               call s_flt2fe(x,yaux,n,m,idif,isp,nsd,phi,beta, 
     +              theta,thetas,ktrue,q,sigfil,indth,n0,tau,sigmadif, 
     +              indfil,rho,cck,iout,yhat,uhat,epshat,st,lscan, 
     +              w,auxm,xy,m1,work2(iw2then),work2(iw2xm),
     +              work2(iw2xp),work2(iw2alfp),work2(iw2r),
     +              work2(iw2v),work2(iw2av),work2(iw2alff),
     +              work2(iw2rv),work2(iw2alfx),work2(iw2xaux), 
     +              work2(iw2v1),work2(iw2v2),work2(iw2ypur), 
     +              idim,work3,idimw3,iwork,idimiw) 
               do i=n0+1,n 
                  uhat(i)=uhat(i)/st(i) 
               enddo 
               call s_calsfe(uhat,n,n0,sigmanew,work2(1),work3) 
               sigmau = sigmanew*sigfil  
               indtipo(indrep)=0 
               t0(indrep)=0 
               do i=1,n 
                  y(i)=yaux(i) 
               enddo 
            endif 
*     
*     This is the case in which the t value of the outlier candidate 
*     is less than the critical value. So, we have finished the search 
*     for outliers 
*     
         else 
            indica=1 
         endif 
      enddo 
      return 
      end
C=======================================================================
      subroutine s_outlfe(x,y,n,m,idif,isp,nsd,k,q,indth,beta,phidif,
     +     theta,thetas,sigmadif,indio,cck,sigfil,critv,nout,indtipo,
     +     t0,wout,lambda,sigmau0,sigmau,idim,work,idimw,iwork,idimiw,
     +     ierror,n0) 
*----------------------------------------------------------------------- 
*     This subroutine assigns locations in the work vector work 
*     to the arrays used in subroutine s_out2fe. Then it calls   
*     subroutine s_out2fe, which detects outliers.  
* 
*         idim: parameter used to compute the dimension of some 
*               arrays. It must be = max0(q+isp+1,k+idif+isp*nsd) 
* 
*     work areas: 
*         work: real vector of dimension idimw 
*            idimw: must be >= n*(18+idim+max(m,1)+1) + 6*idim+6 
*                         + 2*(q+isp)+1 
*                         + 5*idim + 5*idim**2 + 3*(q+isp)+2 + n  
*                         + max(4*idim + 4*idim**2 , n)  
* 
*         iwork: integer work area vector of dimension idimiw 
*            idimiw: must be >= 2*n+idim+1 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      integer q,t0(n),indtipo(n),iwork(idimiw) 
      dimension x(n,m),y(n),beta(m),phidif(idim),theta(idim),wout(n) 
      double precision lambda(n),work(idimw)
*-----------------------------------------------------------------------
      maxqtru=q+isp*indth 
      idim1=maxqtru 
      m1=max0(m,1)
      iw1yhat=1 
      iw1uhat=iw1yhat+n 
      iw1epsh=iw1uhat+n 
      iw1st=iw1epsh+n 
      iw1w=iw1st+n 
      iw1auxm=iw1w+n 
      iw1phi=iw1auxm+n*idim 
      iw1thep=iw1phi+idim+1 
      iw1pols=iw1phi+idim+1 
      iw1podf=iw1pols+idim+1 
      iw1pha1=iw1podf+idim+1 
      iw1pha2=iw1pha1+idim+1 
      iw1resi=iw1pha2+idim+1 
      iw1yaux=iw1resi+n 
      iw1psia=iw1yaux+n 
      iw1psil=iw1psia+n 
      iw1waux=iw1psil+n 
      iw1lamb=iw1waux+3*n 
      iw1rho=iw1lamb+3*n 
      iw1tau=iw1rho+idim+1 
      iw1xy=iw1tau+idim1+1 
      iw1auxi=iw1xy+(m1+1)*n 
      iw1sum=iw1auxi+3*n 
      idimw2=5*idim + 5*idim**2 + 3*(q+isp)+2 + n 
      idimw3=max0(4*idim + 4*idim**2, n)
*-----------------------------------------------------------------------
      call s_out2fe(x,y,n,m,idif,isp,nsd,k,q,indth,beta,phidif,theta,
     +     thetas,sigmadif,indio,cck,sigfil,critv,nout,indtipo,t0, 
     +     wout,lambda,sigmau0,sigmau,work(iw1yhat),work(iw1uhat),
     +     work(iw1epsh),work(iw1st),iwork(1),work(iw1w),work(iw1auxm),
     +     work(iw1phi),work(iw1thep),work(iw1pols),work(iw1podf), 
     +     work(iw1pha1),work(iw1pha2),work(iw1resi),iwork(1+n), 
     +     work(iw1yaux),work(iw1psia),work(iw1psil),work(iw1waux), 
     +     work(iw1lamb),work(iw1rho),work(iw1tau),work(iw1xy),m1, 
     +     work(iw1auxi),idim,work(iw1sum),idimw2, 
     +     work(iw1sum+idimw2),idimw3,iwork(1+2*n),idimiw-2*n,ierror,n0) 
*     
*     We eliminate the repeated outliers 
*     
      nouttru=0  
      do i=1,nout 
         if (t0(i).ne.0) then 
            nouttru=nouttru+1 
            t0(nouttru)=t0(i) 
            indtipo(nouttru)=indtipo(i) 
            wout(nouttru)=wout(i) 
            lambda(nouttru)=lambda(i) 
         endif 
      enddo 
      nout=nouttru 
      sigmau=dsqrt(dble(n-n0)/dble(n-n0-nout))*sigmau
      return
      end 
C=======================================================================
      subroutine s_pindfe(idif,nsd,isp,poldif,ndiff) 
*----------------------------------------------------------------------- 
*     This subroutine computes the composed differences polynomial poldif 
* 
*     Input: 
*           idif: number of ordinary differences 
*           nsd : number of seasonal differences 
*           isp : seasonal period 
* 
*     Output: 
*           poldif: vector containing the coefficients of the composed  
*                   differences polynomial 
*           ndiff : order of poldif 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      dimension poldif(isp*nsd+idif+1) 
      dimension poldr(10),pols(27)  
*----------------------------------------------------------------------- 
      if (idif.eq.0) then 
         poldr(1) = 1.0d0 
      elseif (idif.eq.1) then 
         poldr(1)=1.0d0 
         poldr(2)=-1.0d0 
      else 
         poldr(1)=1.0d0 
         poldr(2)=-2.0d0 
         poldr(3)= 1.0d0 
      endif 
      if (nsd.eq.0) then 
         pols(1)=1.0d0 
      elseif(nsd.eq.1) then 
         pols(1)=1.0d0 
         do i=2,isp 
            pols(i)=0.0d0 
         enddo 
         pols(isp+1)=-1.0d0 
      else 
         pols(1)=1.0d0 
         do i=2,2*isp 
            pols(i)=0.0d0 
         enddo 
         pols(isp+1)=-2.0d0 
         pols(2*isp+1)=1.0d0 
      endif 
      call s_polyfe(poldr,idif,pols,isp*nsd,poldif,ndiff) 
      return 
      end      
C=======================================================================
      subroutine s_polyfe(a,n,b,m,c,nm)                          
*-----------------------------------------------------------------------
*     This subroutine multiplies the polynomials a and b of grades n 
*     and m respectively and it allocates the resulting polynomial  
*     of grade nm (=n+m) in c. 
* 
*     Input: 
*           a : vector containing the coefficients of the first polynomial 
*               a=1+a(2)*x+a(3)*x**2+...+a(n+1)*x**n                           
*           n : grade of a    
*           b : vector containing the coefficients of the second polynomial  
*               b=1+b(2)*x+b(3)*x**2+...+b(m+1)*x**m                           
*           m : grade of b 
* 
*     Output: 
*           c : vector containing the coefficients of the product polynomial 
*               c=1+c(2)*x+c(3)*x**2+...+c(nm+1)*x**(n+m)                      
*           nm: grade of c 
*                                                                     
*     Note: It is assumed that a(n+1) and b(n+1) are not 0 and that  
*           a(1) and b(1) are equal to 1. 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      dimension a(n+1),b(m+1),c(n+m+1)                                 
*-----------------------------------------------------------------------
      nm=n+m                                                         
      c(1)=1.0d0                                                       
      if(nm.lt.1) return                                             
      if (n.eq.0) then                                               
         do i=2,nm+1                                            
            c(i)=b(i)                                                 
         enddo 
         return                                                    
      endif                                                         
      if (m.eq.0) then                                               
         do i=2,nm+1                                            
            c(i)=a(i)                                                 
         enddo 
         return                                                  
      endif                                                   
      do l=2,nm+1                                             
         c(l)=0.0d0
         inic=max0(1,l-m)                                          
         ifin=min0(n+1,l)   
         do i=inic,ifin                                
            c(l)=c(l)+a(i)*b(l+1-i)                                        
         enddo 
      enddo 
      return                                           
      end                     
C=======================================================================
      function s_psiffe(x) 
*-----------------------------------------------------------------------
*     Derivative of s_rhoffe. 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
*-----------------------------------------------------------------------
      g1= -1.944d0 
      g2=  1.728d0 
      g3=  -.312d0 
      g4=   .016d0 
      ax=dabs(x) 
      if (ax.gt.3.0d0) then 
         s_psiffe=0.0d0 
      else if(ax.le.2.d0) then 
         s_psiffe=x 
c     else if(ax.gt.2.d0.and.ax.le.3.d0) then 
      else
         s_psiffe=x*(g4*(x**6)+g3*(x**4)+g2*(x**2)+g1) 
      endif 
      return
      end
C=======================================================================
      subroutine s_rbqrfe(nob,nvar,x,y,eps,itmax,beta,ss,ncons,n,x1,e1, 
     +     f1,e2,s3,res,vv,ww,cof,s2,beta1,y1,auxx,nvarn,nvars)
*-----------------------------------------------------------------------   
*     This subroutine computes a robust estimate for linear regression 
*     using a robust qr method. 
*     
*     Input: 
*             nob      : number of observations  
*             nvar     : number of independent variables including constant 
*             x        : matrix of independent variables. It includes a 
*                        column with 1's when there is an intercept  
*             y        : input series of length nob 
*             eps      : convergence criterium 
*             itmax    : maximum number of iterations 
*             ncons    : 1 if there is an intercept; 0 in other case  
*             n        : number of rows of X in the calling program        
* 
*     Output: 
*             beta     : vector of dimension nvar containing the estimates  
*                        of the regression coefficients  
*             ss       : scale computed by the subroutine 
* 
*     Remark: In case there is an intercept, it must be the first variable. 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z) 
      integer nvarn(nvar),nvars(nvar)
      dimension x(n,nvar),y(nob),beta(nvar),cof(nvar),auxx(3*nob) 
      dimension x1(nob*nvar),y1(nob),res(nob),e1(nvar)
      dimension f1(nvar),e2(nvar),s2(nvar),s3(nvar)
      dimension vv(nvar,nvar),ww(nvar,nvar),beta1(nvar)
      data zero/0.d0/

c initialize variables defined within if
      ss1 = 0.d0

*-----------------------------------------------------------------------   
* 
*     Vectors beta, y1 and res are initialized. 
*     
      iter=0
      do i=1,nvar 
         beta(i)=zero 
      enddo 
      do i=1,nob 
         y1(i)=y(i) 
         res(i)=y(i) 
      enddo 
*     
*     We start the iterative process. In each iteration we try to modify  
*     the regression coefficient of each of the nvar variables. 
* 
 595  iter=iter+1 
* 
*     We copy matrix x on vector x1. 
* 
      iii=0 
      do i=1,nvar 
         do ii=1,nob 
            iii=iii+1 
            x1(iii)=x(ii,i) 
         enddo 
      enddo 
*     
*     We initialize the vector nvarn (initial value : nvarn(i)=i), which 
*     gives the order in which the variables will enter. 
*     
      do i=1,nvar 
         nvarn(i)=i 
      enddo 
*     
*     We begin the process of variable selection. 
*     
      do i=1,nvar 
* 
*     nvar2 is the number of variables to consider.  
* 
         nvar2=nvar-i+1 
* 
*     If i=1 (first variable to consider) and there is intercept,   
*     this varible enters directly, without selection (the intercept is  
*     the first variable 1). 
*     
         if(i.eq.1.and.ncons.eq.1) then 
* 
*     We regress y on the intercept and get the residuals r, the 
*     regression coefficient e1(1) and the residual scale s3(1). 
* 
            call s_vesrfe(x1(1),res,nob,e1(1),s3(1), 
     +           auxx(1),auxx(nob+1),auxx(2*nob+1)) 
* 
*     nn is the chosen variable. 
* 
            nn=1 
*     
*     s2 is the residual scale after the first i variables have entered. 
* 
            s2(i)=s3(1) 
         else 
* 
*     If i is different of 1 or there is not intercept we are going to 
*     regress the residuals on each of the nvar2 variables which have 
*     not been selected yet, and select the one with smallest residual  
*     scale. e1(j) is the regression coefficient and s3(j) the residual  
*     scale. The vector nvarn keeps the variables that have not been  
*     selected yet. 
* 
            do j=1,nvar2  
               call s_vesrfe(x1((nvarn(j)-1)*nob+1),res,nob,e1(j),
     +              s3(j),auxx(1),auxx(nob+1),auxx(2*nob+1))  
            enddo 
* 
*     Here we select the variable witht smallest s3. The order of the  
*     selected variable is nn and the scale after we have selected the  
*     first i variables is s2(i).  
* 
            nn=1  
            s2(i)=s3(1) 
            do k=1,nvar2 
               if (s3(k).lt.s2(i)) then 
                  nn=k 
                  s2(i)=s3(k) 
               endif 
            enddo 
         endif 
* 
*     In res we keep the residual of the selected variable nn: Y-X(nn)*beta. 
*     and in e2(i) the corresponding regression coefficient. 
* 
         do k=1,nob 
            res(k)=res(k)-e1(nn)*x1((nvarn(nn)-1)*nob+k) 
         enddo 
         e2(i)=e1(nn) 
* 
*     vv(i,k) denotes the coeeficient of the i-th variable with respect 
*     to the k-th selected variable, and ww(i,k) the coefficient of the  
*     i-th selected variable in the k-th selected variable. 
*     nvars(i) is the i-th selected variable.  
*     
         if (i.ne.1) then 
            i3=i-1 
            do k=1,i3 
               ww(i,k)=vv(nvarn(nn),k) 
            enddo 
         endif 
         nvars(i)=nvarn(nn) 
* 
*     We redefine the vector  nvarn of the variables which have not 
*     entered yet. 
* 
         if (nn.ne.nvar2) then 
            nvar3=nvar2-1 
            do k=nn,nvar3 
               nvarn(k)=nvarn(k+1) 
            enddo 
         endif 
         nvar1=i 
         if (i.ne.nvar) then 
            nvar3=nvar-i 
* 
*     We start the process of orthogonalizing  the  variables that have 
*     not entered yet with respect to the newly selected variable.  
*     To do this we regress each of these variables with respect to the newly 
*     selected variable and replace them in X1 by the resids. The regression 
*     coefficient  of the nvarn(k) variable on the nvars(i) variable is 
*     kept on vv(nvarn(k),i). 
* 
            do k=1, nvar3 
               call s_vesrfe(x1((nvars(i)-1)*nob+1),
     +              x1((nvarn(k)-1)*nob+1),nob,cof(k),ff,auxx(1),
     +              auxx(nob+1),auxx(2*nob+1)) 
               do l=1,nob 
                  x1((nvarn(k)-1)*nob+l)=x1((nvarn(k)-1)*nob+l)-cof(k)* 
     +                 x1((nvars(i)-1)*nob+l) 
               enddo 
*     
*     In  column i of vv we keep the coefficients cof(k),k=1,nvar-i. 
*  
               vv(nvarn(k),i)=cof(k) 
            enddo 
         endif 
      enddo 
* 
*     Now we choose the optimal number of variables minimizing the  
*     residual scale. 
* 
      nh=1 
      st=s2(1) 
      do i=2,nvar 
         if (s2(i).le.st) then 
            nh=i 
            st=s2(i) 
         endif  
      enddo 
*     
*     We compare this scale with the residual scale of the previous iteration. 
* 
      if (iter.gt.1) then 
         if (st.gt.ss1) then 
            do i=1,nob 
               res(i)=y1(i) 
            enddo 
         endif 
      endif 
      nvar1=nh 
      ss1=st 
      if (nvar1.ne.nvar) then 
         nvar3=nvar1+1 
*     
*     We make 0 the coefficients of the variables which have not entered. 
* 
         do i=nvar3,nvar 
            e2(i)=zero 
         enddo 
      endif 
* 
*     We retransform the variables. 
* 
      do k=1,nvar 
         f1(nvars(k))=e2(k) 
      enddo 
      nvar3=nvar1-1 
      do k=1,nvar3 
         j=nvar1-k 
         do l=1,k 
            f1(nvars(j))=f1(nvars(j))-f1(nvars(j+l))*ww(j+l,j) 
         enddo 
      enddo 
*     
*     We compute the residuals again. 
*     
      do i=1,nob 
         res(i)=y1(i) 
         do j=1,nvar1 
            res(i)=res(i)-f1(nvars(j))*x(i,nvars(j)) 
         enddo 
      enddo 
      ss1=s2(nvar1) 
      do i=1,nvar 
         beta1(i)=beta(i)+f1(i) 
      enddo 
* 
*     If there is perfect fitting the iterations stop. 
*     
      if (ss1.lt..000000001d0) then 
         do i=1,nvar 
            beta(i)=beta1(i) 
         enddo 
         return 
      endif 
* 
*     We compare the scales of the last two iterations to decide if the 
*     iterations should stop. 
* 
      if (iter.eq.1) go to 455 
      aux=dabs(ss1-ss)/dabs(ss) 
      if (aux.lt..02d0) then 
         do i=1,nvar 
            beta(i)=beta1(i) 
         enddo 
         return 
      endif 
* 
*     We compare the regression coefficients to decide if the iterations 
*     should stop.  
*     
      nst=0 
      do j=1,nvar 
         if (dabs(beta(j)).ge..0001d0) then 
            aux=(beta1(j)-beta(j))/beta(j) 
         else 
            aux=beta1(j)-beta(j) 
         endif 
         aux1=dabs(beta1(j)) 
         aux3=dabs(beta(j)) 
         aux1=dmax1(aux1,aux3) 
         aux=dmin1(aux,aux1) 
         if(aux.gt.eps) nst=1 
      enddo
      aux=(ss1-ss)/ss 
      aux=dabs(aux) 
      if (nst.eq.0.and.aux.lt.eps) then 
         do i=1,nvar  
            beta(i)=beta1(i) 
         enddo 
         return 
      endif 
 455  continue 
* 
*     We control if the maximum number of iterations was attained. 
*     
      if (iter.lt.itmax) then 
         do i=1,nvar 
            beta(i)=beta1(i) 
         enddo 
         ss=ss1 
         do i=1,nob 
            y1(i)=res(i) 
         enddo 
         go to 595 
      else 
         do i=1,nob 
            res(i)=y1(i) 
         enddo 
      endif 
      return 
      end 
C=======================================================================
      subroutine s_rcorfe(uhat,st,n,n0,zcor,aux)
*-----------------------------------------------------------------------  
*     This subroutine computes zcor. The acf and pacf of zcor are 
*     robust acf and pcf of uhat. 
*  
*     Input:  
*          uhat    : input series 
*          st      : scales of uhat 
*          n       : number of oservations 
*          n0      : number of initial discarded observations 
* 
*     Output: 
*          zcor    : series whose acf and pacf are robust acf and pacf 
*                    respectively of uhat 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)  
      dimension uhat(n),st(n),zcor(n),aux(2*n)
      data tpf/2.5d0/
*-----------------------------------------------------------------------
      do i=n0+1,n 
         zcor(i)=uhat(i)/st(i) 
      enddo 
*     
*     We compute an M-scale of zcor. 
* 
      call s_calsfe(zcor,n,n0,sigmm,aux(1),aux(n+1)) 
* 
*     We found pseudo observations using a Huber type psi-function. 
*     
      do i=n0+1,n 
         z=uhat(i)/(st(i)*sigmm) 
         if (z .ge. tpf) then 
            zcor(i-n0)=tpf
         elseif (z .le. -tpf) then 
            zcor(i-n0)=-tpf 
         else 
            zcor(i-n0)=z 
         endif 
      enddo
      return 
      end 
C=======================================================================
      subroutine s_regafe(x,y,n,m,idif,isp,nsd,ip,indar,ipfin,iqfin, 
     +     indth,interc,kopt,phiopt,theta,thetas,betaopt,sigmau,bcov,
     +     zcor,zcor1,sigmadif,cck,sigfil,xy,yhat,uhat,epshat,st,
     +     epspred,npred,tauef,infnew,ndim1,ndim2,work1,nw1,iwork1,
     +     niw1,work2,nw2,iwork2,niw2,utol,maxfev,epsmch,dwarf,n0) 
*-----------------------------------------------------------------------  
*     This subroutine distributes the work vectors work1 and iwork1, 
*     of lengths nw1 and niw1 respectively, between different arrays  
*     and vectors to be used by s_rramfe, and then calls it. 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),phiopt(ndim2),theta(ndim1), 
     +     betaopt(m),bcov(m,m),zcor(n),zcor1(n),xy(n,m+1) 
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred), 
     +     epspred(n+npred) 
      dimension work1(nw1),work2(nw2) 
      dimension iwork1(niw1),iwork2(niw2) 
*-----------------------------------------------------------------------  
      n1=1  
      n2=ndim2*ndim2  
      n3=n2+m 
      n4=n3+(ip+1)*(ip+1) 
      n5=n4+ndim2+1 
      n6=n5+ip+1 
      n7=n6+m 
      n8=n7+n*ndim1 
      n9=n8+ndim1 
      n10=n9+n 
      n11=n10+ndim1 
      n12=n11+ndim1 
      n13=n12+ndim1 
      n14=n13+n 
      n15=n14+ndim1 
      n16=n15+ndim2+1 
      n17=n16+ndim2+1 
      n18=n17+ndim2+1 
      n19=n18+ndim1 
      n20=n19+ndim1 
      n21=n20+ndim2+1 
      n22=n21+idif+isp*nsd+1 
      n23=n22+ip 
      n24=n23+ndim2+1 
      n25=n24+ndim2+1 
      n26=n25+iqfin 
      n27=n26+4*n 
      n28=n27+n 
      n29=n28+n+npred
      call s_rramfe(x,y,n,m,ip,idif,isp,nsd,indth,interc,kopt,phiopt, 
     +     theta,thetas,betaopt,indar,ipfin,iqfin, sigmau, 
     +     bcov,zcor,zcor1,sigmadif,cck,sigfil,xy,yhat,uhat, 
     +     epshat,st,epspred,npred,tauef,infnew,ndim1,ndim2, 
     +     work1(n1),work1(n2+1),work1(n3+1),work1(n4+1), 
     +     work1(n5+1),work1(n6+1),work1(n7+1),work1(n8+1), 
     +     work1(n9+1),work1(n10+1),work1(n11+1),work1(n12+1), 
     +     work1(n13+1),work1(n14+1),work1(n15+1),work1(n16+1), 
     +     work1(n17+1),work1(n18+1),work1(n19+1),work1(n20+1), 
     +     work1(n21+1),work1(n22+1),work1(n23+1),work1(n24+1), 
     +     work1(n25+1),work1(n26+1),work1(n27+1),work1(n28+1), 
     +     work1(n29+1),iwork1(1),iwork1(ndim1+1),work2,nw2, 
     +     iwork2,niw2,utol,maxfev,epsmch,dwarf,n0) 
      return 
      end
C=======================================================================
      subroutine s_remvfe(itipo,t0,w,n,ktrue,phi,qtrue,
     +     thetapro,y,yaux,ind,auxil,maxqtru,idim)
*----------------------------------------------------------------------- 
*     If ind=1, this subroutine s_remvfe the effect of an outlier in
*     the observation t0.
*     If ind=0, the subroutine rest the correction made before. 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      dimension y(n),yaux(n),phi(idim),thetapro(maxqtru),auxil(3,n)
      integer qtrue,t,t0
      data zero/0.d0/
*----------------------------------------------------------------------- 
*     
*     For innovation outliers 
*     
      if (itipo.eq.1) then
         do t=1,t0-1
            auxil(1,t)=zero
         enddo
         auxil(1,t0)=w
         do t=t0+1,n
            auxil(1,t)=zero
            do i=1,ktrue
               auxil(1,t)=auxil(1,t)+phi(i)*auxil(1,t-i)
            enddo
            if ((t-t0).le.qtrue) then
               auxil(1,t)=auxil(1,t)-thetapro(t-t0)*w    
            endif       
         enddo
         if (ind.eq.0) then
            do t=1,n
               yaux(t)=y(t)-auxil(1,t)
            enddo
         else
            do t=1,n
               yaux(t)=y(t)+auxil(1,t)
            enddo
         endif
      endif
*     
*     For additive outliers
*     
      if (itipo.eq.2)  then
         do t=1,n
            yaux(t)=y(t)
         enddo
         if (ind.eq.0) then
            yaux(t0)=y(t0)-w
         else
            yaux(t0)=y(t0)+w
         endif 
      endif
*     
*     For level shifts
*     
      if (itipo.eq.3) then
         do t=1,t0-1
            yaux(t)=y(t)
         enddo
         if (ind.eq.0) then
            do t=t0,n
               yaux(t)=y(t)-w
            enddo
         else
            do t=t0,n
               yaux(t)=y(t)+w
            enddo
         endif
      endif
      return
      end
C=======================================================================
      function s_rhoffe(x) 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
*-----------------------------------------------------------------------  
      g1=-1.944d0/2.d0 
      g2= 1.728d0/4.d0 
      g3= -.312d0/6.d0 
      g4=  .016d0/8.d0 
      ax=dabs(x) 
      if (ax.lt.2.d0) then 
         s_rhoffe=x**2/2.d0 
      elseif (ax.gt.3.d0) then 
         s_rhoffe=3.25d0 
      else 
         s_rhoffe=g1*x**2+g2*x**4+g3*x**6+g4*x**8+1.792d0 
      endif 
      return 
      end 
C=======================================================================
      subroutine s_rinife(z,n,np,phi,theta,isp,rho,tau,alfafi,idif,nsd,
     +     p,nq,n0,sigmau,sigma,ndim2,alfa,rhoin,covu,uhat,uuhat, 
     +     w,rhom,b,ipiv)
*----------------------------------------------------------------------- 
*      This subroutine calculates the initial value of the filtered 
*      state vector and its covariance matrix. 
* 
*     Input: 
*           z      : input series (y-x'*beta) 
*           n      : length of the input series 
*           np     : order of the non stationary autoregresive operator 
*           phi    : vector of of the non stationary autoregressive 
*                    coefficients. 
*           theta  : vector of moving average coefficients 
*           isp    : seasonal period 
*           rho    : autocorrelation function of the differenced series 
*           tau    : tau(i) = e(u(t-i+1)*z(t))/e(z(t)^2)  
*           idif   : number of ordinary differences (0 or 1) 
*           nsd    : number of seasonal differences (0 or 1) 
*           nq     : order of the moving average operator 
*           sigmau : scale of the filtered residuals 
*           sigma  : scale of the differenced series 
*           ndim2  : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required 
*                    to dimension several arrays. 
* 
*     Output: 
*           alfafi : initial filtered state vector 
*           p      : covariance matrix of alfafi 
*           n0     : initial observation 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension z(n),phi(ndim2+1),theta(ndim2),rho(0:ndim2), 
     +     tau(0:ndim2),alfafi(ndim2),p(ndim2,ndim2) 
      dimension alfa(ndim2),rhoin(ndim2,ndim2),covu(ndim2,ndim2), 
     +     uhat(ndim2),uuhat(ndim2),w(ndim2),rhom(ndim2,ndim2), 
     +     b(ndim2,ndim2) 
      integer ipiv(ndim2) 
      data zero/0.d0/
*-----------------------------------------------------------------------
      n0=np 
      np1=np 
      do i=1,np 
         alfa(i)=z(np-i+1) 
         w(i)=z(np-i+1) 
      enddo 
* 
*     We construct the initial value of the stationary ARMA(np1,q) process  
*     to use it as predictors of the initial value of the state vector.   
*     We require to do this only if nq>0. 
* 
      if (nq.ne.0) then 
         if (idif.ne.0) then 
            np1=np1-1 
            do i=1,np1 
               w(i)=w(i)-w(i+1) 
            enddo 
         endif 
         if (idif.gt.1) then 
            np1=np1-1 
            do i=1,np1 
               w(i)=w(i)-w(i+1) 
            enddo 
         endif
         if (nsd.ne.0) then 
            np1=np1-isp 
            do i=1,np1 
               w(i)=w(i)-w(i+isp) 
            enddo 
         endif 
      endif 
      if (nsd.gt.1) then 
         np1=np1-isp 
         do i=1,np1 
            w(i)=w(i)-w(i+isp) 
         enddo 
      endif 
      do i=1,np 
         do j=1,np 
            p(i,j)=zero 
         enddo 
      enddo 
*
*     If nq=0 the initial state vector alfafi has as components the  
*     vector (z(np),...,z(1)) or similarly (alfa(1),...,alfa(np)) 
* 
      if (nq.eq.0) then 
         do i=1,np 
            alfafi(i)=alfa(i) 
         enddo 
* 
*     If nq>0, the initial state vector alfafi=(alfafi(1),..., alfafi(np)) 
*     is given by: 
*     alfafi(1)=z(np), and for i>1:  
*     alfafi(i)=phi(i)*z(np-1)+phi(i+1))*z(np-2)+...+phi(np)*z(i-1)) 
*               -theta(i-1)*u(np)-...-theta(q)*u(np+i-1-nq) 
* 
*     We calculate now the first part of the state vector 
*     alfafi(i)=phi(i)*z(np-1)+phi(i+1))*z(np-2)+...+phi(np)*z(i-1)) 
* 
      else 
         alfafi(1)=alfa(1) 
         if (np.gt.1) then 
            do i=2,np 
               alfafi(i)=zero 
               do j=i,np 
                  alfafi(i)=alfafi(i)+phi(j)*alfa(j-i+2) 
               enddo 
            enddo 
         endif 
* 
*     We compute now the covariance matrix of alfafi when there is not 
*     stationary AR operator. In this case alfafi is given by 
*     -theta(i-1)*u(p)-...-theta(q)*u(p+i-1-q) 
* 
         if (np1.eq.0) then 
            do i=2,nq+1 
               do j=2,nq+1 
                  kk=max0(i,j) 
                  do k=1, nq-kk+2 
                     p(i,j)=p(i,j)+theta(i-2+k)*theta(j-2+k) 
                  enddo 
                  p(i,j)=sigmau**2*p(i,j) 
               enddo 
            enddo 
         else 
* 
*     We compute the correlation matrix of the vector w=(w(p),...,w(1)). 
* 
            do i=1, np1
               do j=1,np1 
                  ii=i-j 
                  if (ii.lt.0) ii=-ii 
                  rhom(i,j)=rho(ii) 
               enddo 
            enddo 
            call s_rinvfe(rhom,rhoin,np1,ndim2,b,ipiv) 
* 
*     We calculate the linear predictor of u=(u(p),...,u(1)) based on  
*     w. It is denoted by uuhat.  
*     The formula for uuhat is  
*     uuhat=cov(u,w)*inv(cov(w,w))*w=(cov(u,w)/E(w^2))*rhoin 
*     where rhoin=inv(cov(w,w)/E(w^2)) and tau=cov(u,w)/E(w^2) 
*     (cov is the covariance matrix). 
* 
*     First we obtain uhat=inv(cov(w,w)/E(w^2))*w. 
* 
            do i=1,np1 
               uhat(i)=zero 
               do j=1,np1 
                  uhat(i)=uhat(i)+rhoin(i,j)*w(j) 
               enddo 
            enddo 
* 
*     Now we compute uuhat=(cov(u,w)/E(w^2))*uhat. 
* 
            do i=1,nq 
               uuhat(i)=zero
               iii=min0(i,np1) 
               do j=1,iii 
                  uuhat(i)=uuhat(i)+uhat(j)*tau(i-j) 
               enddo 
            enddo 
* 
*     We calculate the predictor of the state vector 
*     alfafi(i)=alfafi(i)-theta(i-1)*u(np)-...-theta(q)*u(np+i-1-nq) 
*     which is alfafi(i) -theta(i-1)*uuhat(np)-...-theta(q)*uuuhat(np+i-1-nq). 
* 
            do i=2,nq+1 
               do j=i-1,nq 
                  alfafi(i)=alfafi(i)-theta(j)*uuhat(j-i+2) 
               enddo 
            enddo 
*  
*     We calculate covu, the covariance matrix of u-uuhat 
* E(w^2)((I/E(w^2))- (cov(u,w)/E(w^2)))*inv(cov(w,w)/E(w^2))*(cov(w,u)/E(w^2)) 
* 
            do i=1,nq 
               do j=1,nq 
                  if (i.eq.j) then 
                     covu(i,j)=tau(0) 
                  else   
                     covu(i,j)=zero
                  endif 
                  iii=min0(i,np1) 
                  do l=1,iii 
                     jjj=min0(j,np1) 
                     do k=1,jjj 
                        covu(i,j)=covu(i,j)-tau(i-l)*tau(j-k)*rhoin(l,k) 
                     enddo 
                  enddo 
                  covu(i,j)=covu(i,j)*(sigma**2) 
               enddo 
            enddo 
* 
*     We calculate the covariance matrix p of the prediction error of alfafi. 
*     The component i of this vector is: 
*     -theta(i-1)*(u(p)-uuhat(p))-...-theta(q)*(u(p+i-1-q)-uuhat(p+i-1-q)) 
*     and therefore depends of covu. 
* 
            do i=2,nq+1 
               do j=2,nq+1 
                  do k=1,nq-i+2 
                     do l=1,nq-j+2 
                        p(i,j) = p(i,j) +
     +                       theta(i+k-2)*theta(j+l-2)*covu(k,l) 
                     enddo 
                  enddo 
               enddo 
            enddo 
         endif 
      endif 
      return 
      end 
C=======================================================================
      subroutine s_rinvfe(a,c,n,ndim,b,ipiv) 
*-----------------------------------------------------------------------
*     This subroutine computes the inverse of a symmetric matrix a. 
* 
*     Input: 
*           a   : matrix of dimension n x n 
*           n   : number of rows and columns of a 
*           ndim: number of rows of a in the calling program (ndim>=n) 
* 
*     Output: 
*           c: inverse of a 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension a(ndim,ndim),b(ndim,ndim),c(ndim,ndim) 
      integer ipiv(ndim) 
*-----------------------------------------------------------------------
* 
*     We copy matrix a into c and define b=I (identity matrix). 
* 
      do i=1,n 
         do j=1,n 
            c(i,j)=a(i,j) 
            if (i.eq.j) then 
               b(i,j)=1.0d0 
            else 
               b(i,j)=0.0d0 
            endif 
         enddo 
      enddo 
* 
*     We call s_gesvfe in order to solve the system c . x = b 
* 
      call s_gesvfe(n,n,c,ndim,ipiv,b,ndim,ierror) 
* 
*     We transfer the inverse to c. 
* 
      do i=1,n 
         do j=1,n 
            c(i,j)=b(i,j) 
         enddo
      enddo 
      return 
      end 
C=======================================================================
      subroutine s_rqr1fe(nob,nvar,x,y,eps,itmax,beta,ss,ncons,n,work3, 
     +     nw3,iwork3,niw3) 
*-----------------------------------------------------------------------    
*     This subroutine distributes the work vectors work3 and iwork3, 
*     of lengths nw3 and niw3 respectively, between different arrays  
*     and vectors to be used by s_rbqrfe, and then calls it. 
*-----------------------------------------------------------------------    
      implicit double precision (a-h,o-z)    
      dimension x(n,nvar),y(nob),beta(nvar),work3(nw3) 
      integer iwork3(niw3) 
*-----------------------------------------------------------------------    
      n2=nob*nvar 
      n3=n2+nvar 
      n4=n3+nvar 
      n5=n4+nvar 
      n6=n5+nvar 
      n7=n6+nob 
      n8=n7+nvar*nvar 
      n9=n8+nvar*nvar 
      n10=n9+nvar 
      n11=n10+nvar 
      n12=n11+nvar 
      n13=n12+nob 
      call s_rbqrfe(nob,nvar,x,y,eps,itmax,beta,ss,ncons,n,work3(1), 
     +     work3(n2+1),work3(n3+1),work3(n4+1),work3(n5+1), 
     +     work3(n6+1),work3(n7+1),work3(n8+1),work3(n9+1), 
     +     work3(n10+1),work3(n11+1),work3(n12+1),work3(n13+1), 
     +     iwork3(1),iwork3(nvar+1)) 
      return 
      end 
C=======================================================================
      subroutine s_rramfe(x,y,n,m,ip,idif,isp,nsd,indth,interc,kopt, 
     +     phiopt,thetaopt,thetas,bopt,indar,ipfin,iqfin,sigmau,bcov,
     +     zcor,zcor1,sigmadif,cck,sigfil,xy,yhat,uhat,epshat,st,
     +     epspred,npred,tauef,infnew,ndim1,ndim2,phi,beta,phidif,tau,
     +     phiaux,beta0,fjac,qtf,f,wa1,wa2,wa3,wa4,diag,rh,para,para1,
     +     parold,par,phiaux2,poldif,phidif1,phiau2,rho,theta,aux1,
     +     ypure,w,auxm,ipvt,npo,work2,nw2,iwork2,niw2,utol,maxfev,
     +     epsmch,dwarf,n0) 
*-----------------------------------------------------------------------
*     Input: 
*               x       : matrix of independent variables 
*               y       : vector of the input series 
*               n       : number of observations 
*               m       : number of independent variables 
*               ip      : maximum order of the AR model to estimate 
*               idif    : number of ordinary differences 
*               isp     : seasonal period 
*               nsd     : number of seasonal differences 
*               indth   : 1 seasonal MA parameter thetas included 
*               interc  : 1 intercept included 
*                         0 intercept not included 
*               indar   : 1 an AR(p) model will be fit automatically 
*                         0 the ARMA(p,q) model is given  
*               ipfin   : order of the final ordinary AR polynomial 
*               iqfin   : order of the final ordinary MA polynomial 
*               npred   : number of predicted regression errors 
*               ndim1   : ip+iqfin+m+indth+1, required to dimension several  
*                         arrays 
*               ndim2   : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required  
*                         to dimension several arrays 
* 
*     Output: 
*               kopt    : order of the best AR model 
*               phiopt  : estimates of the coefficients of the final AR model 
*               thetaopt: estimates of the coefficients of the final MA model 
*               thetas  : estimate of the seasonal MA coefficient 
*               bopt    : estimates of the regression model 
*               sigmau  : corrected robust estimate of the innovation scale 
*               bcov    : covariance matrix of the regression coefficients 
*                         estimates 
*               zcor    : the acf and pcf of zcor are robust acf and 
*                         pcf of the innovations 
*               zcor1   : the acf and pcf of zcor1 are robust acf and 
*                         pcf of the differenced regression residuals    
*               sigmadif: scale of the differenced regression residuals 
*               cck     : bandwidth of the robust filter 
*               sigfil  : innovation scale 
*               n0      : number of initial non filtered observasions 
*               xy      : n*(m+1) matrix . The first m columns contain 
*                         the filtered x's and the m+1 column contains the 
*                         filtered y. All the columns are obtained using 
*                         the robust filter corresponding to the estimated 
*                         ARIMA model 
*               yhat    : cleaned y-series  
*               uhat    : innovations 
*               epshat  : regression residuals 
*               epspred : predicted regression residuals 
*               st      : scales of the prediction error 
*               tauef   : inverse of the efficiency of the tau estimates 
*               utol    : We make ftol=xtol=gtol=utol in the optimizer
*                         soubroutine s_lmdffe. 
*               maxfev  : Maximum number of calls to the function which 
*                         calculates the pseudo likelihood in s_lmdffe 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension x(n,m),y(n),xy(n,m+1),epspred(n+npred) 
      dimension phiopt(ndim2),thetaopt(iqfin+1),bopt(m),bcov(m,m)
      dimension yhat(n),uhat(n+npred),epshat(n+npred),st(n+npred)
      dimension phi(ndim2,ndim2),beta(m),beta0(m),phidif(ip+1,ip+1), 
     +     tau(0:ndim2),phiaux(ip+1)
      dimension fjac(n,ndim1),qtf(ndim1),f(n),zcor(n),zcor1(n)
      dimension w(n+npred),wa1(ndim1),wa2(ndim1),wa3(ndim1),wa4(n)
      dimension diag(ndim1),rh(ndim2+1)
      dimension para(ndim2+1),para1(ndim2+1),parold(ndim1), 
     +     par(ndim1),phiaux2(ndim2+1),poldif(idif+isp*nsd+1) 
      dimension Phidif1(ip),phiau2(ndim2+1),rho(0:ndim2),theta(iqfin) 
      dimension aux1(4*n),ypure(n),auxm(n+npred,ndim2),work2(nw2) 
      integer ipvt(ndim1),npo(ip+iqfin),iwork2(niw2),iiwork(300)
      data zero/0.d0/
      external s_fc11fe
*-----------------------------------------------------------------------
*
*   We define the size of the different work vectors. 
*
      nw3=9*n+2*m+ndim1*(7+n)+9*ndim2+5*ndim2*ndim2+2*m*m+8+4*ip+ 
     +     idif+isp*nsd+iqfin+(n+npred)*ndim2+npred 
      niw3=max0(ip+iqfin+ndim1,ndim2,m)+1 
      nw4=n*m+7*m+2*m*m+7*n+9*ndim2+4+(n+npred)*ndim2+ 
     +     5*ndim2*ndim2+iqfin+isp+ip 
      niw4=max0(ndim2,2*m)+1 
      nw5=5*ndim2*ndim2+6*ndim2+iqfin+isp+2 
      niw5=ndim2 
      nw6=4*ndim2*ndim2+4*ndim2 
      niw6=ndim2 
* 
*     Constants of s_lmdffe. 
* 
      ftol=zero
      xtol=utol
      gtol=zero
      mode=1 
      nprint=0 
      factor=100.d0 
      ldfjac=n 
*      
*     We compute the differences polynomial poldif whose order is ndiff. 
* 
      call s_pindfe(idif,nsd,isp,poldif,ndiff) 
* 
*     We call s_grd1fe, through s_gd11fe, which estimates beta and phi  
*     for a linear model with errors modelled as AR(1,ndiff). 
*
      call s_gd11fe(x,y,interc,n,m,idif,isp,nsd,iqfin,phi,beta, 
     +     phidif(1,1),sigini,sigmau,akai,cck,parold,ssqold,xy, 
     +     yhat,uhat,epshat,st,epspred,npred,poldif,w,auxm, 
     +     ndim1,ndim2,work2(1),nw3,iwork2(1),niw3,work2(nw3+1), 
     +     nw4,iwork2(niw3+1),niw4,work2(nw3+nw4+1),nw5, 
     +     iwork2(niw3+niw4+1),niw5,work2(nw3+nw4+nw5+1),nw6, 
     +     iwork2(niw3+niw4+niw5+1),niw6,utol,maxfev,epsmch,dwarf)  
      thetas=zero 
* 
*     The optimal k is initialized to 1. 
* 
      kopt=1 
* 
*     The optimal Akaike criterion is initialized to the one corresponding  
*     to the ARI(1,ndiff) found in s_grd1fe. 
* 
      akaiop=akai 
* 
*     The innovation scale is initialized to the one corresponding to the  
*     ARI(1,ndiff) found in s_grd1fe. 
* 
      sopt=sigmau 
* 
*     The optimal filter bandwidth is initialized to the one  
*     corresponding to the ARI(1,ndiff) found in s_grd1fe. 
* 
      cckopt=cck 
* 
*     The regression coefficients are initialized to the ones corresponding  
*     to the ARI(1,ndiff) found in s_grd1fe. 
* 
      do i=1,m 
         beta0(i)=beta(i) 
      enddo 
* 
*     We call s_grdkfe,through s_gdk1fe, which estimates beta and phi for  
*     a linear model with errors modelled as an ARI(k,udif). We do this for 
*     k=2,...,ip, and use an Akaike's criterium to select the "optimal" k. 
*
      do k=2,ip 
         call s_gdk1fe(x,y,n,m,ip,idif,isp,nsd,iqfin,phi,beta,thetas, 
     +        phidif,k,sigini,sigmau,akai,cck,parold,ssqold, 
     +        xy,yhat,uhat,epshat,st,epspred,w,auxm,npred, 
     +        poldif,ndim1,ndim2,work2(1),nw3,iwork2(1),niw3, 
     +        work2(nw3+1),nw4,iwork2(niw3+1),niw4, 
     +        work2(nw3+nw4+1),nw5,iwork2(niw3+niw4+1),niw5, 
     +        work2(nw3+nw4+nw5+1),nw6,iwork2(niw3+niw4+niw5+1),
     +        niw6,utol,maxfev,epsmch,dwarf) 
         if (indar.eq.0.and.iqfin.eq.0) then 
* 
*     In the case of a pure autoregressive model of a fixed order, the 
*     optimal k is the current one. 
* 
            sopt=sigmau 
            kopt=k 
            akaiop=akai 
            cckopt=cck 
            do i=1,m 
               beta0(i)=beta(i) 
            enddo 
* 
*     We determine if the current k is optimal according to the  
*     Akaike criterion. 
* 
         elseif (akai.lt.akaiop) then 
            sopt=sigmau 
            kopt=k 
            akaiop=akai 
            cckopt=cck 
            do i=1,m 
               beta0(i)=beta(i) 
            enddo 
         endif 
      enddo 
      do i=1,m 
         beta(i)=beta0(i) 
      enddo 
      do i=1,kopt 
         phidif1(i)=phidif(kopt,i) 
      enddo  
* 
*     s_gdthfe gives the parameters of a model with the same AR order as 
*     the one found in s_grdkfe and with a MA seasonal parameter. 
* 
      if (indar.eq.1) then 
         ipfin=kopt 
         iqfin=0 
      endif 
      if (kopt.lt.(ipfin+iqfin)) then 
         do i=kopt+1,(ipfin+iqfin) 
            phidif1(i)=zero 
         enddo 
         kopt=ipfin+iqfin 
      endif
      if (indth.eq.0.and.iqfin.eq.0) then 
         do i=1,ipfin 
            phiopt(i)=phidif1(i) 
         enddo 
         do i=1,m 
            bopt(i)=beta(i) 
         enddo 
      elseif (indth.eq.0.and.iqfin.gt.0) then 
         n1=ndim2 
         n2=n1+(ndim2+1)*(ndim2+1) 
         n3=n2+(ndim2+1)
* 
*     We are going to compute the initial ARMA model when there is no 
*     seasonal MA parameter using s_sys3fe. To use s_sys3fe it is 
*     required that s_sys2fe computes the autocorrelations rho's and the  
*     coefficients tau 
* 
         call s_sys2fe(phidif1,theta,thetas,kopt,0,isp,indth,rho,tau, 
     +        work2(1),work2(n1+1),work2(n2+1), 
     +        work2(n3+1),iwork2,ndim2) 
* 
*     Given the rho's and the tau's, s_sys3fe computes the parameters of  
*     the ARMA(p,q) model. 
* 
         call s_sys3fe(phidif1,kopt,phiaux,theta,thetas,ipfin, 
     +        iqfin,isp,indth,rho,tau,ndim2,work2(1), 
     +        work2(ndim2*ndim2+1),iwork2)
*     
*     We check if the AR operator computed by s_sys3fe is stationary. 
*     Otherwise is modified. 
* 
         if (ipfin.gt.0) then 
            call s_yulefe(phiaux,rho,ipfin,work2(1),iiwork,ndim2) 
            call s_durbfe(rho,ipfin,phiau2,ier,work2(1),ndim2) 
            if (ier.eq.1) then 
               do i=1,ipfin 
                  if (phiau2(i).gt. 1.d0) phiau2(i)= .98d0 
                  if (phiau2(i).lt.-1.d0) phiau2(i)=-.98d0 
               enddo 
               call s_invdfe(phiau2,ipfin,phiaux,work2(1),ndim2)   
            endif 
         endif 
* 
*     We check if the MA operator computed by s_sys3fe is invertible. 
*     Otherwise is modified. 
* 
         if (iqfin.gt.0) then 
            call s_yulefe(theta,rho,iqfin,work2(1),iiwork,ndim2) 
            call s_durbfe(rho,iqfin,phiau2,ier,work2(1),ndim2) 
            if (ier.eq.1) then 
               do i=1,iqfin 
                  if (phiau2(i).gt. 1.d0) phiau2(i)= .98d0 
                  if (phiau2(i).lt.-1.d0) phiau2(i)=-.98d0 
               enddo 
               call s_invdfe(phiau2,iqfin,theta,work2(1),ndim2)   
            endif 
         endif 
         do i=1,ipfin 
            phiopt(i)=phiaux(i) 
         enddo 
         do i=1,iqfin 
            thetaopt(i)=theta(i) 
         enddo 
         do i=1,m 
            bopt(i)=beta(i) 
         enddo 
      else 
* 
*     s_gdt1fe compute the initial model when there is a seasonal MA parameter 
* 
         call s_gdt1fe(x,y,n,m,ip,idif,isp,nsd,ipfin,iqfin,beta,thetas, 
     +        phidif1,kopt,sigini,sopt,tau,n0,cckopt,xy,yhat, 
     +        uhat,epshat,st,epspred,w,auxm,npred,poldif, 
     +        theta,ndim2,work2(1),nw3,work2(nw3+1),nw4, 
     +        iwork2(niw3+1),niw4,work2(nw3+nw4+1),nw5, 
     +        iwork2(niw3+niw4+1),niw5,work2(nw3+nw4+nw5+1), 
     +        nw6,iwork2(niw3+niw4+niw5+1),niw6)
         do i=1,ipfin 
            phiopt(i)=phidif1(i) 
         enddo 
         do i=1,m 
            bopt(i)=beta(i) 
         enddo 
         do i=1,iqfin 
            thetaopt(i)=theta(i) 
         enddo 
         thetasop=thetas 
      endif 
* 
*     We will compute a global tau-estimator for the final model. 
*  
*     We set npar to the number of parameters of the optimization  
*     subroutine s_lmdffe and the filter bandwith, par(npar), to 1.0. 
* 
**      npar=ipfin+iqfin+m+indth+1 
      npar=ipfin+iqfin+m+indth
      ccknew=1.d0
* 
*     We call subroutine s_tranfe which computes partial autocorrelations 
*     and transforms them in variables that may take any real value. 
* 
      call s_trasfe(phiopt,thetaopt,thetasop,bopt,ndim2,ipfin,iqfin,
     +     indth,m,para,par,ndim1,rh,work2,nw2,iwork2,niw2,0) 
* 
*     The vector npo contains the powers p of B with non zero coefficients  
*     in the AR and MA operators. In this version 0's are not allowed. 
* 
      do i=1,ipfin+iqfin 
         npo(i)=0 
      enddo 
      do i=1,ipfin 
         npo(i)=i 
      enddo 
      do i=1,iqfin 
         npo(ipfin+i)=i 
      enddo 
      sigmadif=s_xmadfe(x,y,bopt,m,n,aux1(1),aux1(n+1),aux1(2*n+1), 
     +     poldif,ndiff) 
      itte=0 
* 
*     The minimization process starts. We make a loop with three different 
*     starting values of the robust filter bandwidth. 
* 
      do while (itte.le.4) 
         itte=itte+1 
         do i=1,npar-m-1 
            u=1.d0 
            if(par(i).lt.zero) u=-1.d0 
            par(i)=dmin1(dabs(par(i)), 6.d0)*u 
         enddo 
         call s_lmdffe(s_fc11fe,n,npar,par,f,ftol,xtol,gtol,maxfev, 
     +        diag,mode,factor,nprint,info,nfev,fjac,ldfjac, 
     +        ipvt,qtf,wa1,wa2,wa3,wa4,idif,isp,nsd,m, 
     +        ipfin,iqfin,n0,indth,npo,0,x,y,sigman,sigmanew, 
     +        xy,yhat,uhat,epshat,st,epspred,w,auxm,poldif, 
     +        ccknew,ndim1,ndim2,work2(1),nw3,work2(nw3+1),
     +        nw4,iwork2(1),niw4,work2(nw3+nw4+1),nw5, 
     +        iwork2(niw4+1),niw5,epsmch,dwarf)
*  
*     In ssq we compute the value of the goal function for the optimal  
*     solution. 
* 
         ssq=zero
         do jjj=n0+1,n 
            ssq=ssq+f(jjj)*f(jjj) 
         enddo 
*  
*     We call s_fc11fe to compute two estimates of the innovation scale sigman 
*     and sigmanew. The difference between these two estimates is that  
*     sigmanew is computed only using the variance of the differenced  
*     stationary process and sigman corrects this value using the filtered  
*     residuals uhat. 
* 
         call s_fc11fe(n,npar,par,f,iflag,idif,isp,nsd,m,ipfin,iqfin, 
     +        n0,indth,npo,sigman,sigmanew,0,x,y,xy,yhat, 
     +        ccknew,uhat,epshat,st,epspred,w,auxm,poldif,
     +        ndim1,ndim2,work2(1),nw3,work2(nw3+nw4+1),nw5, 
     +        iwork2(niw3+niw4+1),niw5,work2(nw3+nw4+nw5+1),
     +        nw6,iwork2(niw3+niw4+niw5+1),niw6) 
* 
*     We compare the solution with new bandwidth filter with the 
*     current optimal and decide if we should change it. 
*          
         if (itte.eq.1 .or. ssqold.gt.ssq) then 
            ssqold=ssq 
            cck=ccknew
            do i=1,npar 
               parold(i)=par(i) 
            enddo 
            sigmau=sigman 
            sigfil=sigmanew 
            infnew=info 
         endif 
* 
*     We set the new value of the bandwith filter. 
* 
         vv=(0.8d0**itte)*(sopt/sigmanew) 
         ccknew=dmin1(vv,1.d0) 
         if(itte.eq.3) ccknew=1000.d0
         if(itte.eq.4) ccknew=1.d0
      enddo 
*     
*     We antitransform the npar parameters from par, using subroutine 
*     s_tranfe. 
* 
      do i=1,npar 
         par(i)=parold(i) 
      enddo 
      do i=1,ndim2 
         phiopt(i)=zero
      enddo
      call s_tranfe(par,ndim1,ndim2,ipfin,iqfin,indth,0,para,para1, 
     +     work2,phiopt,thetaopt,thetas,bopt) 
* 
*     We allocate the new estimates in row (kopt+idif) of matrix phi. 
*
      para1(1)=1.d0 
      do i=2,kopt+1 
         para1(i)=-phiopt(i-1) 
      enddo 
* 
*     We recover the regression parameters from the optimal solution. 
* 
      do i=1,m 
         bopt(i)=par(ipfin+iqfin+indth+i) 
      enddo
* 
*     In case of automatic selection, the subroutine s_corsfe computes 
*     the series zcor1. The acf of this series is a robust acf of the 
*     differenced residuals. 
* 
      if (indar.eq.1) then 
         call s_corsfe(x,bopt,n,m,idif,isp,nsd,zcor1,yhat,work2(1), 
     +        work2(n+1),work2(2*n+1),work2(nw3+1))
      endif 
* 
*     We determine the final innovation scale in case that the ARIMA 
*     model do not have regular MA parameters. 
* 
      do i=1,ip+1 
         phiaux(i)=zero 
      enddo  
      phiaux(1)=1.0d0 
      do i=1,ipfin 
         phiaux(npo(i)+1)=-phiopt(npo(i)) 
      enddo  
* 
*     We filter with the new values of phi and theta. 
* 
*     We compute the non-stationary autoregressive polynomial operator. 
*  
      call s_polyfe(phiaux,ipfin,poldif,ndiff,phiaux2,k)
      do i=1,k 
         phiaux2(i)=-phiaux2(i+1) 
      enddo 
      do i=k+1,ndim2+1 
         phiaux2(i)=zero 
      enddo 
* 
*     We compute the scale of the stationary AR component of the 
*     regression model, that is, of the regression errors after differencing.  
* 
      sigmadif=s_xmadfe(x,y,bopt,m,n,aux1(1),aux1(n+1),aux1(2*n+1), 
     +     poldif,ndiff) 
* 
*     lfin is the smallest order of the stationary AR operator compatible 
*     with the restriction of the state space represantation of the ARIMA 
*     model. 
* 
      lfin=max0(ipfin+idif+isp*nsd,iqfin+isp*indth+1) 
      lfin=lfin-idif-isp*nsd          
* 
*     s_sys2fe computes the correlation coefficients rho of the 
*     stationary process and the covariances tau between the stationary  
*     ARMA process and the innovations. 
* 
      n1=ndim2 
      n2=n1+(ndim2+1)*(ndim2+1) 
      n3=n2+(ndim2+1) 
      call s_sys2fe(phiopt,thetaopt,thetas,lfin,iqfin,isp,indth, 
     +     rho,tau,work2(1),work2(n1+1),work2(n2+1), 
     +     work2(n3+1),iwork2,ndim2)
* 
*     We filter three times, the first time with indfil=0, the second 
*     time with indfil=1 and the third time with indfil=0. 
*     The first time we filter robustly the residuals getting the values  
*     uhat. The weights used to obtain uhat are stored in w. 
*     The second time we use a linearized version of the filter using  
*     the weights w to filter the x's and the y's obtaining a matrix xy 
*     which will be used to estimate the covariances of the regression  
*     coefficients. The third time we filter again to get uhat.  
*     
      if (m.gt.0) then 
         call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux2,bopt,thetaopt, 
     +        thetas,k,iqfin,sigfil,indth,n0,tau,sigmadif,0, 
     +        rho,cck,0,ypure,xy,yhat,uhat,epshat,st,epspred, 
     +        w,auxm,ndim2,work2(1),nw3,work2(nw3+nw4+nw5+1), 
     +        nw6,iwork2(niw3+niw4+niw5+1),niw6) 
         call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux2,bopt,thetaopt, 
     +        thetas,k,iqfin,sigfil,indth,n0,tau,sigmadif,1, 
     +        rho,cck,0,ypure,xy,yhat,uhat,epshat,st,epspred, 
     +        w,auxm,ndim2,work2(1),nw3,work2(nw3+nw4+nw5+1), 
     +        nw6,iwork2(niw3+niw4+niw5+1),niw6) 
      endif 
      call s_flt1fe(x,y,n,m,idif,isp,nsd,phiaux2,bopt,thetaopt,thetas, 
     +     k,iqfin,sigfil,indth,n0,tau,sigmadif,0,rho,cck, 
     +     npred,ypure,xy,yhat,uhat,epshat,st,epspred,w,auxm, 
     +     ndim2,work2(1),nw3,work2(nw3+nw4+nw5+1),nw6, 
     +     iwork2(niw3+niw4+niw5+1),niw6) 
* 
*     Subroutine s_bdesfe computes the covariance matrix of the regression  
*     coefficients. 
*  
      n1=m*m 
      n2=n1+m*m 
      n3=n2+n 
      n4=n3+n 
      call s_bdesfe(.405d0,n,m,n0,bcov,xy,uhat,st,tauef,work2(1), 
     +     work2(n1+1),work2(n2+1),work2(n3+1),work2(n4+1),iwork2) 
* 
*     The subroutine s_rcorfe computes the series zcor. The acf of  
*     this series is a robust acf of the filtered innovations uhat. 
* 
      call s_rcorfe(uhat,st,n,n0,zcor,work2(1)) 
      return
      end 
C=======================================================================
      subroutine s_sortfe(a,n,iswitch)
*-----------------------------------------------------------------------
*     This subroutine rearranges vector 'a'. 
* 
*     Input: 
*           a       : on input contains the original vector 
*           n       : length of a 
*           iswitch : If iswitch > 0, the subroutine rearranges the  
*                     vector a in increasing way. 
*                     If iswitch <= 0, it rearranges the vector in  
*                     decreasing way.              
* 
*     Output: 
*           a       : on output a contains the ordered vector  
* 
*     Algorithm due to  D. L. Shell (C.A.C.M. July 1959,page 30).              
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)  
      dimension a(n)
*-----------------------------------------------------------------------
      if (n.le.1) go to 999                                           
      m=1                                                             
 106  m=m+m                                                           
      if(m.le.n) go to 106                                            
      m=m-1                                                           
 994  m=m/2                                                          
      if (m.eq.0) go to 999                                           
      kk=n-m                                                          
      j=1                                                             
 992  i=j                                                             
 996  im=i+m                                                          
      if(iswitch .gt. 0) then
      if (a(i).gt.a(im)) go to 110                                  
      else
      if(a(i).lt.a(im)) go to 110                                   
      end if
 995  j=j+1                                                          
      if(j.gt.kk) go to 994                                          
      go to 992                                                      
 110  temp=a(i)                                                      
      a(i)=a(im)                                                     
      a(im)=temp                                                     
      i=i-m                                                         
      if (i.lt.1) go to 995                                          
      go to 996                                                      
 999  return                                                        
      end
C=======================================================================
      subroutine s_sys1fe(phidif1,ip,theta,rho,k,isp,phithe,tau,ierror, 
     +     phiaux,a,b,ipiv,ndim2) 
*----------------------------------------------------------------------- 
*     Consider a time series x(t) following a SARMA model (k,0)*(0,1). 
*     Let phithe be the vector with the AR parameters and theta the 
*     seasonal MA parameter. 
*     This  subroutine computes the AR parameters phithe, given   
*     the autocorrelations rho(i),i=...,k and the parameter theta.  
*     The vector phidif1 is used only as initial value in an iterative  
*     algorithm.   
* 
*     Input: 
*           phidif1     : vector with AR coefficients which are only 
*                         used as initial values of an iterative algorithm 
*           theta       : seasonal MA coefficient 
*           ip          : order of an older adjusted AR model 
*           k           : order of AR model 
*           rho         : vector containing the serial autocorrelations 
*                         rho(i),i=1,k 
*           isp         : seasonal period 
* 
*         
*     Output: 
*           phithe      : vector containing the AR coefficients  
*                         phithe(i), i=1,...,k 
*           ierror      : if ierror=1 there is a numeric problem in 
*                                      subroutine s_gesvfe 
*                         if ierror=2, the iterative algorithm did 
*                                      not converge in maxiter steps. 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      dimension phidif1(ip),rho(0:ndim2),phithe(ndim2),tau(0:ndim2)  
      dimension phiaux(ndim2),a(ip,ip),b(ip) 
      integer ipiv(ip) 
      data tol/0.0000001d0/maxiter/100/
      data zero/0.d0/
*-----------------------------------------------------------------------
      do i=1,ndim2 
         phiaux(i)=zero 
      enddo 
      do i=1,k 
         phiaux(i)=phidif1(i) 
      enddo 
* 
*     We begin the iterative algorithm 
* 
      iter=1 
      xnorm=tol+1.0d0 
      do while(xnorm.gt.tol.and.iter.le.maxiter) 
*   
*     We compute tau(j)=cov(x(t),u(t-j))/var(x(t)), where u(t) is  
*     the innovation series. 
* 
*     We use that tau(0) =  
* 
*   (1 + sum(Phi(i)^2) - 2(sum(phi(i)rho(i))) + 2(sum(phi(i)phi(j)rho(j-i))   
*  -------------------------------------------------------------------------- 
*                            1   +   theta^2        
* 
*   and that 
* 
*         tau(j) = phi(1)*tau(j-1)+...+phi(j)*tau(0),   j=1,...,isp-1 
* 
*     First we compute tau(0) 
* 
         sum1=zero
         sum2=zero 
         sum3=zero 
         do i=1,k 
            sum1=phiaux(i)*phiaux(i)+sum1 
            sum2=phiaux(i)*rho(i)+sum2 
         enddo 
         do i=1,k-1 
            do j=i+1,k 
               sum3=phiaux(i)*phiaux(j)*rho(j-i)+sum3 
            enddo 
         enddo 
         tau(0)=(1.0d0+sum1-2.0d0*sum2+2.0d0*sum3)/(1.0d0+theta*theta) 
* 
*     Here we compute tau(j), for j=1,...,isp-1 
* 
         do j=1,isp-1 
            tau(j)=zero 
            do i=1,j 
               tau(j)=tau(j)+phiaux(i)*tau(j-i) 
            enddo 
         enddo 
* 
*     In order to compute phi(i),i=1,...,k, we have to solve 
*     a linear system or the form  a.x = b. 
*     We call vector x=(phi(1),...,phi(k)) 
* 
*     The elements of the matrix a are: 
*       a(i,i) = 1 ,                       i=1,...,k  
*       a(i,j) = rho(|i-j|) ,              i=1,...,k, j=1,...,k, (i ne j) 
*  
*     The elements of the vector b are: 
*       b(i) = rho(i) + theta*tau(isp-i),  i=1,...,k 
*                    (in this last formula, remember that tau(j)=0 if j<0) 
* 
* 
*      Here we compute the elements of the matrix a and the vector b.  
* 
         do i=1,k 
            do j=1,k 
               if (i.eq.j) then 
                  a(i,j)=1.d0 
               else 
                  lag=abs(j-i) 
                  a(i,j)=rho(lag) 
               endif 
            enddo 
         enddo 
         do i=1,k 
            if (isp.ge.i) then 
               b(i)=rho(i)+theta*tau(isp-i) 
            else 
               b(i)=rho(i) 
            endif 
         enddo 
* 
*     We call subroutine s_gesvfe in order to solve the linear system ax=b. 
* 
         call s_gesvfe(k,1,a,ip,ipiv,b,ip,ierror) 
* 
* 
*     We transfer the results of subroutine s_gesvfe to the vector phithe. 
* 
         do i=1,k 
            phithe(i)=b(i) 
         enddo 
*     
*     We compare the vector phithe with the one obtained in the previous 
*     step of the iterative algotithm 
* 
         xnorm=zero 
         do i=1,k 
            if (phiaux(i).gt.1.d-10) then 
               dife=dabs((phiaux(i)-phithe(i))/phiaux(i))               
            else 
               dife=dabs(phiaux(i)-phithe(i))               
            endif 
            if (dife.gt.xnorm) xnorm=dife 
         enddo 
         do i=1,k 
            phiaux(i)=phithe(i) 
         enddo 
         iter=iter+1 
* 
*       End of the iteration step 
* 
      enddo 
*     
*     If the algorithm does not converges in maxiter iterations, we 
*     set ierror=2 
*           
      if (iter.ge.maxiter) ierror=2    
      return 
      end                                               
C=======================================================================
      subroutine s_sys2fe(phidif,theta,thetas,ip,iq,isp,indth,rho,tau, 
     +     coef,a,b,thetaux,ipiv,ndim2)
*-----------------------------------------------------------------------
*     This subroutine computes the autocorrelations rho(i), i=0,...,ip   
*     and the coefficients tau(j), j=0,...,nqaux, for a time series  
*     following an ARMA model (ip,iq) or a SARMA model (ip,iq)*(0,1). 
*     We call rho(i)=correl(x(t),x(t-i)) and  
*     tau(j)=cov(x(t),u(t-j))/var(x(t)), 
*     where x(t) is the time series and u(t) is the innovation series. 
*     We call nqaux the order of the product of the ordinary and seasonal 
*     MA parts of the model. 
*     The autoregressive and moving average parameters of the model are 
*     inputs of the subroutine. 
* 
*     Input: 
*           phidif      : vector of the AR coefficients 
*           theta       : vector of the ordinary MA coefficients 
*           thetas      : seasonal MA coefficient 
*           ip          : order of the AR model 
*           iq          : order of the ordinary MA model 
*           isp         : seasonal period 
*           indth       : 0 - no seasonal moving average component 
*                         1 - seasonal moving average term included 
*         
*     Output: 
*           rho         : vector containing the autocorrelations   
*                         rho(i), i=0,...,ip 
*           tau         : vector containing the coefficients  
*                         tau(j), j=0,...,iq 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension phidif(ndim2),theta(iq),rho(0:ip),tau(0:iq+indth*isp) 
      dimension coef(ndim2),a(ndim2+1,ndim2+1),b(ndim2+1),thetaux(ndim2) 
      dimension ipiv(ndim2+1)
      data zero,one/0.d0,1.d0/
*-----------------------------------------------------------------------
*     We construct vector thetaux containing the coefficients of the  
*     product of the ordinary and seasonal MA polynomials, and 
*     compute nqaux, the order of the product MA polynomial. 
*-----------------------------------------------------------------------
      nqaux=iq+indth*isp 
      do i=ip+1,nqaux 
         phidif(i)=zero 
      enddo 
      do i=1,ndim2 
         thetaux(i)=zero
      enddo 
      if (iq.gt.0.and.indth.eq.1) then 
         do i=1,iq 
            thetaux(i)=theta(i) 
         enddo 
         thetaux(isp)=thetas 
         do i=1,iq 
            thetaux(isp+i)=-thetas*theta(i) 
         enddo 
      endif 
      if (iq.eq.0.and.indth.eq.1) then 
         thetaux(isp)=thetas 
      endif     
      if (iq.gt.0.and.indth.eq.0) then 
         do i=1,iq 
            thetaux(i)=theta(i) 
         enddo 
      endif 
      rho(0)=one
*     
*     Let coef(j) be the coefficient such that tau(j)=coef(j)*tau(0).  
*     It may be proved that coef(1)=phidif(1)-thetaux(1) and 
*     coef(j)=phidif(j)-theaux(j) + phidif(1)*coef(j-1)+....
*            +phidif(j-1)*coef(1). 
*
*     Here we compute coef(j), j=1,...,nqaux 
*     
      do j=1,nqaux 
         coef(j)=phidif(j)-thetaux(j) 
         do i=1,j-1 
            coef(j)=coef(j)+phidif(i)*coef(j-i) 
         enddo 
      enddo 
*     
*     In order to compute rho(i),i=1,...,ip and tau(0) we have to solve 
*     a linear system or the form  a.x = b. 
*     We call vector x=(rho(1),...,rho(ip),tau(0)) 
* 
*     The elements of the matrix a are: 
*       a(i,i) = -1+phidif(2*i),      i=1,...,ip  
*       a(i,j) = phidif(i-j)+phidif(i+j),  i=1,...,ip, j=1,...,ip, (i ne j)  
*       a(j,ip+1) = -thetaaux(j)-thetaux(j+1)*coef(1)- 
*                    ...-thetaux(nqaux)*coef(nqaux-j) 
*       a(ip+1,j)=phidif(j) 
*       a(p+1,p+1)=1-thetaux(1)*coef(1)-...-thetaux(nqaux)*coef(nqaux) 
* 
*     The elements of the vector b are: 
*        b(i) = -phidif(i), i=1,...ip 
*        b(ip+1) = 1 
* 
*     Here we compute the elements of the matrix a and the vector b.  
*     
      do i=1,ip 
         if (2*i.le.ip) then  
            a(i,i)=-one+phidif(2*i) 
         else 
            a(i,i)=-one 
         endif 
      enddo 
      do i=1,ip 
         do j=1,ip 
            if (i.ne.j) then 
               a(i,j)=zero 
               if ((1.le.i-j).and.(i-j.le.ip)) a(i,j)=a(i,j)+phidif(i-j)  
               if ((1.le.j+i).and.(j+i.le.ip)) a(i,j)=a(i,j)+phidif(j+i)  
            endif 
         enddo 
      enddo 
      do j=1,ip 
         a(ip+1,j)=phidif(j) 
      enddo 
      do j=1,ip 
         a(j,ip+1)=-thetaux(j) 
            do i=j+1,nqaux 
               a(j,ip+1)=a(j,ip+1)-thetaux(i)*coef(i-j) 
            enddo 
      enddo 
      a(ip+1,ip+1)=one 
      do i=1,nqaux 
         a(ip+1,ip+1)=a(ip+1,ip+1)-thetaux(i)*coef(i) 
      enddo 
      do i=1,ip 
         b(i)=-phidif(i) 
      enddo 
      b(ip+1)=one 
      call s_gesvfe(ip+1,1,a,ndim2+1,ipiv,b,ndim2+1,ierror) 
* 
*     We transfer the results of s_gesvfe to rho and tau(0). 
*
      if (ierror.ne.0) ierror=1
      do i=1,ip 
         rho(i)=b(i) 
      enddo 
      tau(0)=b(ip+1) 
* 
*     We compute tau(j),j=1,...,nqaux, using tau(0) and coef. 
* 
      do j=1,nqaux 
         tau(j)=coef(j)*tau(0) 
      enddo 
      return 
      end
C=======================================================================
      subroutine s_sys3fe(phiar,npdif,phidif,theta,thetas,ip,iq,isp, 
     +     indth,rho,tau,ndim2,a,b,ipiv) 
*-----------------------------------------------------------------------
*     This subroutine computes the parameters of a SARMA model (p,q) x (0,1), 
*     given the estimated parameters of an aproximating SAR model  
*     (p*,0) x (0,1). The method is based in matching autocorrelations 
*     of the series and correlations between innovations and the series 
*     of both models. 
* 
*     Input:  
*           phiar: vector of AR coefficients of the initial SAR model 
*           npdif: AR order of the initial SAR model 
*           thetas: seasonal MA coefficient of both models 
*           ip: AR order of the SARMA model 
*           iq: MA order of the SARMA model 
*           isp: seasonal period 
*           indth: 1 if the model contains a seasonal MA operator 
*                  0 if the model does not contain a seasonal MA operator 
*           rho: autocorrelation vector 
*           tau_i=E(u_(t-i)z_t)/E(u_t^2). The subroutine receives as input 
*           tau_i ,i=1,...isp*indth and computes them between isp*indth+1 
*           and isp*indth+iq. 
*           ndim2: max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to  
*                  dimension several arrays 
* 
*     Output: 
*           phidif: vector of AR coefficients of the SARMA model 
*           theta: vector of MA coefficients of the SARMA model 
*----------------------------------------------------------------------- 
      implicit double precision (a-h,o-z) 
      dimension phiar(ndim2),phidif(ip+1),theta(iq),rho(0:ndim2) 
      dimension tau(0:ndim2),a(ndim2,ndim2),b(ndim2) 
      integer ipiv(ndim2) 
      data zero/0.d0/
*-----------------------------------------------------------------------
      k=iq 
      nqaux=isp*indth 
* 
*     Here we compute tau(nqaux+i),i=1,...,iq, using the initial model 
*     with the formula: 
*           tau(i+1)=phi(1)*tau(i-1)+...+phi(npdif)*tau(i-npdif) 
* 
      if (k.gt.0) then 
         do i=1,k 
            tau(nqaux+i)=zero
            npdif1=min0(npdif,nqaux+i) 
            do j=1,npdif1 
               tau(nqaux+i)=tau(nqaux+i)+phiar(j)*tau(nqaux+i-j) 
            enddo 
            if (indth.eq.1)then 
               tau(nqaux+i)=tau(nqaux+i)+thetas*theta(i)*tau(0) 
            else 
               tau(nqaux+i)=tau(nqaux+i)-theta(i)*tau(0) 
            endif 
         enddo 
      endif 
* 
*     We are going to solve the system Az=b where A is (ip+iq)x(ip+iq) and 
*     b is a (ip+iq) column vector. The firs ip elements of z are going to 
*     be the phi's and the last iq the theta's. 
*     We start building the left top part (ipxip) of A=(a(i,j)), 
*      a(i,j)=rho(i-j) 
* 
      do i=1,ip 
         do j=1,ip 
            a(i,j)=rho(abs(i-j))
         enddo 
      enddo 
* 
*     We build the right top(ipxiq) part of A, a(i,ip+j)=0 if i>j  
*     a(i,ip+j)=-tau(j) if i<=j, a(i,ip+j)=a(i,ip+j)+thetas*tau(isp+j-i)  
*     if indth=0 and isp+j-i>=0 
*
      do i=1,ip 
         do j=1,iq 
            a(i,ip+j)=zero 
            k=j-i 
            if (k.ge.0) a(i,ip+j)=-tau(k) 
            if (indth.eq.1) then 
               k=isp+j-i 
               if (k.ge.0) a(i,ip+j)=a(i,ip+j)+thetas*tau(k) 
            endif 
         enddo 
      enddo 
* 
*     We build the left bottom part (iqxip) of A, a(ip+i,j)=0 if i<j,  
*     a(i,ip+j)=-tau(i-j) if i>=j.  
* 
      do i=1,iq 
         do j=1,ip 
            a(ip+i,j)=zero 
            k=i-j 
            if (k.ge.0) a(ip+i,j)=tau(k) 
         enddo 
      enddo 
* 
*     We build the right bottom part (iqxiq) of A, a(ip+i,ip+j)=0 if i<>j,  
*     a(ip+i,ip+j)=-tau(0) if i=j, a(ip+i,ip+j)=thetas*tau(0) if i=isp+j  
*     and indth=1 
* 
      do i=1,iq 
         do j=1,iq 
            a(ip+i,ip+j)=zero 
            if (i.eq.j) a(ip+i,ip+j)=-tau(0) 
            if (indth.eq.1) then 
               if (i.eq.isp+j) a(ip+i,ip+j)=a(ip+i,ip+j)+thetas*tau(0) 
            endif 
         enddo 
      enddo 
* 
*     We build the top part of b (first ip coordinates). 
* 
      do i=1,ip  
         b(i)=rho(i) 
         if (indth.eq.1) then 
            k=isp-i 
            if (k.ge.0) b(i)=b(i)+thetas*tau(k) 
         endif 
      enddo 
* 
*     We build the bottom part of b (last iq coordinates). 
* 
      do i=1,iq 
         b(ip+i)=tau(i) 
         if (indth.eq.1) then 
            if (isp.eq.i) b(ip+i)=b(ip+i)+thetas*tau(0) 
         endif 
      enddo 
* 
*     We solve A*z=b. The solution z is stored in b. 
*         
      call s_gesvfe(ip+iq,1,a,ndim2,ipiv,b,ndim2,ierror) 
* 
*     We extract vector phidif from b. 
* 
      do i=1,ip 
         phidif(i)=b(i) 
      enddo 
* 
*     We extract vector theta from b. 
* 
      do j=1,iq 
         theta(j)=b(ip+j) 
      enddo 
      return 
      end 
C=======================================================================
      subroutine s_tranfe(par,ndim1,ndim2,ip,iq,indth,m,para,para1, 
     +     work,phi,theta,thetas,beta) 
*-----------------------------------------------------------------------    
*     This subroutine antitransform the ndim1 variables from par, using 
*     subroutine s_invdfe. s_invdfe performs the inverse tranformation of the  
*     one which tranforms autoregressive parameters in partial  
*     autocorrelations. 
*  
*     Input:   
*             par      : vector containing the transformed parameters 
*                        of the model 
*             ndim1    : length of par. ndim1=ip+iq+indth+m+1.  
*             ndim2    : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required 
*                        to dimension several arrays. 
*             ip       : order of the AR model  
*             iq       : order of the ordinary MA model  
*             indth    : 0 - no seasonal moving average component 
*                        1 - seasonal moving average term included 
*             m        : number of independent variables 
*
*     Output: 
*             phi      : vector of coefficients of the AR models 
*             theta    : vector of ordinary MA coefficients 
*             thetas   : seasonal moving average coefficient 
*             beta     : vector of coeff. of the independent variables
*-----------------------------------------------------------------------    
      implicit double precision (a-h,o-z)
      dimension par(ndim1),para(ndim2),para1(ndim2),work(ndim2) 
      dimension phi(ndim2),theta(iq),beta(m)
c     data pi/3.1415926927d0/
      data pi/3.1416d0/
*-----------------------------------------------------------------------
      if (ip.gt.0) then       
         do i=1,ip 
            para(i)=par(i) 
         enddo 
         do i=1,ip 
            para1(i)=2.d0*datan(para(i))/pi
         enddo 
         call s_invdfe(para1,ip,phi,work,ndim2) 
      endif  
      if (iq.gt.0) then 
         do i=1,iq 
            para(i)= par(ip+i) 
         enddo 
         do i=1,iq 
            para1(i)=2.d0*datan(para(i))/pi 
         enddo 
         call s_invdfe(para1,iq,theta,work,ndim2) 
      endif 
      if (indth.eq.1) thetas=2.d0*datan(par(ip+iq+1))/pi 
      do i=1,m 
         beta(i) = par(ip+iq+indth+i) 
      enddo 
      return  
      end 
C=======================================================================
      subroutine s_trasfe(phi,theta,thetas,beta,ndim2,ip,iq,indth,m,
     +     para,par,ndim1,rh,work,nw,iwork,niw,irank) 
*-----------------------------------------------------------------------    
*     This subroutine computes partial autocorrelations and then 
*     transforms them in variables that may take any real value. 
* 
*     Input:   
*             phi      : vector of coefficients of the AR models 
*             theta    : vector of ordinary MA coefficients 
*             thetas   : seasonal moving average coefficient 
*             beta     : vector of coeff. of the independent variables         
*             ndim2    : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required 
*                        to dimension several arrays. 
*             ip       : order of the AR model  
*             iq       : order of the ordinary MA model  
*             indth    : 0 - no seasonal moving average component 
*                        1 - seasonal moving average term included 
*             m        : number of independent variables 
*             irank    : 0 - no control is made 
*                        1 - at the output of routine s_durbfe, the values 
*                            of vector para are controlled and corrected if 
*                            they are out of the interval [-1,1] 
*              
*     Output: 
*             par      : vector containing the transformed parameters 
*                        of the model 
*             ndim1    : length of par. ndim1=ip+iq+indth+m+1.  
*-----------------------------------------------------------------------     
      implicit double precision (a-h,o-z) 
      dimension phi(ndim2),theta(iq),beta(m),para(ndim2),par(ndim1) 
      dimension rh(ndim2+1),work(nw)
      integer iwork(niw) 
*-----------------------------------------------------------------------    
      if (ip.gt.0) then 
         call s_yulefe(phi,rh,ip,work(1),iwork(1),ndim2) 
         call s_durbfe(rh,ip,para,ier,work(1),ndim2) 
         if (irank.eq.1) then 
            do i=1,ip 
               if (para(i).ge.1.d0) then 
                  para(i)=.9d0 
               elseif (para(i).le.-1.d0) then 
                  para(i)=-.9d0 
               endif 
            enddo 
         endif 
         do i=1,ip 
            par(i)=dtan(3.1416d0*para(i)/2.d0)        
         enddo 
      endif 
      if (iq.gt.0) then 
         call s_yulefe(theta,rh,iq,work(1),iwork(1),ndim2) 
         call s_durbfe(rh,iq,para,ier,work(1),ndim2) 
         if (irank.eq.1) then 
            do i=1,iq 
               if (para(i).ge.1.d0) then 
                  para(i)=.9d0 
               elseif (para(i).le.-1.d0) then 
                  para(i)=-.9d0 
               endif 
            enddo 
         endif 
         do i=1,iq 
            par(ip+i)=dtan(3.1416d0*para(i)/2.d0) 
         enddo 
      endif 
      if(indth.eq.1) par(ip+iq+1)=dtan(3.1416d0*thetas/2.d0) 
      do i=1,m 
         par(ip+iq+indth+i)=beta(i) 
      enddo 
      return  
      end 
C=======================================================================
      subroutine s_vesrfe(x,y,nob,f1,e2,res,ares1,aux) 
*-----------------------------------------------------------------------  
*     This subroutine estimates the regression coefficient in a model 
*     with only one independent variable through the origin: 
*                    Y = f1 * X + U. 
* 
*     Input:  
*            x    : independent variable 
*            y    : dependent variable 
*            nob  : number of observations 
* 
*     Output: 
*            f1   : estimate of the regression coefficient 
*            e2   : mad of the residuals 
*-----------------------------------------------------------------------  
      implicit double precision (a-h,o-z) 
      dimension x(nob),y(nob),res(nob),ares1(nob),aux(nob)
*-----------------------------------------------------------------------   
      j=0 
*     
*     We compute res=y/x, if it is possible. 
* 
      do i=1,nob 
         if (dabs(x(i)).ge.1.d-15) then 
            j=j+1 
            res(j)=y(i)/x(i) 
         endif 
      enddo 
* 
*     We compute median(res). 
*     
      nob1=j 
      call s_mednfe(res,nob1,f1,aux) 
* 
*     We compute res=y-f1*x (residual corresponding to the estimated 
*     coefficient). 
* 
      do i=1,nob 
         res(i)=y(i)-f1*x(i) 
         ares1(i)=dabs(res(i)) 
      enddo 
* 
*     We compute the mad scale of the residuals. 
* 
      call s_mednfe(ares1,nob,e2,aux) 
      if (e2.ge.1.0D-10) then 
         z=0.d0 
         do i=1,nob 
            w=res(i)/e2 
            if (dabs(w).le. 2.5d0)then 
               z=z+w**2 
            else 
               z=z+6.25d0 
            endif 
         enddo 
         z=z/nob 
         e2=e2*dsqrt(z) 
      endif 
      return 
      end 
C=======================================================================
      function s_xmadfe(x,y,beta,m,n,eps,u,aux,polds,nds) 
*-----------------------------------------------------------------------  
*     This function estimates the scale of the differenced errors of  
*     the regression model. 
* 
*     Input: 
*             x        : matrix of independent variables. 
*             y        : input series  
*             beta     : vector of coeff. of the independent variables         
*             m        : number of independent variables 
*             n        : number of observations 
*             polds    : vector containing the coefficients of the 
*                        differences polynomial 
*             nds      : order of polds 
*-----------------------------------------------------------------------   
      implicit double precision (a-h,o-z)
      dimension x(n,m),y(n),beta(m),polds(nds+1),eps(n),u(n),aux(2*n)
      data zero/0.d0/
*-----------------------------------------------------------------------    
* 
*     We compute the regression residuals. 
* 
      do i=1,n 
         eps(i)=y(i) 
         do j=1,m 
            eps(i)=eps(i)-x(i,j)*beta(j) 
         enddo 
      enddo 
*     
*     The polynomial operator polds is applied to eps, in order to 
*     obtain the differenced series. 
* 
      do i=nds+1,n 
         u(i-nds)=zero
         do j=1,nds+1 
            u(i-nds)=u(i-nds)+polds(j)*eps(i-j+1) 
         enddo 
         u(i-nds)=dabs(u(i-nds)) 
      enddo 
* 
*     We compute an M-scale of the differenced series. 
* 
      call s_calsfe(u,n-nds,0,s_xmadfe,aux(1),aux(n+1)) 
      return 
      end
C=======================================================================
      subroutine s_yulefe(phif,rho,lp,a,ipiv,ndim2) 
*-----------------------------------------------------------------------
*     This subroutine computes lp autocorrelations given the coefficients  
*     of the autoregressive model of order lp, using the Yule-Walker 
*     equations. 
* 
*     Input: 
*          phif   : vector containing the coefficients of the AR model 
*          lp     : length of vector phif 
*          ndim2  : max0(ip+idif+isp*nsd,iqfin+indth*isp+1), required to 
*                   dimension the auxiliary arrays. 
*           
*     Output: 
*          rho  : vector containing lp autocorrelations 
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z) 
      dimension phif(lp),rho(lp),a(ndim2,ndim2) 
      integer ipiv(ndim2) 
*-----------------------------------------------------------------------
      do i=1,lp 
         do j=1,lp 
            a(i,j)=0.0d0 
            if ((i+j).le.lp) a(i,j)=a(i,j)+phif(j+i) 
            if ((i-j).ge.1)  a(i,j)=a(i,j)+phif(i-j) 
            if ((i-j).eq.0)  a(i,j)=a(i,j)-1 
         enddo 
      enddo 
      do i=1,lp 
         rho(i)=-phif(i) 
      enddo 
      call s_gesvfe(lp,1,a,ndim2,ipiv,rho,lp,ierror) 
      return 
      end
*----------------------------------------------------------------------- 
