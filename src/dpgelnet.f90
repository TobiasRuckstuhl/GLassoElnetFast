subroutine dpgelnet(p,SS,la,al,TTh,Wm,UU,maxit,thr,maxit2,thr2,niter,ipen,dlz)
!
!  Arguments
!  =========
!
!  p      (input) integer
!         The dimension of the input matrix SS
!
!  SS     (input) double precision array, dimension p x p
!         The empirical covariance matrix
!
!  la     (input) double precision array, dimension p x p
!         Regularization matrix (symmetric)
!
!  al     (input) double precision
!         Alpha (Alpha=1 for Glasso / Alpha=0 for Rope)
!
!  TTh    (input/output) double precision array, dimension p x p
!         Precision matrix estimate
!
!  Wm     (input/output) double precision array, dimension p x p
!         Covariance matrix estimate
!
!  UU     (input/output) double precision array, dimension p x p
!
!  maxit  (input) integer
!         Maximum number of whole matrix sweeps
!
!  thr    (input) double precision
!         Convergence threshold outer loop
!
!  maxit2 (input) integer
!         Maximum number of inner loop sweeps
!
!  thr2   (input) double precision
!         Convergence threshold inner loop
!
!  niter  (output) integer
!         Actual number of outer loop sweeps
!
!  ipen   (input) logical
!         Should Diagonal be penalized
!
!  dlz    (output) double precision
!         average
!
      implicit double precision (d-h,o-z)
      integer :: p, maxit, maxit2, outer, i,j,k,l,m,n
      double precision SS(p,p), TTh(p,p), Wm(p,p), UU(p,p), la(p,p),l1(p,p),l2(p,p)
      double precision, dimension (:), allocatable :: S,l3,l4,U,Th
      double precision al,thr,thr2,dlz
      integer, allocatable :: i2(:),i3(:)
      integer, allocatable :: i1(:,:)
      logical :: ipen
      l1=al*la
      l2=(1-al)*la
      allocate(i1(2,p))
      allocate(i2(p))
      allocate(i3(p))
!     n5=#length per Components / n6= #number of components
      call dpconnect(p,SS,l1,n6,i1,i2,i3)
      n5=0
10070 do 10080 k1=1,n6
      n5=max(i1(2,k1)-i1(1,k1)+1,n5)
10080 continue
      continue
      n5=n5**2
      allocate(S(1:n5))
      allocate(l3(1:n5))
      allocate(l4(1:n5))
      allocate(U(1:n5))
      allocate(Th(1:n5))
      niter=0
      dlz=0.0
      l=0
10081 do 10089 k1=1,n6
      n=i1(2,k1)-i1(1,k1)+1
      if(n .gt. 1) goto 10082
      k=i2(i1(1,k1))
      Wm(:,k)=0.0
      Wm(k,:)=0.0
      TTh(:,k)=0.0
      TTh(k,:)=0.0
      goto 10089
10082 continue
      k2=i1(1,k1)
      k3=i1(2,k1)
      l=0
      do 10083 k=k2,k3
      i7=i2(k)
      do 10084 j=k2,k3
      i8=i2(j)
      l=l+1
      S(l)=SS(i8,i7)
      l3(l)=l1(i8,i7)
      l4(l)=l2(i8,i7)
      U(l)=UU(i8,i7)
      Th(l)=TTh(i8,i7)
10084 continue
      continue
10083 continue
      continue
      call dpgelnet_loop1(n,S,al,l3,l4,Th,U,maxit,thr,maxit2,thr2,outer,ipen,dly)
      niter=niter+outer
      dlz=dlz+dly
      do 10085 j=k2,k3
      k=i2(j)
      Wm(:,k)=0.0
      Wm(k,:)=0.0
      TTh(:,k)=0.0
      TTh(k,:)=0.0
10085 continue
      continue
      l=0
      do 10088 k=k2,k3
      i7=i2(k)
      do 10087 j=k2,k3
      l=l+1
      TTh(i2(j),i7)=Th(l)
      if (al.eq.0.0) then
         Wm(i2(j),i7)=S(l)+l4(l)*Th(l)
      else
         Wm(i2(j),i7)=S(l)+U(l)+l4(l)*Th(l)
      end if
10086 continue
10087 continue
      continue
10088 continue
      continue
10089 continue
      continue
      dlz=dlz/n6
      do 10093 j=1,p
      if(ipen .eqv. .TRUE.) goto 10090
      Wm(j,j)=SS(j,j)
10090 continue
      if(TTh(j,j).ne.0.0) goto 10093
      if(ipen .eqv. .TRUE.) goto 10092
      TTh(j,j)=1.0/SS(j,j)
      Wm(j,j)=SS(j,j)
      goto 10093
10091 continue
10092 continue
      if(al .eq. 1.0) then
        TTh(j,j)=1/(SS(j,j)+l1(j,j))
      else
        TTh(j,j)=(-SS(j,j)-l1(j,j)+SQRT((SS(j,j)+l1(j,j))**2+4*l2(j,j)))/(2*l2(j,j))
      endif
      Wm(j,j)=1/TTh(j,j)
      continue
10093 continue
      return
      end
! Outer loop of DPGelnet
subroutine dpgelnet_loop1(n,S,al,l3,l4,Th,U,maxit,thr,maxit2,thr2,outer,ipen,dly)
      parameter(eps = 1.0e-7)
      integer n, outer, inner, maxit, maxit2, i, j
      double precision al, shr, thr, thr2, dly, del, objcur, objold, Thd, Ud, a, b, c
      double precision S(n,n), l3(n,n), l4(n,n), Th(n,n), U(n,n)
      double precision t12(n), d(n), l7(n), grad(n)
      logical :: ipen
      do outer = 1,maxit
        dly = 0.0
        do  j = 1,n
          if (al .eq. 0.0) then
            d = l4(:,j)*abs(Th(j,:))
            l7 = 1
          else
            d = (1-al)/al*abs(Th(j,:)) + 1
          end if
          t12(1:n) = 0.0
          do i = 1,n
            if(i .ne. j) then
              t12 = t12 + Th(i,:)*(d(i)*U(j,i)+S(j,i))
            end if
          end do
          t12(j) = 0.0
          objold = DOT_PRODUCT(t12, S(j,:) + U(j,:)*d)
          do inner = 1,maxit2
            dlx = 0.0
            do i = 1,n
              if(i .ne. j) then
                if(d(i) .ne. 0.0) then
                  a = Th(i,i)*U(j,i)*d(i) - t12(i)
                  b = abs(a)/(Th(i,i)*d(i))
                  if(b .lt. l3(i,j)) then
                    c = sign(b,a)
                  else
                    c = sign(l3(i,j),a)
                  end if
                  del = c - U(j,i)
                  if(del .ne. 0.0) then
                    U(j,i) = c
                    !dlx=max(dlx,abs(del))
                    t12 = t12 + del*Th(i,:)*d(i)
                    t12(j) = 0.0
                  end if
                end if
              end if
            end do
            objcur = DOT_PRODUCT(t12, S(j,:) + U(j,:)*d)
            dlx = abs(objcur - objold)/(abs(objold) + eps)
            objold = objcur
            if (dlx .lt. thr2) then
              exit
            end if
          end do
          if(ipen .eqv. .TRUE.) then
            t12 = -t12/(S(j,j) + l3(j,j) + l4(j,j)*Th(j,j))
          else
            t12 = -t12/S(j,j)
          end if
          do i = 1,n
            if(i .ne. j) then
              if(abs(U(j,i)) .lt. l3(i,j)*(1-eps)) t12(i) = 0.0
              if(abs(t12(i)) .lt. eps) t12(i) = 0.0
            end if
          end do
          if(ipen.eqv. .TRUE.) then
            t12(j) = 0.0
            t12(j) = (1 - DOT_PRODUCT(S(j,:) + d*U(j,:), t12))/(S(j,j) + l3(j,j) + l4(j,j)*Th(j,j))
          else
            t12(j) = 0.0
            t12(j) = (1 - DOT_PRODUCT(S(j,:) + d*U(j,:), t12))/S(j,j)
          end if
          !dly=MAXVAL(abs(Told - Th))/MAXVAL(abs(Told))
          !Th(j,j) = t12
          dly = max(dly,MAXVAL(abs((Th(j,:) - t12)/Th(j,:))))
          Th(j,:) = t12
          Th(:,j) = t12
        end do
        if (dly.lt.thr) then
          exit
        end if
      end do
      return
    end subroutine dpgelnet_loop1
! Connected Components
subroutine dpconnect(p,SS,l1,n6,i1,i2,i3)
      implicit double precision(d-h,o-z)
      integer p,n6
      double precision SS(p,p),l1(p,p)
      integer :: i1(2,p),i2(p),i3(p)
      i3=0
      i4=1
      n6=0
10500 do 10501 k=1,p
      if(i3(k).gt.0) goto 10501
      i2(i4)=k
      n6=n6+1
      i3(k)=n6
      i1(1,n6)=i4
      i4=i4+1
      call dprow(n6,1,i2((i4-1):p),p,SS,l1,i3,n2,i2(i4:p))
      if(n2 .ne. 0) goto 10520
      i1(2,n6)=i4-1
      goto 10501
10520 continue
      continue
10521 continue
      n3=n2
      i5=i4
      i6=i5+n3-1
      if(i6.ge.p) goto 10522
      i4=i4+n2
      call dprow(n6,n3,i2(i5:p),p,SS,l1,i3,n2,i2(i4:p))
      if(n2.eq.0) goto 10522
      goto 10521
10522 continue
      i1(2,n6)=i6
10501 continue
      continue
      return
      end
! Needed for Connected Components
      subroutine dprow(n6,n4,jr,p,SS,l1,i3,n2,kr)
      implicit double precision(a-h,o-z)
      integer p
      double precision SS(p,p),l1(p,p)
      integer jr(n4),i3(p),kr(*)
      n2=0
      do 10530 l=1,n4
      k=jr(l)
      do 10531 j=1,p
      if(i3(j).gt.0) goto 10531
      if(j.eq.k)goto 10531
      if(abs(SS(j,k)).le.l1(j,k)) goto 10531
      n2=n2+1
      kr(n2)=j
      i3(j)=n6
10531 continue
      continue
10530 continue
      continue
      return
      end
