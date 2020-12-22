subroutine gelnet(p,SS,la,al,TTh,Wm,TTar,maxit,thr,maxit2,thr2,niter,ipen,iact,dlz)
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
!         Inverse covariance matrix estimate
!
!  Wm     (input/output) double precision array, dimension p x p
!         Covariance matrix estimate
!
!  TTar   (input) double precision array, dimension p x p
!         Target Matrix
!
!  thr    (input) double precision
!         Convergence threshold outer loop
!
!  maxit  (input) integer
!         Maximum number of whole matrix sweeps
!
!  thr2   (input) double precision
!         Convergence threshold inner loop
!
!  maxit2 (input) integer
!         Maximum number of inner loop sweeps
!
!  niter  (output) integer
!         Actual number of outer loop sweeps
!
!  ipen   (input) logical
!         Should Diagonal be penalized
!
!  iact   (input) logical
!         Should Active Set be used
!
!  dlz    (output) double precision
!         average
!
      implicit double precision (a-h, o-z)
      integer :: p, maxit, maxit2, outer, i,j,k,l,m,n
      double precision SS(p,p), TTh(p,p), Wm(p,p), TTar(p,p), la(p,p),l1(p,p),l2(p,p)
      double precision, dimension (:), allocatable :: S,l3,l4,W,Th,Tar
      double precision al,thr,thr2,dlz
      integer, allocatable :: i2(:),i3(:)
      integer, allocatable :: i1(:,:)
      logical :: ipen,iact
      l1=al*la
      l2=(1-al)*la
      allocate(i1(2,p))
      allocate(i2(p))
      allocate(i3(p))
!     n5=#length per Components / n6= #number of components
      call connect(p,SS,l1,n6,i1,i2,i3)
      n5=0
10070 do 10080 k1=1,n6
      n5=max(i1(2,k1)-i1(1,k1)+1,n5)
10080 continue
      continue
      n5=n5**2
      allocate(S(1:n5))
      allocate(l3(1:n5))
      allocate(l4(1:n5))
      allocate(W(1:n5))
      allocate(Th(1:n5))
      allocate(Tar(1:n5))
      niter=0
      dlz=0.0
      l=0
10081 do 10088 k1=1,n6
      n=i1(2,k1)-i1(1,k1)+1
      if(n .gt. 1) goto 10082
      k=i2(i1(1,k1))
      Wm(:,k)=0.0
      Wm(k,:)=0.0
      TTh(:,k)=0.0
      TTh(k,:)=0.0
      goto 10088
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
      W(l)=Wm(i8,i7)
      Th(l)=TTh(i8,i7)
      Tar(l)=TTar(i8,i7)
10084 continue
      continue
10083 continue
      continue
      call gelnet_loop1(n,S,l3,l4,Th,W,Tar,maxit,thr,maxit2,thr2,outer,ipen,iact,dly)
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
      do 10087 k=k2,k3
      i7=i2(k)
      do 10086 j=k2,k3
      l=l+1
      TTh(i2(j),i7)=Th(l)
      Wm(i2(j),i7)=W(l)
10086 continue
      continue
10087 continue
      continue
10088 continue
      continue
      dlz=dlz/n6
      do 10094 j=1,p
      if(Wm(j,j).ne.0.0) goto 10094
      if(ipen .eqv. .TRUE.) goto 10092
      TTh(j,j)=1.0/SS(j,j)
      Wm(j,j)=SS(j,j)
      goto 10094
10091 continue
10092 continue
      if(TTar(j,j) .ne. 0.0) goto 10093
      if(al .eq. 1.0) then
        TTh(j,j)=1/(SS(j,j)+l1(j,j))
      else
        TTh(j,j)=(-SS(j,j)-l1(j,j)+SQRT((SS(j,j)+l1(j,j))**2+4*l2(j,j)))/(2*l2(j,j))
      endif
      Wm(j,j)=1.0/TTh(j,j)
      goto 10094
      continue
10093 continue
      if(TTar(j,j) < 1/(l1(j,j)+SS(j,j))) then
        if(al .eq. 1.0) then
          TTh(j,j)=1/(SS(j,j)+l1(j,j))
        else
          TTh(j,j)=(-SS(j,j)-l1(j,j)+l2(j,j)*TTar(j,j)+SQRT((SS(j,j)+l1(j,j)-l2(j,j)*TTar(j,j))**2+4*l2(j,j)))/(2*l2(j,j))
        endif
      else
        if(l1(j,j) < SS(j,j)) then
          if(TTar(j,j) > 1/(SS(j,j)-l1(j,j))) then
            if(al .eq. 1.0) then
              TTh(j,j) = 1/(SS(j,j)-l1(j,j))
            else
              TTh(j,j)=(-SS(j,j)+l1(j,j)+l2(j,j)*TTar(j,j)+SQRT((SS(j,j)-l1(j,j)-l2(j,j)*TTar(j,j))**2+4*l2(j,j)))/(2*l2(j,j))
            endif
          else
            TTh(j,j)=TTar(j,j)
          endif
        else
          TTh(j,j)=TTar(j,j)
        end if
      end if
      Wm(j,j)=1.0/TTh(j,j)
10094 continue
      return
      end
subroutine gelnet_loop1(n,S,l3,l4,Th,W,Tar,maxit,thr,maxit2,thr2,outer,ipen,iact,dly)
      parameter(eps=1.0e-16)
      integer n, outer, inner, maxit, maxit2, i, j
      double precision shr, thr, thr2, throut, thrin, dly, dlx, test, help, delta, a, b, c, Thd
      double precision S(n,n), l3(n,n), l4(n,n), Th(n,n), W(n,n), Tar(n,n), BoldMat(n,n)
      double precision b12(n), w12(n)
      logical :: ipen, iact
      shr = sum(abs(S))
      do i=1,n
        shr=shr-abs(S(i,i))
      end do
      if (shr .eq. 0.0) then
      ! S is diagonal
        return
      endif
      throut = thr*shr/(n-1)
      thrin = thr2*shr/(n-1)/n
      if(thrin .lt. 2*eps) then
          thrin=2*eps
      endif
      do i = 1,n
        BoldMat(1:n,i) = -Th(1:n,i)/Th(i,i)
        BoldMat(i,i) = 0
      end do
      do outer = 1,maxit
        dly = 0.0
        do j = 1,n
          b12 = BoldMat(:,j)
          b12(j) = 0.0
          w12(1:n) = 0.0
          do i = 1,n
            if(b12(i) .ne. 0.0) then
              w12 = w12 + W(:,i)*b12(i)
            end if
          end do
          do inner = 1,maxit2
            dlx=0.0
            do i = 1,n
              if(i .ne. j) then
                a = S(i,j) - w12(i) + W(i,i)*b12(i)
                b = abs(a) - l3(i,j)
                if(b .gt. 0.0) then
                  c = sign(b,a)/(W(i,i)+l4(i,j)*Th(j,j))
                else
                  c = 0.0
                end if
                delta = c - b12(i)
                !if(delta .ge. thrin/2) then
                if (delta .ne. 0.0) then
                  !Th(i,j) = c*Th(j,j)
                  b12(i) = c
                  w12(1:n) = w12(1:n) + delta*W(:,i)
                  dlx = max(dlx, abs(delta))
                end if
              end if
            end do
            if (dlx.lt.thrin) then
             exit
            end if
          end do
          BoldMat(:,j)=b12
          w12(j)=W(j,j)
          if(ipen .eqv. .TRUE.) then
            if(Tar(j,j) .eq. 0.0) then
              if(l4(j,j) .eq. 0.0 ) then
                  Th(j,j)=1/(S(j,j)+l3(j,j)-DOT_PRODUCT(w12,b12))
                else
                  help= S(j,j)+l3(j,j)-DOT_PRODUCT(w12,b12)
                  Th(j,j)=(-help+SQRT((help)**2+4*l4(j,j)))/(2*l4(j,j))
                endif
                W(j,j) = S(j,j)+l3(j,j)+l4(j,j)*Th(j,j)
            else
              test=1/Tar(j,j)+DOT_PRODUCT(w12,b12)-S(j,j)
              if(test .gt. l3(j,j)) then
                if(l4(j,j) .eq. 0.0 ) then
                  Th(j,j)=1/(S(j,j)+l3(j,j)-DOT_PRODUCT(w12,b12))
                else
                  help= S(j,j)+l3(j,j)-l4(j,j)*Tar(j,j)-DOT_PRODUCT(w12,b12)
                  Th(j,j)=(-help+SQRT((help)**2+4*l4(j,j)))/(2*l4(j,j))
                endif
                W(j,j) = S(j,j)+l3(j,j)+l4(j,j)*max(0.0,Th(j,j)-Tar(j,j))
              elseif(test .lt. -l3(j,j)) then
                if(l4(j,j) .eq. 0.0 ) then
                  Th(j,j)=1/(S(j,j)-l3(j,j)-DOT_PRODUCT(w12,b12))
                else
                  help= S(j,j)-l3(j,j)-l4(j,j)*Tar(j,j)-DOT_PRODUCT(w12,b12)
                  Th(j,j)=(-help+SQRT((help)**2+4*l4(j,j)))/(2*l4(j,j))
                endif
                W(j,j) = S(j,j)-l3(j,j)+l4(j,j)*min(0.0,Th(j,j)-Tar(j,j))
              else
                Th(j,j)=Tar(j,j)
                W(j,j)=S(j,j)+test
              endif
            endif
          else
            Th(j,j)=1/(S(j,j)-DOT_PRODUCT(w12,b12))
          endif
          Thd=Th(j,j)
          w12(j)=W(j,j)
          dly=max(dly,sum(abs(w12(1:n)-W(:,j))))
          Th(j,:) = -Thd*b12(:)
          Th(:,j) = -Thd*b12(:)
          Th(j,j) = Thd
          W(j,:)=w12(1:n)
          W(:,j)=w12(1:n)
        end do
        if (dly.lt.throut) then
          dly=dly/(n-1)
          exit
        end if
      end do
      return
      end subroutine gelnet_loop1
subroutine connect(p,SS,l1,n6,i1,i2,i3)
      implicit double precision(a-h,o-z)
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
      call row(n6,1,i2((i4-1):p),p,SS,l1,i3,n2,i2(i4:p))
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
      call row(n6,n3,i2(i5:p),p,SS,l1,i3,n2,i2(i4:p))
      if(n2.eq.0) goto 10522
      goto 10521
10522 continue
      i1(2,n6)=i6
10501 continue
      continue
      return
      end
! needed for connect
      subroutine row(n6,n4,jr,p,SS,l1,i3,n2,kr)
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
