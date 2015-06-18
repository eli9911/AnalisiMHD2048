!
!------------------------------------------------------
!
!------------------------------------------------------
!
       SUBROUTINE realis(isi,Nx,Ny,lx,ly,f,fr)
!
       IMPLICIT REAL*8 (a-h,o-z)
!
       COMPLEX*16 :: f(lx+1,ly)

       REAL*8 :: fr(Nx,Ny)

       REAL*8, ALLOCATABLE, DIMENSION( :, :) :: ffr !ffr(Nx,Ny)

       COMPLEX*16, ALLOCATABLE, DIMENSION( :, :) :: ff !ff(lx,ly)

       COMPLEX*16, ALLOCATABLE, DIMENSION( : ) :: speq  !speq(ly)
!
!
       dcoeff1=1.d0/(Nx*Ny)    

       IF( isi == -1 ) THEN

         ALLOCATE( ff(lx,ly) )
         ALLOCATE( speq(ly) )
         do my=1,ly
           do mx=1,lx
             ff(mx,my)=f(mx,my) 
           end do
           speq(my)=f(lx+1,my)
         end do
!
         CALL rlft3(ff,speq,Nx,Ny,1,isi)
!     
         DEALLOCATE( speq )
!
         do my=1,ly
           ii=1
           do mx=1,lx
             fr(ii,my)=2*DREAL(ff(mx,my))
             fr(ii+1,my)=2*DIMAG(ff(mx,my))
             ii=ii+2
           end do
         end do
         DEALLOCATE( ff )

       ELSE
!
         ALLOCATE( ffr(Nx,Ny) )
         do iy=1,Ny
           do ix=1,Nx
             ffr(ix,iy)=dcoeff1*fr(ix,iy)
           end do
         end do
!
         ALLOCATE( speq(ly) )
         CALL rlft3(ffr,speq,Nx,Ny,1,isi)
!
         do my=1,ly
           ii=1
           do mx=1,lx
             tmp1=ffr(ii,my)
             tmp2=ffr(ii+1,my)
             f(mx,my)=DCMPLX(tmp1,tmp2)
             ii=ii+2
           end do
         end do
!
         DEALLOCATE ( ffr )
!
         do my=1,ly
           f(lx+1,my)=speq(my)
         end do
!
         DEALLOCATE( speq )

!
       END IF
!    
       END SUBROUTINE realis
!
!--------------------------------------------------------------
!
!--------------------------------------------------------------
!
       SUBROUTINE enforce(lx,ly,f)
!
       IMPLICIT REAL*8 (a-h,o-z)
!
       COMPLEX*16 :: f(lx+1,ly)
!
!       m=  1 2 3 4 5  6  7  8
!      kx = 0,1,2,3,4
!      ky = 0,1,2,3,4,-3,-2,-1    [For ny = 8]
!
!
       mx=1
       mx_c=1
       do my=2,(ly/2)
         my_c=ly-my+2
         f(mx_c,my_c)=dconjg(f(mx,my))
       end do
!
       mx=lx+1
       mx_c=mx
       do my=2,(ly/2)
         my_c=ly-my+2  
         f(mx_c,my_c)=dconjg(f(mx,my))
       end do   
!

       !GIUSTO     
!
       tmp=DREAL(f(1,1))
       f(1,1)=DCMPLX(tmp,0.d0)
!
       tmp=DREAL(f(lx+1,(ly/2)+1))
       f(lx+1,(ly/2)+1)=DCMPLX(tmp,0.d0)
!
       tmp=DREAL(f(lx+1,1)) 
       f(lx+1,1)=DCMPLX(tmp,0.d0)
!
       tmp=DREAL(f(1,(ly/2)+1))
       f(1,(ly/2)+1)=DCMPLX(tmp,0.d0)
!
!
!
       END SUBROUTINE enforce
!
!------------------------------------------------------
!
!------------------------------------------------------
!
      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      DOUBLE PRECISION data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1, &
      k2,n,nprev,nrem,ntot
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=dble(wr)*data(k2)-dble(wi)*data(k2+1)
                tempi=dble(wr)*data(k2+1)+dble(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END SUBROUTINE fourn
!
!
!--------------------------------------------------------------
!
!--------------------------------------------------------------
!
      SUBROUTINE rlft3(data,speq,nn1,nn2,nn3,isign)
      INTEGER isign,nn1,nn2,nn3
      COMPLEX*16 data(nn1/2,nn2,nn3),speq(nn2,nn3)
!CU    USES fourn
      INTEGER i1,i2,i3,j1,j2,j3,nn(3)
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      COMPLEX*16 c1,c2,h1,h2,w
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0*isign)
      theta=6.28318530717959d0/dble(isign*nn1)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      nn(1)=nn1/2
      nn(2)=nn2
      nn(3)=nn3
      if(isign.eq.1)then
        call fourn(data,nn,3,isign)
        do 12 i3=1,nn3
          do 11 i2=1,nn2
            speq(i2,i3)=data(1,i2,i3)
11        continue
12      continue
      endif
      do 15 i3=1,nn3
        j3=1
        if (i3.ne.1) j3=nn3-i3+2
        wr=1.0d0
        wi=0.0d0
        do 14 i1=1,nn1/4+1
          j1=nn1/2-i1+2
          do 13 i2=1,nn2
            j2=1
            if (i2.ne.1) j2=nn2-i2+2
            if(i1.eq.1)then
              h1=c1*(data(1,i2,i3)+conjg(speq(j2,j3)))
              h2=c2*(data(1,i2,i3)-conjg(speq(j2,j3)))
              data(1,i2,i3)=h1+h2
              speq(j2,j3)=conjg(h1-h2)
            else
              h1=c1*(data(i1,i2,i3)+conjg(data(j1,j2,j3)))
              h2=c2*(data(i1,i2,i3)-conjg(data(j1,j2,j3)))
              data(i1,i2,i3)=h1+w*h2
              data(j1,j2,j3)=conjg(h1-w*h2)
            endif
13        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          w=dcmplx(dble(wr),dble(wi))
14      continue
15    continue
      if(isign.eq.-1)then
        call fourn(data,nn,3,isign)
      endif
      return
      END SUBROUTINE rlft3
!
!
!--------------------------------------------------------------
!
!--------------------------------------------------------------
!
