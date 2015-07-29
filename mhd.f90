 program mhd 
 
 use variabili 
 
 implicit none 

 integer:: i,j
 double precision :: E, Hc, Bq, vq, Energy
 double precision, dimension(2048, 2048,1) ::  a_b,  psi_v, bx_prova, by_prova, vx_prova, vy_prova, vxCoer, vyCoer, vxIncoer, vyInc
 double precision, dimension(2048, 2048,1) ::  coeff_fvx, coeff_fvy, vxinv, vyinv,coeff_fo, oinv, coeff_fj, jinv
 double precision, dimension(2048, 2048,1) ::  oinv_coer, coeffcoer_o, coeffincoer_o, oinv_incoer
 
 double precision, dimension(2048, 2048,1) ::  jinv_coer, coeffcoer_j, coeffincoer_j, jinv_incoer
 double precision, dimension(2048, 2048,1) ::  jinv_coer1, coeffcoer_j1, coeffincoer_j1, jinv_incoer1
 double precision, dimension(11) :: Evx, Evy, k, Etot, sigmaEx, sigmaEy, sigmaT, Eo, ko, sigma_o, Ej, kj, sigma_j, Ecoer_o, kco
 double precision, dimension(11) :: sig_coero, Eincoer_o, kinco, sig_incoero, k_A
 integer, dimension(2048, 2048,1):: control_o,  control_j, control_j1, control_A
 double precision :: epsilon_t, epsilon_j, epsilon_j1
 integer :: lx, ly, ii, jj 
 double precision :: media_vort, sigma_vort
 double precision, dimension(2048, 2048, 1) :: lambda, lambda_coer, lambda_incoer 
 double precision, dimension(11) :: sig_coerj, Eincoer_j, Ecoer_j, kinco_j, sig_incoerj, kco_j
 double precision, dimension(2048, 2048,1) :: coeff_fpsi, coeff_fA
 double precision, dimension(2048, 2048,1) ::  psi_coer, psi_incoer, A_coer, A_incoer 
 double precision, dimension(2048, 2048,1) ::  coeff_fvxC, coeff_fvyC
 double precision, dimension(11) :: Evx_Coer, Evy_Coer, Etot_Coer, sigmaExC, sigmaEyC, Evx_Inc, Evy_Inc, Etot_Inc
 double precision, allocatable,  dimension(:, :,:) ::  coeff_fvxI, coeff_fvyI, psi_CSoglia,  psi_ISoglia

 !VARIABILI PER I COEFF DI FOURIER 
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : ) :: u, r, psi
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: vorc, psic, psiCoer_c, psiIncoer_c
 INTEGER*4, ALLOCATABLE, DIMENSION( : ) :: k_x, k_y 
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : ) :: jF, Ab, Bx_F, By_F, Vx_F, Vy_F, VxCoer_F, VyCoer_F,  VxIncoer_F, VyIncoer_F
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: corrc, aa_c
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: bx_c, by_c, vx_c, vy_c, vxCoer_c, vyCoer_c, vxIncoer_c, vyIncoer_c 
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : ) :: ikx, iky, Ekk_psi, Ekk_coer, Ekk_incoer
 integer :: i_n
 double precision :: s
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: oCoer_cmplx, oIncoer_cmplx
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : ) :: O_coer, O_incoer
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , :, : ) :: vx_coer5S,  vx_incoer5S, vy_coer5S,  vy_incoer5S 
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: ECvx_5S, ECvy_5S, EIvx_5S, EIvy_5S, ECoer_5S,  EIncoer_5S

 !PROVA
 integer, dimension(2048, 2048,1):: control_vx, control_vy
 
 double precision, dimension(11) ::  kxi, kyi, sigmaxC, sigmayC, sigmaxI, sigmayI
 
 double precision, dimension(2048, 2048,1) :: vx_coerProva, vy_coerProva, vy_incoerProva, vx_incoerProva
 

 !
 double precision, allocatable,  dimension(:, :,:)  :: vx_coerP, vx_incoerP, vy_coerP, vy_incoerP 
 double precision, allocatable, dimension(:) :: kx, ky, ECvx, ECvy, EIvx, EIvy, ECoer,  EIncoer
 

 !
 integer :: nbins
 double precision, allocatable, dimension(:) :: veloc_x, p, w, erro
 !
 integer :: Mmax, Mmin, nbins1, nbins2
 integer ::  t, Np
 double precision :: dev_standard_Sx

 !
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: term1_c, djx_c, djy_c,  term3_c
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : ) :: term1, djx, djy, term3
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : , : ) :: dissipazione_vort, dissipazione_pot, combinazione, grad_jx, grad_jy
 
 !
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: phi

 double precision, allocatable, dimension(:) :: p_Sx, w_Sx

 !VARIABILI DIVISIONE PSI IN COER ED INCOER CON VX E SOGLIA 3 SIGMA
 DOUBLE PRECISION,ALLOCATABLE, DIMENSION( : , : ) :: psi_c, psi_inc, vy_coer3s, vy_incoer3s, vx_coer3s, vx_incoer3s
 
  DOUBLE PRECISION,ALLOCATABLE, DIMENSION( : , : ) :: A_c, A_inc

 !CORRELAZIONE J O
 double precision ::  correlazione_oj, corr_oj_c, corr_oj_i, correlazione_jA, corr_Aj_i, corr_Aj_c 
 
 double precision, allocatable, dimension(:) :: dev_standard_bx, dev_standard_by
 integer, allocatable, dimension(:) :: num_bins
 double precision, allocatable, dimension(:,:,:) :: coeff_bx, coeff_by, bxinv, byinv
 double precision, allocatable, dimension(:) ::  dev_standard_vx, dev_standard_vy

 double precision, allocatable, dimension(:) :: EC_Bx, EC_By, EI_Bx, EI_By, ECoer_B,  EIncoer_B
 double precision, allocatable, dimension(:,:, :) :: Bx_coer,  Bx_incoer, By_coer,  By_incoer

 double precision, allocatable, dimension(:,:, :) :: j_coer, j_inc
 double precision, allocatable, dimension(:,:) :: j_c, j_i

 double precision, allocatable, dimension(:,:) :: w_c, w_i
 double precision, allocatable, dimension(:,:, :) ::w_coer, w_inc, phi_c, phi_i


 
 !DEFINIZIONE PARAMETRI 
 Nx=2048
 Ny=Nx
 M=11
 ndim2=2**2
 pi=dacos(-1.d0)
 L=2.d0*pi

 !LETTURA CAMPI 
 print*, 'LETTURA FILE'
 call read_file(a, bx, by, dj, vx, vy, o)
 print*, 'FINE LETTURA FILE'

 open(unit=24, status='UNKNOWN', file='modulovelocit.dat',recl=10000)
 do j= 1,Ny
   do i = 1,Nx 
      write(24,*), i, j, sqrt(vx(i,j,1)*vx(i,j,1)+vy(i,j,1)*vy(i,j,1))
    enddo
     write(24,*), ''
 enddo
 close(24)
 open(unit=88, status='UNKNOWN', file='vort.txt',recl=10000)
 open(unit=22, status='UNKNOWN', file='vorticita.dat',recl=10000)
 do j= 1,Ny
   do i = 1,Nx 
      write(22,*), i, j, o(i,j,1)
      write(88,*), o(i,j,1)
    enddo
     write(22,*), ''
 enddo

 open(unit=23, status='UNKNOWN', file='corrente.dat',recl=10000)
 do j= 1,Ny
   do i = 1,Nx 
      write(23,*), i, j, dj(i,j,1)
    enddo
     write(23,*), ''
 enddo

 media_vort = sum(o)/(1.d0*Nx*Ny)
 sigma_vort = dsqrt(sum((o-media_vort)**2.d0)/(1.d0*Nx*Ny))
 !print*, 'La media del campo di vorticità è: ', media_vort
 !print*, 'Il sigma della vorticità è: ', sigma_vort
 !CALCOLO DEI COEFF WAVELET DEI CAMPI 

 call init_haar

 allocate(coeff_bx(2048,20148,1), coeff_by(2048,20148,1))
 allocate( bxinv(2048,20148,1), byinv(2048,20148,1))

 coeff_fvx = 0.d0
 coeff_fvy = 0.d0
 coeff_fo = 0.d0
 coeff_fj = 0.d0 
 coeff_fpsi= 0.d0
 coeff_fA=0.d0
 sig_coero=0.d0
 coeff_bx = 0.d0
 coeff_by = 0.d0

 !COEFF vx
 call directhaar(vx,coeff_fvx)
 call inverse(coeff_fvx,vxinv)
 !print*, maxval(dabs(vx-vxinv))
 
 !COEFF vy
 call directhaar(vy,coeff_fvy)
 call inverse(coeff_fvy,vyinv)
 !print*, maxval(dabs(vy-vyinv))


 !COEFF bx
 call directhaar(bx,coeff_bx)
 call inverse(coeff_bx,bxinv)
 !print*, maxval(dabs(bx-bxinv))

 !COEFF by
 call directhaar(by,coeff_by)
 call inverse(coeff_by,byinv)
 !print*, maxval(dabs(by-byinv))
 
 !COEFF VORTICITA' 
 call directhaar(o,coeff_fo)
 call inverse(coeff_fo,oinv)
 !print*, maxval(o-oinv)

 !COEFF CORRENTE 
 call directhaar(dj,coeff_fj)
 call inverse(coeff_fj,jinv)
 !print*, maxval(dabs(dj-jinv))


 !CALCOLO DELLO SPETTRO WAVELET PER LA VELOCITA'
 call spettro(coeff_fvx,Evx,k, sigmaEx)
 call spettro(coeff_fvy,Evy,k, sigmaEy)

 open(unit=1, status='UNKNOWN', file='spWvelocita.dat',recl=10000)
 open(unit=35, status='UNKNOWN', file='spdualePv.dat',recl=10000)
 open(unit=39, status='UNKNOWN', file='spdualeMv.dat',recl=10000)

 do i=1,M
   Etot(i)=Evx(i)+Evy(i)
   sigmaT(i)= sigmaEx(i)+ sigmaEy(i)
   write(1 ,*), k(i), Etot(i) 
   write(35,*), k(i), Etot(i)+sigmaT(i)
   write(39,*), k(i), Etot(i)-sigmaT(i)
 enddo 
 close(1)
 close(35)
 close(39)

 !CALCOLO DELLO SPETTRO WAVELET PER LA VORTICITA' 
 call spettro(coeff_fo, Eo, ko, sigma_o)
 print*, 'HA CALCOLATO SPETTRO VORTICITA '
 open(unit=40, status='UNKNOWN', file='spWvorticita.dat',recl=10000)
 open(unit=41, status='UNKNOWN', file='spdualePo.dat',recl=10000)
 open(unit=42, status='UNKNOWN', file='spdualeMo.dat',recl=10000)
 do i=1,M
   write(40,*), ko(i), Eo(i)
   write(41,*), ko(i), Eo(i)+sigma_o(i)
   write(42,*), ko(i), dabs(Eo(i)-sigma_o(i))
 enddo
 close(40)
 close(41)
 close(42)

 !CALCOLO DELLO SPETTRO WAVELET PER LA CORRENTE 
 call spettro(coeff_fj, Ej, kj, sigma_j)
 print*, 'HA CALCOLATO SPETTRO CORRENTE '
 open(unit=50, status='UNKNOWN', file='spWcorrente.dat',recl=10000)
 open(unit=51, status='UNKNOWN', file='spdualePj.dat',recl=10000)
 open(unit=53, status='UNKNOWN', file='spdualeMj.dat',recl=10000)
 do i=1,M
   write(50,*), kj(i), Ej(i)
   write(51,*), kj(i), Ej(i)+sigma_j(i)
   write(53,*), kj(i), Ej(i)-sigma_j(i)
 enddo
 close(50)
 close(51)
 close(53)

 !print*,'Teorema di Parseval per la velocità lungo x '
 !print*,'Spazio reale:   ', sum(vx*vx)
 !print*,'Spazio wavelet: ', sum(coeff_fvx*coeff_fvx)

 !print*,'Teorema di Parseval per la velocità lungo y '
 !print*,'Spazio reale:   ', sum(vy*vy)
 !print*,'Spazio wavelet: ', sum(coeff_fvy*coeff_fvy)

 !print*,'Teorema di Parseval per la vorticità '
 !print*,'Spazio reale:   ', sum(o*o)
 !print*,'Spazio wavelet: ', sum(coeff_fo*coeff_fo)

 !print*,'Teorema di Parseval per la corrente '
 !print*,'Spazio reale:   ', sum(dj*dj)
 !print*,'Spazio wavelet: ', sum(coeff_fj*coeff_fj)  

 
 !CALCOLO COEFF INVERSI E CAMPO INVERSO PER LA VORTICITA'
 !DIVISIONE DEL CAMPO DI VORTICITA' IN COERENTE E INCOERENTE 
 
 epsilon_t = (2.d0 * sum(o**2.d0)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0

 call coeff_coerenti(coeff_fo, o, control_o, epsilon_t)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_o(i,j,1) = control_o(i,j,1)*coeff_fo(i,j,1)
      coeffincoer_o(i,j,1)= coeff_fo(i,j,1)*(1-control_o(i,j,1))
    enddo
 enddo

 call inverse(coeffcoer_o, oinv_coer)

 open(unit=17, status='UNKNOWN', file='oinv_coer.dat',recl=10000)
 open(unit=701, status='UNKNOWN', file='ocoer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(17,*), i, j, oinv_coer(i,j,1)
       write(701,*),  oinv_coer(i,j,1)
   enddo
  write(17,*), ''
 enddo
 close(17)


 call inverse(coeffincoer_o, oinv_incoer)
 
 open(unit=71, status='UNKNOWN', file='oinv_incoer.dat',recl=10000)
 open(unit=700, status='UNKNOWN', file='oincoer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(71,*), i, j, oinv_incoer(i,j,1)
        write(700,*), oinv_incoer(i,j,1)
   enddo
  write(71,*), ''
 enddo
 
  
   
 !print*, 'La media del campo di vorticità inc è: ', sum(oinv_incoer)/(1.d0*Nx*Ny)
 !print*, 'Il sigma della vorticità inc è: ', dsqrt(sum((oinv_incoer-(sum(oinv_incoer)/(1.d0*Nx*Ny)))**2.d0)/(1.d0*Nx*Ny))
 
 
 !CALCOLO COEFF INVERSI E CAMPO INVERSO PER LA CORRENTE 
 !DIVISIONE DEL CAMPO DI CORRENTE IN COERENTE E INCOERENTE 

 epsilon_j = (2.d0 * sum(dj**2.d0)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0

 call coeff_coerenti(coeff_fj, dj, control_j, epsilon_j)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_j(i,j,1) = control_j(i,j,1)*coeff_fj(i,j,1)
      coeffincoer_j(i,j,1)= coeff_fj(i,j,1)*(1-control_j(i,j,1))
    enddo
 enddo

 call inverse(coeffcoer_j, jinv_coer)

  open(unit=67, status='UNKNOWN', file='jinv_coer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(67,*), i, j, jinv_coer(i,j,1)
   enddo
  write(67,*), ''
 enddo
 close(67)

 call inverse(coeffincoer_j, jinv_incoer)
 
 open(unit=76, status='UNKNOWN', file='jinv_incoer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(76,*), i, j, jinv_incoer(i,j,1)
   enddo
  write(76,*), ''
 enddo

 epsilon_j1 = (2.d0 * sum(jinv_incoer**2.d0)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0

 call coeff_coerenti(coeff_fj, dj, control_j1, epsilon_j1)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_j1(i,j,1) = control_j1(i,j,1)*coeff_fj(i,j,1)
      coeffincoer_j1(i,j,1)= coeff_fj(i,j,1)*(1-control_j1(i,j,1))
    enddo
 enddo

 call inverse(coeffincoer_j1, jinv_incoer1)
 
 open(unit=44, status='UNKNOWN', file='jinv_incoer1.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(44,*), i, j, jinv_incoer1(i,j,1)
   enddo
  write(44,*), ''
 enddo

 !CALCOLO DELLO SPETTRO DELLA PARTE COERENTE E INCOERENTE DELLA VORTICITA'

 call spettro(coeffcoer_o, Ecoer_o, kco, sig_coero)
 call spettro(coeffincoer_o, Eincoer_o, kinco, sig_incoero)

 open(unit=900, status='UNKNOWN', file='spvort_coer.dat',recl=10000)
 open(unit=700, status='UNKNOWN', file='spvort_incoer.dat',recl=10000)
 do j=1,M 
      write(900,*), kco(j), Ecoer_o(j)
      write(700,*), kinco(j), Eincoer_o(j)
 enddo
 close(900)
 close(700)

 !CALCOLO DELLO SPETTRO DELLA PARTE COERENTE E INCOERENTE DELLA CORRENTE
 call spettro(coeffcoer_j, Ecoer_j, kco_j, sig_coerj)
 call spettro(coeffincoer_j, Eincoer_j, kinco_j, sig_incoerj)

 open(unit=90, status='UNKNOWN', file='spcoerr_coer.dat',recl=10000)
 open(unit=72, status='UNKNOWN', file='spcorr_incoer.dat',recl=10000)
 do j=1,M 
      write(90,*), kco_j(j), Ecoer_j(j)
      write(72,*), kinco_j(j), Eincoer_j(j)
 enddo
 close(90)
 close(72)

 
 !CALCOLO DEI COEFF DI FOURIER 
 
 lx = Nx/2 !Punti nello spazio fisico
 ly = Ny
!~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  
! K-vectors definition
!           1 2 3 4 5  6  7  8     
!      kx = 0,1,2..nx/2  
!      ky = 0,1,2,3,4,-3,-2,-1    [For ny = 8]

 ALLOCATE( k_x(lx+1), k_y(ly) )
 k_x = 0
 k_y = 0

 DO ii = 1, ly
    IF( ii <= (ly/2) )      k_y(ii) = ii - 1
    IF( ii ==  ((ly/2)+1) ) k_y(ii) = (ly/2)
    IF( ii > ((ly/2)+1) )   k_y(ii) = ii - ly - 1
 END DO

 DO ii = 1, lx + 1
     k_x(ii) = ii - 1
 END DO
!~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  

 !CALCOLO DELLA STREAM FUNCTION 

 ALLOCATE( u(Nx,Ny), vorc(lx+1,ly), r(lx+1,ly) )
 allocate(Ekk_psi(lx+1))

 psi_v= 0.d0
 u = 0.d0
 do j=1,Ny
   do i=1,Nx
      u(i,j)=o(i,j,1)
    enddo
 enddo 
 vorc = DCMPLX(0.d0,0.d0)

 !isi = 1 (trasformata diretta), -1 (inversa)
 CALL realis(1,Nx,Ny,lx,ly,vorc,u)  

 ALLOCATE(psic(lx+1,ly), psi(Nx,Ny))


 do j=1,ly
   do i=1,lx+1
     if((k_x(i).NE. 0).or.(k_y(j).NE. 0)) psic(i,j)= vorc(i,j)/(k_x(i)*k_x(i)+k_y(j)*k_y(j))
   enddo
 enddo

 open(unit=400, status='UNKNOWN', file='prova.dat',recl=10000)
 do j=1,ly
   do i=1,lx+1
      write(400,*), i, j, psic(i,j)
   enddo
  write(400,*), ''
 enddo
 close(400)
 
 !SPETTRO FOURIER
 do i_n=1,lx+1
   Ekk_psi(i_n)=0.d0
   do j=1,ly
    do i=1,lx+1
       s = dsqrt( i*i* 1.d0 +j*j* 1.d0)
      if((s.ge.i_n*1.d0).and.(s < i_n*1.d0 + 1.d0)) then
        Ekk_psi(i_n)= (k_x(i)*k_x(i)+k_y(j)*k_y(j)) * psic(i,j)* CONJG(psic(i,j)) + Ekk_psi(i_n)
      endif
    enddo
   enddo
 enddo

 open(unit=2346, status='UNKNOWN', file='Spettro_Fourier_psi.dat',recl=10000)
 do i_n=1,lx+1
   write(2346,*), i_n, DREAL(Ekk_psi(i_n)) 
 enddo
 close(2346)


 CALL realis(-1,Nx,Ny,lx,ly,psic,psi) 
 
 do j=1,Ny
   do i=1,Nx
      psi_v(i,j,1)= psi(i,j)
    enddo
 enddo 
 
 open(unit=100, status='UNKNOWN', file='psi.dat',recl=10000)
 do j=1,Ny
   do i=1,Nx
      write(100,*), i, j, psi_v(i,j,1)
   enddo
  write(100,*), ''
 enddo

 close(100)


 !CALCOLO SPETTRI CON FOURIER 
 ALLOCATE( O_coer(Nx,Ny), O_incoer(Nx,Ny), oCoer_cmplx(lx+1,ly), oIncoer_cmplx(lx+1,ly) )
 allocate(Ekk_coer(lx+1), Ekk_incoer(lx+1))

 do j=1,Ny
   do i=1,Nx
      O_coer(i,j)= oinv_coer(i,j,1)
      O_incoer(i,j)= oinv_incoer(i,j,1)
    enddo
 enddo 
 oCoer_cmplx = DCMPLX(0.d0,0.d0)
 oIncoer_cmplx = DCMPLX(0.d0,0.d0)

 !isi = 1 (trasformata diretta), -1 (inversa)
 CALL realis(1,Nx,Ny,lx,ly,oCoer_cmplx , O_coer) 
 CALL realis(1,Nx,Ny,lx,ly,oIncoer_cmplx , O_incoer) 

 do i_n=1,lx+1
   Ekk_coer(i_n)=DCMPLX(0.d0,0.d0)
   Ekk_incoer(i_n) = DCMPLX(0.d0,0.d0)
   do j=1,ly
    do i=1,lx+1
       s = dsqrt( i*i* 1.d0 +j*j* 1.d0)
      if((s.ge.i_n*1.d0).and.(s < i_n*1.d0 + 1.d0)) then
      Ekk_coer(i_n)=  oCoer_cmplx(i,j)* CONJG(oCoer_cmplx(i,j)) / (k_x(i)*k_x(i)+k_y(j)*k_y(j)) + Ekk_coer(i_n)
    Ekk_incoer(i_n)= oIncoer_cmplx(i,j)* CONJG(oIncoer_cmplx(i,j)) / (k_x(i)*k_x(i)+k_y(j)*k_y(j))  + Ekk_incoer(i_n)
      endif
    enddo
   enddo
 enddo

 open(unit=2846, status='UNKNOWN', file='Spettro_Fourier_diviso.dat',recl=10000)
 do i_n=1,lx+1
   write(2846,*), i_n, DREAL(Ekk_coer(i_n)), DREAL(Ekk_incoer(i_n)) 
 enddo
 close(2846)
 

 !CALCOLO DEL POTENZIALE VETTORE CON FOURIER 

 ALLOCATE( jF(Nx,Ny), corrc(lx+1,ly) )

 a_b= 0.d0
 jF = 0.d0
 do j=1,Ny
   do i=1,Nx
     jF(i,j)=dj(i,j,1)
    enddo
 enddo 
 corrc = DCMPLX(0.d0,0.d0)

 !isi = 1 (trasformata diretta), -1 (inversa)
 CALL realis(1,Nx,Ny,lx,ly,corrc,jF) 

 ALLOCATE(aa_c(lx+1,ly), Ab(Nx,Ny))


 do j=1,ly
   do i=1,lx+1
     if((k_x(i).NE. 0).or.(k_y(j).NE. 0)) aa_c(i,j)= corrc(i,j)/(k_x(i)*k_x(i)+k_y(j)*k_y(j))
   enddo
 enddo

! open(unit=4200, status='UNKNOWN', file='prova.dat',recl=10000)
! do j=1,ly
!   do i=1,lx+1
!      write(400,*), i, j, psic(i,j)
!  enddo
!  write(400,*), ''
! enddo

 CALL realis(-1,Nx,Ny,lx,ly, aa_c, Ab) 
 
 do j=1,Ny
   do i=1,Nx
      a_b(i,j,1)= Ab(i,j)
    enddo
 enddo  

 open(unit=600, status='UNKNOWN', file='potenziale_vettoreEl.dat',recl=10000)
 do j=1,Ny
   do i=1,Nx
      write(600,*), i, j, a_b(i,j,1), a(i,j,1)-a_b(i,j,1)
   enddo
  write(600,*), ''
 enddo
 close(600)

 !print*, 'Errore sul calcolo di fourier del potenziale vettore ', maxval(dabs(a-a_b))

 !CALCOLO DEL CAMPO MAGNETICO A PARTIRE DAL POTENZIALE VETTORE 

 ALLOCATE(ikx(lx+1), iky(ly))
 ALLOCATE(bx_c(lx+1,ly), by_c(lx+1,ly), Bx_F(Nx,Ny), By_F(Nx,Ny))

 do j=1,ly
   do i=1,lx+1
      iky(j)=DCMPLX(0.d0,-k_y(j))
      ikx(i)=DCMPLX(0.d0, k_x(i))
      bx_c(i,j) = iky(j)* aa_c(i,j)
      by_c(i,j) = ikx(i)* aa_c(i,j)
    enddo
 enddo 

 CALL realis(-1,Nx,Ny,lx,ly, bx_c, Bx_F) 
 CALL realis(-1,Nx,Ny,lx,ly, by_c, By_F) 

  do j=1,Ny
   do i=1,Nx
      bx_prova(i,j,1)= Bx_F(i,j)
      by_prova(i,j,1)= By_F(i,j)
    enddo
 enddo  

 !print*,'Errore sul calcolo di fourier del campo magnetico lungo x',  maxval(abs(bx_prova-bx))
 !print*,'Errore sul calcolo di fourier del campo magnetico lungo y',  maxval(abs(by_prova-by))

 !CALCOLO DEL CAMPO DI VELOCITA' A PARTIRE DALLA STREAM FUNCTION 

 ALLOCATE(vx_c(lx+1,ly), vy_c(lx+1,ly), Vx_F(Nx,Ny), Vy_F(Nx,Ny))

 do j=1,ly
   do i=1,lx+1
      iky(j)=DCMPLX(0.d0,-k_y(j))
      ikx(i)=DCMPLX(0.d0, k_x(i))
      vx_c(i,j) = iky(j)* psic(i,j)
      vy_c(i,j) = ikx(i)* psic(i,j)
    enddo
 enddo 

 CALL realis(-1,Nx,Ny,lx,ly, vx_c, Vx_F) 
 CALL realis(-1,Nx,Ny,lx,ly, vy_c, Vy_F) 

  do j=1,Ny
   do i=1,Nx
      vx_prova(i,j,1)= Vx_F(i,j)
      vy_prova(i,j,1)= Vy_F(i,j)
    enddo
 enddo  

 !print*,'Errore sul calcolo di fourier del campo di velocità lungo x',  maxval(abs(vx_prova-vx))
 !print*,'Errore sul calcolo di fourier del campo di velocità lungo y',  maxval(abs(vy_prova-vy))

 !PDF vx, vy, bx, by ALLE DIVERSE SCALE

 allocate(dev_standard_vx(11), num_bins(11), dev_standard_vy(11))
 allocate(dev_standard_bx(11), dev_standard_by(11))

 num_bins(1)= 100
 num_bins(2)= 100
 num_bins(3)= 80
 num_bins(4)= 50
 num_bins(5)= 40
 num_bins(6)= 40
 num_bins(7)= 40
 num_bins(8)= 20
 num_bins(9)= 6
 num_bins(10)=2
 num_bins(11)=1

 do i=1,11
    call calcolo_dev_standard(i, num_bins(i), coeff_fvx, dev_standard_vx(i), 'vx')
    call calcolo_dev_standard(i, num_bins(i), coeff_fvy, dev_standard_vy(i), 'vy')
    call calcolo_dev_standard(i, num_bins(i), coeff_bx, dev_standard_bx(i), 'bx')
    call calcolo_dev_standard(i, num_bins(i), coeff_by, dev_standard_by(i), 'by')
 enddo 


 !CALCOLO DI v COERENTE ED INCOERENTE CON LA SOGLIA 3*SIGMA

 allocate( kx(11), ky(11), ECvx(11), ECvy(11), EIvx(11), EIvy(11))
 allocate( ECoer(11),  EIncoer(11))
 allocate(vx_coerP (2048,2048,1),  vx_incoerP (2048,2048,1) )
 allocate(vy_coerP (2048,2048,1),  vy_incoerP (2048,2048,1) )

 call Divisione_campo(coeff_fvx, dev_standard_vx, vx_coerP, vx_incoerP, EIvx, ECvx, kx, 3.d0, .false.)
 
 call Divisione_campo(coeff_fvy, dev_standard_vy, vy_coerP, vy_incoerP, EIvy, ECvy, ky, 3.d0, .true.)
 
 open(unit=8834, status='UNKNOWN', file='spettro_velocitasigma.dat',recl=10000)
 
 do i=1,M
   ECoer(i)=ECvx(i)+ECvy(i)
   EIncoer(i)=EIvx(i)+EIvy(i)
   write(8834 ,*), kx(i), ECoer(i),  EIncoer(i)
 enddo 
 close(8834)
 
 !print*,'controllo vx',  maxval(abs(vx- vx_coerP- vx_incoerP))
 !print*, 'controllo vy', maxval(abs(vy- vy_coerP- vy_incoerP))

 !COEFF WAVELET STREAM FUNCTION E POTENZIALE VETTORE 
 call directhaar(psi_v,coeff_fpsi)
 call directhaar(a,coeff_fA) 

 !DIVISIONE DI PSI IN PARTE COERENTE E INCOERENTE CON LA SOGLIA DEI 3 SIGMA
 
 ALLOCATE( psi_c(Nx,Ny), psi_inc(Nx,Ny))
 
 call divisione_psi_A(vy_coerP, vy_incoerP, vx_coerP, vx_incoerP, psi_c, psi_inc) 
 
 open(unit=58, status='UNKNOWN', file='psi_coer.dat',recl=10000)
 open(unit=85, status='UNKNOWN', file='psi_incoer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       psi_coer(i,j,1)= psi_c(i,j)
       psi_incoer(i,j,1)= psi_inc(i,j)
       write(58,*), i, j, psi_coer(i,j,1) 
       write(85,*), i, j, psi_incoer(i,j,1)
   enddo
   write(58,*), ''
   write(85,*), ''
 enddo
 close(58)
 close(85)

 !print*, 'prova sulla divisione di psi: ', maxval(abs(psi_v- psi_coer - psi_incoer))
 !print*, 'media psi_coer:', sum(psi_coer)/(Nx*Ny), 'media psi_incoer:', sum(psi_incoer)/(Nx*Ny)
 !print*, 'maxval psi_coer:', maxval(psi_coer), 'minval psi_coer:', minval(psi_coer)
 !print*, 'maxval psi_incoer:', maxval(psi_incoer), 'minval psi_incoer:', minval(psi_incoer)

 open(unit=5899, status='UNKNOWN', file='psi_coersoglia.dat',recl=10000)
 allocate(psi_CSoglia(Nx,Ny,1))
 psi_CSoglia=0.d0
 do j= 1,Ny
  do i = 1,Nx 
      if(dabs(psi_coer(i,j,1)) .ge. 0.1d0) then
           psi_CSoglia(i,j,1)=1.d0
           write(5899,*), i, j, psi_CSoglia(i,j,1) 
       endif
   enddo
   !write(5899,*), ''
 enddo
 close(5899)
      

 open(unit=5998, status='UNKNOWN', file='psi_incoersoglia.dat',recl=10000)
 allocate(psi_ISoglia(Nx,Ny,1))
 psi_ISoglia=0.d0
 do j= 1,Ny
  do i = 1,Nx 
      if(dabs(psi_incoer(i,j,1)) .ge. 0.1d0) then
           psi_ISoglia(i,j,1)=1.d0
           write(5998,*), i, j, psi_ISoglia(i,j,1) 
       endif
   enddo
   !write(5899,*), ''
 enddo
 close(5998)
      

 !CALCOLO DELLA VORTICITA' COER ED INCOER USANDO V CON LA SOGLIA 3 SIGMA 

 allocate(w_c(Nx,Ny), w_i(Nx,Ny))
 allocate(w_coer(Nx,Ny,1), w_inc(Nx,Ny,1))

 call divisione_omega_j(vy_coerP, vy_incoerP, vx_coerP, vx_incoerP, w_c, w_i)

 open(unit=7117, status='UNKNOWN', file='w_coer.dat',recl=10000) 
 open(unit=8118, status='UNKNOWN', file='w_incoer.dat',recl=10000)

  do j= 1,Ny
   do i = 1,Nx 
       w_coer(i,j,1)= w_c(i,j)
       w_inc(i,j,1)= w_i(i,j)
       write(7117,*), i, j, w_coer(i,j,1)
       write(8118,*), i, j, w_inc(i,j,1)
   enddo
   write(7117,*), ''
   write(8118,*), ''
 enddo
 close(8118)
 close(7117)


 !CALCOLO DI B COERENTE ED INCOERENTE CON LA SOGLIA 3*SIGMA

 allocate( EC_Bx(11), EC_By(11), EI_Bx(11), EI_By(11))
 allocate( ECoer_B(11),  EIncoer_B(11))
 allocate(Bx_coer (2048,2048,1),  Bx_incoer (2048,2048,1) )
 allocate(By_coer (2048,2048,1),  By_incoer (2048,2048,1) )

 call Divisione_campo(coeff_bx, dev_standard_bx, bx_coer, bx_incoer, EI_bx, EC_bx, kx, 3.d0, .false.)
 
 call Divisione_campo(coeff_by, dev_standard_by, by_coer, by_incoer, EI_by, EC_by, ky, 3.d0, .true.)
 
 open(unit=2000, status='UNKNOWN', file='spettro_Bsigma.dat',recl=10000)
 
 do i=1,M
   ECoer_B(i)=EC_Bx(i)+EC_By(i)
   EIncoer_B(i)=EI_Bx(i)+EI_By(i)
   write(2000 ,*), kx(i), ECoer_B(i),  EIncoer_B(i)
 enddo 
 close(2000)
 
 print*,'controllo bx',  maxval(abs(bx- bx_coer- bx_incoer))
 print*, 'controllo by', maxval(abs(by- by_coer- by_incoer))

 
 !DIVISIONE DEL POTENZIALE VETTORE IN PARTE COERENTE ED INCOERENTE USANDO BX E BY DIVISI CON LA SOGLIA 3 SIGMA 

 ALLOCATE( A_c(Nx,Ny), A_inc(Nx,Ny))
 
 call divisione_psi_A(by_coer, by_incoer, bx_coer, bx_incoer, A_c, A_inc) 
 
 open(unit=87, status='UNKNOWN', file='A_coer.dat',recl=10000) 
 open(unit=78, status='UNKNOWN', file='A_incoer.dat',recl=10000)

  do j= 1,Ny
   do i = 1,Nx 
       A_coer(i,j,1)= A_c(i,j)
       A_incoer(i,j,1)= A_inc(i,j)
       write(87,*), i, j, A_coer(i,j,1)
       write(78,*), i, j, A_incoer(i,j,1)
   enddo
   write(87,*), ''
   write(78,*), ''
 enddo
 close(78)
 close(87)


 open(unit=8531, status='UNKNOWN', file='controllo_A.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(8531,*), i, j, A_incoer(i,j,1) + A_coer(i,j,1), a(i,j,1)
        write(8531,*),  abs(a(i,j,1)-A_incoer(i,j,1)-A_coer(i,j,1))
   enddo
 enddo
 close(8531)

 !CALCOLO DI J COER ED INCOER USANDO B CON LA SOGLIA 3 SIGMA 

 allocate(j_c(Nx,Ny), j_i(Nx,Ny))
 allocate(j_coer(Nx,Ny,1), j_inc(Nx,Ny,1))

 call divisione_omega_j(by_coer, by_incoer, bx_coer, bx_incoer, j_c, j_i)

 open(unit=8117, status='UNKNOWN', file='j_coer.dat',recl=10000) 
 open(unit=7118, status='UNKNOWN', file='j_incoer.dat',recl=10000)

  do j= 1,Ny
   do i = 1,Nx 
       j_coer(i,j,1)= j_c(i,j)
       j_inc(i,j,1)= j_i(i,j)
       write(8117,*), i, j, j_coer(i,j,1)
       write(7118,*), i, j, j_inc(i,j,1)
   enddo
   write(8117,*), ''
   write(7118,*), ''
 enddo
 close(7118)
 close(8117)




 !CALCOLO DELLA VELOCITA' COERENTE ED INCOERENTE MEDIANTE PSI 
 !QUI PSI DIVISO E' QUELLO OTTENUTO PRIMA DALLA DIVISIONE CON V E LA SOGLIA DI SIGMA 
 !DOVEREBBE ESSERE UN CONTROLLO

 ALLOCATE( psiCoer_c (lx+1,ly) )
 psiCoer_c = DCMPLX(0.d0,0.d0)

 CALL realis(1,Nx,Ny,lx,ly, psiCoer_c, psi_coer) 

 ALLOCATE( psiIncoer_c (lx+1,ly) )
 psiIncoer_c = DCMPLX(0.d0,0.d0)

 CALL realis(1,Nx,Ny,lx,ly, psiIncoer_c, psi_incoer) 
 
 ALLOCATE(vxCoer_c(lx+1,ly), vyCoer_c(lx+1,ly), VxCoer_F(Nx,Ny), VyCoer_F(Nx,Ny))
 ALLOCATE(vxIncoer_c(lx+1,ly), vyIncoer_c(lx+1,ly), VxIncoer_F(Nx,Ny), VyIncoer_F(Nx,Ny))

 do j=1,ly
   do i=1,lx+1
      iky(j)=DCMPLX(0.d0,-k_y(j))
      ikx(i)=DCMPLX(0.d0, k_x(i))
      vxCoer_c(i,j) = iky(j)* psiCoer_c(i,j)
      vyCoer_c(i,j) = ikx(i)* psiCoer_c(i,j)
      vxIncoer_c(i,j) = iky(j)* psiIncoer_c(i,j)
      vyIncoer_c(i,j) = ikx(i)* psiIncoer_c(i,j)
    enddo
 enddo 

 CALL realis(-1,Nx,Ny,lx,ly, vxCoer_c, VxCoer_F) 
 CALL realis(-1,Nx,Ny,lx,ly, vyCoer_c, VyCoer_F) 

 CALL realis(-1,Nx,Ny,lx,ly, vxIncoer_c, VxIncoer_F) 
 CALL realis(-1,Nx,Ny,lx,ly, vyIncoer_c, VyIncoer_F) 

  do j=1,Ny
   do i=1,Nx
      vxCoer(i,j,1)= VxCoer_F(i,j)
      vyCoer(i,j,1)= VyCoer_F(i,j)
      vxIncoer(i,j,1)= VxIncoer_F(i,j)
      vyInc(i,j,1) = VyIncoer_F(i,j) 
      !vyIncoer(i,j,1)= VyIncoer_F(i,j)
    enddo
 enddo 

 !CONTROLLO SULLA DIVISIONE
 open(unit=531, status='UNKNOWN', file='controllo_v.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(531,*), 'vx'
       write(531,*), i, j, abs(- vxIncoer(i,j,1) - vxCoer(i,j,1)+ vx(i,j,1))
       write(531,*), 'vy'
       write(531,*), i, j, abs(- vyInc(i,j,1) - vyCoer(i,j,1)+ vy(i,j,1))
   enddo
 enddo
 close(531)

 !CALCOLO DELLO SPETTRO WAVELET PER LA VELOCITA' COERENTE ED INCOERENTE 

 allocate(coeff_fvxI(2048,2048,1), coeff_fvyI(2048,2048,1))
 
 k=0.d0
 call directhaar(vxCoer,coeff_fvxC)
 call directhaar(vyCoer,coeff_fvyC)
 call directhaar(vxIncoer,coeff_fvxI)
 call directhaar(vyInc,coeff_fvyI)

 call spettro(coeff_fvxC,Evx_Coer,k, sigmaExC)
 call spettro(coeff_fvyC,Evy_Coer,k, sigmaEyC)
 call spettro(coeff_fvxI,Evx_Inc,k, sigmaExC)
 call spettro(coeff_fvyI,Evy_Inc,k, sigmaEyC)


 open(unit=991, status='UNKNOWN', file='spettro_velocitaCoerInc.dat',recl=10000)
 
 do i=1,M
   Etot_Coer(i)=Evx_Coer(i)+Evy_Coer(i)
   Etot_Inc(i)=Evx_Inc(i)+Evy_Inc(i)
   !sigmaT(i)= sigmaEx(i)+ sigmaEx(i)
   write(991 ,*), k(i), Etot_Coer(i), Etot_Inc(i), Etot(i)
 enddo 
 close(991)
 

 

 !CALCOLO DI LAMBDA 
 !CALCOLO DI LAMBDA COER ED INCOER CON LA SOGLIA DEI 3 SIGMA SU V
 open(unit=333, status='UNKNOWN', file='lambda.dat',recl=10000)
 open(unit=555, status='UNKNOWN', file='lambdacoer.dat',recl=10000)
 open(unit=566, status='UNKNOWN', file='lambdaincoer.dat',recl=10000)

 do j=1,Ny
   do i= 1,Nx
      lambda(i,j,1)= o(i,j,1)/dj(i,j,1)
      lambda_coer(i,j,1)= w_coer(i,j,1)/j_coer(i,j,1)
      lambda_incoer(i,j,1)= w_inc(i,j,1)/j_inc(i,j,1) 
     write(333,*), i, j, lambda(i,j,1)
     write(555,*), i, j, lambda_coer(i,j,1)
     write(566,*), i, j, lambda_incoer(i,j,1)
    enddo
   write(333,*), ''
   write(555,*), ''
   write(566,*), ''
 enddo 
 close(333)
 close(555)
 close(566)

 open(unit=31, status='UNKNOWN', file='controllo_lambda.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(31,*), i, j, abs(lambda(i,j,1)-lambda_incoer(i,j,1)-lambda_coer(i,j,1))
   enddo
        write(31,*),  
 enddo
 close(31)


 !CALCOLO PDF SIGMA DI V
 
 allocate(p_Sx(20), w_Sx(20) )
 call pdfeli(sigmaEx, M, 20, maxval(sigmaEx), minval(sigmaEx), w_Sx, p_Sx)
 dev_standard_Sx= dsqrt( sum((w_Sx- sum(w_Sx)/(20.d0))**2.d0) /(20.d0))

 open(unit=5521, status='unknown', file='pdf_Sigma_vx.dat', recl=100000)
 do i=1,20
     write(5521,*), i, w_Sx(i)/dev_standard_Sx, p_Sx(i)
 enddo
 close(5521)

 !CALCOLO PDF VELOCITA' VX

 nbins = 100
 allocate(veloc_x(Nx*Nx))
 allocate(p(nbins),w(nbins),erro(nbins))

 open(unit=834, status='UNKNOWN', file='confrontov.dat',recl=10000)
 do j=1,Nx
  do i=1,Nx
     veloc_x((j-1)*Nx+i)=vx(i,j,1)
     write(834,*), vx(i,j,1),  veloc_x((j-1)*Nx+i)
  enddo
 enddo
 close(834)
 
 call pdfeli(veloc_x, Nx*Nx, nbins,maxval(veloc_x), minval(veloc_x), w, p)

 open(unit=521, status='unknown', file='pdf_vx.dat', recl=100000)
 do i=1,nbins
     write(521,*), i, w(i), p(i)
 enddo
 close(521)


 !CALCOLO DERIVATA TOTALE PER LA VORTICITA' E IL POTENZIALE VETTORE

 ALLOCATE(term1_c(lx+1,ly), djx_c(lx+1,ly), djy_c(lx+1,ly), term3_c(lx+1,ly), term1(Nx,Ny), djx(Nx,Ny), djy(Nx,Ny), term3(Nx,Ny))
 allocate(dissipazione_vort(Nx,Ny,1), dissipazione_pot(Nx,Ny,1), combinazione(Nx,Ny,1), grad_jx(Nx,Ny,1), grad_jy(Nx,Ny,1) )
 
 do j=1,ly
   do i=1,lx+1
     term1_c(i,j)= vorc(i,j)*(k_x(i)*k_x(i)+k_y(j)*k_y(j))
   enddo
 enddo
 

 CALL realis(-1,Nx,Ny,lx,ly, term1_c, term1) 
 
 open(unit=63, status='unknown', file='dissipazione_vort.dat', recl=100000)
 do j=1,Ny
   do i=1,Nx
      dissipazione_vort(i,j,1)= 0.0003d0* term1(i,j)
      write(63,*), i, j, dissipazione_vort(i,j,1)
    enddo
    write(63,*), ''
 enddo
 close(63)  
 

 do j=1,ly
   do i=1,lx+1
    term3_c(i,j)= aa_c(i,j)*(k_x(i)*k_x(i)+k_y(j)*k_y(j))
   enddo
 enddo
 
 
  CALL realis(-1,Nx,Ny,lx,ly, term3_c, term3) 
 
 open(unit=64, status='unknown', file='dissipazione_pot.dat', recl=100000)
  do j=1,Ny
   do i=1,Nx
      dissipazione_pot(i,j,1)= 0.0003d0* term3(i,j)
      write(64,*), i, j, dissipazione_pot(i,j,1)
    enddo
    write(64,*), ''
 enddo 
 close(64)
 
 
 iky=DCMPLX(0.d0, 0.d0)
 ikx=DCMPLX(0.d0, 0.d0) 


 do j=1,ly
  do i=1,lx+1
     ikx(i)=DCMPLX(0.d0, -k_x(i))
     iky(j)=DCMPLX(0.d0, -k_y(j))
     djx_c(i,j) = ikx(i)*corrc(i,j)
     djy_c(i,j) = iky(j)*corrc(i,j)
  enddo
 enddo 


  CALL realis(-1,Nx,Ny,lx,ly, djx_c, djx) 
  CALL realis(-1,Nx,Ny,lx,ly, djy_c, djy) 

 do j=1,Ny
   do i=1,Nx
    grad_jx(i,j,1)=djx(i,j)
    grad_jy(i,j,1)=djy(i,j)
   enddo
 enddo


 open(unit=65, status='unknown', file='combinazione.dat', recl=100000)
 do j=1,Ny
   do i=1,Nx
     combinazione(i,j,1)= bx(i,j,1)*grad_jx(i,j,1)+ by(i,j,1)*grad_jy(i,j,1)
     write(65,*), i, j, combinazione(i,j,1)
   enddo
   write(65,*), ''
 enddo 
 close(65) 


 open(unit=69, status='unknown', file='somma_term.dat', recl=100000)
 do j=1,Ny
   do i=1,Nx
     write(69,*), i, j, dissipazione_vort(i,j,1) + combinazione(i,j,1)
   enddo
   write(69,*), ''
 enddo 
 close(69) 


 !CALCOLO PHI 

 allocate(phi(Nx,Nx,1))
 
 phi=0.d0
 allocate(phi_c(Nx,Ny,1), phi_i(Nx,Ny,1))

 open(unit=6900, status='unknown', file='phi.dat', recl=100000)
 open(unit=69110, status='unknown', file='phicoer.dat', recl=100000)
 open(unit=69220, status='unknown', file='phiinc.dat', recl=100000)
 do j=1,Nx
   do i =1, Nx
     if(a(i,j,1) .NE. 0.d0 ) then
      phi(i,j,1)=  0.5d0*(lambda(i,j,1)*o(i,j,1)- dj(i,j,1) )/ a(i,j,1)
      phi_c(i,j,1)=0.5d0*(lambda_coer(i,j,1)*w_coer(i,j,1)- j_coer(i,j,1) )/A_coer(i,j,1)
      phi_i(i,j,1)=0.5d0*(lambda_incoer(i,j,1)*w_inc(i,j,1)- j_inc(i,j,1) )/A_incoer(i,j,1)
    endif 
    write(6900,*), i, j, phi(i,j,1)
    write(69110,*), i, j, phi_c(i,j,1)
    write(69220,*), i, j, phi_i(i,j,1)
    enddo
    write(6900,*), ''
    write(69110,*), ''
    write(69220,*), ''
 enddo
 close(6900)
 close(69110)
 close(69220)

 !CALCOLO CORRELAZIONE FRA OMEGA E J
 call correlazione(dj, o, correlazione_oj)
 call correlazione(dj, a, correlazione_jA)
 call correlazione(j_coer, w_coer, corr_oj_c)
 call correlazione(j_inc, w_inc, corr_oj_i)
 call correlazione(j_inc, A_incoer, corr_Aj_i)
 call correlazione(j_coer, A_coer, corr_Aj_c)

 open(unit=69400, status='unknown', file='correlazione.dat', recl=100000)
 write(69400,*), 'correlazione tra vorticità e corrente:', correlazione_oj 
 write(69400,*), 'correlazione tra vorticità e corrente coerente:', corr_oj_c 
 write(69400,*), 'correlazione tra vorticità e corrente incoerente:', corr_oj_i 
 write(69400,*), 'correlazione tra potenziale vettore e corrente:', correlazione_jA 
 write(69400,*), 'correlazione tra potenziale vettore e corrente coer:', corr_Aj_c 
 write(69400,*), 'correlazione tra potenziale vettore e corrente incoer:', corr_Aj_i 
 close(69400)

 end program mhd 


 INCLUDE "S_FFT2D_JH.f90"

!----------------------------------------------------------
 subroutine read_file(aa, B, BB, JJ, VV, V, Omega)
!----------------------------------------------------------
  implicit none 
  integer:: i, j, Nx, Ny 
  double precision, dimension(2048,2048,1) :: aa, BB, B, JJ, VV, V, Omega 
 
  Nx=2048
  Ny=Nx
 
  open(unit=152, status='unknown', file='CAMPI2D_t1.2.dat', recl=100000)
  do j=1,Nx
    do i=1,Nx
      read(152,*), aa(i,j,1), B(i,j,1), BB(i,j,1), JJ(i,j,1), VV(i,j,1), V(i,j,1), Omega(i,j,1)
    enddo
  enddo 
 close(152)
 end subroutine read_file


!--------------------------------------
  subroutine init_haar
!--------------------------------------
  use variabili
  implicit none
  integer :: N2 , ip, jp
   
  allocate(Ni(M+1),Nf(M+1))
  dx=L/dfloat(Nx)
  dy=dx
  dz=dx
  N2=nx/2
  Ni(1)=1
  Nf(1)=N2
  do ip=2,M
    Ni(ip)=N2+Ni(ip-1)
    N2=N2/2
    Nf(ip)=N2+Nf(ip-1)
  end do
  Ni(M+1)=Nx
  Nf(M+1)=Nx

  allocate(Lscale(M))
  do ip=1,M
    Lscale(ip)=2**ip*dx
  enddo

  allocate (sigq(ndim2,ndim2))
  allocate (Isigq(ndim2,ndim2))
 
  sigq=0.d0
  Isigq=0.d0
 
 !Riempimento della matrice sigq
  sigq(1,1)=  1.d0
  sigq(1,2)= -1.d0
  sigq(1,3)= -1.d0
  sigq(1,4)=  1.d0

  sigq(2,1)=  1.d0
  sigq(2,2)=  1.d0
  sigq(2,3)= -1.d0
  sigq(2,4)= -1.d0

  sigq(3,1)=  1.d0
  sigq(3,2)= -1.d0
  sigq(3,3)=  1.d0
  sigq(3,4)= -1.d0

  sigq(4,1:4)=  1.d0
  
  do ip=1,ndim2
    Isigq(ip,1:ndim2) = sigq(1:ndim2,ip)
  enddo

  end subroutine init_haar


!-----------------------------------------------
 subroutine directhaar(field,coeff_f)
!-----------------------------------------------
  use variabili
  implicit none

  integer :: ip,Np
  integer :: iq,iq1,q1,q2,iav1,iav2
  integer :: k1,k2,n1,n2

  double precision :: facto, mean
  double precision, dimension(:), allocatable :: wm1,Sm
  double precision, dimension(:,:,:), allocatable :: work
  double precision, dimension(Nx,Ny,1) :: coeff_f,field

  mean=sum(field)/dfloat(Nx*Ny)

  allocate (work(Nx,Ny,1))
  allocate(Sm(ndim2),wm1(ndim2))
  
  work=field

  facto=(0.5d0)**(2.d0/2.d0) 

   do ip=1,M
     Np=2**(M-ip)
     do k2=1,Np
       do k1=1,Np
         Sm(1)=work(2*k1-1,2*k2-1,1)
         Sm(2)=work(2*k1-1,2*k2  ,1)
         Sm(3)=work(2*k1  ,2*k2-1,1)
         Sm(4)=work(2*k1  ,2*k2  ,1)
         iav1=k1
         iav2=k2

         do iq=1,ndim2
           wm1(iq)=facto*sum(sigq(iq,:)*Sm(:))
         enddo
         
         do iq=1,ndim2-1
           iq1=iq-1
           q1=iq1/2;  q2=(iq1-q1*2)
           n1=k1 + q1*Np + Ni(ip) -1
           n2=k2 + q2*Np + Ni(ip) -1

           coeff_f(n1,n2,1)=wm1(iq)
         enddo
         work(iav1,iav2,1) = wm1(ndim2)
      enddo
    enddo
  enddo
  coeff_f(nx,ny,1) = wm1(ndim2)
 

  deallocate(work) 
  deallocate(Sm,wm1)
  end subroutine directhaar 

!--------------------------------------------
 subroutine  inverse(coeff_f, fieldinv)
!--------------------------------------------
 use variabili 
 implicit none

  integer :: i, j, ip, Np
  integer :: ix, iy
  integer :: n1, n2,  q1, q2, q3, iq1, iq 
  integer :: k1, k2
  double precision :: mean, facto
  double precision, dimension(:), allocatable :: wm1
  double precision, dimension(:,:,:), allocatable ::  SS1
  double precision, dimension(Nx,Ny,1) :: coeff_f, fieldinv 

  allocate (SS1(Nx,Ny,1))
  allocate(wm1(ndim2))
  fieldinv=0.d0
  facto=(0.5d0)**(2.d0/2.d0)
 
  do ip= M,1,-1
   Np= 2**(M-ip)
   do k2=1,Np
    do k1=1,Np
     if (ip.eq. M) then
      do iq=1,ndim2
       iq1=iq-1
       q1=iq1/2;  q2=(iq1-q1*2)
       n1=k1 + q1*Np + Ni(ip) -1
       n2=k2 + q2*Np + Ni(ip) -1
       wm1(iq)=coeff_f(n1,n2,1)
      enddo

      fieldinv(2*k1-1,2*k2-1,1)=facto*sum(Isigq(1,:)*wm1(:))
      fieldinv(2*k1-1,2*k2,1)=facto*sum(Isigq(2,:)*wm1(:))
      fieldinv(2*k1,2*k2-1,1)=facto*sum(Isigq(3,:)*wm1(:))
      fieldinv(2*k1,2*k2,1)=facto*sum(Isigq(4,:)*wm1(:))
      
    

     else
    
      do iq=1,ndim2-1
       iq1=iq-1
       q1=iq1/2;  q2=(iq1-q1*2)
       n1=k1 + q1*Np + Ni(ip) -1
       n2=k2 + q2*Np + Ni(ip) -1
       wm1(iq)=coeff_f(n1,n2,1)
      enddo
      wm1(4)=fieldinv(k1,k2,1)
      
      SS1(2*k1-1,2*k2-1,1)=facto*sum(Isigq(1,:)*wm1(:))
      SS1(2*k1-1,2*k2,1)=facto*sum(Isigq(2,:)*wm1(:))
      SS1(2*k1,2*k2-1,1)=facto*sum(Isigq(3,:)*wm1(:))
      SS1(2*k1,2*k2,1)=facto*sum(Isigq(4,:)*wm1(:))

     endif
    enddo
   enddo
  if (ip.lt. M) then
     fieldinv=SS1
   endif
  enddo

  deallocate(SS1)
  deallocate(wm1)
 
 end subroutine inverse

!--------------------------------------------
 subroutine spettro(coeff_f, E, k, sigmaE)
!--------------------------------------------
  use variabili 
  implicit none 
  
 integer :: ip, Np, iq
 integer :: n1, n2 
 integer :: k1, k2 
 integer :: iq1,q1,q2
 double precision :: S, kk,cost
 double precision, dimension(M) :: E, k, sigmaE
 double precision, dimension(:), allocatable :: wm1 
 double precision, dimension(:,:), allocatable :: em 
 double precision, dimension(Nx,Ny,1) :: coeff_f 

 allocate (wm1(ndim2-1))
 allocate (em(2**(M-1),2**(M-1)))

 E=0.d0
 k=0.d0
 
 do ip=1,M
   em=0.d0
   Np = 2**(M-ip)
   k(ip)= 2.d0*pi*2.d0**(-ip)*(dx*dy)**(-0.5d0)
  
   cost=2.d0**(-2.d0*ip)*(dx*dy)**(0.5d0)/(2.d0*pi*dlog(2.d0))
     do k2=1,Np
       do k1=1,Np
           wm1=0.d0
           do iq=1,ndim2-1
             iq1=iq-1
             q1=iq1/2;  q2=(iq1-q1*2) 
             n1=k1 + q1*Np + Ni(ip) -1
             n2=k2 + q2*Np + Ni(ip) -1
             wm1(iq)= coeff_f(n1,n2,1)
           enddo
           S=0.d0
           do iq=1,ndim2-1
             S=S+(wm1(iq))**2.d0
	   enddo
          !em(k1,k2,1)= 0.5d0*S 
          em(k1,k2)= S 
       enddo 
     enddo
  E(ip)=0.d0 
  do k2=1,Np
    do k1=1,Np
       E(ip)=E(ip)+em(k1,k2)
     enddo
   enddo
   !E(ip)=4.d0*pi*pi*ip*E(ip)/(dk*dk*1024*1024)
   !E(ip)= E(ip)/(Np*1024*dlog(2))**2.d0
    E(ip)= E(ip)/(2048*dlog(2.d0)*2.d0**(ip)*Np*Np)
   
   
    !sigmaE(ip) = cost* (sum(em**2.d0)/(Np*Np)-(sum(em)/(Np*Np))**2.d0)**0.5d0
    sigmaE(ip)=  (sum(em**2.d0)/(Np*Np)-(sum(em)/(Np*Np))**2.d0)**0.5d0 / (2048*dlog(2.d0)*2.d0**(ip))
    
 enddo

 end subroutine spettro


!----------------------------------------------------------------------------------------
 subroutine coeff_coerenti(coeff_f, field, coeff_coer,espilon_t)
!----------------------------------------------------------------------------------------
 
 use variabili 
 
 double precision, dimension(Nx,Ny,1) :: coeff_f, field 
 integer, dimension(Nx,Ny,1) :: coeff_coer
 !double precision, allocatable, dimension(:,:,:) :: coeff_coer 
 integer :: i,  j, cont 
 double precision :: epsilon_t 


! epsilon_t = (2.d0 * sum(field**2.d0)/(Nx*Ny) *dlog(Nx*1.d0) )**0.5d0
 
 !print*, epsilon_t

 cont=0
 coeff_coer = 0 
 
 open(unit=80, status='UNKNOWN', file='coeffcoerenti.dat',recl=10000) 
 
 do j= 1,Ny
   do i = 1,Nx
     if ( coeff_f(i,j,1) > epsilon_t ) then
          coeff_coer(i,j,1) = 1
          write(80,*), i,j, coeff_f(i,j,1)
          cont=cont+1
     endif 
   enddo
 enddo 
 write(80,*), ''
 write(80,*), cont 
 close(80)


 end subroutine coeff_coerenti

!----------------------------------------------------------------------------------------
 subroutine Divisione_campo(coeff_f,  sigma, fieldcoer, fieldinco, EI, EC, k, s, B )
!----------------------------------------------------------------------------------------
 
 use variabili 

  implicit none 
  
 integer :: ip, Np, iq, i, j
 integer :: n1, n2 
 integer :: k1, k2 
 integer :: iq1,q1,q2
 integer, dimension(Nx,Ny,1) :: control 
 double precision, dimension(Nx,Ny,1) :: coeff_f,  coeff_coer,  coeff_incoer, fieldcoer, fieldinco, prova
 double precision, dimension(M) :: sigma 
 double precision, dimension(M) :: EI, EC, k, sigmaC, sigmaI
 double precision :: s
 logical :: B

 open(unit=7773, status='unknown', file='punti_coer_vy.dat', recl=100000)
 M=11
 prova=1.d0
 control=0
 do ip=M,1, -1
   Np = 2**(M-ip)
  !if(s.eq. 5.d0) print*,'Sigma',ip,  sigma(ip) 
   do k2=1,Np
      do k1=1,Np
         do iq=1,ndim2-1
           iq1=iq-1
            q1=iq1/2;  q2=(iq1-q1*2) 
            n1=k1 + q1*Np + Ni(ip) -1
            n2=k2 + q2*Np + Ni(ip) -1
            if ( dabs(coeff_f(n1,n2,1)) > s*sigma(ip) ) then
               control(n1,n2,1)=1
               if (B.eqv. (.true.)) write(7773,*), n1, n2, 1, ip
            endif           
          enddo
        enddo 
     enddo
 enddo
  do j=1,Ny
    do i =1,Nx
      coeff_coer(i,j,1) = control(i,j,1)*coeff_f(i,j,1)
      coeff_incoer(i,j,1)= coeff_f(i,j,1)*(1-control(i,j,1))
    enddo
 enddo
 close(7773)
 

 !print*, 'massimo valore coeff coer', maxval(coeff_coer), 'minimo valore coeff coer', minval(coeff_coer)
! print*, 'massimo valore ', maxval(sigma), 'minimo valore coeff ', minval(sigma)
 !print*, 'massimo valore coeff ', maxval(coeff_f), 'minimo valore coeff ', minval(coeff_f)

 call inverse(coeff_coer, fieldcoer)
 call inverse(coeff_incoer, fieldinco)

 !CALCOLO DELLO SPETTRO

 call spettro(coeff_coer,EC, k, sigmaC)


 call spettro(coeff_incoer,EI, k, sigmaI)
  


 end subroutine Divisione_campo




!-------------------------------------------------------------------------------
 subroutine pdfeli(x, n, nbins, xmax, xmin, w, p)
!-------------------------------------------------------------------------------
 
 use variabili 
 implicit none
 
 integer :: n, i, j, nbins, cont
 double precision :: xmin, xmax, bin
 double precision, dimension(n) :: x, y 
 double precision, dimension(nbins) ::  w, p


 cont=0
 open(unit=5253, status='unknown', file='y.dat', recl=100000)

 do i=1,n
   if((x(i) >= xmin) .and. (x(i) <= xmax)) then
     cont = cont + 1 
     y(cont)=x(i)
     write(5253,*), cont, y(cont),x(i)
   endif 
  enddo
 close(5253)
    
 bin = (xmax- xmin)/(nbins*1.d0)

 p=0.d0
 open(unit=5251, status='unknown', file='p.dat', recl=100000)

 do i=1,cont 
   j= (y(i)-xmin)/bin + 1
   p(j)=p(j)+1.d0
   write(5251,*), i, j, p(j)
 enddo
 close(5251)

 do i=1,nbins
   w(i)=xmin+(i*1.d0 + 0.5d0)*bin
   p(i)=p(i)/(n)
 enddo
 
	   
 end subroutine pdfeli 



 !-----------------------------------------------------------------------------------
  subroutine pdf_coeff(M_s, nbins, coeff_f, w_Mm, p_M, dev_standard, media)
 !-----------------------------------------------------------------------------------

 use variabili 
 implicit none

  integer :: M_s, Np, nbins, t, i, j 
  double precision :: dev_standard, media
  double precision, dimension(Nx,Ny,1) :: coeff_f
  double precision,  dimension(nbins) :: w_Mm, p_M
  double precision, dimension(2**(2*(M-M_s))*3) :: coeff_M


 Np = 2**(M-M_s)

 coeff_M= 0.d0

 t=1

 do j=Ni(M_s), Np+ Ni(M_s)-1
  do i=Ni(M_s), Np+ Ni(M_s)-1
      coeff_M(t)=coeff_f(i,j,1)
     t=t+1
   enddo
 enddo


 do j=Ni(M_s), Np+ Ni(M_s)-1
  do i= Np+ Ni(M_s), 2**M
      coeff_M(t)=coeff_f(i,j,1)
     t=t+1
   enddo
 enddo
 

 do j= Np+ Ni(M_s), 2**M
  do i= Ni(M_s), Np+ Ni(M_s)-1
     coeff_M(t)=coeff_f(i,j,1)
     t=t+1
   enddo
 enddo

 !dev_standard= dsqrt( sum((coeff_M- sum(coeff_M)/(Np*Np*3.d0))**2.d0) /(Np*Np*3.d0))
 media = sum(coeff_M)/(Np*Np*3.d0)
 !dev_standard = sqrt(sum((coeff_M-media)**2.d0)/(Np*Np*3.d0))
 dev_standard = dsqrt(sum((coeff_M)**2.d0)/ (Np*Np*3.d0) - media*media)

 open(unit=52541, status='unknown', file='dev_standard_vy.dat', recl=100000)
 write(52541,*), M_s,  dev_standard
 close(52541)

 call pdfeli(coeff_M, Np*Np*3 , nbins ,maxval(coeff_M), minval(coeff_M), w_Mm, p_M)
 !print*, 'area', M_s, ':', sum(p_M)
 
 end subroutine pdf_coeff

 !----------------------------------------------------------------------------------------------
 subroutine divisione_psi_A(Fy_c, Fy_i, Fx_c, Fx_i, Field_c, Field_inc)
 !------------------------------------------------------------------------------------------------
 !CALCOLA LA STREAM FUNCTION (O IL POTENZIALE VETTORE) DIVISA IN PARTE COERENTE ED INCOERENTE  A 
 !PARTIRE DAL CAMPO DI VELOCITA' (MAGNETICO) DIVISO USANDO LA SOGLIA DI SIGMA.
 !

 
 use variabili 
 implicit none
 
 integer :: lx, ly, i ,j, ii 
 INTEGER*4, DIMENSION( Nx/2+1 ) :: k_x
 INTEGER*4, DIMENSION(Ny) ::  k_y 
 DOUBLE COMPLEX,  DIMENSION(Nx/2+1) :: ikx
 DOUBLE COMPLEX,  DIMENSION(Ny) :: iky
 DOUBLE COMPLEX, DIMENSION( Nx/2+1 , Ny ) :: Fieldcomplex_c, Fieldcomplex_inc, Fycomplex_c, Fycomplex_inc
 DOUBLE COMPLEX, DIMENSION( Nx/2+1 , Ny ) :: Fxcomplex_c, Fxcomplex_inc
 DOUBLE PRECISION, DIMENSION(Nx, Ny) :: Field_c, Field_inc, Fy_coer, Fy_incoer, Fx_coer, Fx_incoer
 DOUBLE PRECISION, DIMENSION(Nx, Ny, 1) :: Fy_c, Fy_i, Fx_c, Fx_i

 lx = Nx/2 !Punti nello spazio fisico
 ly = Ny
!~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  
! K-vectors definition
!           1 2 3 4 5  6  7  8     
!      kx = 0,1,2..nx/2  
!      ky = 0,1,2,3,4,-3,-2,-1    [For ny = 8]

 
 k_x = 0
 k_y = 0

 DO ii = 1, ly
    IF( ii <= (ly/2) )      k_y(ii) = ii - 1
    IF( ii ==  ((ly/2)+1) ) k_y(ii) = (ly/2)
    IF( ii > ((ly/2)+1) )   k_y(ii) = ii - ly - 1
 END DO

 DO ii = 1, lx + 1
     k_x(ii) = ii - 1
 END DO


 do j=1,Ny
   do i=1,Nx
     Fy_coer(i,j)  =Fy_c(i,j,1)
     Fy_incoer(i,j)=Fy_i(i,j,1)
     Fx_coer(i,j)  =Fx_c(i,j,1)
     Fx_incoer(i,j)=Fx_i(i,j,1)
    enddo
 enddo 

 Fycomplex_c   = DCMPLX(0.d0,0.d0)
 Fycomplex_inc = DCMPLX(0.d0,0.d0)
 Fxcomplex_c   = DCMPLX(0.d0,0.d0)
 Fxcomplex_inc = DCMPLX(0.d0,0.d0)

 Fieldcomplex_c(i,j)   = DCMPLX(0.d0,0.d0)
 Fieldcomplex_inc(i,j) = DCMPLX(0.d0,0.d0)

 CALL realis(1,Nx,Ny,lx,ly,Fycomplex_c,Fy_coer)

 CALL realis(1,Nx,Ny,lx,ly,Fycomplex_inc,Fy_incoer)

 CALL realis(1,Nx,Ny,lx,ly,Fxcomplex_c,Fx_coer)

 CALL realis(1,Nx,Ny,lx,ly,Fxcomplex_inc,Fx_incoer)

 do j=1,ly
   do i=1,lx+1
      ikx(i)=DCMPLX(0.d0,k_x(i))
      iky(j)=DCMPLX(0.d0,-k_y(j))
      if((k_y(j).eq.0).and. (k_x(i).NE.0))then
            Fieldcomplex_c(i,j)   = Fycomplex_c(i,j)   /ikx(i)
            Fieldcomplex_inc(i,j) = Fycomplex_inc(i,j) /ikx(i) 
     
       elseif((k_y(j).eq.0).and. (k_x(i).eq.0))then
           Fieldcomplex_c(i,j)  = DCMPLX(0.d0,0.d0)
           Fieldcomplex_inc(i,j)= DCMPLX(0.d0,0.d0)
       else
           Fieldcomplex_c(i,j)   = Fxcomplex_c(i,j)   /iky(j)
           Fieldcomplex_inc(i,j) = Fxcomplex_inc(i,j) /iky(j) 
      endif
    enddo
 enddo 

 CALL realis(-1,Nx,Ny,lx,ly, Fieldcomplex_c,   Field_c) 
 CALL realis(-1,Nx,Ny,lx,ly, Fieldcomplex_inc, Field_inc) 

 end subroutine divisione_psi_A



!-------------------------------------------------------------------------
  subroutine calcolo_dev_standard(M, nbins, coeff_w, dev_standard, campo)
!-------------------------------------------------------------------------

 implicit none 
 integer :: M, nbins, i, fid
 double precision, dimension(nbins) :: p, w
 double precision, dimension(2048, 2048, 1) :: coeff_w
 double precision :: media 
 double precision :: dev_standard
 character(len=32) :: filename
 character(len=5) :: char_i
 character(len=2) ::  campo


 call  pdf_coeff(M, nbins, coeff_w, w, p, dev_standard, media )

 fid=1000+M
 write(char_i, '(I5)'), M
 filename= 'pdf_' // campo // 'M' //  trim(adjustl(char_i)) // '.dat'
 open(fid, file=filename)

 do i=1, nbins
      write(fid,*), (w(i)-media)/dev_standard, p(i),  w(i)/dev_standard
 enddo

 close(fid)
 end subroutine calcolo_dev_standard

!----------------------------------------------------------------------------------------------
 subroutine divisione_omega_j(Fy_c, Fy_i, Fx_c, Fx_i, Field_c, Field_inc)
 !------------------------------------------------------------------------------------------------
 !CALCOLA LA VORICITA' (O LA CORRENTE) DIVISA IN PARTE COERENTE ED INCOERENTE  A 
 !PARTIRE DAL CAMPO DI VELOCITA' (MAGNETICO) DIVISO USANDO LA SOGLIA DI SIGMA.
 !

 use variabili 
 implicit none
 
 integer :: lx, ly, i ,j, ii 
 INTEGER*4, DIMENSION( Nx/2+1 ) :: k_x
 INTEGER*4, DIMENSION(Ny) ::  k_y 
 DOUBLE COMPLEX,  DIMENSION(Nx/2+1) :: ikx
 DOUBLE COMPLEX,  DIMENSION(Ny) :: iky
 DOUBLE COMPLEX, DIMENSION( Nx/2+1 , Ny ) :: Fieldcomplex_c, Fieldcomplex_inc, Fycomplex_c, Fycomplex_inc
 DOUBLE COMPLEX, DIMENSION( Nx/2+1 , Ny ) :: Fxcomplex_c, Fxcomplex_inc
 DOUBLE PRECISION, DIMENSION(Nx, Ny) :: Field_c, Field_inc, Fy_coer, Fy_incoer, Fx_coer, Fx_incoer
 DOUBLE PRECISION, DIMENSION(Nx, Ny, 1) :: Fy_c, Fy_i, Fx_c, Fx_i

 lx = Nx/2 !Punti nello spazio fisico
 ly = Ny
!~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  
! K-vectors definition
!           1 2 3 4 5  6  7  8     
!      kx = 0,1,2..nx/2  
!      ky = 0,1,2,3,4,-3,-2,-1    [For ny = 8]

 
 k_x = 0
 k_y = 0

 DO ii = 1, ly
    IF( ii <= (ly/2) )      k_y(ii) = ii - 1
    IF( ii ==  ((ly/2)+1) ) k_y(ii) = (ly/2)
    IF( ii > ((ly/2)+1) )   k_y(ii) = ii - ly - 1
 END DO

 DO ii = 1, lx + 1
     k_x(ii) = ii - 1
 END DO


 do j=1,Ny
   do i=1,Nx
     Fy_coer(i,j)  =Fy_c(i,j,1)
     Fy_incoer(i,j)=Fy_i(i,j,1)
     Fx_coer(i,j)  =Fx_c(i,j,1)
     Fx_incoer(i,j)=Fx_i(i,j,1)
    enddo
 enddo 

 Fycomplex_c   = DCMPLX(0.d0,0.d0)
 Fycomplex_inc = DCMPLX(0.d0,0.d0)
 Fxcomplex_c   = DCMPLX(0.d0,0.d0)
 Fxcomplex_inc = DCMPLX(0.d0,0.d0)

 Fieldcomplex_c(i,j)   = DCMPLX(0.d0,0.d0)
 Fieldcomplex_inc(i,j) = DCMPLX(0.d0,0.d0)

 CALL realis(1,Nx,Ny,lx,ly,Fycomplex_c,Fy_coer)

 CALL realis(1,Nx,Ny,lx,ly,Fycomplex_inc,Fy_incoer)

 CALL realis(1,Nx,Ny,lx,ly,Fxcomplex_c,Fx_coer)

 CALL realis(1,Nx,Ny,lx,ly,Fxcomplex_inc,Fx_incoer)

 do j=1,ly
   do i=1,lx+1
      ikx(i)=DCMPLX(0.d0,k_x(i))
      iky(j)=DCMPLX(0.d0,-k_y(j))
      Fieldcomplex_c(i,j)   = ikx(i)*Fycomplex_c(i,j)- iky(j)*Fxcomplex_c(i,j)  
      Fieldcomplex_inc(i,j) = ikx(i)*Fycomplex_inc(i,j)- iky(j)*Fxcomplex_inc(i,j)    

   enddo
 enddo 

 CALL realis(-1,Nx,Ny,lx,ly, Fieldcomplex_c,   Field_c) 
 CALL realis(-1,Nx,Ny,lx,ly, Fieldcomplex_inc, Field_inc) 

 end subroutine divisione_omega_j


!--------------------------------------------------------------
 subroutine correlazione(field1, field2, corr)
!--------------------------------------------------------------

  use variabili 
 implicit none
 
 double precision :: media1, media2, sigma, dev_standard1, dev_standard2, corr
 double precision, dimension(Nx,Ny) :: field1, field2
 
  
 media1= sum(field1) / (Nx*Ny)
 media2= sum(field2) / (Nx*Ny)
 sigma = sum((field1-media1)*(field2-media2))

 dev_standard1 = dsqrt(sum((field1-media1)**2.d0))
 dev_standard2 = dsqrt(sum((field2-media2)**2.d0))

 corr = sigma / (dev_standard1* dev_standard2)
 
 end subroutine correlazione









