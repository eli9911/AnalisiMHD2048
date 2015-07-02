 program mhd 
 
 use variabili 
 
 implicit none 

 integer:: i,j
 double precision :: E, Hc, Bq, vq, Energy
 double precision, dimension(2048, 2048,1) ::  a_b,  psi_v, bx_prova, by_prova, vx_prova, vy_prova, vxCoer, vyCoer, vxIncoer, vyInc
 double precision, dimension(2048, 2048,1) ::  coeff_fvx, coeff_fvy, vxinv, vyinv,coeff_fo, oinv, coeff_fj, jinv
 double precision, dimension(2048, 2048,1) ::  oinv_coer, coeffcoer_o, coeffincoer_o, oinv_incoer
 double precision, dimension(2048, 2048,1) ::  oinv_coer1, coeffcoer_o1, coeffincoer_o1, oinv_incoer1
 double precision, dimension(2048, 2048,1) ::  jinv_coer, coeffcoer_j, coeffincoer_j, jinv_incoer
 double precision, dimension(2048, 2048,1) ::  jinv_coer1, coeffcoer_j1, coeffincoer_j1, jinv_incoer1
 double precision, dimension(11) :: Evx, Evy, k, Etot, sigmaEx, sigmaEy, sigmaT, Eo, ko, sigma_o, Ej, kj, sigma_j, Ecoer_o, kco
 double precision, dimension(11) :: sig_coero, Eincoer_o, kinco, sig_incoero
 integer, dimension(2048, 2048,1):: control_o, control_o1,  control_j, control_j1, control_A
 double precision :: epsilon_t, epsilon_t1, epsilon_j, epsilon_j1, epsilon_A
 integer :: lx, ly, ii, jj 
 double precision :: media_vort, sigma_vort
 double precision, dimension(2048, 2048, 1) :: lambda, lambda_coer, lambda_incoer 
 double precision, dimension(11) :: sig_coerj, Eincoer_j, Ecoer_j, kinco_j, sig_incoerj, kco_j
 double precision, dimension(2048, 2048,1) :: coeff_fpsi, coeff_fA, coeffcoer_A, coeffincoer_A
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
 double precision :: epsilon_vx, epsilon_vy 
 double precision, dimension(11) :: EvxC, EvyC,  EvxI, EvyI, kxc, kyc, kxi, kyi, sigmaxC, sigmayC, sigmaxI, sigmayI
 double precision, allocatable,  dimension(:, :,:) ::  coeffcoer_vx, coeffcoer_vy, coeffincoer_vx, coeffincoer_vy 
 double precision, dimension(2048, 2048,1) :: vx_coerProva, vy_coerProva, vy_incoerProva, vx_incoerProva
 double precision, dimension(11) :: E_Coer, E_Incoer

 !
 double precision, allocatable,  dimension(:, :,:)  :: vx_coerP, vx_incoerP, vy_coerP, vy_incoerP 
 double precision, allocatable, dimension(:) :: kx, ky, ECvx, ECvy, EIvx, EIvy, ECoer,  EIncoer
 

 !
 integer :: nbins
 double precision, allocatable, dimension(:) :: veloc_x, p, w, erro
 !
 integer :: Mmax, Mmin, nbins1, nbins2
 integer ::  t, Np
 double precision :: dev_standard_Sx, media_vyM2, media_vyM1, media_vyM3, media_vyM4, media_vyM5, media_vyM6, media_vyM7
 double precision :: media_vyM8, media_vyM9, media_vyM10, media_vyM11

 double precision :: media_vxM2, media_vxM1, media_vxM3, media_vxM4, media_vxM5, media_vxM6, media_vxM7
 double precision :: media_vxM8, media_vxM9, media_vxM10, media_vxM11

 double precision, allocatable, dimension(:) ::  dev_standard_vx, dev_standard_vy
 double precision, allocatable, dimension(:) ::  p_vxM9, w_vxM9,  p_vxM7, w_vxM7,  p_vxM6, w_vxM6
 double precision, allocatable, dimension(:) ::   p_vxM5, w_vxM5, p_vxM4, w_vxM4, p_vxM3, w_vxM3,  p_vxM1, w_vxM1
 double precision, allocatable, dimension(:) ::  p_vxM11, w_vxM11,  p_vxM10, w_vxM10, p_vxM2, w_vxM2, p_vxM8, w_vxM8

 double precision, allocatable, dimension(:) ::  p_vyM9, w_vyM9,  p_vyM7, w_vyM7,  p_vyM6, w_vyM6
 double precision, allocatable, dimension(:) ::   p_vyM5, w_vyM5, p_vyM4, w_vyM4, p_vyM3, w_vyM3,  p_vyM1, w_vyM1
 double precision, allocatable, dimension(:) ::  p_vyM11, w_vyM11,  p_vyM10, w_vyM10, p_vyM2, w_vyM2, p_vyM8, w_vyM8

 !
 DOUBLE COMPLEX, ALLOCATABLE, DIMENSION( : , : ) :: term1_c, djx_c, djy_c,  term3_c
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : ) :: term1, djx, djy, term3
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : , : , : ) :: dissipazione_vort, dissipazione_pot, combinazione, grad_jx, grad_jy
 
 !
 DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: phi
 double precision, allocatable, dimension(:) ::EI_A, EC_A, k_A, sigma_A, E_A 
 double precision, allocatable, dimension(:,:,:) :: coeff_lambda
 double precision, allocatable, dimension(:) :: EI_lambda, EC_lambda, k_lambda, sigma_lambda, E_lambda
 double precision, allocatable, dimension(:) :: EvxC_prova,  sigma_provaC
 double precision, allocatable, dimension(:) :: p_Sx, w_Sx

 !VARIABILI DIVISIONE PSI IN COER ED INCOER CON VX E SOGLIA 3 SIGMA
 DOUBLE PRECISION,ALLOCATABLE, DIMENSION( : , : ) :: psi_c, psi_inc, vy_coer3s, vy_incoer3s, vx_coer3s, vx_incoer3s

 !CORRELAZIONE J O
 double precision ::  media_j, sigma_oj, dev_standard_j, dev_standard_o, correlazione_oj 
 
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
 print*, 'La media del campo di vorticità è: ', media_vort
 print*, 'Il sigma della vorticità è: ', sigma_vort
 !CALCOLO DEI COEFF WAVELET DEI CAMPI 

 call init_haar
 coeff_fvx = 0.d0
 coeff_fvy = 0.d0
 coeff_fo = 0.d0
 coeff_fj = 0.d0 
 coeff_fpsi= 0.d0
 coeff_fA=0.d0
 sig_coero=0.d0

 !COEFF vx
 call directhaar(vx,coeff_fvx)
 call inverse(coeff_fvx,vxinv)
 print*, maxval(dabs(vx-vxinv))
 
 !COEFF vy
 call directhaar(vy,coeff_fvy)
 call inverse(coeff_fvy,vyinv)
 print*, maxval(dabs(vy-vyinv))


 !COEFF VORTICITA' 
 call directhaar(o,coeff_fo)
 call inverse(coeff_fo,oinv)
 print*, maxval(o-oinv)

 !COEFF CORRENTE 
 call directhaar(dj,coeff_fj)
 call inverse(coeff_fj,jinv)
 print*, maxval(dabs(dj-jinv))


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

 print*,'Teorema di Parseval per la velocità lungo x '
 print*,'Spazio reale:   ', sum(vx*vx)
 print*,'Spazio wavelet: ', sum(coeff_fvx*coeff_fvx)

 print*,'Teorema di Parseval per la velocità lungo y '
 print*,'Spazio reale:   ', sum(vy*vy)
 print*,'Spazio wavelet: ', sum(coeff_fvy*coeff_fvy)

 print*,'Teorema di Parseval per la vorticità '
 print*,'Spazio reale:   ', sum(o*o)
 print*,'Spazio wavelet: ', sum(coeff_fo*coeff_fo)

 print*,'Teorema di Parseval per la corrente '
 print*,'Spazio reale:   ', sum(dj*dj)
 print*,'Spazio wavelet: ', sum(coeff_fj*coeff_fj)  

 
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
 
 epsilon_t1=(2.d0* sum(oinv_incoer**2.d0)/(1.d0*Nx*Ny) * dlog(1.d0*Nx))**0.5d0
 
 call coeff_coerenti(coeff_fo, o, control_o1, epsilon_t1)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_o1(i,j,1) = control_o1(i,j,1)*coeff_fo(i,j,1)
      coeffincoer_o1(i,j,1)= coeff_fo(i,j,1)*(1-control_o1(i,j,1))
    enddo
 enddo

 call inverse(coeffcoer_o1, oinv_coer1)

  open(unit=171, status='UNKNOWN', file='oinv_coer1.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(171,*), i, j, oinv_coer1(i,j,1)
   enddo
  write(171,*), ''
 enddo
 close(171)

 call inverse(coeffincoer_o1, oinv_incoer1)
 
 open(unit=711, status='UNKNOWN', file='oinv_incoer1.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(711,*), i, j, oinv_incoer1(i,j,1)
   enddo
  write(711,*), ''
 enddo
 
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

 print*, 'Errore sul calcolo di fourier del potenziale vettore ', maxval(dabs(a-a_b))

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

 print*,'Errore sul calcolo di fourier del campo magnetico lungo x',  maxval(abs(bx_prova-bx))
 print*,'Errore sul calcolo di fourier del campo magnetico lungo y',  maxval(abs(by_prova-by))

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

 print*,'Errore sul calcolo di fourier del campo di velocità lungo x',  maxval(abs(vx_prova-vx))
 print*,'Errore sul calcolo di fourier del campo di velocità lungo y',  maxval(abs(vy_prova-vy))

 !PDF vx ALLE DIVERSE SCALE

 allocate(dev_standard_vx(M))

!----------------------------------------------
 !M=11
  
 allocate(p_vxM11(1), w_vxM11(1))

 call  pdf_coeff(11, 1, coeff_fvx, w_vxM11, p_vxM11, dev_standard_vx(11), media_vxM11)

 
 open(unit=41140, status='unknown', file='pdf_vxM11.dat', recl=100000)
 write(41140,*),  w_vxM11(1)/dev_standard_vx(11), p_vxM11(1)
 close(41140)

 !----------------------------------------------
 !M=10
  
 allocate(p_vxM10(2), w_vxM10(2))

 call  pdf_coeff(10, 2, coeff_fvx, w_vxM10, p_vxM10, dev_standard_vx(10), media_vxM10)

 open(unit=41040, status='unknown', file='pdf_vxM10.dat', recl=100000)
 do i=1,2
      write(41040,*),  (w_vxM10(i)-media_vxM10)/dev_standard_vx(10), p_vxM10(i)
 enddo
 close(41040)
 

 

 !-----------------------------------------------
 ! M=8

 allocate(p_vxM8(10), w_vxM8(10))

 call  pdf_coeff(8, 10, coeff_fvx, w_vxM8, p_vxM8, dev_standard_vx(8), media_vxM8)

 open(unit=4444, status='unknown', file='pdf_vxM8.dat', recl=100000)
 do i=1,10
      write(4444,*), (w_vxM8(i)-media_vxM8)/dev_standard_vx(8), p_vxM8(i),  w_vxM8(i)/dev_standard_vx(8)
 enddo
 close(4444)

 

 !----------------------------------------------
 !M=9
 
 allocate(p_vxM9(6), w_vxM9(6))

 call  pdf_coeff(9, 6, coeff_fvx, w_vxM9, p_vxM9, dev_standard_vx(9), media_vxM9) 
 
 open(unit=4440, status='unknown', file='pdf_vxM9.dat', recl=100000)
 do i=1,6
      write(4440,*), (w_vxM9(i)-media_vxM9)/dev_standard_vx(9), p_vxM9(i), w_vxM9(i)/dev_standard_vx(9)
 enddo
 close(4440)


 !----------------------------------------------
 !M=7
  
 allocate(p_vxM7(10), w_vxM7(10))
 
 call  pdf_coeff(7, 10, coeff_fvx, w_vxM7, p_vxM7, dev_standard_vx(7), media_vxM7) 
 
 open(unit=4740, status='unknown', file='pdf_vxM7.dat', recl=100000)
 do i=1,10
      write(4740,*), (w_vxM7(i)- media_vxM7)/dev_standard_vx(7), p_vxM7(i),  w_vxM7(i)/dev_standard_vx(7)
 enddo
 close(4740)

 !----------------------------------------------
 !M=6
  
 allocate(p_vxM6(26), w_vxM6(26))
 
 call  pdf_coeff(6, 26, coeff_fvx, w_vxM6, p_vxM6, dev_standard_vx(6), media_vxM6) 
 
 open(unit=4640, status='unknown', file='pdf_vxM6.dat', recl=100000)
 do i=1,26
      write(4640,*), (w_vxM6(i)-media_vxM6)/dev_standard_vx(6), p_vxM6(i), w_vxM6(i)/dev_standard_vx(6)
 enddo
 close(4640)


 !----------------------------------------------
 !M=5
 
 allocate(p_vxM5(26), w_vxM5(26))
 call  pdf_coeff(5, 26, coeff_fvx, w_vxM5, p_vxM5, dev_standard_vx(5), media_vxM5) 
 
 open(unit=4540, status='unknown', file='pdf_vxM5.dat', recl=100000)
 do i=1,26
      write(4540,*), (w_vxM5(i)-media_vxM5)/dev_standard_vx(5), p_vxM5(i), w_vxM5(i)/dev_standard_vx(5)
 enddo
 close(4540)


  !----------------------------------------------
 !M=4
  
 allocate(p_vxM4(30), w_vxM4(30))
 
 call  pdf_coeff(4, 30, coeff_fvx, w_vxM4, p_vxM4, dev_standard_vx(4), media_vxM4) 

 open(unit=44440, status='unknown', file='pdf_vxM4.dat', recl=100000)
 do i=1,30
      write(44440,*), (w_vxM4(i)-media_vxM3)/dev_standard_vx(4), p_vxM4(i), w_vxM4(i)/dev_standard_vx(4)
 enddo
 close(44440)


 !----------------------------------------------
 !M=3
  
 allocate(p_vxM3(50), w_vxM3(50))

 call  pdf_coeff(3, 50, coeff_fvx, w_vxM3, p_vxM3, dev_standard_vx(3), media_vxM3)

 open(unit=4340, status='unknown', file='pdf_vxM3.dat', recl=100000)
 do i=1,50
      write(4340,*), (w_vxM3(i)- media_vxM3)/dev_standard_vx(3), p_vxM3(i), w_vxM3(i)/dev_standard_vx(3)
 enddo
 close(4340)
 

 !--------------------------------------------------------------
 !M=2

 allocate(p_vxM2(100), w_vxM2(100))

 call  pdf_coeff(2, 100, coeff_fvx, w_vxM2, p_vxM2, dev_standard_vx(2), media_vxM2)

 open(unit=2222, status='unknown', file='pdf_vxM2.dat', recl=100000)
 do i=1,100
     write(2222,*), (w_vxM2(i)-media_vxM2)/dev_standard_vx(2), p_vxM2(i),  w_vxM2(i)/dev_standard_vx(2)
 enddo
 close(2222)


!----------------------------------------------
 !M=1
  
 allocate(p_vxM1(100), w_vxM1(100))

 call  pdf_coeff(1, 100, coeff_fvx, w_vxM1, p_vxM1, dev_standard_vx(1), media_vxM1)

 open(unit=4140, status='unknown', file='pdf_vxM1.dat', recl=100000)
 do i=1,100
      write(4140,*), (w_vxM1(i)-media_vxM1)/dev_standard_vx(1), p_vxM1(i), w_vxM1(i)/dev_standard_vx(1)
 enddo
 close(4140)
 

 !PDF vy ALLE DIVERSE SCALE

 allocate(dev_standard_vy(M))

!----------------------------------------------
 !M=11
  
 allocate(p_vyM11(1), w_vyM11(1))

 call  pdf_coeff(11, 1, coeff_fvy, w_vyM11, p_vyM11, dev_standard_vy(11), media_vyM11)
 
 open(unit=4114, status='unknown', file='pdf_vyM11.dat', recl=100000)
 write(4114,*), w_vyM11(1)/dev_standard_vy(11), p_vyM11(1)
 close(4114)

 !----------------------------------------------
 !M=10
  
 allocate(p_vyM10(2), w_vyM10(2))

 call  pdf_coeff(10, 2, coeff_fvy, w_vyM10, p_vyM10, dev_standard_vy(10), media_vyM10)

 open(unit=4104, status='unknown', file='pdf_vyM10.dat', recl=100000)
 do i=1,2
      write(4104,*), (w_vyM10(i)- media_vyM10)/dev_standard_vy(10), p_vyM10(i)
 enddo
 close(4104)
 

 

 !-----------------------------------------------
 ! M=8

 allocate(p_vyM8(30), w_vyM8(30))

 call  pdf_coeff(8, 30, coeff_fvy, w_vyM8, p_vyM8, dev_standard_vy(8), media_vyM8)
 
 open(unit=444, status='unknown', file='pdf_vyM8.dat', recl=100000)
 do i=1,30
      write(444,*), (w_vyM8(i)- media_vyM8)/dev_standard_vy(8), p_vyM8(i), w_vyM8(i)/dev_standard_vy(8)
 enddo
 close(444)

 

 !----------------------------------------------
 !M=9
 
 allocate(p_vyM9(6), w_vyM9(6))

 call  pdf_coeff(9, 6, coeff_fvy, w_vyM9, p_vyM9, dev_standard_vy(9), media_vyM9) 
 
 open(unit=440, status='unknown', file='pdf_vyM9.dat', recl=100000)
 do i=1,6
      write(440,*), (w_vyM9(i) - media_vyM9)/dev_standard_vy(9), p_vyM9(i), w_vyM9(i)/dev_standard_vy(9)
 enddo
 close(440)


 !----------------------------------------------
 !M=7
  
 allocate(p_vyM7(40), w_vyM7(40))
 
 call  pdf_coeff(7, 40, coeff_fvy, w_vyM7, p_vyM7, dev_standard_vy(7), media_vyM7) 
 
 open(unit=740, status='unknown', file='pdf_vyM7.dat', recl=100000)
 do i=1,40
      write(740,*), (w_vyM7(i)-media_vyM7)/dev_standard_vy(7), p_vyM7(i), w_vyM7(i)/dev_standard_vy(7)
 enddo
 close(740)

 !----------------------------------------------
 !M=6
  
 allocate(p_vyM6(60), w_vyM6(60))
 
 call  pdf_coeff(6, 60, coeff_fvy, w_vyM6, p_vyM6, dev_standard_vy(6), media_vyM6) 
 
 open(unit=640, status='unknown', file='pdf_vyM6.dat', recl=100000)
 do i=1,60
      write(640,*), (w_vyM6(i)-media_vyM6)/dev_standard_vy(6), p_vyM6(i), w_vyM6(i)/dev_standard_vy(6)
 enddo
 close(640)


 !----------------------------------------------
 !M=5
 
 allocate(p_vyM5(70), w_vyM5(70))
 call  pdf_coeff(5, 70, coeff_fvy, w_vyM5, p_vyM5, dev_standard_vy(5), media_vyM5) 

 open(unit=540, status='unknown', file='pdf_vyM5.dat', recl=100000)
 do i=1,70
      write(540,*), (w_vyM5(i)-media_vyM5)/dev_standard_vy(5), p_vyM5(i),  w_vyM5(i)/dev_standard_vy(5)
 enddo
 close(540)


  !----------------------------------------------
 !M=4
  
 allocate(p_vyM4(80), w_vyM4(80))
 
 call  pdf_coeff(4, 80, coeff_fvy, w_vyM4, p_vyM4, dev_standard_vy(4), media_vyM4) 

 open(unit=444340, status='unknown', file='pdf_vyM4.dat', recl=100000)
 do i=1,80
      write(444340,*), (w_vyM4(i)-media_vyM4)/dev_standard_vy(4), p_vyM4(i), w_vyM4(i)/dev_standard_vy(4)
 enddo
 close(444340)


 !----------------------------------------------
 !M=3
  
 allocate(p_vyM3(90), w_vyM3(90))

 call  pdf_coeff(3, 90, coeff_fvy, w_vyM3, p_vyM3, dev_standard_vy(3), media_vyM3 )

 !dev_standard_vy(3)= dsqrt( sum((w_vyM3- sum(w_vyM3)/(90.d0))**2.d0) /(90.d0))
 
 open(unit=43340, status='unknown', file='pdf_vyM3.dat', recl=100000)
 do i=1,90
      write(43340,*), (w_vyM3(i)-media_vyM3)/dev_standard_vy(3), p_vyM3(i),  w_vyM3(i)/dev_standard_vy(3)
 enddo
 close(43340)
 

 !--------------------------------------------------------------
 !M=2

 allocate(p_vyM2(100), w_vyM2(100))

 call  pdf_coeff(2, 100, coeff_fvy, w_vyM2, p_vyM2, dev_standard_vy(2), media_vyM2)

 !dev_standard_vy(2)= dsqrt( sum(p_vyM2*(w_vyM2- sum(w_vyM2)/(100.d0))**2.d0) /(100.d0))
 
 open(unit=222, status='unknown', file='pdf_vyM2.dat', recl=100000)
 do i=1,100
     write(222,*), (w_vyM2(i)-media_vyM2)/dev_standard_vy(2), p_vyM2(i), w_vyM2(i)/dev_standard_vy(2)
 enddo
 close(222)


!----------------------------------------------
 !M=1
  
 allocate(p_vyM1(150), w_vyM1(150))

 call  pdf_coeff(1, 150, coeff_fvy, w_vyM1, p_vyM1, dev_standard_vy(1), media_vyM1)

!dev_standard_vy(1)= dsqrt( sum((p_vyM1*w_vyM1- sum(w_vyM1*p_vyM1))**2.d0) )
 
 open(unit=41240, status='unknown', file='pdf_vyM1.dat', recl=100000)
 do i=1,150
      write(41240,*), (w_vyM1(i)-media_vyM1)/dev_standard_vy(1), p_vyM1(i), w_vyM1(i)/dev_standard_vy(1)
 enddo
 close(41240)


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
 
 print*,'controllo vx',  maxval(abs(vx- vx_coerP- vx_incoerP))
 print*, 'controllo vy', maxval(abs(vy- vy_coerP- vy_incoerP))

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

 print*, 'prova sulla divisione di psi: ', maxval(abs(psi_v- psi_coer - psi_incoer))
 print*, 'media psi_coer:', sum(psi_coer)/(Nx*Ny), 'media psi_incoer:', sum(psi_incoer)/(Nx*Ny)
 print*, 'maxval psi_coer:', maxval(psi_coer), 'minval psi_coer:', minval(psi_coer)
 print*, 'maxval psi_incoer:', maxval(psi_incoer), 'minval psi_incoer:', minval(psi_incoer)

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
      
 
 !DIVISIONE DEL POTENZIALE VETTORE IN PARTE COERENTE ED INCOERENTE 
 
 epsilon_A = (2.d0 * sum(a*a)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0

 call coeff_coerenti(coeff_fA, a, control_A, epsilon_A)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_A(i,j,1) = control_A(i,j,1)*coeff_fA(i,j,1)
      coeffincoer_A(i,j,1)= coeff_fA(i,j,1)*(1-control_A(i,j,1))
    enddo
 enddo

 !USO DELLA SOGLIA 3 SIGMA 
 allocate(EI_A(11), EC_A(11), k_A(11), sigma_A(11), E_A(11)) 

 call spettro(coeff_fA,E_A, k_A, sigma_A)
 call Divisione_campo(coeff_fA,  sigma_A, A_coer, A_incoer, EI_A, EC_A, k_A, 3.d0, .false.)

 ! call inverse(coeffcoer_A, A_coer)

 open(unit=87, status='UNKNOWN', file='A_coer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(87,*), i, j, A_coer(i,j,1)
   enddo
  write(87,*), ''
 enddo
 close(87)

 !call inverse(coeffincoer_A, A_incoer)
 
 open(unit=78, status='UNKNOWN', file='A_incoer.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(78,*), i, j, A_incoer(i,j,1)
   enddo
  write(78,*), ''
 enddo
 close(78)

 open(unit=2031, status='UNKNOWN', file='spettroW_A.dat',recl=10000)
 do j=1,11
     write(2031,*), k_A(j), EI_A(j), EC_A(j), E_A(j) 
 enddo
 close(2031)

 open(unit=8531, status='UNKNOWN', file='controllo_A.dat',recl=10000)
  do j= 1,Ny
   do i = 1,Nx 
       write(8531,*), i, j, A_incoer(i,j,1) + A_coer(i,j,1), a(i,j,1)
        write(8531,*),  abs(a(i,j,1)-A_incoer(i,j,1)-A_coer(i,j,1))
   enddo
 enddo
 close(8531)

 !CALCOLO DELLA VELOCITA' COERENTE ED INCOERENTE MEDIANTE PSI 

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
 allocate(EvxC_prova(11),  sigma_provaC(11) )
 k=0.d0
 call directhaar(vxCoer,coeff_fvxC)
 call directhaar(vyCoer,coeff_fvyC)
 call directhaar(vxIncoer,coeff_fvxI)
 call directhaar(vyInc,coeff_fvyI)

 call spettro(coeff_fvxC,Evx_Coer,k, sigmaExC)
 call spettro(coeff_fvyC,Evy_Coer,k, sigmaEyC)
 call spettro(coeff_fvxI,Evx_Inc,k, sigmaExC)
 call spettro(coeff_fvyI,Evy_Inc,k, sigmaEyC)
 call spettro(coeff_fvxC+coeff_fvyC, EvxC_prova, k, sigma_provaC)

 open(unit=991, status='UNKNOWN', file='spettro_velocitaCoerInc.dat',recl=10000)
 
 do i=1,M
   Etot_Coer(i)=Evx_Coer(i)+Evy_Coer(i)
   Etot_Inc(i)=Evx_Inc(i)+Evy_Inc(i)
   !sigmaT(i)= sigmaEx(i)+ sigmaEx(i)
   write(991 ,*), k(i), Etot_Coer(i), Etot_Inc(i),  EvxC_prova(i), Etot(i)
 enddo 
 close(991)
 
 

 !PROVA PER CALCOLARE LA VELOCITA' COERENTE ED INCOERENTE IN MODO DIRETTO

 allocate(coeffcoer_vx(2048,2048,1), coeffcoer_vy(2048,2048,1), coeffincoer_vx(2048,2048,1), coeffincoer_vy(2048,2048,1))

 epsilon_vx = (2.d0 * sum(vx**2.d0)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0
 epsilon_vy = (2.d0 * sum(vy**2.d0)/(1.d0*Nx*Ny) *dlog(Nx*1.d0) )**0.5d0

 call coeff_coerenti(coeff_fvx, vx, control_vx, epsilon_vx)

 call coeff_coerenti(coeff_fvy, vy, control_vy, epsilon_vy)

 do j=1,Ny
   do i =1,Nx
      coeffcoer_vx(i,j,1) = control_vx(i,j,1)*coeff_fvx(i,j,1)
      coeffincoer_vx(i,j,1)= coeff_fvx(i,j,1)*(1-control_vx(i,j,1))

      coeffcoer_vy(i,j,1) = control_vy(i,j,1)*coeff_fvy(i,j,1)
      coeffincoer_vy(i,j,1)= coeff_fvy(i,j,1)*(1-control_vy(i,j,1))
    enddo
 enddo

 call inverse(coeffcoer_vx, vx_coerProva)
 call inverse(coeffcoer_vy, vy_coerProva)


  call inverse(coeffincoer_vx, vx_incoerProva)
  call inverse(coeffincoer_vy, vy_incoerProva)
 
 
 
  call spettro(coeffcoer_vx,EvxC, k, sigmaxC)
  call spettro(coeffcoer_vy,EvyC, k, sigmayC)

  call spettro(coeffincoer_vx,EvxI, kxi, sigmaxI)
  call spettro(coeffincoer_vy,EvyI, kyi, sigmayI)

 open(unit=883, status='UNKNOWN', file='spettro_velocitaProva.dat',recl=10000)
 
 do i=1,M
   E_Coer(i)=EvxC(i)+EvyC(i)
   E_InCoer(i)=EvxI(i)+EvyI(i)
   sigmaT(i)= sigmaEx(i)+ sigmaEx(i)
   write(883 ,*), k(i), E_Coer(i), kyi(i), E_InCoer(i)
 enddo 
 close(883)


 !CONTROLLO SULLE PARTI COERENTI E INCOERENTI DELLA VELOCITA'

 open(unit=445, status='UNKNOWN', file='confronto_vxC.dat', recl=10000)

 do j=1,Ny
   do i=1,Nx
     write(445,*), dabs(coeff_fvxc(i,j,1)-coeffcoer_vx(i,j,1)), dabs(coeff_fvyc(i,j,1)-coeffcoer_vy(i,j,1))
    enddo
 enddo
 close(445)
 

 !CALCOLO DI LAMBDA E DELLA  PARTE COERENTE E INCOERENTE 

 open(unit=333, status='UNKNOWN', file='lambda.dat',recl=10000)

 do j=1,Ny
   do i= 1,Nx
      lambda(i,j,1)= o(i,j,1)/dj(i,j,1)
      write(333,*), i, j, lambda(i,j,1)
    enddo
   write(333,*), ''
 enddo 
 close(333)
 

 !USO DEI TRE SIGMA !QUESTO E' SBAGLIATO

 allocate(coeff_lambda(2048,2048,1))
 allocate(EI_lambda(11), EC_lambda(11), k_lambda(11), sigma_lambda(11), E_lambda(11)) 

 call directhaar(lambda,coeff_lambda)
 call spettro(coeff_lambda,E_lambda, k_A, sigma_lambda)
 call Divisione_campo(coeff_lambda,  sigma_lambda, lambda_coer,lambda_incoer,EI_lambda,EC_lambda,k_lambda,3.d0, .false.)

 open(unit=555, status='UNKNOWN', file='lambdacoer.dat',recl=10000)
 open(unit=566, status='UNKNOWN', file='lambdaincoer.dat',recl=10000)

 do j=1,Ny
   do i= 1,Nx
      write(555,*), i, j, lambda_coer(i,j,1)
      write(566,*), i, j, lambda_incoer(i,j,1)
    enddo
   write(555,*), ''
   write(566,*), ''
 enddo 
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
 
 open(unit=6900, status='unknown', file='phi.dat', recl=100000)
 do j=1,Nx
   do i =1, Nx
     if(a(i,j,1) .NE. 0.d0 ) phi(i,j,1)=  0.5d0*(lambda(i,j,1)*o(i,j,1)- dj(i,j,1) )/ a(i,j,1)
     write(6900,*), i, j, phi(i,j,1)
    enddo
    write(6900,*), ''
 enddo



 !CALCOLO CORRELAZIONE FRA OMEGA E J
 media_j = sum(dj) / (Nx*Ny)
 sigma_oj = sum(dj*o)/(Nx*Ny) - media_j * media_vort

 dev_standard_j = dsqrt(sum((dj)**2.d0)/ (Nx*Ny) - media_j * media_j)
 dev_standard_o = dsqrt(sum((o)**2.d0)/ (Nx*Ny) - media_vort * media_vort)

 correlazione_oj = sigma_oj / (dev_standard_j* dev_standard_o)

 open(unit=69400, status='unknown', file='correlazione.dat', recl=100000)
 write(69400,*), 'correlazione tra vorticità e corrente:', correlazione_oj 
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
 print*, 'area', M_s, ':', sum(p_M)
 
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















