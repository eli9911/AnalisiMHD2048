 module variabili 

 implicit none 
 
 double precision, dimension(2048,2048,1) :: a,  bx, by, dj, vx, vy, o
 integer :: Nx, Ny , M ,  ndim2 
 double precision:: L , dx, dy, dz, pi
 double precision, allocatable , dimension(:) :: Lscale
 integer, allocatable, dimension(:) :: Ni, Nf

 double precision, allocatable, dimension(:,:) :: sigq, Isigq


 end module variabili
