module globals
 integer::nelem,nelemV,nelemH,nnod,nrows
 integer::option1,numt
 integer,dimension(:),allocatable::BC_type,convert
 integer,dimension(:),allocatable::o_lud
 integer,dimension(:,:),allocatable::elem_mat
 double precision::deltaV,deltaH,initial,deltat,time,sumlim
 double precision::Dx,Dy,velx,vely,tf,writetime
 double precision,dimension(:),allocatable::q_lud
 double precision,dimension(:),allocatable::BC_value,RHS_full,RHS_reduced
 double precision,dimension(:,:),allocatable::stiff_full,stiff_reduced
 double precision,dimension(:,:),allocatable::stress_full,stress_reduced
 double precision,dimension(:,:),allocatable::xy_coord,phi
 character(500)::str
 save
 
end module globals