module legendre
use linsolver
integer,dimension(2)::node,dn
contains

function quad(node1,node2,dn1,dn2,xy)
 integer::ii,jj,kk,node1,node2,dn1,dn2
 double precision,dimension(5)::x,w
 double precision,dimension(2,2)::J,invJ
 double precision,dimension(2,4)::P
 double precision,dimension(4,2)::xy

 x(1) = 0.0D0
 x(2) = 1.0D0/3.0D0*sqrt(5.0D0-2.0D0*sqrt(1.0D1/7.0D0))
 x(3) = -1.0D0/3.0D0*sqrt(5.0D0-2.0D0*sqrt(1.0D1/7.0D0))
 x(4) = 1.0D0/3.0D0*sqrt(5.0D0+2.0D0*sqrt(1.0D1/7.0D0))
 x(5) = -1.0D0/3.0D0*sqrt(5.0D0+2.0D0*sqrt(1.0D1/7.0D0))
 w(1) = 1.28D2/2.25D2
 w(2) = (3.22D2+1.3D1*sqrt(7.0D1))/9.0D2
 w(3) = (3.22D2+1.3D1*sqrt(7.0D1))/9.0D2
 w(4) = (3.22D2-1.3D1*sqrt(7.0D1))/9.0D2
 w(5) = (3.22D2-1.3D1*sqrt(7.0D1))/9.0D2
 
 node(1) = node1
 node(2) = node2
 dn(1) = dn1
 dn(2) = dn2
 
 quad = 0.0D0
 
 do ii = 1,5
  do jj = 1,5
   quad = quad + f(x(ii),x(jj),xy)*w(ii)*w(jj)
  end do
 end do
 
end function quad

function f(zeta,eta,xy)
 integer::ii,jj,kk,m
 double precision::zeta,eta,f
 double precision,dimension(2,2)::J,invJ,tranJ,invtranJ
 double precision,dimension(2,4)::P
 double precision,dimension(4,2)::xy
 double precision,dimension(2)::N
 
 N = 0D0
 
 do jj = 1,4
  P(1,jj) = quad_der_zeta(eta,jj)
  P(2,jj) = quad_der_eta(zeta,jj)
 end do
 
! write(*,*) dn1,dn(1),dn2,dn(2)
! write(*,*)
 
! write(*,'(4(e13.6,1x))') P(1:2,1:4)
! write(*,*)
 
 J = matmul(P,xy)
 tranJ = transpose(J)
 detJ = det(J)
 call inverse(J,invJ,2)
 call inverse(tranJ,invtranJ,2)
 
! write(*,'(2(e13.6,1x))') J(1:2,1:2)
! write(*,*)
! write(*,'(2(e13.6,1x))') tranJ(1:2,1:2)
! write(*,*)
! write(*,*) detJ
! write(*,*)
! write(*,'(2(e13.6,1x))') invJ(1:2,1:2)
! write(*,*)
! write(*,'(2(e13.6,1x))') invtranJ(1:2,1:2)
! write(*,*)
 
 do kk = 1,2
  if(dn(kk).eq.0) then
    N(kk) = quad_bas(zeta,eta,node(kk))
  else
    N(kk) = invJ(dn(kk),1)*quad_der_eta(zeta,node(kk)) + &
          & invJ(dn(kk),2)*quad_der_zeta(eta,node(kk))
  end if
 end do
 
 f = N(1)*N(2)*detJ
 
end function f

function det(a)
 double precision::det
 double precision,dimension(2,2)::a
 
 det = a(1,1)*a(2,2)-a(2,1)*a(1,2)
 if(det.eq.0D0) then
  write(*,'(a)') "Jacobian is Singular, Determinant is equal to 0"
  write(*,'(a)') "Not Good!"
  stop
 end if
 
end function det

function quad_bas(zeta,eta,coord)
integer::coord
double precision::zeta,eta
double precision::quad_bas

if(coord.eq.1) quad_bas = 2.5D-1*(1.0D0-zeta)*(1.0D0-eta)
if(coord.eq.2) quad_bas = 2.5D-1*(1.0D0+zeta)*(1.0D0-eta)
if(coord.eq.3) quad_bas = 2.5D-1*(1.0D0+zeta)*(1.0D0+eta)
if(coord.eq.4) quad_bas = 2.5D-1*(1.0D0-zeta)*(1.0D0+eta)
 
end function quad_bas

function quad_der_zeta(eta,coord)
integer::coord
double precision::eta
double precision::quad_der_zeta

if(coord.eq.1) quad_der_zeta = -2.5D-1*(1.0D0-eta)
if(coord.eq.2) quad_der_zeta = 2.5D-1*(1.0D0-eta)
if(coord.eq.3) quad_der_zeta = 2.5D-1*(1.0D0+eta)
if(coord.eq.4) quad_der_zeta = -2.5D-1*(1.0D0+eta)

end function quad_der_zeta

function quad_der_eta(zeta,coord)
integer::coord
double precision::zeta
double precision::quad_der_eta

if(coord.eq.1) quad_der_eta = -2.5D-1*(1.0D0-zeta)
if(coord.eq.2) quad_der_eta = -2.5D-1*(1.0D0+zeta)
if(coord.eq.3) quad_der_eta = 2.5D-1*(1.0D0+zeta)
if(coord.eq.4) quad_der_eta = 2.5D-1*(1.0D0-zeta)

end function quad_der_eta

end module legendre