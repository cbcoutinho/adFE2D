module linsolver
use globals
integer::errflag
integer,dimension(:),allocatable::o
double precision,parameter::epsi1=1D-7
double precision,dimension(:),allocatable::q,z,y
double precision,dimension(:,:),allocatable::L,diag,decomp
logical::matcalc
save
contains

subroutine lud(a,b,x,n)
    integer::n
    double precision,dimension(n)::b
    double precision,dimension(n)::x
    double precision,dimension(n,n)::a,decomp1
    
    errflag2 = 0
    allocate(o(n),q(n))
    
    if(matcalc) then
        call decompose(a,n)
        decomp = a
        o_lud = o
        q_lud = q
        matcalc = .false.
    end if
    decomp1 = decomp
    o = o_lud
    q = q_lud
    if(errflag2 .ne. -1) then
        call substitute(decomp1,b,x,n)
    end if
    
    deallocate(o,q)

end subroutine lud

subroutine decompose(a,n)
    integer::n,i,j,k
    double precision::factor
    double precision,dimension(n,n)::a
    
    do i=1,n
        o(i)=i
        q(i)=abs(a(i,1))
        do j=2,n
            if(abs(a(i,j)) .gt. q(i)) q(i) = abs(a(i,j))
        end do
    end do

    do k = 1,n-1
        call pivot(a,k,n)
        if(abs(a(o(k),k)/q(o(k))) .lt. epsi1) then
            errflag2=-1
!            write(unitnum(9),'(a)')"Error Flag in Decomposition Routine"
!            write(unitnum(9),'(a,f10.3)') "Something? (check book): ",a(o(k),k)/q(o(k))
            exit
        end if
        do i = k+1,n
            factor = a(o(i),k)/a(o(k),k)
            a(o(i),k) = factor
            do j = k+1,n
                a(o(i),j) = a(o(i),j)-factor*a(o(k),j)
            end do
        end do
    end do
 
    if(abs(a(o(k),k)/q(o(k))) .lt. epsi1) then
        errflag2=-1
        write(2,*) a(o(k),k)/q(o(k))
    end if
 
end subroutine decompose

subroutine pivot(a,k,n)
    integer::n,k,p,ii,dummyint
    double precision::big,dummyreal
    double precision,dimension(n,n)::a
    
    p = k
    big = abs(a(o(k),k)/q(o(k)))
    do ii = k+1,n
        dummyreal = abs(a(o(ii),k)/q(o(ii)))
        if(dummyreal.gt.big) then
            big = dummyreal
            p = ii
        end if
    end do
    
    dummyint = o(p)
    o(p) = o(k)
    o(k) = dummyint
 
end subroutine pivot

subroutine substitute(a,b,x,n)
    integer::n,i,j
    double precision::sumb
    double precision,dimension(:),intent(out)::x,b
    double precision,dimension(:,:),intent(in)::a
    
    do i = 2,n
        sumb = b(o(i))
        do j = 1,i-1
            sumb = sumb - a(o(i),j)*b(o(j))
        end do
        b(o(i)) = sumb
    end do
    
    x(n) = b(o(n))/a(o(n),n)
    
    do i = n-1,1,-1
        sumb = 0
        do j = i+1,n
            sumb = sumb+a(o(i),j)*x(j)
        end do
        x(i) = (b(o(i)) - sumb)/a(o(i),i)
    end do

end subroutine substitute

subroutine LDLT(a,b,x,n)
    integer::n
    double precision,dimension(n)::b
    double precision,dimension(n)::x
    double precision,dimension(n,n)::a
    
    errflag2 = 0
    allocate(L(n,n),diag(n,n))
    allocate(y(n),z(n))
    
    L = 0.0D0
    diag = 0.0D0
    y = 0.0D0
    z = 0.0D0
    x = 0.0D0
    
    call assign(a,n)
    call forward_sub(b,n)
    call diag_scale(n)
    call back_sub(x,n)
    
    deallocate(L,diag,y,z)
    
end subroutine LDLT

subroutine assign(a,n)
    integer::i,j,k,n
    double precision,dimension(n,n)::a
    
    do i = 1,n
        do j = 1,i
            if(i.eq.1) then
                diag(i,j) = a(i,j)
                L(i,j) = 1.0D0
            else if(j.eq.1) then
                L(i,j) = a(i,j)/diag(j,j)
            else if(i.eq.j) then
                L(i,j) = 1.0D0
                do k = 1,j-1
                    diag(j,j) = diag(j,j) + (L(j,k)**2.0D0*diag(k,k))
                end do
                diag(j,j) = a(j,j) - diag(j,j)
            else
                do k = 1,j-1
                    L(i,j) = L(i,j) + (L(i,k)*diag(k,k)*L(j,k))
                end do
                L(i,j) = a(i,j) - L(i,j)
                L(i,j) = L(i,j)/diag(j,j)
            end if
        end do
    end do
    
end subroutine assign

subroutine forward_sub(b,n)
    integer::i,k,n
    double precision,dimension(n)::b
    
    do i = 1,n
        if(i.eq.1) then
            z(i) = b(i)
        else
            do k = 1,i-1
                z(i) = L(i,k)*z(k)
            end do
            z(i) = b(i) - z(i)
        end if
    end do
    
end subroutine forward_sub

subroutine diag_scale(n)
    integer::i,n
    
    do i = 1,n
        y(i) = z(i)/diag(i,i)
    end do

end subroutine diag_scale

subroutine back_sub(x,n)
    integer::i,k,n
    double precision,dimension(n)::x
    
    do i=n,1,-1
        if(i.eq.n) then
            x(i) = y(i)
        else
            do k = i+1,n
                x(i) = x(i) + L(k,i)*x(k)
            end do
            x(i) = y(i) - x(i)
        end if
    end do
    
end subroutine back_sub

subroutine inverse(a,c,n)
!============================================================
!
! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
!
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================

    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
    L=0.0D0
    U=0.0D0
    b=0.0D0

! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.0D0
    end do
! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

! Step 3: compute columns of the inverse matrix C
    do k=1,n
        b(k)=1.0D0
        d(1) = b(1)
        
! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do
        
! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
        end do
        
! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
            c(i,k) = x(i)
        end do
        b(k)=0.0D0
    end do
end subroutine inverse

end module linsolver