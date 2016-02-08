module construct
use globals
use legendre
contains

subroutine build
    integer::i,j,k,elem,dum1,dum2
    double precision::test
    
    stiff_full = 0D0
    stress_full = 0D0
    RHS_full = 0D0
    stiff_reduced = 0D0
    stress_reduced = 0D0
    RHS_reduced = 0D0
    
    do i = 1,nelem
        do j = 1,4
            do k = 1,4
                stiff_full(elem_mat(i,j),elem_mat(i,k)) = &
                    & stiff_full(elem_mat(i,j),elem_mat(i,k)) - &
                    & Dx*quad(j,k,1,1,xy_coord(elem_mat(i,1:4),1:2)) - &
                    & Dy*quad(j,k,2,2,xy_coord(elem_mat(i,1:4),1:2)) - &
                    & velx*quad(j,k,0,1,xy_coord(elem_mat(i,1:4),1:2)) - &
                    & vely*quad(j,k,0,2,xy_coord(elem_mat(i,1:4),1:2))
                stress_full(elem_mat(i,j),elem_mat(i,k)) = &
                    & stress_full(elem_mat(i,j),elem_mat(i,k)) - &
                    & quad(j,k,0,0,xy_coord(elem_mat(i,1:4),1:2))
            end do
        end do
    end do
    
    stress_full = stress_full/deltat
    
    do i = 1,nnod
        if(BC_type(i).eq.0) then
            do j = 1,nnod
                if(i.ne.j) then
                    RHS_full(j) = RHS_full(j) + BC_value(i)*stiff_full(j,i)
                end if
            end do
        else if(BC_type(i).eq.1) then
            RHS_full(i) = RHS_full(i) + BC_value(i)
        end if
    end do
 
!    do i = 1,nnod
!        write(103,str) BC_value(i)
!    end do
!    write(103,*)
     
!    do i = 1,nnod
!        write(103,str) stiff_full(i,:),RHS_full(i)
!    end do
!    write(103,*)
    
    dum1 = 1
    dum2 = nnod
    do j = 1,nnod    
        if(BC_type(j).eq.0) then
            stiff_full(dum1:nnod-1,:) = stiff_full(dum1+1:nnod,:)
            stiff_full(:,dum1:nnod-1) = stiff_full(:,dum1+1:nnod)
            stress_full(dum1:nnod-1,:) = stress_full(dum1+1:nnod,:)
            stress_full(:,dum1:nnod-1) = stress_full(:,dum1+1:nnod)
            RHS_full(dum1:nnod-1) = RHS_full(dum1+1:nnod)
            dum2 = dum2 - 1
        else
            dum1 = dum1 + 1
        end if
    end do
    
    stiff_reduced(1:nrows,1:nrows) = stiff_full(1:nrows,1:nrows)
    stress_reduced(1:nrows,1:nrows) = stress_full(1:nrows,1:nrows)
    RHS_reduced(1:nrows) = RHS_full(1:nrows)
    
!    do i = 1,nrows
!        write(103,str) stiff_reduced(i,:),RHS_reduced(i)
!    end do
!    write(103,*)
    
    100 format(2(e13.6,1x))
    101 format(4(i3,1x))
    
end subroutine build

end module construct