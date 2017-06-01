module construct
use globals
use legendre
use auxiliary
contains

subroutine build
    
    stiff_full = 0D0
    stress_full = 0D0
    RHS_full = 0D0
    stiff_reduced = 0D0
    stress_reduced = 0D0
    RHS_reduced = 0D0
    
    call build_Stiff_Full
    call build_Stress_Full
    call build_RHS_Full
    call build_Full2Reduced
    
end subroutine build

subroutine build_Stiff_Full
    integer::i,j,k
    
    do i = 1,nelem
        do j = 1,numNodes(elemType(i))
            do k = 1,numNodes(elemType(i))
                
! Construct the full stiffness matrix, term by term
! The terms of the Stiffness matrix are:
!   1. The stiffness matrix value itself
!   2. The diffusive term in the X direction
!   3. The diffusive term in the Y direction
!   4. The advective term in the X direction
!   5. The advective term in the Y direction
                stiff_full(elem_mat(i,j),elem_mat(i,k)) = &
                    & stiff_full(elem_mat(i,j),elem_mat(i,k))                                                       - &
                    & Dx*quad(j,k,1,1,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))      - &
                    & Dy*quad(j,k,2,2,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))      - &
                    & velx*quad(j,k,0,2,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))    - &
                    & vely*quad(j,k,0,1,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))
    
            end do
        end do
    end do

end subroutine build_Stiff_Full

subroutine build_Stress_Full
    integer::i,j,k
    
    do i = 1,nelem
        do j = 1,numNodes(elemType(i))
            do k = 1,numNodes(elemType(i))

! Construct the full stress matrix, term by term
! The terms of the stress matrix are:
!   1. The stress matrix value itself
!   2. The time dependent term
                stress_full(elem_mat(i,j),elem_mat(i,k)) = &
                    & stress_full(elem_mat(i,j),elem_mat(i,k))          - &
                    & quad(j,k,0,0,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))
    
            end do
        end do
    end do

end subroutine build_Stress_Full

subroutine build_RHS_Full
    integer::i,j,k
    double precision::Q=-1d-1
    
    do i = 1,nelem
        do j = 1,numNodes(elemType(i))
            do k = 1,numNodes(elemType(i))
                    
! Add the source term to the RHS vector
                RHS_full(i) = RHS_full(i) + &
                    & Q*quad(j,k,0,0,xy_coord(elem_mat(i,1:numNodes(elemType(i))),1:2),numNodes(elemType(i)))
    
            end do
        end do
    end do
    
! Add Boundary conditions to the RHS_full vector
! Possible Boundary conditions are:
!   1. Dirchlet Boundary Conditions (BC_type = 0)
!   2. Neumann Boundary Conditions (BC_type = 1)
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

end subroutine build_RHS_Full

subroutine build_Full2Reduced
    integer::i,dum1,dum2
    
    dum1 = 1
    dum2 = nnod
    do i = 1,nnod
        if(BC_type(i).eq.0) then
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


end subroutine build_Full2Reduced

end module construct
