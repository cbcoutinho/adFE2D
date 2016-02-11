module auxiliary
! This module provides some auxiliary functions needed by other modules

contains
function numNodes(elemType) result(num)

! This function returns the number of nodes in an element
! The input element type is based on a GMSH created file, which uses
! the GMSH element type numbering system. For more details:
!   `http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format`

    integer, intent(in)::   elemType
    integer::  num
    
    select case (elemType)
!        case(1)     
!            num = 2    ! 1-D Linear line element with 2 nodes
        case(2)     
            num = 3    ! 2-D Linear triangle element with 3 nodes
        case(3)    
            num = 4    ! 2-D Linear quadrilateral element with 4 nodes
!        case(4)    
!            num = 4    ! 3-D Linear tetrahedral element with 4 nodes
        case default
            write(*,*) "ERROR: Detected unknown element type"
            write(*,*) "Element type = `", elemType, "` is not supported"
            stop
    end select
    
end function numNodes

end module auxiliary
