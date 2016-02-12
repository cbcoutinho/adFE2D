module input
use globals
use auxiliary
contains

subroutine input_all
!    integer::i,num,dum
    
    call read_params
    call read_gmsh
    
    allocate(BC_type(nnod),BC_value(nnod),convert(nnod))
    allocate(stiff_full(nnod,nnod),stress_full(nnod,nnod))
    allocate(RHS_full(nnod))
    
    call read_bc
    
    allocate(stiff_reduced(nrows,nrows),stress_reduced(nrows,nrows))
    allocate(RHS_reduced(nrows))
 
end subroutine input_all

subroutine read_params
!    double precision::length,width
    
    open(100,file='params.in',status='old')
    read(100,*) ! Read blank line
    
    read(100,*) deltat
    read(100,*) tf
    read(100,*) writetime
    read(100,*) Dx
    read(100,*) Dy
    read(100,*) velx
    read(100,*) vely
    read(100,*) sumlim
    
    close(100)
    
end subroutine read_params

subroutine read_gmsh
! This function reads a GMSH .msh file and returns:
! 1. The number of nodes
! 2. The number of elements
! 3. The xy coordinates of the elemental nodes
! 4. The element connectivitiy information (i.e. node ordering)

    integer::i,j,dum,nElemOrig,elemTypeInt,numTag
    character(len=32)::filename
    
    
!~~~~~~~~~~~~~ Read filename from input arguement
!~~~
!~~~
!~~~
    if(iargc().lt.1) then
        write(*,*) "No .msh file detected. Input file name to adFE2D"
        stop
    end if
    
    do i = 1,iargc()
        call getarg(i, filename)
        write(*,*) "The gmsh file `", trim(filename), "` is being used"
    end do
!~~~
!~~~
!~~~
!~~~~~~~~~~~~~ End input arguement reading
    
    
    
!~~~~~~~~~~~~~ Read header files in .msh file
!~~~
!~~~
!~~~
    open(101,file=filename,status='old')
    do i=1,4
        read(101,*)
    end do
!~~~
!~~~
!~~~
!~~~~~~~~~~~~~ Done reading header lines
    
    
    
!~~~~~~~~~~~~~ Read node xy coordinates
!~~~
!~~~
!~~~
    read(101,*) nnod
    allocate(xy_coord(nnod,2))
    
    do i=1,nnod
        read(101,*) dum,(xy_coord(i,j),j=1,2)   ! Read each line in .msh file
    end do
!~~~
!~~~
!~~~
!~~~~~~~~~~~~~ Done Reading node xy coordinates

    
!~~~~~~~~~~~~~ Read element connectivitiy information
! Number of elements in .msh files also contain more information than needed
! We only need number of quadrilateral elements, thus node and line elements
! are not needed. For more details, see:
!   `http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format`
!
! The idea is to count number of elements using nElemOrig, read lines until
! element type (column 2) equals 3. The 'i' variable is used to keep track
! of how many non-quadrilateral elements there are, and then subtracts that
! number from nElemOrig to get nelem. Then the rest of the data is read
! normally. This method needs to be changed for future .msh files that contain
! mesh types other than quads.

! Read 2 blank lines
    do i=1,2
        read(101,*)
    end do
    
    read(101,*) nElemOrig
    do while (elemTypeInt.ne.3)
        read(101,*) dum, elemTypeInt
    end do
    backspace(101)
    
    nelem = nElemOrig-(dum-1)
    allocate(elem_mat(nelem,4), elemType(nelem))
    
    do i = 1,nelem
        read(101,*) dum, &                                  ! Original element number
            & elemType(i), &                                ! Element type
            & numTag, &                                     ! Number of tags
            & (dum,j=1,numTag), &                           ! Read tags into dummy variable
            & (elem_mat(i,j),j = 1,numNodes(elemType(i)))   ! Element connectivitiy information
!        write(*,*) (elem_mat(i,j),j = 1,numNodes(elemType(i)))
!        write(*,*) i, elemType(i), numNodes(elemType(i))
    end do
    
    close(101)
    
end subroutine read_gmsh

subroutine read_bc
    integer::i, j, dum, num, coord, dum_bc_type, indx
    integer,dimension(:),allocatable::indx_subset
    double precision::dum_coord_value, dum_bc_value
    double precision, parameter::eps=1d-8
    logical,dimension(:),allocatable::query
    
    BC_type = 1             ! Initializes all node boundaries as natural boundaries
    BC_value = 0D0          ! Sets gradient to zero  (zero neumann)
    
    open(103,file='bc.in',status='old')
    
    read(103,*)             ! Read the header line
    read(103,*) num         ! Number of Unique boundaries
    read(103,*)             ! Read the header line
    
    if (num .ge. 1) then    ! If there is at least one boundary condition, then read
        do i = 1,num
            read(103,*) dum,BC_type(dum),BC_value(dum)
        end do
    end if
    
    read(103,*)             ! Read the header line
    read(103,*) num         ! Number of Coordinate-based boundaries
    read(103,*)             ! Read the header line
    
    if (num .ge. 1) then
        write(*,*) "Number of Coordinate-based boundaries =", num
        
        allocate(query(nnod))
        
        do i = 1,num
            query = .false.
            
            read(103,*) coord, dum_coord_value, dum_bc_type, dum_bc_value
            
            write(*,*) coord, dum_coord_value, dum_bc_type, dum_bc_value
            write(*,*)
            write(*,*)

!~~~~~~~~~~~~~~ Some error checking
!~~~
!~~~
!~~~
            if ((coord.lt.1) .or. (coord.gt.2)) then
                write(*,*) "Coordinate =", coord, " is not supported"
                write(*,*) "Currently only accept X=1 and Y=2 coordinates"
                stop
            end if
            
            if ((dum_bc_type.lt.0) .or. (dum_bc_type.gt.1)) then
                write(*,*) "Boundary Condition type =", coord, " is not supported"
                write(*,*) "Currently only accept Dirchlet=0 and Neumann=1 BCs"
                stop
            end if
            
            if ((dum_coord_value .gt. maxval(xy_coord(:,coord))+eps) .or. &
                & (dum_coord_value .lt. minval(xy_coord(:,coord))-eps)) then
                write(*,*) "Boundary Value =", dum_coord_value, " exceeds minimum or"
                write(*,*) "maximum value of coordinate. Pay attention to geometry bounds"
                stop
            end if
!~~~
!~~~
!~~~
!~~~~~~~~~~~~~ Finish error checking
            
            do j=1,nnod
                query(j) = dum_coord_value > (xy_coord(j,coord) - eps) .and. &
                    & dum_coord_value < (xy_coord(j,coord) + eps)
            end do
            
            write(*,*) pack([(indx,indx=1,nnod)],query)
            write(*,*) 
            write(*,*) 
            
            allocate(indx_subset(count(query)))
            
            indx_subset = pack([(indx,indx=1,nnod)],query)
            
            do j=1,count(query)
                BC_type(indx_subset(j)) = dum_bc_type
                BC_value(indx_subset(j)) = dum_bc_value
            end do
            
            deallocate(indx_subset)
            
        end do
        
        deallocate(query)
    end if
    
    
    
!    value = minval(xy_coord(:,1))               ! Minimum X coordinate
!    write(*,*) pack([(indx,indx=1,nnod)],query)
!    write(*,*) minval(xy_coord,1), maxval(xy_coord,1)
!    deallocate(query)
            
    close(103)

    stop
    
    convert = 0
    nrows = 0
    do i = 1,nnod
        if(BC_type(i).ne.0) then
            nrows = nrows + 1
            convert(i) = nrows
!            write(*,*) i,convert(i)
        end if
    end do
 
    write(str,100) nnod+1
!    write(*,*) trim(str)
    
    100 format('(',i5,'(es17.10,1x))')
    
end subroutine read_bc

subroutine read_ic
    integer::i, dum
    
    open(104,file='ic.in',status='old')
    
    read(104,*) ! Read blank/header file
    read(104,*) num_ic ! Number of Initial Conditions
    read(104,*) ! Read blank/header file
    
    if(num_ic.lt.1) then
        write(*,*) "Must have at least 1 Initial Condition"
        stop
    end if
    
    read(104,*) dum,initial ! First Value should always have dum==0
    
    if(dum.ne.0) then
        write(*,*) "1st Initial Condition must initialize all nodes"
        write(*,*) "Set first `Node #` to zero, to initialize all nodes"
        stop
    end if
    
    x(:,:) = initial
    if(num_ic.gt.1) then
        do i = 2,num_ic
            read(104,*) dum,x(dum,2)
        end do
    end if
    
    read(104,*) ! Read blank/header file
    read(104,*) num_src ! Number of Sources/Sinks
    read(104,*) ! Read blank/header file
    if(num_src.gt.0) then
     occur = 0
        do i = 1,num_src
            read(104,*) dum,xplay(dum),xsrc(dum,1),xsrc(dum,2)
            occur(dum) = 1
        end do
    end if
    
    close(104)
    
end subroutine read_ic

end module input