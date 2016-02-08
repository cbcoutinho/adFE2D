module input
use globals
contains

subroutine input_all
    integer::i,num,dum
    
    call read_params
    call read_gmsh
    call read_bc
 
end subroutine input_all

subroutine read_params
    double precision::length,width,numWrite
    
    open(100,file='params.in',status='old')
    read(100,*) ! Read blank line
    
    read(100,*) deltat
    read(100,*) tf
    read(100,*) numWrite
    read(100,*) Dx
    read(100,*) Dy
    read(100,*) velx
    read(100,*) vely
    read(100,*) sumlim
    
    writetime = tf/dble(numWrite)
    
    close(100)
 
end subroutine read_params

subroutine read_gmsh
    integer::i,j,dum,nElemOrig,elemType,numTag
    character(len=32)::filename
    
! Read filename from input arguement
    if(iargc().lt.1) then
        write(*,*) "No .msh file detected. Input file name to adFE2D"
        stop
    end if
    
    do i = 1,iargc()
        call getarg(i, filename)
        write(*,*) "The gmsh file `", trim(filename), "` is being used"
    end do
! End input arguement reading

! Read header files in .msh file
    open(101,file=filename,status='old')
    do i=1,4
        read(101,*)
    end do
! Done reading header lines
    
! Read node coordinates
    read(101,*) nnod
    allocate(xy_coord(nnod,2))
    
    do i=1,nnod
        read(101,*) dum,(xy_coord(i,j),j=1,2) ! Read each line in .msh file
    end do
! Done Reading node coordinates

! Read 2 blank lines
    do i=1,2
        read(101,*)
    end do
    
! Read element connectivitiy information
! Number of elements in .msh files also contain more information than needed
! We only need number of quadrilateral elements, thus node and line elements
! are not needed. For more details, see:
!   'http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format'
!
! The idea is to count number of elements using nElemOrig, read lines until
! element type (column 2) equals 3. The 'i' variable is used to keep track
! of how many non-quadrilateral elements there are, and then subtracts that
! number from nElemOrig to get nelem. Then the rest of the data is read
! normally. This method needs to be changed for future .msh files that contain
! mesh types other than quads.

    read(101,*) nElemOrig
    do while (elemType.ne.3)
        read(101,*) dum, elemType
    end do
    backspace(101)
    
    nelem = nElemOrig-(dum-1)
    allocate(elem_mat(nelem,4))
    
    do i = 1,nelem
        read(101,*) dum, elemType, numTag, (dum,j=1,numTag), (elem_mat(i,j),j = 1,4)
!        write(*,*) (elem_mat(i,j),j = 1,4)
    end do
    
    close(101)
    
end subroutine read_gmsh

subroutine read_bc
    integer::i,j,dum,num
    
    allocate(BC_type(nnod),BC_value(nnod),convert(nnod))
    allocate(stiff_full(nnod,nnod),stress_full(nnod,nnod))
    allocate(RHS_full(nnod))
    
    BC_type = 1     ! Initializes all node boundaries as natural boundaries
    BC_value = 0D0  ! Sets gradient to zero  (zero neumann)
    
    open(103,file='bc.in',status='old')
    read(103,*) ! Read the header line
 
    read(103,*) num
    read(103,*) ! Read the header line
    do i = 1,num
     read(103,*) dum,BC_type(dum),BC_value(dum)
    end do
    
    close(103)
    
    convert = 0
    nrows = 0
    do i = 1,nnod
        if(BC_type(i).ne.0) then
            nrows = nrows + 1
            convert(i) = nrows
!            write(*,*) i,convert(i)
        end if
    end do
    
    allocate(stiff_reduced(nrows,nrows),stress_reduced(nrows,nrows))
    allocate(RHS_reduced(nrows))
 
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