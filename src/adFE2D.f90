program adFE2D
use input
use globals
use construct
use linsolver
    integer::i,j,t,num1,num2,dum
    integer,dimension(:),allocatable::occur
    double precision::delt,delt_out,dum1,dum2
    double precision,dimension(:),allocatable::xplay,RHS_save
    double precision,dimension(:,:),allocatable::x,inv,xsrc
    
    
    call input_all
    call build
    
    open(99,file='data.out',status='replace')
    open(909,file='scratch.out',status='replace')
    write(99,*) nnod
    
    allocate(x(nrows,2))
    allocate(inv(nrows,nrows))
    allocate(xsrc(nrows,2),xplay(nrows),occur(nrows))
    allocate(decomp(nrows,nrows),o_lud(n),q_lud(n))
    allocate(RHS_save(nrows))
    
    open(104,file='ic.in',status='old')
    
    read(104,*) ! Read blank/header file
    read(104,*) num1 ! Number of Initial Conditions
    read(104,*) ! Read blank/header file
    
    if(num1.lt.1) then
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
    if(num1.gt.1) then
        do i = 2,num1
            read(104,*) dum,x(dum,2)
        end do
    end if
    
    read(104,*) ! Read blank/header file
    read(104,*) num2 ! Number of Sources/Sinks
    read(104,*) ! Read blank/header file
    if(num2.gt.0) then
     occur = 0
        do i = 1,num2
            read(104,*) dum,xplay(dum),xsrc(dum,1),xsrc(dum,2)
            occur(dum) = 1
        end do
    end if
    
    close(104)
    
!        do i = 1,nnod
!            write(99,*) (xy_coord(i,j),j=1,2)
!        end do
!        write(99,*) 
!        do i = 1,nelem
!            write(99,*) (elem_mat(i,j),j = 1,4)
!        end do
!        stop
    
    
!        write(*,*) "Calling inverse"
!        call inverse(stiff_reduced+stress_reduced,inv,nrows)
!        write(*,*) "Finish inverse"


!!!!!!!!!!!!!! Attempt at doing Steady State
!
!    call lud(stiff_reduced, &
!            & -RHS_reduced, &
!            & x(:,2),nrows)
!    
!    
!    stop
!    
!    sumsq = 0d0
!    
!    do i = 1,nnod
!        if(BC_type(i).eq.0) then
!            write(99,100) xy_coord(i,:),BC_value(i)
!        else
!            write(99,100) xy_coord(i,:),x(convert(i),2)
!            sumsq = sumsq + (x(convert(i),2)-x(convert(i),1))**2D0
!        end if
!    end do
!    
!    do i=1,nrows
!        write(909,*) (x(i,j),j=1,2)
!    end do
!
!!!!!!!!!!!!!!!!!!!
    
    numt = 0
    numt_out = 0
    time = 0d0
    delt_out = 0d0
    delt = deltat
    RHS_save = RHS_reduced
    matcalc = .true.
    
    sumsq = 1d0
    
    do while (sumsq.gt.sumlim.and.time.lt.tf)
        RHS_reduced = RHS_save
        numt = numt + 1
        time = time + delt
        delt_out = delt_out + delt
        x(:,1) = x(:,2)
        
        if(num2.gt.0) then
            do i = 1,nrows
                do j = 1,num2
                    if(time.gt.xsrc(i,1)) then
                        if(time.lt.(xsrc(i,1)+xsrc(i,2)).or.occur(i).eq.1) then
                            occur(i) = 0
                            RHS_reduced(i) = RHS_reduced(i)+xplay(i)
                        end if
                    end if
                end do
            end do
        end if
        
!        write(*,*) "Calling LUD"
        call lud(stiff_reduced+stress_reduced, &
            & -(RHS_reduced-matmul(stress_reduced,x(:,1))), &
            & x(:,2),nrows)
!        write(*,*) "Finish LUD"
        
        do i=1,nrows
            write(909,*) (x(i,j),j=1,2)
        end do
            
!        write(*,*) "Calling matmul"
!        x(:,2) = matmul(inv,-(RHS_reduced-matmul(stress_reduced,x(:,1))))
     
        if(delt_out.ge.writetime) then
            sumsq = 0D0
            delt_out = 0d0
            numt_out = numt_out + 1
            
            do i = 1,nnod
                if(BC_type(i).eq.0) then
                    write(99,100) i, xy_coord(i,:), BC_value(i)
                else
                    write(99,100) i, xy_coord(i,:),x(convert(i),2)
                    sumsq = sumsq + (x(convert(i),2)-x(convert(i),1))**2D0
                end if
            end do
            
            sumsq = sqrt(sumsq/nnod)
            write(*,*) numt,numt_out,time,sumsq
            write(909,*) time,sumsq
        end if
     
!        if(sumsq.le.sumlim.or.time.ge.tf) exit
!        if(time.ge.tf) exit
!    if(numt_out.ge.10) stop
    end do

    
    write(99,*) numt_out,nnod
    
    100 format (i4,3(es13.6,1x))
    
    deallocate(BC_type,BC_value)
    deallocate(stiff_full,stiff_reduced)
    deallocate(RHS_full,RHS_reduced,x,inv)
    deallocate(xy_coord,elem_mat,convert)
    deallocate(xsrc,xplay,occur,decomp,o_lud,q_lud)
 
end program adFE2D