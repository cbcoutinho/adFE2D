module input
use globals
contains

subroutine ingest
 integer::i,num,dum
 
 call read_in1
 
 if(option1.eq.1) call inputgen
 
! open(999,file='test.msh',status='old')
! 
! read(999,*) nnod
! allocate(xy_coord(nnod,2))
! do i = 1,nnod
! read(999,*) dum,(xy_coord(i,j),j=1,2)
! end do
! 
! read(999,*) nelem
! allocate(elem_mat(nelem,4))
! do i = 1,nelem
!  read(999,*) dum,(elem_mat(i,j),j = 1,4)
! end do
! 
! close(999)
 
 allocate(BC_type(nnod),BC_value(nnod),convert(nnod))
 allocate(stiff_full(nnod,nnod),stress_full(nnod,nnod))
 allocate(RHS_full(nnod))
 
 BC_type = 1
 BC_value = 0D0
 
 call read_in2
 
 convert = 0
 
 nrows = 0
 do i = 1,nnod
  if(BC_type(i).ne.0) then
   nrows = nrows + 1
   convert(i) = nrows
!   write(*,*) i,convert(i)
  end if
 end do
 
 allocate(stiff_reduced(nrows,nrows),stress_reduced(nrows,nrows))
 allocate(RHS_reduced(nrows))
 
 write(str,100) nnod+1
 
 100 format('(',i5,'(es17.10,1x))')
 
end subroutine ingest

subroutine read_in1
 double precision::length,width,num

 open(100,file='params.in',status='old')
 
 read(100,*) nelemV
 read(100,*) nelemH
 read(100,*) length
 read(100,*) width
 read(100,*) nnod
 read(100,*) nelem
 read(100,*) option1
 read(100,*) deltat
 read(100,*) tf
 read(100,*) num
 read(100,*) Dx
 read(100,*) Dy
 read(100,*) velx
 read(100,*) vely
 read(100,*) sumlim
 
 deltaV = length/nelemV
 deltaH = width/nelemH
! nelem = nelemV*nelemH
 writetime = tf/dble(num)
 
 close(100)
 
end subroutine read_in1

subroutine read_in2
 integer::i,j,dum,num

 open(101,file='elem.in',status='old')
 open(102,file='xycoord.in',status='old')
 open(103,file='bc.in',status='old')
 
 read(101,*) ! Read the header line
 read(102,*) ! Read the header line
 read(103,*) ! Read the header line
 
 read(103,*) num
 read(103,*)
 do i = 1,num
  read(103,*) dum,BC_type(dum),BC_value(dum)
 end do
 
! write(*,*) "Blue Footed Booby"
! stop
 
! do i = 1,nelem
!  read(101,*) dum,(elem_mat(i,j),j = 1,4)
! end do
 
! do i = 1,nnod
!  read(102,*) dum,(xy_coord(i,j),j=1,2)
! end do
 
 close(101)
 close(102)
 close(103)
 
end subroutine read_in2

subroutine inputgen
 integer::i,j,num
 
 open(101,file='elem.in',status='replace')
 open(102,file='xycoord.in',status='replace')
 
 write(101,100)
 write(102,102)
 
 num = 1
 i = 1
 do
  do j = i,nelemV+i-1
   write(101,101) num,j,j+nelemV+1,j+nelemV+1+1,j+1
   num = num+1
   nnod = j+nelemV+1+1
   if(num.gt.nelem)exit
  end do
   if(num.gt.nelem)exit
   i = i+nelemV+1
 end do
 
 num = 0
 do i = 1,nelemH+1
  do j = 1,nelemV+1
   num = num+1
   write(102,103) num,(i-1)*deltaH,(j-1)*deltaV
  end do
 end do
 
 nnod = num
 
 close(101)
 close(102)
 write(*,104)
 stop
 
 100 format ("Elem #",t9,"Node 1",t17,"Node 2",t25,"Node 3",t33,"Node 4")
 101 format (i5,t9,i5,t17,i5,t25,i5,t33,i5)
 102 format ("Node #",t10,"X-coord",t21,"Y-coord")
 103 format (i5,t9,2(es17.10,1x),i5,t41,2(es17.10,1x))
 104 format (/15x,"Input file completed"/ &
            &7x,"Please alter it, assign option1 to 0"/ &
            &16x,"Then Rerun program"/)
 
end subroutine inputgen

end module input