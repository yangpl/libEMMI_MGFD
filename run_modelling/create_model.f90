program create_model
  implicit none

  integer :: i, j, k, kk
  integer :: n1, n2, n3
  real :: d1, d2, d3
  real :: o1, o2, o3
  real, dimension(:), allocatable :: x1, x2, x3
  real, dimension(:,:,:), allocatable :: rho
  real, dimension(:,:), allocatable :: bathy
  real :: x, y, z, val
  
  n1 = 200
  n2 = 200
  n3 = 100
  d1 = 100
  d2 = 100
  d3 = 40
  
  allocate(x1(n1+1))
  allocate(x2(n2+1))
  allocate(x3(n3+1))
  allocate(rho(n1, n2, n3))
  allocate(bathy(n1+1, n2+1))
  
  o1 = -n1*d1/2
  o2 = -n2*d2/2
  o3 = 0
  do i=1,n1+1
     x1(i) = (i-1)*d1 + o1
  end do
  do j=1,n2+1
     x2(j) = (j-1)*d2 + o2
  end do
  do k=1,n3+1
     x3(k) = (k-1)*d3 + o3
  end do

  do k=1,n3
     z = o3 + (k-0.5)*d3
     do j=1,n2
        y = o2 + (j-0.5)*d2
        do i=1,n1
           x = o1 + (i-0.5)*d1
           
           if(z<1000) then
              val = 0.3
           else if(z<1500) then
              val = 1.5
           else if(z<1600) then
              val = 100.
           else
              val =  2.5
           end if

           rho(i, j, k) = val
        enddo
     enddo
  enddo


  open(30,file='frho',access='direct',recl=4*n1*n2*n3,status='replace')
  write(30,rec=1) rho
  close(30)

  bathy(:,:) = 1000
  open(30,file='fbathy',access='direct',recl=4*(n1+1)*(n2+1),status='replace')
  write(30,rec=1) bathy
  close(30)
  
  !here we create a resistivity model for emg3d whose z is positive upward
  !and we include one layer of air on the top
  o3 = -40
  do k=n3+1,1,-1 !we start from bottom to the top 
     z = o3 + (k-0.5)*d3
     print *, z
     do j=1,n2
        y = o2 + (j-0.5)*d2
        do i=1,n1
           x = o1 + (i-0.5)*d1

           if(z<0) then
              val = 1e8
           else if(z<1000) then
              val = 0.3
           else if(z<1500) then
              val = 1.5
           else if(z<1600) then
              val = 100.
           else
              val =  2.5
           end if

           write(7,*) val!output txt file as fort.7, as input for emg3d_test.py
           if(z<0) print *, val
        enddo
     enddo
  enddo
  
  deallocate(x1)
  deallocate(x2)
  deallocate(x3)
  deallocate(rho)
  deallocate(bathy)
end program create_model
