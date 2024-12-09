!Project: Advection Equation by Finite Difference and Finite Volume
!Name: Yaman Sanghavi and Pradhyumna Parthasarathy

module conditions
        implicit none
        save

!u is the velocity of the propogation
!dx is the discrete space gap
!C is CFL number which should be smaller than 1 to keep things stable
!Nx is the number of zones in x axis
!Nt is the number of edges in discretized t axis

        real :: u = 1.0
        real :: dx = 1.0
        real :: C = 0.7
        Integer :: Nx = 100
        Integer :: Nt = 100
        
end module

!------------------------------------------------------------

program one

      use conditions

      implicit none

      integer :: i,  j , s

      !-------------------------------------------------------

      ! Initial Conditions is a
      real :: x0
      real :: minmod

      real , dimension(:) , allocatable :: a
      real :: dt
      
      !at is the wavefunction in space and time
      !Dat is the space derivative of the wavefunction
      
      real , dimension(:,:),allocatable :: at , Dat

      real :: d1 , d2

      allocate(a(Nx))
      allocate(at(Nt,Nx))
      allocate(Dat(Nt,Nx))
      
      dt = C*dx/u
      
      write(*,*) "For initial conditions: enter 1 for Top Hat function or 2 for Gaussian function"
      read(*,*) s


      !Calling for initial conditions 
        if (s==1) then

        !Calling top hat conditions

                call top_hat(a)
      
        else if (s==2) then
      
        !calling gaussian condition

                call gaussian(a)

        else 
                write(*,*) "Please enter either 1 or 2"
        end if

      s = 0;

      !Initial conditions assigned to our main function array at
      do i = 1 , Nx
      at(1,i) = a(i)
      end do

      write(*,*) "Choose among the following methods:"
      write(*,*) "Enter 1 for Finite Difference method with FCTS method"
      write(*,*) "Enter 2 for Finite Difference method with Upwinding method"
      write(*,*) "Enter 3 for Finite Volume method without any limiter"
      write(*,*) "Enter 4 for Finite Volume method with minmod limiter"

      read(*,*) s

      !-------------------------------------------------------------

      if (s==1) then
      !Finite difference method with FCTS
      do i = 1 , Nt-1     

                do j = 2 , Nx-1
                at(i+1,j) = at(i,j) - 0.5*C*(at(i,j+1)-at(i,j-1))
                end do
      
                !Periodic Boundary Conditons at at x=1
                j = 1
                at(i+1,j) = at(i,j) - 0.5*C*(at(i,j+1)-at(i,Nx))

                !Periodic Boundary Condition at x = Nx
                j = Nx
                at(i+1,j) = at(i,j) - 0.5*C*(at(i,1)-at(i,j-1))
      end do

      call store(at)

      !-------------------------------------------------------------

      else if (s==2) then

      !Finite difference method with upwinding
      do i = 1 , Nt-1

                 do j = 2 , Nx
                 at(i+1,j) = at(i,j) - C*(at(i,j)-at(i,j-1))
                 end do

                 !Periodic Boundary Conditions at x = 1
                 j = 1
                 at(i+1,j) = at(i,j) - C*(at(i,j)-at(i,Nx))
      end do

      call store(at)
      !-------------------------------------------------------------

      else if (s==3) then

      !Finite volume method without any limiter
      do i = 1 , Nt-1

      !Making an array for the space derivative by the name Dat
                do j = 2 , Nx-1
                Dat(i,j) = 0.5*( at(i,j+1)-at(i,j-1) )/dx
                end do

                !Periodic Boundary Conditions at x = 1 and x = Nx
                Dat(i,Nx) = 0.5*( at(i,1)-at(i,Nx-1) )/dx
                Dat(i,1) = 0.5*( at(i,2)-at(i,Nx) )/dx

       !After finding the derivative, we solve the function at using the derivative Dat
                do j = 2 , Nx
                at(i+1,j) = at(i,j)  -   C*(  at(i,j)-at(i,j-1)  +   0.5*dx*(1 - C)*(Dat(i,j)-Dat(i,j-1))  )
                end do

                j = 1
                at(i+1,j) = at(i,j)  -   C*(  at(i,j)-at(i,Nx)  +   0.5*dx*(1 - C)*(Dat(i,j)-Dat(i,Nx))  )

      end do

      call store(at)

      !---------------------------------------------------------------

      else if (s==4) then

      !Finite volume with minmod limiter
       do i = 1 , Nt-1

      !Making an array for the space derivative by the name Dat
                do j = 2 , Nx-1
                !Dat(i,j) = 0.5*( at(i,j+1)-at(i,j-1) )/dx
                d1  = ( at(i,j) - at(i,j-1) )/dx
                d2  = ( at(i,j+1) - at(i,j) )/dx
                Dat(i,j) = minmod(d1,d2)
                end do

                !Periodic boundary conditions at x = Nx
                j = Nx
                d1  = ( at(i,j) - at(i,j-1) )/dx
                d2  = ( at(i,1) - at(i,j) )/dx
                Dat(i,j) = minmod(d1,d2)

                !Periodic boundary condition at x = 1
                j = 1
                d1  = ( at(i,j) - at(i,Nx) )/dx
                d2  = ( at(i,j+1) - at(i,j) )/dx
                Dat(i,j) = minmod(d1,d2)

       !After finding the derivative, we solve the function at using the derivative Dat
                do j = 2 , Nx
                at(i+1,j) = at(i,j)  -   C*(  at(i,j)-at(i,j-1)  +   0.5*dx*(1 - C)*(Dat(i,j)-Dat(i,j-1))  )
                end do

                !Periodic Boundary Conditions
                j = 1
                at(i+1,j) = at(i,j)  -   C*(  at(i,j)-at(i,Nx)  +   0.5*dx*(1 - C)*(Dat(i,j)-Dat(i,Nx))  )

      end do

      call store(at)
!------------------------------------------------------------------------

      else 
              write(*,*) "Please re-run the program and enter an integer from 1 to 4"

      end if
!------------------------------------------------------------------------

!Deallocating the arrays

      deallocate(a)
      deallocate(at)
      deallocate(Dat)

      !writing something for debugging
      write(*,*) "---------------"
      !write(*,*) minmod(3.4,4.5)
      
      stop 0
      end program

!------------------------------------------------------


function minmod(p,q) result(z)
        use conditions
        implicit none

        real, intent(out) :: p , q
        real :: z

        if (p*q>0) then
                if (p*p>q*q) then
            !            return q
            z = q    
                else
             !           return p
            z = p
                end if
        else 
              !  return 0.0 
              z = 0.0
        end if

end function minmod

!------------------------------------------------------

subroutine top_hat(b)

      use conditions
      implicit none

      real, intent(out) :: b(Nx)
      integer :: i

      do i = 1 , Nx/3 - 1
      b(i) = 0.0
      end do

      do i = Nx/3 , 2*Nx/3
      b(i) = 1.0
      end do

      do i = 2*Nx/3 + 1  , 100
      b(i) = 0.0
      end do
 
end subroutine top_hat

!---------------------------------------------------------

subroutine gaussian(g)
        use conditions
        implicit none

        real, intent(out) :: g(Nx)
        integer :: i

        do i = 1 , Nx
        g(i) = EXP(-(i-(0.5*Nx))*(i-(0.5*Nx))/(Nx))
        end do

        end subroutine gaussian

!----------------------------------------------------------
!Subroutine for storing the values in .dat files

subroutine store(f)
        use conditions
        implicit none

        real, dimension(Nt,Nx), intent(in) :: f

        integer :: i, L

        !Printing the mid time solution in a file

      
      open(newunit=L,file='out12.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(12*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out11.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(11*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out10.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(Nt/2,i)
      end do
      close(unit=L)
  
      open(newunit=L,file='out9.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(9*Nt/20,i)
      end do
      close(unit=L)
  
      open(newunit=L,file='out8.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(8*Nt/20,i)
      end do
      close(unit=L)
       
      open(newunit=L,file='out7.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(7*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out6.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(6*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out5.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(5*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out4.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(4*Nt/20,i)
      end do
      close(unit=L)
        
      open(newunit=L,file='out3.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(3*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out2.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i , f(2*Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out1.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i ,f(Nt/20,i)
      end do
      close(unit=L)

      open(newunit=L,file='out0.dat',status='REPLACE')
      do i = 1 , Nx
      write(L,*) i ,f(1,i)
      end do
      close(unit=L)

write(*,*) "The outputs for t=0.0 , t=0.1 , ... to t=1.2 are saved in the files by the name out0.dat, out1.dat, ... ,out12.dat"

end subroutine store

!-------------------------------------------------------------------------------------------



