Module system
    Use pub
    Use lapack95   
    Implicit None
    Private
    Integer,parameter :: DP=Selected_Real_kind(50,14),n=2
    Complex(kind=DP) :: H(n,n),v_x(n,n),v_y(n,n),v_z(n,n)
    Real(kind=DP) :: k_x,k_y,k_z,m,c,l,omega
    Real(kind=DP) :: d1=l*k_x*k_z,d2=-l*k_y*k_z,d3=m-c*(k_x**2+k_y**2+k_z**2),kT
        Call Pauli()
        Subroutine H(k_x,k_y,k_z)
        H=d1*sx+d2*sy+d3*sz
        End Subroutine
        Subroutine v_x(k_x,k_y,k_z)
        v_x=l*k_z*sx-2.*c*k_x**2*sz
        End Subroutine
        Subroutine v_y(k_x,k_y,k_z)
        v_y=-l*kz*sy-2.*c*k_y**2*sz
        End Subroutine
        Subroutine v_z(k_x,k_y,k_z)
        v_z=l*k_x*sx-l*k_y*sy-2.*c*k_z*sz
        End Subroutine
        Subroutine optical(k_x,k_y,k_z)
        !declaring variable
        Integer :: i,j
        Real(kind=DP) :: wr(n),wi(n),vl(n,n),vr(n,n),E0,eta,Fermi(n)
        Complex(kind=DP) :: VxV(n,n),VyV(n,n),VzV(n,n),&
            &f_xx,f_yx,f_zx,f_zz,R_xx,R_yx,R_zx,R_zz,dos
        !call functions
        Call H(k_x,k_y,k_z)
        Call v_x(k_x,k_y,k_z)
        Call v_y(k_x,k_y,k_z)
        Call v_z(k_x,k_y,k_z)
        Call geev(H,wr,wi,vl,vr)
        Fermi=fermi(wr,E0)
        VxV=matmul(conjg(transpose(vr)),matmul(v_x,vr)))
        VyV=matmul(conjg(transpose(vr)),matmul(v_y,vr)))
        VzV=matmul(conjg(transpose(vr)),matmul(v_z,vr)))
        f_xx=0.0_DP
        f_yx=0.0_DP
        f_zx=0.0_DP
        f_zz=0.0_DP
        do i=1,n
            do j=1,n
                if (m /=n) then
                    R_xx=(Fermi(i)-Fermi(j))*VxV(i,j)*VxV(j,i)/(omega+Fermi(i)-Fermi(j)+cmplx(0,1)*eta)
                    R_yx=(Fermi(i)-Fermi(j))*VyV(i,j)*VxV(j,i)/(omega+Fermi(i)-Fermi(j)+cmplx(0,1)*eta)
                    R_zx=(Fermi(i)-Fermi(j))*VzV(i,j)*VxV(j,i)/(omega+Fermi(i)-Fermi(j)+cmplx(0,1)*eta)
                    R_zz=(Fermi(i)-Fermi(j))*VzV(i,j)*VzV(j,i)/(omega+Fermi(i)-Fermi(j)+cmplx(0,1)*eta)
                    f_xx=f_xx+imag(R_xx)
                    f_yx=f_yx+imag(R_yx)
                    f_zx=f_zx+imag(R_zx)
                    f_zz=f_zz+imag(R_zz)
                    dos=dos+eta/(omega+Fermi(i)-Fermi(j)+cmplx(0,1)*eta)
                end do
            end do
        end do      
        End Subroutine        
        Subroutine Gauss_Integral(k_x,k_y,k_z)
        Real(kind=DP) :: 
        
        
        
        
        
        
        
        End Subroutine    
        Real Function fermi(E,E0)
        Real(kind=DP),allocate,intent(in) :: E(:),E0(:)
        Real(kind=DP),intent(out) :: Fermi(size(E))
        s=(E-E0)/kT
        if (s<60.0) then
            Fermi=1/(1.0+exp(s/kT))
        elseif
            Fermi=exp(-s/kT)/(exp(-s/kT)+1)
        endif
        End Function fermi
	    SUBROUTINE gauleg(x1,x2,x,w)
        USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror 
	    REAL(SP), INTENT(IN) :: x1,x2
	    REAL(SP), DIMENSION(:), INTENT(OUT) :: x,w
	    REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	    INTEGER(I4B) :: its,j,m,n
	    INTEGER(I4B), PARAMETER :: MAXIT=10
	    REAL(DP) :: xl,xm
	    REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	    LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	    n=assert_eq(size(x),size(w),'gauleg')
	    m=(n+1)/2
	    xm=0.5_dp*(x2+x1)
	    xl=0.5_dp*(x2-x1)
	    z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	    unfinished=.true.
	    do its=1,MAXIT
		    where (unfinished)
			    p1=1.0
			    p2=0.0
		    end where
		    do j=1,n
			    where (unfinished)
				    p3=p2
				    p2=p1
				    p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			    end where
		    end do
		    where (unfinished)
			    pp=n*(z*p1-p2)/(z*z-1.0_dp)
			    z1=z
			    z=z1-p1/pp
			    unfinished=(abs(z-z1) > EPS)
		    end where
		    if (.not. any(unfinished)) exit
	    end do
	    if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	    x(1:m)=xm-xl*z
	    x(n:n-m+1:-1)=xm+xl*z
	    w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	    w(n:n-m+1:-1)=w(1:m)
	    END SUBROUTINE gauleg
    
    End Module system

Module pub
    Implicit none
    Complex,parameter::im = (0.0,1.0) 
    Real(kind=DP),parameter::pi = 3.14159265359
    Real(kind=DP) :: m=0.2,c=0.1,l=0.1,kT=8.617333262e-4
    Complex(kind=DP) :: sx(2,2),sy(2,2),sz(2,2),s0(2,2)   
    End Module pub 
    
    
Program main
    !-------------Call libaray---------------------------------------------------------------------
    Include 'mpif.h'
    Use lapack95
    Use system
    !----------------------------------------------------------------------------------------------
    Implicit none  
    !---------------------------Declaring variable-------------------------------------------------
    Integer,parameter :: FL=Selected_Int_kind(50,14)
    Integer :: my_rank, ierr,np,time_start,time_end,N
    Real(kind=FL) :: 
    Complex(kind=FL):: D(2,2)
    !----------------------------------------------------------------------------------------------
    
    
    call system_clock(time_start)
        
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    
    
    write(*,*)a
    
    
    
    
    

    call system_clock(time_end)
    if(myid.eq.0.and.i_write.eq.1)write(*,*)'Calculation time (s):'
     &            ,dble(time_end-time_start)/1000.d0/10.d0

    call mpi_finalize(ierr) 
    End Program main
    !===============================================================================================
    Subroutine Pauli()
    use pub
    sx(1,2) = 1
    sx(2,1) = 1
    !----
    sy(1,2) = -im
    sy(2,1) = im
    !-----
    sz(1,1) = 1
    sz(2,2) = -1
    !-----
    s0(1,1) = 1
    s0(2,2) = 1
    End Subroutine Pauli
    