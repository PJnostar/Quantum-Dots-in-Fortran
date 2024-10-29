MODULE mod_zdiagonalize
 IMPLICIT  NONE
 SAVE
CONTAINS

!-----------------------------------------------------------------------------
!--------Diagonalizacja complex-----------------------------------------------
!-----------------------------------------------------------------------------
SUBROUTINE zdiagonalize( cnz, nz, jnz, nrst, rp, n_ions, nst, nsz, cust, enes )
  REAL*8, INTENT(IN) :: rp(nsz, 4)
  REAL*8, INTENT(INOUT) :: enes(nst)
  COMPLEX*16, INTENT(IN) :: cnz(jnz)
  COMPLEX*16, INTENT(INOUT) :: cust(n_ions,nst)
  INTEGER, INTENT(IN) :: nrst(n_ions,0:3), nz(jnz,2)
  INTEGER, INTENT(IN) :: n_ions, nsz, nst, jnz

  INTEGER :: ido, n, nev, ncv, ldv, lworkl, info, iparam(1:11), ipntr(1:11)
  INTEGER :: maxn, maxnev, maxncv
  INTEGER :: ishfts, maxitr, mode, ierr
  REAL*8 :: tol, sigma
  REAL*8, ALLOCATABLE :: rwork(:)
  COMPLEX*16, ALLOCATABLE :: resid(:), v(:,:), d(:,:), workd(:), workl(:), workev(:)
  CHARACTER :: bmat*1, which*2
  LOGICAL :: rvec
  LOGICAL, ALLOCATABLE :: select(:)
  
  REAL*8 :: temp_A_up, temp_B_up
  INTEGER :: i, j, deb, petla
  
  petla=0
  write(*,*) 'szukam ', nst, ' wartosci wlasnych'
  
!     %--------------------------------------------------%
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%
!

  ! write(*,*) n,nsi,jnz,nst
  maxn=165000
  maxnev=420
  maxncv=400
  ldv=maxn
  n=n_ions
  nev=nst
  ncv=nev*2+1
  lworkl  = 3*ncv**2+5*ncv
  tol    = 0.00000001
  ido    = 0
  info   = 0
  ALLOCATE(resid(maxn))
  ALLOCATE(v(maxn,maxncv))
  ALLOCATE(d(maxncv,2))
  ALLOCATE(workd(3*maxn))
  ALLOCATE(workl(lworkl))
  ALLOCATE(rwork(ncv))
  ALLOCATE(workev(3*maxn))
  ALLOCATE(select(maxncv))
  resid=0
  v=0
  d=0
  workd=0
  workl=0
  rwork=0
  workev=0
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
	  
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 500000
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

!     %------------------------------------------------%
!     | M A I N   L O O P (Reverse communication loop) |
!     %------------------------------------------------%
!
 10   continue 
		  ! petla=petla+1
		  ! print*, petla
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
		! write(*,*) 'calling dsaupd'
		! read(*,*) deb
         call znaupd ( ido, bmat, n, which, nev, tol, resid, &
     &                    ncv, v, ldv, iparam, ipntr, workd, &
     &                    workl, lworkl, rwork, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix-vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix-vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
            call av(n, workd(ipntr(1)), workd(ipntr(2)), cnz, nz, jnz, nrst, nsz)
!
!           %-----------------------------------------% 
!           | L O O P   B A C K to call DSAUPD again. | 
!           %-----------------------------------------% 
!        
            go to 10
!        
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        |                                           |
!        | The routine DSEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!           
          rvec = .true.
          call zneupd ( rvec, 'A', select, d, v, ldv, sigma, &
     &         workev, bmat, n, which, nev, tol, resid, ncv, v, &
     &         ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr )
!
!         %----------------------------------------------%
!         | Eigenvalues are returned in the first column |
!         | of the two dimensional array D and the       |
!         | corresponding eigenvectors are returned in   |
!         | the first NCONV (=IPARAM(5)) columns of the  |
!         | two dimensional array V if requested.        |
!         | Otherwise, an orthogonal basis for the       |
!         | invariant subspace corresponding to the      |
!         | eigenvalues in D is returned in V.           |
!         %----------------------------------------------%
		  ! write(*,*) 'rozmiary v', size(v,1), size(v,2)
		  do i=1,nst
		    temp_A_up = 0
			temp_B_up = 0
			enes(i) = d(i,1)
			do j=1,n_ions
			  cust(j,i)=v(j,i)
			  temp_A_up = temp_A_up+zabs(cust(j,i))**2
			end do
			! do j=1,n_ions
			  ! cust(j,i)=cust(j,i)/sqrt(temp_A_up)
			! end do
		  end do
		  ! write(*,*) 'temp_A_up = ', temp_A_up
		  
		  !POSORTOWAC FUNKCJE I ENERGIE

!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
       if ( info .eq. 1) then
           print *, ' '
           print *, ' Maximum number of iterations reached.'
           print *, ' '
       else if ( info .eq. 3) then
           print *, ' '
           print *, ' No shifts could be applied during implicit &
   &                  Arnoldi update, try increasing NCV.'
           print *, ' '
       end if
!
!        print *, ' '

!        print *, ' '
!        print *, ' Size of the matrix is ', n
!        print *, ' The number of Ritz values requested is ', nev
!        print *, ' The number of Arnoldi vectors generated',
!    &            ' (NCV) is ', ncv
!        print *, ' What portion of the spectrum: ', which
!        print *, ' The number of converged Ritz values is ',
!    &              nconv
       print *, ' The number of Implicit Arnoldi update', &
   &            ' iterations taken is ', iparam(3)
!        print *, ' The number of OP*x is ', iparam(9)
!        print *, ' The convergence criterion is ', tol
!        print *, ' '
	  end if
9000 continue
  write(*,*) 'dseupd skonczony'
  
  
  OPEN(22, FILE="OUT/d.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1, nst
    write(22,*) enes(i)*1e3*27.2116
  end do
  CLOSE(22)
		 
  DEALLOCATE(resid)
  DEALLOCATE(v)
  DEALLOCATE(d)
  DEALLOCATE(workd)
  DEALLOCATE(workl)
  DEALLOCATE(rwork)
  DEALLOCATE(workev)
  DEALLOCATE(select)
  
END SUBROUTINE zdiagonalize


!------------------------------------------------------------------------------------
!----------------------------AV mnozenie macierzy dla arpacka------------------------
!------------------------------------------------------------------------------------
SUBROUTINE av(n, v, w, cnz, nz, jnz, nrst, nsz)
  COMPLEX*16, INTENT(INOUT) :: v(n), w(n)
  COMPLEX*16, INTENT(IN) :: cnz(jnz)
  INTEGER, INTENT(IN) :: nrst(n,0:3), nz(jnz,2)
  INTEGER, INTENT(IN) :: n, nsz, jnz

  INTEGER :: i, j, k, is
  
  ! Computes w <--- OP*v,
  
  w=0.
  do k=1,jnz
    i=nz(k,1)
	j=nz(k,2)
	if ( i .eq. j ) then
	  w(i)=w(i)+cnz(k)*v(j)
	else
	  w(i)=w(i)+cnz(k)*v(j)
	  w(j)=w(j)+conjg(cnz(k))*v(i)
	end if
  end do
  
END SUBROUTINE av




END MODULE mod_zdiagonalize