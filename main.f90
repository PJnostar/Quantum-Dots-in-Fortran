! compile with:  gfortran mod_subr.f90 mod_zdiagonalize.f90 main2.f90 ~3jurkowski/Downloads/ARPACK/libarpack_SUN4.a ../../../../../../../../usr/lib/lapack/liblapack.so.3.7.0 && ./a.out


program testt
  USE mod_subr
  USE mod_zdiagonalize
  USE mod_codc
  USE mod_multi_el
  USE mod_elmnt
  Use mod_time
  implicit none
  integer, parameter :: nsz=600000
  
  real* 8 :: a1(3), a2(3), tau(3), rp(nsz,4)
  real* 8 :: a, d, dxl
  real* 8 :: rmax, rxmax, rymax
  real* 8 :: w_u, w_d, w, Rpp, d_dot_center, x_shift, y_shift, d_asym, w_asym
  real* 8 :: t, t1, t2, t3, Fz
  real* 8 :: Bz, Bz0
  real* 8 :: en_shift		!for changing the part of energy spectrum calculated by arpack
  real* 8 :: F_ac			!F_ac is the amplitute of the ac electric field
  real* 8 :: au2T = 2.35e5
  real* 8 :: au2ns = 2.418884e-8
  integer, ALLOCATABLE :: nz(:,:), nrst(:,:), nnrst(:,:,:)
  real* 8, ALLOCATABLE :: enes(:), enes_no(:), Vp(:)
  complex*16, ALLOCATABLE :: cust(:,:), cust_no(:,:), cnz(:), elmnt(:,:)
  integer :: nst = 50
  integer :: no = 30
  integer :: n_ions, n_1st, jnz
  integer :: sf, ii, i_ddc, n_ddc, i, j, k, l, iw


  !---------Flagi------------------------
  logical :: flag_generate_Ham = .FALSE.	!false = size of cnz is unknown, changed to true after 1st iteration
  

  real* 8 :: deb			!useful for debugging

	sf=1
	a=.38936/.05292*sf
	! a - stala sieci
	d=.2248/.05292*sf
	! d NN w plaszczyznie
	dxl=.046/.05292
	! dxl=2*l , czyli odleglosc w 'z'
	a1(1)=0.5*a
	a1(2)=a*sqrt(3.0)/2
	a1(3)=0
	a2(1)=a
	a2(2)=0
	a2(3)=0
	! wektory sieci
	tau(1)=0
	tau(2)=d
	tau(3)=dxl

	en_shift = 0e-3/27.21139
	w = 200e-3/27.21139
	! w_u = 200e-3/27.21139
	! w_d = 1.03*w_u
	Rpp = 12e-9/5.291772e-11
	d_dot_center = 24e-9/5.291772e-11
	x_shift = 0.e-9/5.291772e-11
	y_shift = 0.e-9/5.291772e-11

	t   = 1.6/27.21139/sf
	t1  = 7./15.*1.e-3/27.21139
	t2  = 0.75e-3/27.21139
	t3  = 1.e-5/27.21139							!ueV -> a.u.
	Fz  = 17.*(1.e-4*5.292772/27.21139)			!meV/A -> a.u.
	! Bz  = 1./au2T

	deb=XX1

  
	!--------------------------Generowanie siatki-----------------------------
	write(*,*) 'poczatek generowania siatki'
	!rp(n, 4): n-nr atomu, 1-x, 2-y, 3-spin, 4-z
	rmax = 25.0/.05291772
	call generate_grid_sdot( a1, a2, tau, rp, nsz, a, d, dxl, rmax, n_ions )
	! rxmax = 20.0/.05291772
	! rymax = 46.3/.05291772
	! call generate_grid_ddot( a1, a2, tau, rp, nsz, a, d, dxl, rxmax, rymax &
	  ! & , n_ions )
	ALLOCATE(nrst(n_ions,0:3))
	ALLOCATE(nnrst(n_ions,0:12,2))
	ALLOCATE(Vp(n_ions))
	write(*,*) 'koniec generowania siatki'
	if( size(rp,1) .lt. n_ions) then 
		write(*,*) '--------------------------------------'
		write(*,*) ' PRZYDZIELONO ZBYT MALO PAMIECI NA Rp'
		write(*,*) '--------------------------------------'
		go to 10000
	end if

	!--------------------------Generowanie siatki sasiadow--------------------
	write(*,*) 'poczatek generowania siatki sasiadow'
	call generate_nearest( rp, nsz, n_ions, nrst, nnrst, d, a )
	write(*,*) 'koniec generowania siatki sasiadow'

	!--------------------------Generowanie potencjalu-------------------------
	print*, ' '
	print*, '======================================'
	print*, 'Poczatek generowania potencjalu'
	! print*, 'w = ', w/1.e-3*27.2116, 'meV'
	! call generate_Vp( rp, nsz, d_dot_center, n_ions, w_u, w_d, Vp, Rpp)
	! call generate_Vp_sdot( rp, nsz, n_ions, w, Vp, Rpp )
	d_asym = 0.9*Rpp
	w_asym = -0.1*w
	call generate_Vp_sdot_asymmetric( rp, nsz, n_ions, w, Vp, Rpp, d_asym, w_asym )
	! call generate_Vp_ddot_min( rp, nsz, d_dot_center, n_ions, w_u, w_d, Vp, Rpp )
	print*, 'Koniec generowania potencjalu'
	! go to 10000

	!--------------------------Tworzenie Hamiltonianu-------------------------
	write(*,*) 'poczatek skladania Hamiltonianu'
	Bz = XX1/au2T
	flag_generate_Ham = .FALSE.
	call generate_Ham( rp, Vp, nsz, n_ions, nrst, nnrst, cnz, nz, jnz, a, d, t, t1, t2 &
	& , t3, w, dxl, Fz , en_shift, Bz, flag_generate_Ham )
	write(*,*) 'koniec skladania Hamiltonianu'
  
	!--------------------------Rozwiazanie problemu---------------------------
	write(*,*) 'poczatek diagonalizacji'
	allocate(cust(n_ions,nst))
	allocate(enes(nst))
	call zdiagonalize( cnz, nz, jnz, nrst, rp, n_ions, nst, nsz, cust, enes )  
	write(*,*) 'koniec diagonalizacji'

	!--------------------------Sortowanie energii/funkcji wlasnych------------
	!--------------------------znalezienie pierwszej nieujemnej energii-------
	!--------------------------i zapisanie energii/funkcji do pliku----------
	call sort_eigenstates( nst, n_ions, n_1st, enes, cust )
	write(*,*) 'n_1st = ', n_1st

	i=0
	allocate(cust_no(n_ions,no))
	allocate(enes_no(no))
	do j=n_1st, n_1st+no-1
		i=i+1
		enes_no(i) = enes(j)
		cust_no(:,i) = cust(:,j)
	enddo
	do ii=1,no
		write(*,*) (enes_no(ii)+en_shift)*27211.39, 'meV'
	enddo
	open(111, FILE="../E_B/enes_XX2.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
	write(111,'(100F20.15)') Bz*au2T, ((enes_no(i)+en_shift)*27211.39, i=1,no)
	close(111)
	
	! call fun_which_dott( nsz, no, n_ions, en_shift, rp, cust_no, enes_no, Rpp, d_dot_center )
	if (deb==1.) then
		call print_wavefunctions( nsz, no, n_ions, en_shift, rp, cust_no, enes_no )
	endif
	go to 10000
	
	open(111, FILE="READ/enes_no.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
	do i=1,no
		write(111,*) enes_no(i)+en_shift
	enddo
	close(111)
	open(222, FILE="READ/cust_no.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
	do i=1,n_ions
		do j=1,no
			write(222,*) cust_no(i,j)
		enddo
	enddo
	close(222)
	
	
	!--------------------------Liczenie elementow macierzowych----------------
	call generate_elmnt_sigmay( no, n_ions, nsz, elmnt, cust_no, rp )
	open(333, FILE="OUT/elmnt_sigmay.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
	do i=1,no
		write(333,'(100E15.6)') cdabs(elmnt(i,:))
	enddo
	close(333)
	open(444, FILE="READ/elmnt_sigmay.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
	do i=1,no
		do j=1,no
			write(444,*) elmnt(i,j)
		enddo
	enddo
	close(444)
	
	!----------------------------Koniec obliczen------------------------------
	DEALLOCATE(cnz)
	DEALLOCATE(nz)
	DEALLOCATE(cust)
	DEALLOCATE(enes)
	DEALLOCATE(cust_no)
	DEALLOCATE(enes_no)
	DEALLOCATE(elmnt)
  ! DEALLOCATE(which)
  ! DEALLOCATE(time_half)
  ! DEALLOCATE(time_98)
  ! DEALLOCATE(c_max)
  
10000 continue

END program



