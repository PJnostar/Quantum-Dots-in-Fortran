MODULE mod_subr
 IMPLICIT  NONE
 SAVE
CONTAINS


!-----------------------------------------------------------------------------
!--------------------------------Generowanie Siatki Single Dot----------------
!-----------------------------------------------------------------------------
SUBROUTINE generate_grid_sdot( a1, a2, tau, rp, nsz, a, d, dxl, rmax, n_ions )
  integer :: nsz
  REAL* 8, INTENT(IN) :: a, d, dxl
  REAL* 8, INTENT(IN) :: a1(3), a2(3), tau(3), rmax
  REAL* 8, INTENT(INOUT) :: rp(nsz, 4)
  integer, INTENT(INOUT) :: n_ions
  
  REAL* 8 :: r(4), y1, y2, y3, y4
  integer :: ia1, ia2, i, j
  
  n_ions=0
  ! OPEN(11, FILE="OUT/st.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do ia1=-100, 100
    do ia2=-100, 100
	!------------Sieć A; spin up----------------
      r(1)=ia1*a1(1)+ia2*a2(1)
      r(2)=ia1*a1(2)+ia2*a2(2)+d/2*2
      r(4)=-dxl/2
      r(3)=+1
	  y1 =  r(1)/sqrt(3.0)+rmax
      y2 = -r(1)/sqrt(3.0)-rmax
      y3 =  r(1)/sqrt(3.0)-rmax
      y4 = -r(1)/sqrt(3.0)+rmax
      if( abs(r(1)).lt.rmax*sqrt(3.0)/2.0 .AND. &
        & r(2).lt.y1 .AND. r(2).gt.y2 .AND. r(2).gt.y3 .AND. r(2).lt.y4 ) then
        n_ions=n_ions+1
        do i=1,4
          rp(n_ions,i)=r(i)
        end do
        ! WRITE(11,*) n_ions, (rp(n_ions,i)*.05292, i=1,4)
       	! WRITE(11,*) n_ions, rp(n_ions,1)*.05292, rp(n_ions,2)*.05292, rp(n_ions,4)*.05292, rp(n_ions,3)
      end if
      !------------Sieć B; spin up----------------
      r(1)=ia1*a1(1)+ia2*a2(1)+tau(1)
      r(2)=ia1*a1(2)+ia2*a2(2)+tau(2)+d/2*2
      r(4)=-dxl/2+tau(3)
      r(3)=+1 
      y1 =  r(1)/sqrt(3.0)+rmax
      y2 = -r(1)/sqrt(3.0)-rmax
      y3 =  r(1)/sqrt(3.0)-rmax
      y4 = -r(1)/sqrt(3.0)+rmax
      if( abs(r(1)).lt.rmax*sqrt(3.0)/2.0 .AND. &
        & r(2).lt.y1 .AND. r(2).gt.y2 .AND. r(2).gt.y3 .AND. r(2).lt.y4 ) then
        n_ions=n_ions+1
        do i=1,4
          rp(n_ions,i)=r(i)
        end do
        ! WRITE(11,*) n_ions, (rp(n_ions,i)*.05292, i=1,4)
        ! WRITE(11,*) n_ions, rp(n_ions,1)*.05292, rp(n_ions,2)*.05292, rp(n_ions,4)*.05292, rp(n_ions,3)
      end if
      !------------Sieć A, B; spin down-------------
      do i=1,n_ions
        do j=1,4
        rp(i+n_ions,j)=rp(i,j)
        end do
        rp(i+n_ions,3)=-1
      end do
    end do
  end do
  ! CLOSE(11)
  n_ions=2*n_ions			!poniewaz sa jeszcze spiny down
  write(*,*) 'znaleziono ', n_ions, 'atomow podsieci A i B, spin up i down'
END SUBROUTINE generate_grid_sdot



!-----------------------------------------------------------------------------
!--------------------------------Generowanie Siatki Double Dot----------------
!-----------------------------------------------------------------------------
SUBROUTINE generate_grid_ddot( a1, a2, tau, rp, nsz, a, d, dxl, rxmax, rymax &
      & , n_ions )
  integer :: nsz
  REAL* 8, INTENT(IN) :: a, d, dxl
  REAL* 8, INTENT(IN) :: a1(3), a2(3), tau(3), rxmax, rymax
  REAL* 8, INTENT(INOUT) :: rp(nsz,4)
  integer, INTENT(INOUT) :: n_ions
  
  REAL* 8 :: r(4), y1, y2, y3, y4
  integer :: ia1, ia2, i, j

  rp=0
  n_ions=0
  OPEN(11, FILE="OUT/st.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do ia1=-150, 150
    do ia2=-150, 150
	!------------Sieć A; spin up----------------
      r(1)=ia1*a1(1)+ia2*a2(1)
      r(2)=ia1*a1(2)+ia2*a2(2)+d/2*2
      r(4)=-dxl/2
      r(3)=+1
	  y1 =  r(1)/sqrt(3.0)+rymax*0.5+rxmax/sqrt(3.)
	  y2 = -r(1)/sqrt(3.0)-rymax*0.5-rxmax/sqrt(3.)
	  y3 =  r(1)/sqrt(3.0)-rymax*0.5-rxmax/sqrt(3.)
	  y4 = -r(1)/sqrt(3.0)+rymax*0.5+rxmax/sqrt(3.)
	  if( abs(r(1)).lt.rxmax .AND. &
	    & r(2).lt.y1 .AND. r(2).gt.y2 .AND. r(2).gt.y3 .AND. r(2).lt.y4 ) then
	    n_ions=n_ions+1
	    do i=1,4
          rp(n_ions,i)=r(i)
	    end do
	    ! WRITE(11,*) n_ions, (rp(n_ions,i)*.05292, i=1,4)
		WRITE(11,*) n_ions, rp(n_ions,1)*.05292, rp(n_ions,2)*.05292, rp(n_ions,4)*.05292, rp(n_ions,3)
	  end if
	  !------------Sieć B; spin up----------------
      r(1)=ia1*a1(1)+ia2*a2(1)+tau(1)
      r(2)=ia1*a1(2)+ia2*a2(2)+tau(2)+d/2*2
      r(4)=-dxl/2+tau(3)
      r(3)=+1 
	  y1 =  r(1)/sqrt(3.0)+rymax*0.5+rxmax/sqrt(3.)
	  y2 = -r(1)/sqrt(3.0)-rymax*0.5-rxmax/sqrt(3.)
	  y3 =  r(1)/sqrt(3.0)-rymax*0.5-rxmax/sqrt(3.)
	  y4 = -r(1)/sqrt(3.0)+rymax*0.5+rxmax/sqrt(3.)
	  if( abs(r(1)).lt.rxmax .AND. &
	    & r(2).lt.y1 .AND. r(2).gt.y2 .AND. r(2).gt.y3 .AND. r(2).lt.y4 ) then
	    n_ions=n_ions+1
	    do i=1,4
          rp(n_ions,i)=r(i)
	    end do
	    ! WRITE(11,*) n_ions, (rp(n_ions,i)*.05292, i=1,4)
		WRITE(11,*) n_ions, rp(n_ions,1)*.05292, rp(n_ions,2)*.05292, rp(n_ions,4)*.05292, rp(n_ions,3)
	  end if
	  !------------Sieć A, B; spin down-------------
	  do i=1,n_ions
        do j=1,4
        rp(i+n_ions,j)=rp(i,j)
        end do
        rp(i+n_ions,3)=-1
      end do
    end do
  end do
  CLOSE(11)
  n_ions=2*n_ions			!poniewaz sa jeszcze spiny down
  write(*,*) 'znaleziono ', n_ions, 'atomow podsieci A i B, spin up i down'
  
END SUBROUTINE generate_grid_ddot


!-----------------------------------------------------------------------------
!--------------------------------Generowanie potencjalu-----------------------
!--------------------------------------double dot-----------------------------
SUBROUTINE generate_Vp( rp, nsz, d_dot_center, n_ions, w_u, w_d, Vp, Rpp )
  integer :: nsz
  REAL* 8, INTENT(IN) :: w_u, w_d, Rpp, d_dot_center, rp(nsz, 4)
  REAL* 8, INTENT(INOUT) :: Vp(n_ions)
  integer, INTENT(IN) :: n_ions
  
  REAL* 8 :: r(4), y1, y2, y3, y4, v_u, v_d
  integer :: ia1, ia2, i, j
  
  OPEN(11, FILE="OUT/vp_A.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  OPEN(22, FILE="OUT/vp_B.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1, n_ions
    v_u = - w_u*EXP(-( rp(i,1)**2 + (rp(i,2)-d_dot_center/2)**2 )/Rpp**2) 
	v_d = - w_d*EXP(-( rp(i,1)**2 + (rp(i,2)+d_dot_center/2)**2 )/Rpp**2) 
    if (rp(i,4).gt.0.) then					!sublattice A
	  ! Vp(i) = w_u + MIN(v_u, v_d)
	  Vp(i) = w_u + v_u + v_d
	  WRITE(11,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	else									!sublattice B
	  ! Vp(i) = -w_u + MIN(v_u, v_d)
	  Vp(i) = -w_u + v_u + v_d
	  WRITE(22,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	end if
  end do
  CLOSE(11)
  CLOSE(22)
  
END SUBROUTINE generate_Vp



!-----------------------------------------------------------------------------
!--------------------------------Generowanie potencjalu-----------------------
!--------------------------------------double dot-----------------------------
SUBROUTINE generate_Vp_ddot_min( rp, nsz, d_dot_center, n_ions, w_u, w_d, Vp, Rpp )
  integer :: nsz
  REAL* 8, INTENT(IN) :: w_u, w_d, Rpp, d_dot_center, rp(nsz, 4)
  REAL* 8, INTENT(INOUT) :: Vp(n_ions)
  integer, INTENT(IN) :: n_ions
  
  REAL* 8 :: r(4), y1, y2, y3, y4, v_u, v_d
  integer :: ia1, ia2, i, j
  
  OPEN(11, FILE="OUT/vp_A.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  OPEN(22, FILE="OUT/vp_B.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1, n_ions
    v_u = - w_u*EXP(-( rp(i,1)**2 + (rp(i,2)-d_dot_center/2)**2 )/Rpp**2) 
	v_d = - w_d*EXP(-( rp(i,1)**2 + (rp(i,2)+d_dot_center/2)**2 )/Rpp**2) 
	if (rp(i,4).gt.0.) then					!sublattice A
	  Vp(i) = w_u + MIN(v_u, v_d)
	  ! Vp(i) = w_u + v_u + v_d
	  WRITE(11,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	else									!sublattice B
	  Vp(i) = -w_u + MIN(v_u, v_d)
	  ! Vp(i) = -w_u + v_u + v_d
	  WRITE(22,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	end if
  end do
  CLOSE(11)
  CLOSE(22)
  
END SUBROUTINE generate_Vp_ddot_min



!-----------------------------------------------------------------------------
!--------------------------------Generowanie potencjalu-----------------------
!--------------------------------------single dot-----------------------------
SUBROUTINE generate_Vp_sdot( rp, nsz, n_ions, w, Vp, Rpp )
  integer :: nsz
  REAL* 8, INTENT(IN) :: w, Rpp, rp(nsz, 4)
  REAL* 8, INTENT(INOUT) :: Vp(n_ions)
  integer, INTENT(IN) :: n_ions
  
  REAL* 8 :: r(4), y1, y2, y3, y4, v
  integer :: ia1, ia2, i, j
  
  OPEN(11, FILE="OUT/vp_A.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  OPEN(22, FILE="OUT/vp_B.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1, n_ions
    v = - w*EXP(-( rp(i,1)**2 + (rp(i,2))**2 )/Rpp**2) 
    if (rp(i,4).gt.0.) then					!sublattice A
	  ! Vp(i) = w_u + MIN(v_u, v_d)
	  Vp(i) = w + v
	  WRITE(11,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	else									!sublattice B
	  ! Vp(i) = -w_u + MIN(v_u, v_d)
	  Vp(i) = -w + v
	  WRITE(22,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	end if
  end do
  CLOSE(11)
  CLOSE(22)
  
END SUBROUTINE generate_Vp_sdot



!-----------------------------------------------------------------------------
!--------------------------------Generowanie potencjalu-----------------------
!--------------------------------single dot assymmetric-----------------------
SUBROUTINE generate_Vp_sdot_asymmetric( rp, nsz, n_ions, w, Vp, Rpp, d_asym, w_asym )
  integer :: nsz
  REAL* 8, INTENT(IN) :: w, Rpp, rp(nsz, 4), d_asym, w_asym
  REAL* 8, INTENT(INOUT) :: Vp(n_ions)
  integer, INTENT(IN) :: n_ions
  
  REAL* 8 :: v
  integer :: ia1, ia2, i, j
  
  OPEN(11, FILE="OUT/vp_A.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  OPEN(22, FILE="OUT/vp_B.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1, n_ions
    v = - w*EXP(-( rp(i,1)**2 + (rp(i,2))**2 )/Rpp**2) &
	  & - w_asym*EXP(-( (rp(i,1)-d_asym)**2 + (rp(i,2))**2 )/Rpp**2)
    if (rp(i,4).gt.0.) then					!sublattice A
	  ! Vp(i) = w_u + MIN(v_u, v_d)
	  Vp(i) = w + v
	  WRITE(11,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	else									!sublattice B
	  ! Vp(i) = -w_u + MIN(v_u, v_d)
	  Vp(i) = -w + v
	  WRITE(22,*) rp(i,1)*.05292, rp(i,2)*.05292, Vp(i)*27.2116*1.e3
	end if
  end do
  CLOSE(11)
  CLOSE(22)
  
END SUBROUTINE generate_Vp_sdot_asymmetric



!-----------------------------------------------------------------------------
!--------Budowanie siatki najblizszych sasiadow i kolejnych-------------------
!-----------------------------------------------------------------------------
SUBROUTINE generate_nearest( rp, nsz, n_ions, nrst, nnrst, d, a )
  integer :: nsz
  REAL* 8, INTENT(INOUT) :: rp(nsz,4)
  integer, INTENT(INOUT) :: nrst(n_ions,0:3), nnrst(n_ions,0:12,2)
  ! nnrst (numer,nsasiad,1) - sasiad drugi
  ! nnrst (numer,nsasiad,2) - wspolny najblizszy
  REAL* 8, INTENT(IN) :: d, a
  integer, INTENT(IN) :: n_ions
  
  integer :: i, j, k, is
  REAL*8 :: dd, ddki, ddkj
  
    do 111 i=1,n_ions
	!----------Nearest neighbour----------
	is=0
	do 222 j=1,n_ions
	dd = (rp(i,1)-rp(j,1))**2+(rp(i,2)-rp(j,2))**2-d**2
	if ( i.ne.j .AND. rp(i,3)*rp(j,3).gt.0 .AND. &
	  & (rp(i,1)-rp(j,1))**2 + (rp(i,2)-rp(j,2))**2 .lt. d**2*1.05 ) then
	  is=is+1
      nrst(i,is)=j
      nrst(i,0)=is
	end if
222 continue
    ! write(*,*) is
	!----------Next nearest neighbour----------
	is=0
	do 333 j=1,n_ions
	dd = (rp(i,1)-rp(j,1))**2+(rp(i,2)-rp(j,2))**2-a**2
	if(rp(i,3)*rp(j,3).gt.0 .AND. abs(dd).lt.a**2/1000) then
	do 444 k=1,n_ions
	ddki=(rp(i,1)-rp(k,1))**2+(rp(i,2)-rp(k,2))**2-d**2
    ddkj=(rp(j,1)-rp(k,1))**2+(rp(j,2)-rp(k,2))**2-d**2
    if ( rp(i,3)*rp(k,3).gt.0 .AND. abs(ddki).lt.d**2/1000 &
       & .AND. abs(ddkj).lt.d**2/1000 ) then
    is=is+1
    nnrst(i,is,1)=j
    nnrst(i,0,1)=is
    nnrst(i,is,2)=k
    end if
444 continue
	end if
333 continue
111 continue
    
	! OPEN(11, FILE="OUT/nrst.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
	! OPEN(22, FILE="OUT/nnrst_1.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
	! OPEN(33, FILE="OUT/nnrst_2.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
	! do 666 i=1,n_ions
	! write(11,'(I4,4X,I4,4X,I4,4X,I4)') (INT(nrst(i,j)), j=0,3)
	! write(22,*) (nnrst(i,j,1), j=0,12) 
	! write(33,*) (nnrst(i,j,2), j=0,12) 
! 666 continue	
	! CLOSE(11)
	! CLOSE(22)
	! CLOSE(33)
  
END SUBROUTINE generate_nearest


!-----------------------------------------------------------------------------
!--------Tworzenie Hamiltonianu-----------------------------------------------
!-----------------------------------------------------------------------------
SUBROUTINE generate_Ham( rp, Vp, nsz, n_ions, nrst, nnrst, cnz, nz, jnz, a, d, t, t1, t2 &
   & , t3, w_u, dxl, Fz, en_shift, Bz, flag_generate_Ham )
  integer :: nsz
  REAL* 8, INTENT(IN) :: rp(nsz,4), Vp(n_ions), a, d, t, t1, t2, t3, w_u, dxl,Fz, en_shift, Bz
  COMPLEX*16, ALLOCATABLE, INTENT(INOUT) :: cnz(:)
  integer, INTENT(INOUT) :: jnz
  integer, ALLOCATABLE, INTENT(INOUT) :: nz(:,:)
  integer, INTENT(IN) :: nrst(n_ions,0:3), nnrst(n_ions,0:12,2), n_ions
  logical, intent(IN) :: flag_generate_Ham
  
  integer :: size_cnz, ip, il, ik, is, i, j
  REAL*8 :: dd
  REAL*8 :: x_lk, y_lk, x_kp, y_kp, znak, d_lp(1:3)
  REAL*8 :: g=2								!Lande Factor
  REAL*8 :: mi_b=0.5						!Bohr magneton in au
  COMPLEX*16 :: ci=(0.,1.)
  
  ! REAL*8 :: dkj(1:2), dik(1:2)
  ! COMPLEX*16 :: cm
	  
  ! OPEN(111, FILE="OUT/debug1.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
  if (flag_generate_Ham .eqv. .FALSE.) then
  size_cnz = 0										!-------Number of matrix elements
  do 11 i=1,n_ions
  do 12 il=1,INT(nrst(i,0))
  if (nrst(i,il).gt.i) then						!-------Hopping elements
    size_cnz=size_cnz+1
  end if
12 continue
  do 13 il=1,INT(nnrst(i,0,1))
  if (nnrst(i,il,1).gt.i) then					!-------Intrinsic SO elements
    size_cnz=size_cnz+1
  end if
13 continue
  do 14 j=i+1,n_ions							!-------Built-in Rashba elements
	dd=(rp(i,1)-rp(j,1))**2+(rp(i,2)-rp(j,2))**2
	if( dd.lt.a**2*1.05 .AND. dd.gt.a**2*.95 &
	  & .AND. rp(i,3)*rp(j,3).lt.0 ) then
	  size_cnz=size_cnz+1
	endif 
	if( dd.lt.d**2*1.05 .AND. dd.gt.d**2*.95 &	!-------Extrinsic Rashba elements
	  & .AND. rp(i,3)*rp(j,3).lt.0 ) then
	  size_cnz=size_cnz+1
	endif 
14 continue
11 continue
  size_cnz = size_cnz+n_ions					!-------Potential & Zeeman elements
  ALLOCATE(cnz(size_cnz))
  ALLOCATE(nz(size_cnz,2))
  cnz=0
  nz=0
  else
  cnz=0
  nz=0
  end if
  ! close(111)
  ! OPEN(11, FILE="OUT/cnz.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
	! 10 Format(' ',F7.4,SP,F7.4,2X,$)
  ! OPEN(11, FILE="OUT/debug.dat",ACTION="WRITE", &
	! & STATUS="REPLACE", FORM="FORMATTED")
  jnz=0
    do 111 ip=1,n_ions
	jnz=jnz+1
	cnz(jnz) = Vp(ip)-en_shift &			!--------potential------
	  & + 0.5*g*mi_b*Bz*rp(ip,3)			!--------Zeeman---------
	nz(jnz,1) = ip
	nz(jnz,2) = ip
	! write(11,*) nz(jnz,1), nz(jnz,2), cnz(jnz)
111 continue
    ! write(11,*) 'break'
	do 222 ip=1,n_ions						!--------hopping--------
	do 333 il=1,INT(nrst(ip,0))
	is = nrst(ip,il)
	if (is.gt.ip) then
	jnz=jnz+1
	cnz(jnz) = -t*exp(0.5*ci*Bz*(rp(is,1)*rp(ip,2) - rp(ip,1)*rp(is,2)))
	nz(jnz,1) = ip
	nz(jnz,2) = is
	! write(11,*) nz(jnz,1), nz(jnz,2), cnz(jnz)
	end if
333 continue
222 continue
    ! print*, 'Intrinsic SO'
    ! write(11,*) 'break'
    do 555 il=1,n_ions						!--------Intrinsic SO & built-in Rashba---
	do 556 is=1,INT(nnrst(il,0,1))			!---Intrinsic SO
	ip = nnrst(il,is,1)
    ik = nnrst(il,is,2)
	if (ip.gt.il) then
	jnz=jnz+1
    nz(jnz,1) = il
    nz(jnz,2) = ip
	x_lk = rp(ik,1) - rp(il,1)
	y_lk = rp(ik,2) - rp(il,2)
	x_kp = rp(ip,1) - rp(ik,1)
	y_kp = rp(ip,2) - rp(ik,2)
	znak = x_lk*y_kp - x_kp*y_lk
	if (znak.gt.0) then
      cnz(jnz) = ci*t2*rp(ip,3)*exp( 0.5*ci*Bz	 &
			  & *(rp(ip,1)*rp(il,2) - rp(il,1)*rp(ip,2)) )
	else
	  cnz(jnz) = -ci*t2*rp(ip,3)*exp(0.5*ci*Bz	 &
			  & *(rp(ip,1)*rp(il,2) - rp(il,1)*rp(ip,2)) )
	end if
	! write(11,*) nz(jnz,1), nz(jnz,2), AIMAG(cnz(jnz))
	end if
556 continue
    do 557 ip=il+1,n_ions					!---Built-in Rashba
	dd=(rp(il,1)-rp(ip,1))**2+(rp(il,2)-rp(ip,2))**2
	d_lp(1) = rp(ip,1) - rp(il,1)
	d_lp(2) = rp(ip,2) - rp(il,2)
	d_lp(3) = sqrt( d_lp(1)**2 + d_lp(2)**2 )
	d_lp(1) = d_lp(1)/d_lp(3)
	d_lp(2) = d_lp(2)/d_lp(3)
	if( dd.lt.a**2*1.05 .AND. dd.gt.a**2*.95 &
	  & .AND. rp(il,3).gt.0 .AND. rp(ip,3).lt.0 ) then
	  znak = -rp(il,4)/abs(rp(il,4))
	  jnz=jnz+1
      nz(jnz,1) = il
      nz(jnz,2) = ip
	  cnz(jnz) = -ci*t1*znak*(d_lp(2)+ci*d_lp(1))  &
	          & *exp( 0.5*ci*Bz*(rp(ip,1)*rp(il,2) - rp(il,1)*rp(ip,2)) )
	end if
	if( dd.lt.d**2*1.05 .AND. dd.gt.d**2*.95 &	!---Extrinsic Rashba 
	  & .AND. rp(il,3).gt.0 .AND. rp(ip,3).lt.0 ) then
	  znak = -rp(il,4)/abs(rp(il,4))
	  jnz=jnz+1
      nz(jnz,1) = il
      nz(jnz,2) = ip
	  cnz(jnz) = ci*t3*(2*w_u/dxl)/Fz*(d_lp(2)+ci*d_lp(1))  &
	          & *exp( 0.5*ci*Bz*(rp(ip,1)*rp(il,2) - rp(il,1)*rp(ip,2)) )
	  ! write(11,*) il, ip, cnz(jnz)
	end if
557 continue
555 continue

  print*, 'jnz', jnz, size(cnz,1)
  if (jnz.ne.size(cnz,1)) then
    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
	print*, 'Nie zgadzaja sie wymiary macierzy Hamiltonianu'
	print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if
  ! CLOSE(11)
END SUBROUTINE generate_Ham



!-----------------------------------------------------------------------------
!----------------------------Sortowanie energii/funkcji wlasnych--------------
!-------------------------i znalezienie pierwszej nieujemnej energii----------
!-----------------------------------------------------------------------------
SUBROUTINE sort_eigenstates( nst, n_ions, n_1st, enes, cust )
  INTEGER, INTENT(IN) :: nst, n_ions
  INTEGER, INTENT(INOUT) :: n_1st			!number of the first state with E>0
  REAL* 8, INTENT(INOUT) :: enes(nst)
  COMPLEX*16, INTENT(INOUT) :: cust(n_ions,nst)
  
  REAL* 8 :: emin, e0
  INTEGER :: ii, jj, imin, ix
  COMPLEX*16 :: c0
  
    do 111 ii=1,nst
    emin=enes(ii)
    imin=ii
    do 112 jj=ii+1,nst
    if(enes(jj).lt.emin) then
    emin=enes(jj)
    imin=jj
    endif
112 continue
    e0=enes(ii)
    enes(ii)=emin
    enes(imin)=e0
    do 113 ix=1,n_ions
    c0=cust(ix,ii)
    cust(ix,ii)=cust(ix,imin)  
    cust(ix,imin)=c0
113 continue
111 continue

   n_1st=1
   do 222 ii=2,nst
   if (ii.gt.1 .and. enes(ii).gt.0 .and. enes(ii-1).lt.0) then
     n_1st=ii
   end if
222 continue

END SUBROUTINE sort_eigenstates


!-----------------------------------------------------------------------------
!----------------------------Obliczenie w ktorej kropce jest elektron---------
!-----------------------------------------------------------------------------
SUBROUTINE fun_which_dott( nsz, nst, n_ions, en_shift, rp, cust, enes, Rpp, d_dot_center )
	integer :: nsz
	INTEGER, INTENT(IN) :: nst, n_ions
	REAL*8, INTENT(IN) :: rp(nsz,4), enes(nst), en_shift, Rpp, d_dot_center
	COMPLEX*16, INTENT(IN) :: cust(n_ions,nst)

	INTEGER :: i, j, k
	REAL*8 :: sum_top, sum_bot, sum_spin
	
	do i=1,nst
		sum_top=0.
		sum_bot=0.
		do j=1,n_ions
			! TOP dot
			if ( (rp(j,1))**2+(rp(j,2)-d_dot_center/2)**2 < (1.1*Rpp)**2 ) then
				sum_top = sum_top + cdabs(cust(j,i))**2
			! BOTTOM dot
			elseif ( (rp(j,1))**2+(rp(j,2)+d_dot_center/2)**2 < (1.1*Rpp)**2 ) then
				sum_bot = sum_bot + cdabs(cust(j,i))**2
			endif
			! do k=1,n_ions
				! if ( rp(j,3)*rp(k,3)>0 ) then
					! sum_spin = sum_spin + cdabs(cust(k,i))*cust(j,i)*rp(j,3)
				! endif
			! enddo
			sum_spin = sum_spin + cdabs(cust(j,i))**2*rp(j,3)
		enddo
		
		if (sum_top>sum_bot) then
			write(*,*) 'TOP', sum_top, sum_bot, 'spin:', sum_spin, &
				&'energia:', enes(i)*27211.39
		else
			write(*,*) 'BOT', sum_top, sum_bot, 'spin:', sum_spin, &
				&'energia:', enes(i)*27211.39
		endif
	enddo

END SUBROUTINE fun_which_dott



!-----------------------------------------------------------------------------
!----------------------------Wypisanie funkcji falowych do pliku--------------
!-----------------------------------------------------------------------------
SUBROUTINE print_wavefunctions( nsz, nst, n_ions, en_shift, rp, cust, enes )
  integer :: nsz
  INTEGER, INTENT(IN) :: nst, n_ions
  REAL*8, INTENT(IN) :: rp(nsz,4), enes(nst), en_shift
  COMPLEX*16, INTENT(IN) :: cust(n_ions,nst)
  
  INTEGER :: i, j
  REAL*8 :: sum_au, sum_ad, sum_bu, sum_bd
  character(len=8) :: fmt
  character(len=2) :: x1
  
  fmt = '(I2.2)'
  
   do 111 i=1,nst
   write(x1,fmt) i
   OPEN(11, FILE="OUT/cust/"//trim(x1)//"_Au.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
   OPEN(22, FILE="OUT/cust/"//trim(x1)//"_Bu.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
   OPEN(33, FILE="OUT/cust/"//trim(x1)//"_Ad.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
   OPEN(44, FILE="OUT/cust/"//trim(x1)//"_Bd.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
   do 222 j=1,n_ions
   if (rp(j,4).lt.0. .AND. rp(j,3).gt.0) then		!sublattice A, spin up
      write(11,*) rp(j,1)*.05292, rp(j,2)*.05292, zabs(cust(j,i))**2, rp(j,4)*.05292
	else if(rp(j,4).gt.0. .AND. rp(j,3).gt.0) then	!sublattice B, spin up
	  write(22,*) rp(j,1)*.05292, rp(j,2)*.05292, zabs(cust(j,i))**2, rp(j,4)*.05292
	else if (rp(j,4).lt.0. .AND. rp(j,3).lt.0) then	!sublattice A, spin down
      write(33,*) rp(j,1)*.05292, rp(j,2)*.05292, zabs(cust(j,i))**2, rp(j,4)*.05292
	else if(rp(j,4).gt.0. .AND. rp(j,3).lt.0) then	!sublattice B, spin down
	  write(44,*) rp(j,1)*.05292, rp(j,2)*.05292, zabs(cust(j,i))**2, rp(j,4)*.05292
	end if
222 continue
    CLOSE(11)
    CLOSE(22)
    CLOSE(33)
    CLOSE(44)
111 continue

  ! OPEN(11, FILE="OUT/wave_fun.dat", ACTION="WRITE", Status='unknown', FORM="FORMATTED")
  ! do 333 i=1,nst		!enes(i), cust(:,i)
  ! sum_au = 0
  ! sum_ad = 0
  ! sum_bu = 0
  ! sum_bd = 0
  ! do 444 j=i,n_ions
  ! if 	  (rp(j,4).lt.0 .AND. rp(j,3).gt.0 ) then	!A up
  ! sum_au = sum_au + cdabs(cust(j,i))
  ! else if (rp(j,4).lt.0 .AND. rp(j,3).lt.0 ) then	!A down
  ! sum_ad = sum_ad + cdabs(cust(j,i))
  ! else if (rp(j,4).gt.0 .AND. rp(j,3).gt.0 ) then	!B up
  ! sum_bu = sum_bu + cdabs(cust(j,i))
  ! else if (rp(j,4).gt.0 .AND. rp(j,3).lt.0 ) then	!B down
  ! sum_bd = sum_bd + cdabs(cust(j,i))
  ! end if
! 444 continue
  ! write(11,'(5F10.2)') (enes(i)+en_shift)*1e3*27.2116, sum_au, sum_ad, sum_bu, sum_bd
! 333 continue
  
  CLOSE(11)  

END SUBROUTINE print_wavefunctions





!-----------------------------------------------------------------------------
!--------------------------------READ DATA FROM FILE--------------------------
!-----------------------------------------------------------------------------
SUBROUTINE READ_ENES_Y( nst, enes, y_ij )
  INTEGER, INTENT(IN) :: nst
  REAL* 8, INTENT(INOUT) :: enes(nst)
  COMPLEX*16, INTENT(INOUT) :: y_ij(nst,nst)
  
  INTEGER :: i
  enes = 0
  y_ij = 0
  
  OPEN(111, FILE="OUT/enes.dat", ACTION="READ", Status='unknown', FORM="FORMATTED")
  READ(111,*) enes
  CLOSE(111)
  
  OPEN(111, FILE="OUT/y_ij.dat", ACTION="READ", Status='unknown', FORM="FORMATTED")
  do i=1,nst
  READ(111,*) y_ij(i,:)
  end do
  close(111)
  
END SUBROUTINE READ_ENES_Y



END MODULE mod_subr