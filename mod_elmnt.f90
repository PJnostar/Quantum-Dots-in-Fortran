MODULE mod_elmnt
  USE mod_multi_el
 IMPLICIT  NONE
 SAVE
CONTAINS

! generate_elmnt
! generate_elmnt_sigmay

!-----------------------------------------------------------------------------
!------------------------------- Multi electrons <psi_i|y_ij|psi_j> ----------
!-----------------------------------------------------------------------------
SUBROUTINE generate_elmnt( elmnt, cust, rp, cih, n_el, no, numerod, n_ions, nsz, & 
	& iperm )
	integer, intent(in) :: n_el, no, n_ions, nsz
	integer, intent(inout) :: numerod
	integer, intent(inout) :: iperm(720,0:6)
	real*8, intent(in) :: rp(nsz,4)
	complex*16, intent(in) :: cust(n_ions,no), cih(numerod,numerod)
	complex*16, allocatable, intent(inout) :: elmnt(:,:,:)
	
	integer :: i, j, k, ist1, ist2, ie1, ie2, states_check, ip, ip_sign, states_max, ip_max
	integer :: ndlt, l, ie, ix
	integer, allocatable :: komb(:,:)
	real*8 :: lperm
	complex*16 :: y_ij(no,no), cppy
	! complex*16 :: csxau, csyau, cszau, csxbu, csybu, cszbu, csxad, csyad, cszad &
				! & , csxbd, csybd, cszbd, csx, csy, csz

	write(*,*) 'wszedl do gen_elmnt'
	
	call wk(n_el,no,komb,numerod)
	call perm(n_el,iperm)
	call ddgamma( n_el, lperm )	!lperm - liczba permutacji aka n_el!
	allocate(elmnt(numerod,numerod,3))
	elmnt = 0.

	write(*,*) lperm
	! do i=1,lperm
		! write(*,*) i, iperm(i,0), iperm(i,1), iperm(i,2)
	! enddo
	
	! do i=1,numerod
		! write(*,*) i, komb(i,1), komb(i,2)
		! do ist1=1,no
		! enddo
	! enddo
	
	y_ij = 0.
	do i=1,no
		do j=1,no
			do k=1,n_ions
				y_ij(i,j) = y_ij(i,j) + dconjg(cust(k,i))*cust(k,j)*rp(k,2)
			enddo
			write(11,'(2I6.3,F15.6)') i,j, cdabs(y_ij(i,j))
		enddo
	enddo
	
	! open(11, FILE="OUT/cppy_4.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
	! write(11,*) 'ist1 ist2 k l ip ie ix komb(k,ie) komb(l,iperm(ip,ie)) cppy cih(k,ist1) cih(l,ist2) y_ij(komb(k,ie) komb(l,iperm(ip,ie)))'
	do ist1=1,numerod
		do ist2=1,numerod
			cppy=0
			do k=1,numerod
				do l=1,numerod
					do ip=1,lperm
						do ie=1,n_el
							ndlt=0
							do ix=1,n_el
								if (ix.ne.ie.and.komb(k,ix).eq.komb(l,iperm(ip,ix))) then
									ndlt=ndlt+1
								endif
							enddo
							if (ndlt.eq.n_el-1) then
								cppy = cppy+(-1.)**iperm(ip,0) 	& 
								& *dconjg(cih(k,ist1))*cih(l,ist2)	& 
								& * y_ij(komb(k,ie),komb(l,iperm(ip,ie)))
								! write(11,'(9I5.3,20E20.10)') ist1, ist2, k, l, ip, ie, ix, komb(k,ie), komb(l,iperm(ip,ie)), cppy, cih(k,ist1), cih(l,ist2), y_ij(komb(k,ie), komb(l,iperm(ip,ie)))								
							endif
						enddo
					enddo
				enddo
			enddo
			elmnt(ist1,ist2,2)=cppy
		enddo
	enddo
	! close(11)

	
	open(11, FILE="OUT/elmnt.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
	open(22, FILE="OUT/elmnt_matrix.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
	do ist1=1,numerod
		do ist2=1,numerod
			write(11,*) elmnt(ist1,ist2,2)
			write(22,'(F10.4)', advance='no') cdabs(elmnt(ist1,ist2,2))
		enddo
		write(22,*) ' '
	enddo
	close(11)
	close(22)

	
END SUBROUTINE generate_elmnt




!-----------------------------------------------------------------------------
!------------------------------- Single electron <psi_i|s_y|psi_j> -----------
!-----------------------------------------------------------------------------
SUBROUTINE generate_elmnt_sigmay( no, n_ions, nsz, elmnt, cust, rp )
	integer, intent(in) :: no, n_ions, nsz
	real*8, intent(in) :: rp(nsz,4)
	complex*16, intent(in) :: cust(n_ions,no)
	complex*16, allocatable, intent(inout) :: elmnt(:,:)
	
	integer :: i, j, ist, ist1, ist2
	complex*16 :: ci = (0.,1.)
	
	allocate(elmnt(no,no))
	elmnt = 0.

	do ist1=1,no				!----------licze elementy macierzowe-----------
		do ist2=1,no
			do i=1,n_ions/2		!i - spin up, j - spin down
				j=i+n_ions/2
				elmnt(ist1,ist2) = elmnt(ist1,ist2) - ci*dconjg(cust(i,ist1))*cust(j,ist2)
			enddo
			do i=n_ions/2+1,n_ions	!i - spin down, j - spin up
				j=i-n_ions/2
				elmnt(ist1,ist2) = elmnt(ist1,ist2) + ci*dconjg(cust(i,ist1))*cust(j,ist2)
			enddo
		enddo
	enddo
	
END SUBROUTINE generate_elmnt_sigmay


	
END MODULE mod_elmnt