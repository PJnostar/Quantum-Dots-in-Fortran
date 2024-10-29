MODULE mod_codc
 IMPLICIT  NONE
 SAVE
CONTAINS

!-----------------------------------------------------------------------------
!------------------------------- Multi electrons interaction -----------------
!-----------------------------------------------------------------------------
SUBROUTINE generate_codc( nsz, n_ions, jnz, no, rp, enes, cust, codc)
  
  INTEGER, INTENT(IN) :: nsz, n_ions, jnz, no
  REAL*8, INTENT(IN) :: rp(nsz,4), enes(no)
  COMPLEX*16, INTENT(IN) :: cust(n_ions,no)
  COMPLEX*16, INTENT(INOUT) :: codc(no,no,no,no)
  
  INTEGER :: i, j, k, l, ir1, ir2
  REAL*8 :: Zeffe, Zupar, eps, r12, odl
  REAL*8 :: au2T = 2.35e5
  COMPLEX*16 :: ci, cro(n_ions), cpot(n_ions)

  
  codc=0
  eps=9.1
  Zeffe=4.15 
  Zupar= 3577./46080.*Zeffe
  ci = (0.0,1.0)
 
  write(*,*) no
  write(*,*) n_ions
  write(*,*) enes(:)

  do i=1,no
  ! i=XX1
	print*, i
	do k=1,no
		do ir1=1,n_ions
		cro(ir1) = dconjg(cust(ir1,i))*cust(ir1,k)
		end do
		cpot=0
		do ir1=1,n_ions
			do ir2=1,n_ions
				r12=sqrt((rp(ir1,1)-rp(ir2,1))**2 &
				&    +(rp(ir1,2)-rp(ir2,2))**2 &
				&    +(rp(ir1,4)-rp(ir2,4))**2)     
				odl=1/eps/(r12+1/(Zupar)*exp(-r12*10))
				cpot(ir2)=cpot(ir2)+odl*cro(ir1)
			end do
		end do

		do j=1,no
			do l=1,no
				do ir2=1,n_ions
					codc(i,j,k,l)=codc(i,j,k,l)+cpot(ir2)* &
					& dconjg(cust(ir2,j))* &
					& cust(ir2,l)
				end do
			end do
		end do
	end do
  end do


  OPEN(11, FILE="OUT/codc.dat",ACTION="WRITE", Status='unknown', FORM="FORMATTED")
  do i=1,no
	do j=1,no
		do k=1,no
			do l=1,no
				! write(11,'(4I5,2ES30.18)') XX1, j, k, l, codc(XX1,j,k,l)
				write(11,'(4I5,2E30.18)') i, j, k, l, codc(i,j,k,l)
			enddo
		enddo
	enddo
  enddo
  close(11)
  
END SUBROUTINE generate_codc


END MODULE mod_codc