MODULE mod_multi_el
 IMPLICIT  NONE
 SAVE
CONTAINS

!-----------------------------------------------------------------------------
!------------------------------- Multi electrons interaction -----------------
!-----------------------------------------------------------------------------
SUBROUTINE multi_el( nsz, n_ions, jnz, nst, no, n_el, iperm, komb, cust, enes, &
	& rp, cnz, cnz_B0, codc, cih, enci, nz, d, Bz, flag_read_codc )
	
  USE mod_subr
  ! USE mod_zdiagonalize
  INTEGER, INTENT(IN) :: nsz, n_ions, jnz, nst, no, n_el, nz(jnz,2)
  INTEGER, INTENT(INOUT) :: iperm(720,0:6)
  integer, ALLOCATABLE, INTENT(INOUT) :: komb(:,:)
  REAL*8, INTENT(IN) :: rp(nsz,4), enes(no), d, Bz
  REAL*8, ALLOCATABLE, INTENT(INOUT) :: enci(:)
  COMPLEX*16, INTENT(IN) :: cust(n_ions,no), cnz_B0(jnz)
  COMPLEX*16, INTENT(INOUT) :: cnz(jnz), codc(no,no,no,no)
  COMPLEX*16, ALLOCATABLE, INTENT(INOUT) :: cih(:,:)
  LOGICAL, INTENT(IN) :: flag_read_codc

  INTEGER :: i, j, k, l, i1, i2, ie, il, ik, ip, ix, ir1, ir2, ii, jj, kk, ll, llt &
		&, numerod, iaa, ibb, icc, idd, info, ndlt, ist, ist1, ist2, ik1, ik2, ik3, ik4
  INTEGER, ALLOCATABLE :: nz_cih(:,:)
  REAL*8 :: Zeffe, Zupar, eps, r12, odl, lperm, ene, d2, spi2(720,3), wnu
  REAL*8 :: au2T = 2.35e5
  COMPLEX*16 :: csuma, c1e, ccod, cro(n_ions), cpot(n_ions) &
		&, cdif(no,no), cppx, cppy, cppz, elspin(nst,nst,3), elmnt(nst,nst,3) &
		&, csxau, csyau, cszau, csxbu, csybu, cszbu, csxad, csyad, cszad, csxbd    &
		&, csybd, cszbd, cswx, cswy, cswz, csx, csy, csz, ci, xyz(6,6,3)
  COMPLEX*16, ALLOCATABLE :: Cih0(:,:), CWORK(:), RWORK(:)
  CHARACTER(LEN=100) :: fmt_len
  CHARACTER(LEN=100) :: fmt1

  INTEGER :: it, j1, j2
  REAL*8 :: dt, t, tpdt, tmax, wmax, w, F_ac, pi, xmax(no,2), el, ek, au2ns = 2.418884e-8
  COMPLEX*16 :: cm1(no,no), cm2(no,no), co(no), cp(no), ipiv(no)
  !n_el - ilosc elektronow w kropkach
  !no - ilosc badanych stanow (tych z obliczen jednoelektronowych)
  
  eps=9.1
  Zeffe=4.15 
  Zupar= 3577./46080.*Zeffe
  print*, '====================================================',Bz*au2T
  ci = (0.0,1.0)
 
  if ( flag_read_codc .eqv. .FALSE. ) then   
  ! OPEN(11, FILE="OUT/codc.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  ! OPEN(22, FILE="OUT/codc2.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
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
  else
  ! OPEN(11, FILE="calc_codc/OUT/codc.dat",ACTION="READ", STATUS="UNKNOWN", FORM="FORMATTED")
  ! do 444 i=1,nst
  ! do 444 k=1,nst
  ! do 444 j=1,nst
  ! do 444 l=1,nst
  ! read(11,*) ii,kk,jj,ll, codc(i,j,k,l)
! 444 continue
  ! close(11)
  ! close(22)

  end if
  
! tutaj konczy sie rachunek  elementow macierzowych
! oddzialywania. od tego miejsca budujemy macierz do oddzialywania elektron-elektron
  numerod=0
  call perm(n_el,iperm)
  call ddgamma( n_el, lperm )	!lperm - liczba permutacji aka n_el!
  
  call wk(n_el,no,komb,numerod)
  ALLOCATE(cih(numerod,numerod))
  ALLOCATE(Cih0(numerod,numerod))
  cih=0
  Cih0=0
  
  ! go to 9999
  do il=1,jnz
	cnz(il)=cnz(il)-cnz_B0(il)
  enddo
  cdif=0
  do i1=1,no
	do i2=1,no
		csuma=0
		do ik=1,jnz
			il=nz(ik,1)
			ip=nz(ik,2)
			! csuma = csuma + cnz(ik)*dconjg(cust(il,i1))*cust(ip,i2)
			csuma = csuma + cnz(ik)*dconjg(cust(il,i1))*cust(ip,i2)
			if(il.ne.ip) then
				! csuma = csuma + dconjg(cnz(ik))*dconjg(cust(ip,i1))*cust(il,i2)
				csuma = csuma + dconjg(cnz(ik))*dconjg(cust(ip,i1))*cust(il,i2)
			endif
		enddo
		cdif(i1,i2)=csuma
	enddo
  enddo
  
  do i=1,numerod
	ene=0
	do ie=1,n_el
		ene=ene+enes(komb(i,ie))
	enddo
	cih(i,i)=ene
	print*, ene

	do j=1,numerod
		c1e=0
		do i1=1,n_el
			do ip=1,lperm
				llt=0 
				do ie=1,n_el
					if ( ie.ne.i1.and.komb(i,ie).eq.komb(j,iperm(ip,ie)) ) then  
						llt=llt+1
					endif
				enddo
					if ( llt.eq.n_el-1 ) then
						iaa=komb(i,i1)
						icc=komb(j,iperm(ip,i1))
						c1e=c1e+(-1)**iperm(ip,0)*cdif(iaa,icc) 
					endif 
			enddo
		enddo
		cih(i,j)=cih(i,j)+c1e
	enddo      


	do j=i,numerod
		! petla po parach elektronow
		ccod=0
		do i1=1,n_el
			do i2=i1+1,n_el
			! petla po permutacjach elementu l
				do ip=1,lperm
					! write(*,*) i1,i2,ip
					! sprawdzenie delt dla elektronow
					! o numerach innych niz i oraz j
					llt=0
					do ie=1,n_el
						if (ie.ne.i1.and.ie.ne.i2.and.komb(i,ie).eq.komb(j,iperm(ip,ie))) then
							llt=llt+1
						endif
					enddo
					! write(*,*) n_el,llt ,'lut'
					if(llt.eq.n_el-2) then
						! write(*,*) i,i1,komb(i,i1)
						! write(*,*) i,i2,komb(i,i2)
						iaa=komb(i,i1)
						ibb=komb(i,i2)
						icc=komb(j,iperm(ip,i1))
						idd=komb(j,iperm(ip,i2))  
						ccod=ccod+codc(iaa,ibb,icc,idd)*(-1)**iperm(ip,0)
					endif
				enddo
			enddo
		enddo
		cih(i,j)=cih(i,j)+ccod
		! write(*,*) i,j,dreal(ccod)*27211.6
	enddo
  enddo

  do k=1,numerod
	do l=k,numerod
		cih(l,k)=dconjg(cih(k,l))
	enddo
  enddo
  Cih0=Cih

  print*, 'numerodnumerodnumerodnumerodnumerodnumerodnumerodnumerod',numerod
  ALLOCATE(CWORK(2*numerod-1))
  ALLOCATE(RWORK(3*numerod-2))
  ALLOCATE(enci(numerod))
  ! ALLOCATE(CWORK(6000))
  ! ALLOCATE(RWORK(9000-2))
  call ZHEEV('v','u',Numerod,CiH,Numerod,ENCI,CWORK,2*Numerod-1, RWORK,INFO )
  ! call ZHEEV('v','u',Numerod,CiH,3000,ENCI,CWORK,6000, RWORK,INFO )
  
  OPEN(11, FILE="OUT/cih.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1,numerod
	do j=1,numerod
		! write(*,'(2I5.1,2E30.18)') i, j, cih(i,j)
		write(11,*) cih(i,j)
	enddo
  enddo
  close(11)
  

  write(*,*) 'Wypisuje energie'
  OPEN(11, FILE="OUT/enci_meV.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  write(fmt_len, '(I3.3)') numerod
  fmt1 = "(F7.3," // trim(fmt_len) // "F15.6/)"
  write(11,trim(fmt1), advance='no') Bz*au2T,enci*27211.39
  close(11)
  OPEN(11, FILE="OUT/enci.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
  do i=1,numerod
	write(11,*) enci(i)
  enddo
  close(11)
  
  ene=0
  do 991 i=1,numerod
  do 991 j=1,numerod
  ene=ene+dconjg(cih(i,1))*cih0(i,j)*cih(j,1)
991   continue
  write(*,*) 'ene = ', ene
    
  
END SUBROUTINE multi_el




!-----------------------------------------------------------------------------
!------------------------------- Permutation ---------------------------------
!-----------------------------------------------------------------------------
SUBROUTINE perm( n, iperm )
! generacja permutacji 
! iperm(i,j) - funkcja falowa obsadzona przez j-ty elektron w i-tej iteracji
! n - ilosc elektronow
! iperm(i,0) - znak i-tej permutacji
!               j= 123   123
! przyklad: n=3 -> fgh - gfh + ...
!				i=  1     2
! iperm(2,2) = f
  INTEGER, INTENT(INOUT) :: iperm(720,0:6)
  INTEGER, INTENT(IN) :: n
  INTEGER :: ip, is, np, ib, i, j, k, ki
  INTEGER :: ipe(0:13), ip2(0:13)
  REAL*8 :: icontr
  iperm = 0
  iperm(1,0)=0

!  diagonala
  do i=1,n
  iperm(1,i)=i
  enddo
  
  if(n.eq.1) then
  go to 1972		!exits the subroutine
  else 
  
  ip=1			!bierzaca permutowana ip
  is=1			!znak permutacji
  np=1   		!liczba znalezionych permutacji  
  
1     continue

! kopiowanie bierzacego
  do i=1,n
  ipe(i)=iperm(ip,i)
  enddo
  
! petla po probnych zmianach
  do 3 i=2,n
  do 333 ki=0,12
  ip2(ki)=ipe(ki)
333   continue
  ip2(0)=is
  ip2(1)=ipe(i)
  ip2(i)=ipe(1)     

! sprawdzenie czy taka permutacja juz byla
  do 4 j=1,np
  ib=0
  do 5 k=1,n
  ib=abs(ip2(k)-iperm(j,k))+ib
5     continue
  if(ib.eq.0) then
  go to 3		!znaleziona permutacja juz byla, idz szukac innej
  end if
4     continue

! dopisanie nowej
  np=np+1
  do 6 j=0,n
  iperm(np,j)=ip2(j)
6     continue

! write(*,122) (iperm(np,j),j=1,n)

  call ddgamma(n, icontr)
  if(np.eq.icontr) then
  go to 1972
  end if
3     continue
  ip=ip+1
  is=mod(iperm(ip,0)+1,2)
  go to 1
	  
  end if
  
1972  continue
      
END SUBROUTINE perm



!-----------------------------------------------------------------------------
!------------------------------- Factorial n!---------------------------------
!-----------------------------------------------------------------------------
SUBROUTINE ddgamma( i, dgamma )
  REAL*8, INTENT(INOUT) :: dgamma
  INTEGER, INTENT(IN) :: i
  INTEGER :: j
  dgamma=1.0
  do j=1,i
  dgamma=dgamma*j
  enddo
END SUBROUTINE ddgamma


!-----------------------------------------------------------------------------
!------------------------------- Kombinacje M po N  --------------------------
!------------------------------- M elementowy zbior --------------------------
!------------------------------- N elementow kombinacji ----------------------
SUBROUTINE wk( N, M, KNZM, nmb )
  INTEGER, INTENT(IN) :: N, M
  INTEGER, INTENT(INOUT) :: nmb
  integer, ALLOCATABLE, INTENT(INOUT) :: KNZM(:,:)
  
  INTEGER :: i, i1, i2, i3, i4, i5, in
! w tablicy KNZM zapisuje M po N kombinacji N elementow
! z M-elementowego zbioru aka (M N) = M!/(M-N)!N!
! nmb - liczba kombinacji
  nmb = 1
  
  do i=1,N
    NMB = NMB * (M+1-i)/i
  enddo

  ALLOCATE(KNZM(nmb,6))
  KNZM = 0
  
  if( N>M ) then			!test rozmiarow
  write(*,*) 'N>M w KNZM'
  stop
  endif

! kombinacja startowa        
  if( N==1 ) then
  do i=1,M
  KNZM(i,1)=i
  enddo
  endif
  if( N==2 ) then
  in=0
  do 2 i1=1,M
  do 2 i2=i1+1,M
  in=in+1
  KNZM(in,1)=i1
2     KNZM(in,2)=i2   
  endif
  if( N==3 ) then
  in=0
  do 3 i1=1,M
  do 3 i2=i1+1,M
  do 3 i3=i2+1,M
  in=in+1
  KNZM(in,1)=i1
  KNZM(in,2)=i2
  KNZM(in,3)=i3
3     continue
  endif
  if( N==4 ) then
  in=0
  do 4 i1=1,M
  do 4 i2=i1+1,M
  do 4 i3=i2+1,M
  do 4 i4=i3+1,M
  in=in+1
  KNZM(in,1)=i1
  KNZM(in,2)=i2
  KNZM(in,3)=i3
  KNZM(in,4)=i4
4     continue
  endif
  if( N==5 ) then
  in=0
  do 5 i1=1,M
  do 5 i2=i1+1,M
  do 5 i3=i2+1,M
  do 5 i4=i3+1,M
  do 5 i5=i4+1,M

  in=in+1
  KNZM(in,1)=i1
  KNZM(in,2)=i2
  KNZM(in,3)=i3
  KNZM(in,4)=i4
  KNZM(in,5)=i5
5     continue
  endif
  
END SUBROUTINE wk


END MODULE mod_multi_el