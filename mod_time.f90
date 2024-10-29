MODULE mod_time
 IMPLICIT  NONE
 SAVE
CONTAINS

! function_time_AC_1el
! function_time_AC_multi_el


!-----------------------------------------------------------------------------
!----------------------------Rachunek z czasem--------------------------------
!-------------------------------Askar-Cakmak----------------------------------
!-----------------------------------1 el--------------------------------------
SUBROUTINE function_time_AC_1el( nsz, nst, enes, mat_ij, Ampl, au2ns, flag_time_mod,&
           & state, which, time_half, time_max, flag_time_half, flag_time_max, c_max )
	integer, INTENT(IN) :: nst, nsz, state, flag_time_mod
	integer, INTENT(INOUT) :: which(1:nst)
	REAL* 8, INTENT(IN) :: enes(nst), Ampl, au2ns
	REAL* 8, INTENT(INOUT) :: time_half(1:nst), time_max(1:nst), c_max(1:nst)
	COMPLEX*16, INTENT(IN) :: mat_ij(1:nst,1:nst)			!matrix elements between states i&j
	logical, INTENT(INOUT) :: flag_time_half, flag_time_max

	REAL*8 :: t, w, dw, w_min, w_max, tmax, dt, pi
	integer :: i, j, k, it, iw, deb
	complex*16 :: ci, c_p(1:nst), c_c(1:nst), c_n(1:nst)
	! logical :: flag_time_half = .FALSE.	!duration of Ampl pulse upon which the probability of finding
									!the electron in another state exceeds 50%
									!.FALSE. -> time was not found yet
									!.TRUE. -> time was found

	ci   = (0,1)
	pi   = atan(1.0) *4
	dt   = 1.e-5/au2ns							!ns -> a.u.

	! open(11, FILE="OUT/time.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
	if (flag_time_mod .eq. 0) then
		write(*,*) 'Sprawdzam konkretne czestotliwosci'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Access='append', Status='old', FORM="FORMATTED")
		tmax = 3.74/au2ns								!ns -> a.u.
		w_min = (enes(state)-enes(1)) 										!E=hbar*w -> w=E/hbar 
		w_max = (enes(state)-enes(1)) +0.2e-3/27.21139
		dw    = 3.e-3/27.21139			!doesnt matter, it just has to be bigger than 0.2e-3/27.21139
	else if (flag_time_mod .eq. 1) then
		write(*,*) 'Sprawdzam waskie spektrum czestotliwosci z malym krokiem'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Access='append',Status='old', FORM="FORMATTED")
		! tmax = time_max(state)*1.1/au2ns
		tmax = 3.74/au2ns								!ns -> a.u.
		w_min = (enes(state)-enes(1)) -0.2e-3/27.21139						!E=hbar*w -> w=E/hbar 
		w_max = (enes(state)-enes(1)) +0.2e-3/27.21139
		dw    = 0.005e-3/27.21139
	else if (flag_time_mod .eq. 2) then
		write(*,*) 'Sprawdzam szerokie spektrum czestotliwosci z duzym krokiem'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Status='unknown', FORM="FORMATTED")
		tmax = 3.74/au2ns
		w_min = 6.9/27211.39			!E=hbar*w -> w=E/hbar 
		! w_max = enes(nst)-enes(1) + 1.e-3/27.21139
		w_max = 7./27211.39
		dw    = 0.25e-3/27.21139
	else 
		write(*,*) 'Podano zla wartosc flag_time_mod:', flag_time_mod
	end if
	w     = w_min
	print*, 'w_min:', w_min*1.e3*27.21139
	print*, 'w_max:', w_max*1.e3*27.21139
	print*, 'tmax', tmax
	
	Frequency_Loop: do while (w .lt. w_max)
		c_max  = 0
		c_p    = 0		!c in previous step (t-dt)
		c_p(1) = 1
		c_c    = 0		!c in current step  (t)
		c_c(1) = 1
		c_n    = 0		!c in next step     (t+dt)
		t=2*dt
		it=0
		Time_Loop: do while(t .lt. tmax)
			it=it+1
			t = t+dt
			! C_l(t+dt)=C_l(t-dt)+2dt/i/\hbar \sum_k H'_lk(t) c_k(t) \exp(-i(El-Ek)/hbar*t)
			C_i_next_Loop: do i=1,nst
				c_n(i) = c_p(i)
				H_ij_sum_Loop: do j=1,nst
					c_n(i) = c_n(i) + 2*dt/ci*Ampl*sin(w*t)*mat_ij(i,j)*exp(ci*(enes(i)-enes(j))*t)*c_c(j)
				end do H_ij_sum_Loop
			end do C_i_next_Loop
			c_p = c_c
			c_c = c_n
			c_n = 0
			do i=1,nst
				if (cdabs(c_c(i))**2 .gt. c_max(i)) then
					c_max(i) = cdabs(c_c(i))**2
					time_max(i) = t*au2ns
				end if
				if (flag_time_mod .eq. 2) then
					if (c_max(i) .gt. 0.05) then
						which(i) = 1
					end if
				end if
			enddo

			if (flag_time_mod .eq. 0) then
				if (flag_time_half .eqv. .FALSE.) then
					if (c_max(state) .ge. 0.5) then
						flag_time_half = .TRUE.
						time_half(state) = t*au2ns
						write(*,*) 'czas po ktorym... 50%:', time_half(state)
					end if
				end if
				if (flag_time_max .eqv. .FALSE.) then
					if (c_max(state) .ge. 0.1) then
						flag_time_max = .TRUE.
						! time_max(state) = t*au2ns
						write(*,*) 'czas po ktorym... 10%:', time_max(state) 
						exit
					end if
				end if
			end if

			! if (mod(it,100) .eq. 0) then
			! write(11,*) w*1.e3*27.21139, t*au2ns, (cdabs(c_c(i))**2,i=1,nst)
			! end if
		end do Time_Loop

		write(22,'(F20.12,100E20.12)') w*1.e3*27.21139, c_max
		print*,  w*1.e3*27.21139, c_max(state)

		w=w+dw
	end do Frequency_Loop

	if (flag_time_mod .eq. 0) then
		do i=2,nst
			if (which(i) .eq. 1) then
				write(*,'(I3)',advance="no") i
			end if
		enddo
		write(*,*) ' '
	end if
	! close(11)
	close(22)
  
END SUBROUTINE function_time_AC_1el



!-----------------------------------------------------------------------------
!----------------------------Rachunek z czasem--------------------------------
!-------------------------------Askar-Cakmak----------------------------------
!---------------------------------multi el------------------------------------
SUBROUTINE function_time_AC_multi_el( nsz, numerod, enci, elmnt, F_ac, au2ns, &
           & flag_time_mod, state, which, time_half, time_max, flag_time_half, &
		   & flag_time_max, c_max )
  integer, INTENT(IN) :: nsz, numerod, state, flag_time_mod
  integer, INTENT(INOUT) :: which(1:numerod)
  REAL* 8, INTENT(IN) :: enci(numerod), F_ac, au2ns
  REAL* 8, INTENT(INOUT) :: time_half(1:numerod), time_max(1:numerod), c_max(1:numerod)
  COMPLEX*16, INTENT(IN) :: elmnt(numerod,numerod,3)
  logical, INTENT(INOUT) :: flag_time_half, flag_time_max
 
  REAL*8 :: t, w, dw, w_min, w_max, tmax, dt, pi
  integer :: i, j, k, it, iw, deb
  complex*16 :: ci, c_p(1:numerod), c_c(1:numerod), c_n(1:numerod)
  ! logical :: flag_time_half = .FALSE.	!duration of F_ac pulse upon which the probability of finding
									!the electron in another state exceeds 50%
									!.FALSE. -> time was not found yet
									!.TRUE. -> time was found
  
  ci   = (0,1)
  pi   = atan(1.0) *4
  dt   = 1.e-5/au2ns							!ns -> a.u.
  
  ! open(11, FILE="OUT/time.dat",ACTION="WRITE", STATUS="UNKNOWN", FORM="FORMATTED")
	if (flag_time_mod .eq. 0) then
		write(*,*) 'Sprawdzam konkretne czestotliwosci'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Access='append', Status='old', FORM="FORMATTED")
		! tmax = 10./au2ns								!ns -> a.u.
		tmax = 3.74/au2ns
		w_min = (enci(state)-enci(1)) 										!E=hbar*w -> w=E/hbar 
		w_max = (enci(state)-enci(1)) +0.2e-3/27.21139
		dw    = 3.e-3/27.21139			!doesnt matter, it just has to be bigger than 0.2e-3/27.21139
	else if (flag_time_mod .eq. 1) then
		write(*,*) 'Sprawdzam waskie spektrum czestotliwosci z malym krokiem'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Access='append',Status='old', FORM="FORMATTED")
		! tmax = time_max(state)*1.1/au2ns
		tmax = 3.74/au2ns
		w_min = (enci(state)-enci(1)) -0.2e-3/27.21139						!E=hbar*w -> w=E/hbar 
		w_max = (enci(state)-enci(1)) +0.2e-3/27.21139
		dw    = 0.005e-3/27.21139
	else if (flag_time_mod .eq. 2) then
		write(*,*) 'Sprawdzam szerokie spektrum czestotliwosci z duzym krokiem'
		open(22, FILE="OUT/prob_w.dat",ACTION="WRITE", Status='unknown', FORM="FORMATTED")
		! tmax = 3./au2ns
		tmax = 3.74/au2ns
		w_min = 0			!E=hbar*w -> w=E/hbar 
		! w_max = enci(8)-enci(1) + 1.e-3/27.21139
		w_max = 13./27211.39
		dw    = 0.5e-3/27.21139
	else 
		write(*,*) 'Podano zla wartosc flag_time_mod:', flag_time_mod
	end if
	
	w     = w_min
	print*, 'w_min:', w_min*1.e3*27.21139
	print*, 'w_max:', w_max*1.e3*27.21139
	print*, 'tmax', tmax
	Frequency_Loop: do while (w .lt. w_max)
		c_max  = 0
		c_p    = 0		!c in previous step (t-dt)
		c_p(1) = 1
		c_c    = 0		!c in current step  (t)
		c_c(1) = 1
		c_n    = 0		!c in next step     (t+dt)
		t=2*dt
		it=0
		Time_Loop: do while(t .lt. tmax)
			it=it+1
			t = t+dt
			! C_l(t+dt)=C_l(t-dt)+2dt/i/\hbar \sum_k H'_lk(t) c_k(t) \exp(-i(El-Ek)/hbar*t)
			C_i_next_Loop: do i=1,numerod
				c_n(i) = c_p(i)
				H_ij_sum_Loop: do j=1,numerod
					c_n(i) = c_n(i) + 2*dt/ci*F_ac*sin(w*t)*elmnt(i,j,2)*exp(ci*(enci(i)-enci(j))*t)*c_c(j)
				end do H_ij_sum_Loop
			end do C_i_next_Loop
			c_p   = c_c
			c_c   = c_n
			c_n   = 0
			do i=1,numerod
				if (cdabs(c_c(i))**2 .gt. c_max(i)) then
					c_max(i) = cdabs(c_c(i))**2
					time_max(i) = t*au2ns
				end if
				if (flag_time_mod .eq. 2) then
					if (c_max(i) .gt. 0.05) then
						which(i) = 1
					end if
				end if
			end do

			if (flag_time_mod .eq. 0) then
				if (flag_time_half .eqv. .FALSE.) then
					if (c_max(state) .ge. 0.5) then
						flag_time_half = .TRUE.
						time_half(state) = t*au2ns
						write(*,*) 'czas po ktorym... 50%:', time_half(state)
					end if
				end if
				if (flag_time_max .eqv. .FALSE.) then
					if (c_max(state) .ge. 0.98) then
						flag_time_max = .TRUE.
						! time_max(state) = t*au2ns
						write(*,*) 'czas po ktorym... 98%:', time_max(state) 
						exit
					end if
				end if
			end if

			! if (mod(it,100) .eq. 0) then
			! write(11,*) w*1.e3*27.21139, t*au2ns, (cdabs(c_c(i))**2,i=1,numerod)
			! end if
		end do Time_Loop

		write(22,'(100F20.12)') w*1.e3*27.21139, c_max
		print*,  w*1.e3*27.21139, c_max(state)

		w=w+dw
	end do Frequency_Loop

	if (flag_time_mod .eq. 0) then
		do i=2,numerod
			if (which(i) .eq. 1) then
				write(*,'(I3)',advance="no") i
			end if
		enddo
		write(*,*) ' '
	end if
	! close(11)
	close(22)
  
END SUBROUTINE function_time_AC_multi_el


END MODULE mod_time