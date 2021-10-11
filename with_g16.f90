subroutine chrg_server( code, natm, naqm, atmn, cord, ener, qfit, grad, hess )
	implicit none
	integer, intent( in ) :: code, natm, naqm
	real*8, dimension(1:natm), intent( in ) :: atmn
	real*8, dimension(1:3*natm), intent( in ) :: cord
	real*8, intent( inout ) :: ener
	real*8, dimension(1:naqm), intent( inout ) :: qfit
	real*8, dimension(1:*), intent( inout ) :: grad
	real*8, dimension(1:*), intent( inout ) :: hess

	integer :: i, j, k, l, n3, nh, ierr
	real*8 :: smm
	character( len=256 ) :: str

	real*8, parameter :: bohr = 0.529177249d0
	character( len=2 ), dimension(1:109), parameter :: smb = (/ &
			"H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", &
            "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", &
            "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", &
            "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", &
            "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", &
            "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", &
            "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", &
            "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt" /)


! The unique thing you must change is the charge and the QM method, Charge 1 and MP2/6-31+G(d,p) in this case
	open( unit=999, file = "calc.com", action = "write", form = "formatted" )
	write( 999, "(a/a/a/a)" ) "%Chk=calc", &
		"%Nproc=16", &
		"%Mem=4gb", &
		"#M062X/6-31+g(d,p) charge scf=(direct,conver=6) nosymm"
!               "#B3LYP/6-31+g(d) charge scf=(direct,conver=8) nosymm"
!               "#B3LYP/6-31g(d) charge scf=(direct,conver=8) nosymm"
!               "#B3LYP/3-21g charge scf=(direct,conver=6) nosymm"
!                "#PM7 charge scf=(direct,conver=8) nosymm"
!                "#RHF/3-21g charge scf=(direct,conver=8) nosymm"
	if( code == 1 ) then
		write(*,*) "cab_initio gradient evaluation"
		write( 999, "(a/)" ) "force pop=(CHelpG,ReadRadii)"
	else if( code == 2 ) then
		write(*,*) "cab_initio hessian evaluation"
		write( 999, "(a/)" ) "freq=noraman pop=(CHelpG,ReadRadii)"
	else 
		write(*,*) "cab_initio energy evaluation"
		write( 999, "(a/)" ) "pop=(CHelpG,ReadRadii)"
	end if
	write( 999, "(a//2i4)" ) "- slave -", -1, 1
	do i = 1, naqm
		j = ( i - 1 ) * 3
		write( 999, "(a8,3f20.10)" ) smb(int(atmn(i))), cord(j+1), cord(j+2), cord(j+3)
	end do
	write( 999, "(a)" ) ""
	smm = 0.0d0
	if( naqm == natm ) then
		write( 999, "(f28.10,2f20.10,f8.3)" ) 999.0d0, 999.0d0, 999.0d0, 0.0d0
	else
		do i = naqm + 1, natm
			j = ( i - 1 ) * 3
			write( 999, "(f28.10,2f20.10,f8.3)" ) cord(j+1), cord(j+2), cord(j+3), atmn(i)
			do k = i + 1, natm
				l = ( k - 1 ) * 3
				smm = smm + atmn(i) * atmn(k) / dsqrt( sum( ( cord(j+1:j+3) - cord(l+1:l+3) ) ** 2 ) ) * bohr
			end do
		end do
	end if
	write( 999, "(a)" ) ""
	write( 999, "(a)" ) "35 1.85"
	write( 999, "(//)" )
	close( 999 )

!	write(*,"(a,f20.10)") "self_charges:", smm

	call system( ". $g16root/g16/bsd/g16.profile; g16 calc; formchk calc.chk calc.fchk" )

!	write(*,"(a)") "Gaussian job has run"

	ener = 0.0d0
	qfit = 0.0d0
	open( unit = 999, file = "calc.fchk", action = "read", form = "formatted" ) 
	read( 999, "(a)", end = 999 ) str
	do while( .true. )
		if( str(1:12) == "Total Energy" ) read( str(50:71), "(f22.15)", end = 999 ) ener

		if( str(1:11) == "ESP Charges" ) then
			k = 1
			do i = 1, int( naqm / 5 )
				read( 999, "(a)", end = 999 ) str
				do j = 1, 5
					read( str((j-1)*16+1:j*16), "(f16.8)" ) qfit(k)
					k = k + 1
				end do
			end do
			i = mod( naqm, 5 )
			if( i > 0 ) then
				read( 999, "(a)", end = 999 ) str
				do j = 1, i
					read( str((j-1)*16+1:j*16), "(f16.8)" ) qfit(k)
					k = k + 1
				end do
			end if
		end if

		if( str(1:18) == "Cartesian Gradient" .and. ( code == 1 .or. code == 2 ) ) then
			n3 = 3 * naqm
			grad(1:n3) = 0.0d0
			k = 1
			do i = 1, int( n3 / 5 )
				read( 999, "(a)", end = 999 ) str
				do j = 1, 5
					read( str((j-1)*16+1:j*16), "(f16.8)" ) grad(k)
					k = k + 1
				end do
			end do
			i = mod( n3, 5 )
			if( i > 0 ) then
				read( 999, "(a)", end = 999 ) str
				do j = 1, i
					read( str((j-1)*16+1:j*16), "(f16.8)" ) grad(k)
					k = k + 1
				end do
			end if
		end if

		if( str(1:25) == "Cartesian Force Constants" .and. code == 2 ) then
			nh = 3 * naqm * ( 3 * naqm + 1 ) / 2
			hess(1:nh) = 0.0d0
			k = 1
			do i = 1, int( nh / 5 )
				read( 999, "(a)", end = 999 ) str
				do j = 1, 5
					read( str((j-1)*16+1:j*16), "(f16.8)" ) hess(k)
					k = k + 1
				end do
			end do
			i = mod( nh, 5 )
			if( i > 0 ) then
				read( 999, "(a)", end = 999 ) str
				do j = 1, i
					read( str((j-1)*16+1:j*16), "(f16.8)" ) hess(k)
					k = k + 1
				end do
			end if
		end if

		read( 999, "(a)", end = 999 ) str
	end do
999     continue
	close( 999 )

	ener = ener - smm
!	write(*,*) "ener =",ener
end subroutine
