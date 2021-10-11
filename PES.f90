program potential_energy_scan
! include one or more water molecules in QM region
        use dynamo
        implicit none
        logical, dimension(:), allocatable      :: acs, flg
        character( len=6 )                      :: si, sj
        integer                                 :: a1, a2, a3, a4, a5,a6, i, j, k, my_random
        real( kind=dp ), dimension(1:2)         :: coef
        call dynamo_header


        call mm_file_process ( 'borra', 'opls' )
        call mm_system_construct ( 'borra', 'seq1' )
! edit the next line to specify the initial coordinates
        call coordinates_read ( "jack" )
        allocate( acs(1:natoms), flg(1:natoms) )
        acs = .false.
        acs = atom_selection( &
                                subsystem      = (/ "C" /), &
                                residue_number = (/ 1 /), &
                                residue_name   = (/ "MTA" /) )
        ! change 999 to the number of the water molecule to include
        acs = acs .or. atom_selection( &
                                subsystem      = (/ "A" /), &
                                residue_number = (/ 197 /), &
                                residue_name   = (/ "ASP" /) )
        ! add more lines as necessary to include more waters in QM region
        acs = acs .or. atom_selection( &
                                subsystem      = (/ "A" /), &
                                residue_number = (/ 151 /), &
                                residue_name   = (/ "PHE" /) )
        ! add more lines as necessary to include more waters in QM region
        acs = acs .or. atom_selection( &
                                subsystem      = (/ "A" /), &
                                residue_number = (/ 152 /), &
                                residue_name   = (/ "VAL" /) )
        ! add more lines as necessary to include more waters in QM region
        acs = acs .or. atom_selection( &
                                subsystem      = (/ "A" /), &
                                residue_number = (/ 76 /), &
                                residue_name   = (/ "SER" /) )
        ! add more lines as necessary to include more waters in QM region
        acs = acs .or. atom_selection( &
                                subsystem      = (/ "A" /), &
                                residue_number = (/ 174 /), &
                                residue_name   = (/ "GLU" /) )


        flg = .false.
        call my_sele_qmnb( flg )
!        flg = flg .and. .not. acs
        call atoms_fix ( .not. flg )
        write(*,*) "Flexible atoms", count( flg ), "Frozen Atoms", count ( .not. flg )

        call mopac_setup ( &
                method    = "AM1", &
                charge    = -1, &
                selection = acs )
        call mopac_scf_options ( iterations = 500000 )
        call energy_initialize
        call energy_non_bonding_options ( &
                list_cutoff   = 15.5_dp, &
                outer_cutoff  = 13.0_dp, &
                inner_cutoff  = 12.5_dp, &
                minimum_image = .true. )
 

! Define atoms that make up the constraints

        a1 = atom_number( subsystem = 'A', residue_number = 197, atom_name = 'OD2' )
        a2 = atom_number( subsystem = 'C', residue_number = 1, atom_name = 'H32' )
        a3 = atom_number( subsystem = 'C', residue_number = 1, atom_name = 'H32' )
        a4 = atom_number( subsystem = 'C', residue_number = 1, atom_name = 'N14' )
        a5 = atom_number( subsystem = 'C', residue_number = 1, atom_name = 'C11' )
        a6 = atom_number( subsystem = 'C', residue_number = 1, atom_name = 'N12' )

! Define weights of multiple distance

        coef = (/ 1._dp, -1._dp /)
        
! Run first constraint loop

do i = 0, 12

        ! Convett i into a string

        call encode_integer( i, si, "(i3)" )
        
        ! Run second constraint loop

        do j = 0, 16

                ! Convett i into a string

                call encode_integer( j, sj, "(i3)" )

                call constraint_initialize

                call constraint_point_define( atom_selection( atom_number = (/ a1 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a2 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a3 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a4 /) ) )
                call constraint_define( &
                        type = 'MULTIPLE_DISTANCE', &
                        eq   = 1.28_dp - 0.2_dp*real(i - 2, kind=dp), &
                        fc   = 3500._dp, &
                        weights = coef )

                call constraint_point_define( atom_selection( atom_number = (/ a5 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a6 /) ) )
                call constraint_define( &
                        type = 'DISTANCE', &
                        eq   = 1.460_dp + 0.1_dp*real(j - 2, kind=dp), &
                        fc   = 3500._dp )


                call gradient
                do k = 1, 100
                        if( grms < 10._dp ) cycle
                        atmcrd = atmcrd - atmder / dsqrt( sum( atmder ** 2 ) ) * 0.1_dp
                        call gradient
                end do

                call optimize_lbfgsb( &
                        step_number = 100, &
                        gradient_tolerance = 1._dp, &
                        print_frequency = 1 )

                call coordinates_write( "pes." // trim( si ) // "." // trim( sj ) )
                

                write( 900, "(3f20.10)" ) &
                        geometry_distance( atmcrd, a1, a2) - geometry_distance( atmcrd, a3, a4 ), &
                        geometry_distance( atmcrd, a5, a6), &
                        etotal
                call flush( 900 )

        end do

                call coordinates_read( "pes."  // trim( si ) // ".0" )
                write( 900, "(/)" )

end do
call dynamo_footer

end program
include "in15.f90"
