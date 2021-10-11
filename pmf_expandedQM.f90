program potential_mean_force
! include one or more water molecules in QM region
        use dynamo
        implicit none
        logical, dimension(:), allocatable      :: acs, flg
        character( len=6 )                      :: si, sj
        integer                                 :: a1, a2, a3, a4, a5, a6, i, j, my_random
        call dynamo_header
        call mm_file_process ( 'borra', 'opls' )
        call mm_system_construct ( 'borra', 'seq1' )
! edit the next line to specify the initial coordinates
        read(*,*) i,j
                call encode_integer ( i, si, '(i3)' )
                call encode_integer ( j, sj, '(i3)' )
        call coordinates_read ( "pes." // trim( si ) // '.' // trim( sj ) )
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
!        flg = .false.
 !       flg = atom_selection( &
  !                              subsystem      = (/ "D" /), &
   !                             residue_number = (/ 1 /), &
    !                            residue_name   = (/ "MTA" /) ) 

        flg = .false.
        call my_sele_qmnb( flg )
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
        a1 = atom_number ( &
                atom_name      = "OD2", &
                residue_number = 197, &
                subsystem      = "A" )
        a2 = atom_number ( &
                atom_name      = "H32", &
                residue_number = 1, &
                subsystem      = "C" )
        a3 = atom_number ( &
                atom_name      = "H32", &
                residue_number = 1, &
                subsystem      = "C" )
        a4 = atom_number ( &
                atom_name      = "N14", &
                residue_number = 1, &
                subsystem      = "C" )
        a5 = atom_number ( &
                atom_name      = "C11", &
                residue_number = 1, &
                subsystem      = "C" )
        a6 = atom_number ( &
                atom_name      = "N12", &
                residue_number = 1, &
                subsystem      = "C" )
                call constraint_initialize
                call constraint_point_define( atom_selection( atom_number = (/ a1 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a2 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a3 /) ) )
                call constraint_point_define( atom_selection( atom_number = (/ a4 /) ) )
                call constraint_define( type = "MULTIPLE_DISTANCE", fc = 2800.0_dp, &
                                          eq = 1.66_dp - ( i )* 0.20_dp, weights = (/ 1.0_dp, -1.0_dp /), &
                                        file = "dat_OHN." // trim( si ) // "." // trim(sj))

                call constraint_point_define ( atom_selection ( atom_number = (/ a5 /) ) )
                call constraint_point_define ( atom_selection ( atom_number = (/ a6 /) ) )
                call constraint_define ( &
                        type        = "DISTANCE", &
                        fc          = 2800.0_dp, &
                        eq          = 1.370_dp + ( j - 1 )* 0.10_dp, &
                        file        = "dat_CN." // trim ( si ) // "." // trim(sj) )
                call gradient
! Relax MD
                call random_initialize ( my_random() + i )
                call velocity_assign ( 300._dp, .false. )
                call dynamics_options ( &
                        time_step       = .001_dp, &
                        print_frequency = 100, &
                        steps           = 5000 )
                call langevin_verlet_dynamics ( 300._dp, 100._dp )
! Production MD
                call dynamics_options ( &
                        time_step       = .001_dp, &
                        print_frequency = 500, &
                        steps           = 10000 )
                call constraint_writing_start
                call langevin_verlet_dynamics ( 300._dp, 100._dp )
                call constraint_writing_stop
                call coordinates_write ( 'm_' // trim( si ) // '.' // trim( sj ) // ".crd" )
                call velocity_write ( 'm_' // trim( si ) // '.' // trim( sj ) // ".vel" )
                call pdb_write( 'm_' // trim( si ) // '.' // trim( sj ) // '.pdb' )
        deallocate( acs )
                call dynamo_footer
end program
include "in15.f90"
