program potential_mean_force
! include one or more water molecules in QM region
        use dynamo
        implicit none
        logical, dimension(:), allocatable      :: acs, flg
        logical                                 :: skip_abinitio
        character( len=6 )                      :: si, sj
        integer                                 :: a1, a2, a3, a4, i, j, my_random
        double precision                        :: eq
        character(len=20)                       :: opls, res, atm_1, atm_2, atm_3, atm_4
       
        read(*,*) i,j 
        call encode_integer ( i, si, '(i3)' )
        call encode_integer ( j, sj, '(i3)' )
        call dynamo_header
        call mm_file_process ( 'borra', 'opls' )
        call mm_system_construct ( 'borra', 'seq1' )
! edit the next line to specify the initial coordinates
        call coordinates_read ( "pes." // trim(si) // "." // trim(sj) )
        allocate( acs(1:natoms), flg(1:natoms) )
        acs = .false.
        flg = .false.
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

!QM Setup
        call cabinitio_setup ( acs )
        write(*,*) 'Number of quantum atoms', count ( acs )
        skip_abinitio = .false.
        call energy_initialize
        call energy_non_bonding_options ( &
                list_cutoff   = 15.5_dp, &
                outer_cutoff  = 13.0_dp, &
                inner_cutoff  = 12.5_dp, &
                minimum_image = .true. )

        call energy

        deallocate( acs )
        call dynamo_footer
end program
include "with_g16.f90"
include "in15.f90"
