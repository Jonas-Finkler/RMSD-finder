! Copyright (C) 2020 Jonas A. Finkler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


program rmsdFinder
    use precision
    use atomicStructure
    use quaternions
    use rmsds
    use random
    use argumentParser
    use quaternionLists
    implicit none

    integer :: nquats
    type(atStruct) :: atsA, atsB
    integer, allocatable :: assignment(:)
    real(dp), allocatable :: quatList(:,:)
    real(dp) :: rmsd, q(4)
    logical :: inv
    real(dp) :: cmA(3), cmB(3)
    integer :: i
    real(dp) :: rMat(3,3), angle, axis(3)

    type(argumentList) :: argList
    type(intArgument) :: argnQuats
    type(stringArgument) :: argStructA
    type(stringArgument) :: argStructB
    type(logicalArgument) :: argInversion
    type(stringArgument) :: argOutFile
    type(logicalArgument) :: argMass
    type(stringArgument) :: argRotFormat

    character(len=*), parameter :: nl = new_line(' ')
    character(len=*), parameter :: programName = 'rmsdFinder'
    character(len=*), parameter :: version = '1.0'
    character(len=*), parameter :: helpText = &
      programName // ' (verions=' // version // '), written by Jonas A. Finkler' // nl // &
      ' '// nl // &
      ' The source code for this program can be found here: https://github.com/Jonas-Finkler/RMSD-finder'// nl // &
      ' '// nl // &
      ' If you use this code, please cite:'// nl // &
      '   The Journal of Chemical Physics 152.16 (2020): 164106.   ' // nl // &
      '   Finkler, Jonas A., and Stefan Goedecker. '// nl // &
      '   "Funnel hopping Monte Carlo: An efficient method to overcome broken ergodicity."'// nl // &
      '   https://doi.org/10.1063/5.0004106'// nl // &
      ' '// nl // &
      ' Usage: ' // programName //' -A<structure-A> -B<structure-B> [-n<num-rotations>] [-i<use-inversion>]'// nl // &
      '       [-o<out-structure>] [-m<use-mass-weighted-RMSD>] [-r<rotation-format>]'// nl // &
      '   <structure-A> is translated, rotated and permuted to minimize the RMSD to <structure-B>'// nl // &
      '       All structures must be provided in xyz coordinates.'// nl // &
      '   <num-rotations> defines the number of rotations that are tried. The default is 2000. '// nl // &
      '       Possible values are 100, 200, ..., 1900, 2000.'// nl // &
      '       A higher number increases the computation time but also the probability'// nl // &
      '       that the true minimal RMSD is found.'// nl // &
      '   <use-inversion> can either be T (true) or F (false).'// nl // &
      '       If true, inverted structures are also considered when searching the minimal RMSD.'// nl // &
      '   <out-structure> is the file to which the transformed copy of <structure-A> is saved'// nl // &
      '   <use-mass-weighted-RMSD> can either be T (true) or F (false).'// nl // &
      '       If true, the mass weighted RMSD definition is used: RMSD = sqrt((sum_i  m_i |R_i - r_i|^2) / N).'// nl // &
      '       Otherwise RMSD = sqrt((sum_i |R_i - r_i|^2) / N) is used.'// nl // &
      '   <rotation-format> specifies the format in which the rotation is printed.'// nl // &
      '       Can be either Q for quaternion (default), M for rotation matrix or A for angle-axis representation.'// nl // &
      ' '// nl // &
      ' Assignment convention: Atom number assignment(i) of structure A has been matched to atom i in structure B.'// nl // &
      ' '
      !' '// nl // &

    argStructA%name = 'A'
    argStructA%required = .true.

    argStructB%name = 'B'
    argStructB%required = .true.

    argnQuats%name = 'n'
    argnQuats%required = .false.
    argnQuats%default = 2000

    argInversion%name = 'i'
    argInversion%required = .false.
    argInversion%default = .false.

    argOutFile%name = 'o'
    argOutFile%required = .false.
    argOutFile%default = 'out.xyz'

    argMass%name = 'm'
    argMass%required = .false.
    argMass%default = .false.

    argRotFormat%name = 'r'
    argRotFormat%required = .false.
    argRotFormat%default = 'Q'

    call initArgumentParser(argList, programName, helpText, version)
    call addArgument(argList, argnQuats)
    call addArgument(argList, argStructA)
    call addArgument(argList, argStructB)
    call addArgument(argList, argInversion)
    call addArgument(argList, argOutFile)
    call addArgument(argList, argMass)
    call addArgument(argList, argRotFormat)

    call parseArguments(argList)

    nquats = argnQuats%value
    allocate(quatList(4,nquats))
    call getQuaternionList(nquats, quatList)

    call as_readXYZ(argStructA%value, atsA)
    call as_readXYZ(argStructB%value, atsB)

    if (.not. as_compare(atsA, atsB)) then
        print*, 'Error: Both structures must have same number of atoms for each element.'
        stop
    end if

    if (argMass%value) then
        call as_getCenterOfMass(atsA, cmA)
        call as_getCenterOfMass(atsB, cmB)
        call as_centerMass(atsA)
        call as_centerMass(atsB)
    else
        call as_getCenter(atsA, cmA)
        call as_getCenter(atsB, cmB)
        call as_center(atsA)
        call as_center(atsB)
    end if



    allocate(assignment(atsA%nat))
    if (argInversion%value) then
        call as_findRmsd(atsA, atsB, nquats, quatList, 0._dp, argMass%value, q, assignment, rmsd, inv)
    else
        call as_findRmsd(atsA, atsB, nquats, quatList, 0._dp, argMass%value, q, assignment, rmsd)
    end if

    print*, 'RMSD: ', rmsd
    print*, 'Translation: ', cmB - cmA
    select case (argRotFormat%value)
        case ('Q')
            print*, 'Rotation [quaternion]: ', q
        case ('M')
            rMat = buildRotationMatrix(q)
            print*, 'Rotation [Matrix]: '
            print*, '    ', rMat(1,:)
            print*, '    ', rMat(2,:)
            print*, '    ', rMat(3,:)
        case ('A')
            call quatToVectorAndAngle(q, axis, angle)
            print*, 'Rotation [angle, axis]: ', angle, ',   ', axis
        case default
            print*, 'Error: Unknown rotation format "' // argRotFormat%value // '"'
    end select
    print*, 'Assignment: ', assignment



    if (argOutFile%wasFound) then
        call as_reassign(atsA, assignment)
        call rotatePointSet(atsA%nat, atsA%ats, q)
        do i=1,atsA%nat
            atsA%ats(:,i) = atsA%ats(:,i) - cmA + cmB
        end do
        call as_writeXYZ(argOutFile%value, atsA)
        print*, 'Transformed structure (' // argStructA%value // ') has been written to: ' // argOutFile%value
    end if



end program rmsdFinder
