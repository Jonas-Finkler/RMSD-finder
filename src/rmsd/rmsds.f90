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



module rmsds
    use precision
    use quaternions
    use assignmentProblem
    use atomicStructure
    use constants
    use util
    use random
    use linalg
    implicit none

contains

    ! Minimises the RMSD ->         RMSD = sqrt((sum_i      |R_i - r_i|^2/N)
    ! If massWeighted == .true. ->  RMSD = sqrt((sum_i  m_i |R_i - r_i|^2)/N)

    ! convention: atsA (first argument) is the one being rotated and permutated, while the second one is kept fixed
    ! convention: atom assignment(i) of structure A is matched to atom i of structure B

    ! IMPORTANT: structures have to be centered
    subroutine as_findRmsd(atsA, atsB, nquats, quats, eps, massWeighted, q, assignment, rmsd, inversion)
        type(atStruct), intent(in) :: atsA, atsB
        integer, intent(in) :: nquats
        real(dp), intent(in) :: quats(4, nquats)
        real(dp), intent(in) :: eps ! stop the search if an rmsd < eps has been found
        logical, intent(in) :: massWeighted
        real(dp), intent(out) :: q(4)
        integer, intent(out) :: assignment(atsA%nat)
        real(dp), intent(out) :: rmsd
        logical, intent(out), optional :: inversion

        type(atStruct) :: atsAA, atsBB, atsBBInv
        logical :: inv
        integer :: quatOrder(nquats)
        integer :: i
        real(dp) :: tmpRmsd, tmpQ(4)
        integer :: tmpAssignment(atsA%nat)

        if (.not. as_compare(atsA, atsB)) then
            stop 'can only find RMSD alignment for structures with same number of atoms of each element'
        end if

        do i=1,nquats
            quatOrder(i) = i
        end do

        call shuffle(nquats, quatOrder)

        rmsd = huge(rmsd)

        call as_copy(atsA, atsAA)
        call as_copy(atsB, atsBB)

        if (massWeighted) then
            do i=1,atsAA%nat
                atsAA%ats(:,i) = atsAA%ats(:,i) * sqrt(getAtomicMass(atsAA%el(i)))
                atsBB%ats(:,i) = atsBB%ats(:,i) * sqrt(getAtomicMass(atsBB%el(i)))
            end do
        end if

        if (present(inversion)) then
            call as_copy(atsBB, atsBBInv)
            atsBBInv%ats = -1._dp * atsBBInv%ats
        end if

        do i=1,nquats

            tmpQ(:) = quats(:,quatOrder(i))
            call as_iHung(atsAA, atsBB, tmpQ, tmpAssignment, tmpRmsd)
            if(tmpRmsd < rmsd) then
                rmsd = tmpRmsd
                assignment = tmpAssignment
                q = tmpQ
                inv = .false.
                if(rmsd < eps) exit
            end if
            if(present(inversion)) then
                call as_iHung(atsAA, atsBBInv, tmpQ, tmpAssignment, tmpRmsd) !todo: recalculation of -1._dp * atsB is unnecessary
                if(tmpRmsd < rmsd) then
                    rmsd = tmpRmsd
                    assignment = tmpAssignment
                    q = tmpQ
                    inv = .true.
                    if(rmsd < eps) exit
                end if
            end if
        end do

        if(present(inversion)) inversion = inv

    end subroutine

    subroutine as_iHung(atsA, atsB, q, assignment, rmsd)
        type(atStruct), intent(in) :: atsA, atsB
        real(dp), intent(inout) :: q(4)
        real(dp), intent(out) :: rmsd
        integer, intent(out) :: assignment(atsA%nat)
        integer :: lastAssignment(atsA%nat), i
        real(dp) :: qTmp(4)
        type(atStruct) :: atsTmp
        rmsd = huge(rmsd)
        assignment = 1
        lastAssignment = 0

        call as_copy(atsA, atsTmp)

        i = 0
        do while (.not. arrayCompare(assignment, lastAssignment))
            lastAssignment = assignment
            atsTmp%ats = atsA%ats
            atsTmp%el = atsA%el
            call rotatePointSet(atsTmp%nat, atsTmp%ats, q)
            call as_getOptimalAssignment(atsTmp, atsB, assignment, rmsd)
            call as_reassign(atsTmp, assignment)
            call as_quaternionFromAssignment(atsTmp, atsB, qTmp)
            q = quatMult(qTmp, q)
            i = i + 1
        end do
        !write(*,*) "iHung iterations: ", i

    end subroutine as_iHung

    subroutine as_getOptimalAssignment(atsA, atsB, assignment, rmsd)
        type(atStruct), intent(in) :: atsA, atsB
        integer, intent(out) :: assignment(atsA%nat)
        real(dp), intent(out) :: rmsd
        !real(dp) :: cost(atsA%nat,atsA%nat)
        integer :: i, j, ic, jc
        integer :: nel
        integer, allocatable :: els(:), elcount(:)
        real(dp), allocatable :: cost(:,:)
        integer :: iel
        real(dp) :: tmprmsd
        integer, allocatable :: tmpassignment(:), elatsA(:), elAtsB(:)

        rmsd = 0._dp
        assignment = -1

        if (atsA%nat /= atsB%nat) stop 'Error: Not same number of atoms in getOptimalAssignment'

        call getDistinctElements(atsA%nat, atsA%el, nel, els, elcount)

        ! we solve on assignment problem for each element
        do iel=1,nel
            allocate(cost(elcount(iel), elcount(iel)))
            allocate(tmpassignment(elcount(iel)))
            allocate(elatsA(elcount(iel)))
            allocate(elatsB(elcount(iel)))

            ! find which atoms are of that element
            ic = 0
            jc = 0
            do i=1,atsA%nat
                if (atsA%el(i) == els(iel)) then
                    ic = ic + 1
                    elatsA(ic) = i
                end if
                if (atsB%el(i) == els(iel)) then
                    jc = jc + 1
                    elatsB(jc) = i
                end if
            end do

            ! build cost matrix
            do j = 1,elcount(iel)
                do i = 1, elcount(iel)
                    cost(i,j) = sum((atsB%ats(:,elatsB(i)) - atsA%ats(:,elatsA(j)))**2)
                end do
            end do

            ! solve assignment problem
            call solveAP(elcount(iel), cost, tmpassignment, tmprmsd)

            rmsd = rmsd + tmprmsd

            ! put result back to structure wide assignment array
            do i = 1,elcount(iel)
                assignment(elatsB(i)) = elAtsA(tmpassignment(i))
            end do

            deallocate(cost)
            deallocate(tmpassignment)
            deallocate(elAtsA)
            deallocate(elAtsB)
        end do

        rmsd = sqrt(rmsd / atsA%nat)

    end subroutine




end module rmsds