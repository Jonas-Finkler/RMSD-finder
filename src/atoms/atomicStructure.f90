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



! convention: start subroutines acting on atomic structures with 'as_'
! this is done because the many subroutines have to be rewritten to use the new datatype
! backward compatibility has to be kept
! todo: maybe rethink or refactor this in the future
module atomicStructure
    use precision
    use constants
    use util

    implicit none

    ! type to store atomic structures
    type atStruct
        integer :: nat
        logical :: periodic = .false.
        real(dp)  :: lat(3,3)
        real(dp), allocatable :: ats(:,:)
        integer, allocatable :: el(:)
        ! only used with volumetric data
        real(dp), allocatable :: nValenceElectrons(:)
        ! only used with density fitting
        real(dp), allocatable :: rnuc(:)
    end type atStruct

contains

    subroutine as_center(ats)
        type(atStruct) :: ats
        real(dp) :: cm(3)
        integer :: i

        cm = sum(ats%ats, 2) / ats%nat
        do i = 1, ats%nat
            ats%ats(:, i) = ats%ats(:, i) - cm(:)
        end do

    end subroutine

    subroutine as_getCenter(ats, c)
        type(atStruct) :: ats
        real(dp), intent(out) :: c(3)
        c = sum(ats%ats, 2) / ats%nat
    end subroutine

    subroutine as_getCenterOfMass(ats, c)
        type(atStruct) :: ats
        real(dp), intent(out) :: c(3)
        integer :: i
        real(dp) :: mtot, m

        c(:) = 0._dp
        mtot = 0._dp
        do i=1,ats%nat
            m = getAtomicMass(ats%el(i))
            c(:) = c(:) + ats%ats(:,i) * m
            mtot = mtot + m
        end do

        c(:) = c(:) / mtot
    end subroutine

    subroutine as_centerMass(ats)
        type(atStruct) :: ats
        real(dp) :: c(3)

        call as_getCenterOfMass(ats, c)

        ats%ats(1,:) = ats%ats(1,:) - c(1)
        ats%ats(2,:) = ats%ats(2,:) - c(2)
        ats%ats(3,:) = ats%ats(3,:) - c(3)

    end subroutine

    subroutine as_copy(from, to)
        type(atStruct), intent(in) :: from
        type(atStruct), intent(out) :: to

        to%nat = from%nat
        if (allocated(to%ats)) then
            deallocate(to%ats)
        end if
        if (allocated(from%ats)) then
            allocate(to%ats(3,to%nat))
            to%ats(:,:) = from%ats(:,:)
        end if

        if (allocated(to%el)) then
            deallocate(to%el)
        end if
        if (allocated(from%el)) then
            allocate(to%el(to%nat))
            to%el(:) = from%el(:)
        end if

        if (allocated(to%nValenceElectrons)) then
            deallocate(to%nValenceElectrons)
        end if
        if (allocated(from%nValenceElectrons)) then
            allocate(to%nValenceElectrons(to%nat))
            to%nValenceElectrons(:) = from%nValenceElectrons(:)
        end if

        to%periodic = from%periodic
        to%lat = from%lat

    end subroutine

    subroutine as_readXYZ(filename, ats)
        character(len = *), intent(in) :: filename
        type(atStruct), intent(out) :: ats
        integer :: i, filestat
        character(len=2) :: el

        open(41, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat == 0) then
            !write(*,*) "sucessfully opened file"
            read(41, *) ats%nat
            read(41, *)
            allocate(ats%ats(3, ats%nat))
            allocate(ats%el(ats%nat))
            do i = 1, ats%nat
                read(41, *) el, ats%ats(:, i)
                ats%el(i) = elemSymToNum(el)
            end do
            ats%periodic = .false.
            close(41)
        else
            print*, 'Error: Could not open file (' // filename // ')'
            stop
        end if
    end subroutine as_readXYZ

    subroutine as_writeXYZ(filename, ats)
        character(len = *), intent(in) :: filename
        type(atStruct), intent(in) :: ats

        integer :: i, filestat

        open(41, file = filename, iostat = filestat, action = "write")
        if(filestat == 0) then
            write(41, *) ats%nat
            write(41, *)
            do i = 1, ats%nat
                ! todo: test format statement
                write(41, '(A2,3E15.6)') elemNumToSym(ats%el(i)), ats%ats(:, i)
            end do
            close(41)
        else
            stop "could not open file"
        end if

    end subroutine as_writeXYZ

    subroutine as_elemCount(ats, nelem)
        type(atStruct), intent(in) :: ats
        integer, intent(out) :: nelem(maxElemNum)
        integer :: i

        nelem(:) = 0

        do i=1,ats%nat
            nelem(ats%el(i)) = nelem(ats%el(i)) + 1
        end do

    end subroutine

    subroutine as_reassign(ats, assignment)
        type(atStruct), intent(inout) :: ats
        integer, intent(in) :: assignment(ats%nat)
        real(dp) :: tmp(3, ats%nat)
        integer :: tmpEl(ats%nat)
        integer :: i

        tmp = ats%ats
        tmpEl = ats%el

        do i=1,ats%nat
            ats%ats(:,i) = tmp(:,assignment(i))
            ats%el(i) = tmpEl(assignment(i))
        end do

    end subroutine as_reassign

    subroutine as_reassignInv(ats, assignment)
        type(atStruct), intent(inout) :: ats
        integer, intent(in) :: assignment(ats%nat)
        real(dp) :: tmp(3, ats%nat)
        integer :: tmpEl(ats%nat)
        integer :: i

        tmp = ats%ats
        tmpEl = ats%el

        do i=1,ats%nat
            ats%ats(:,assignment(i)) = tmp(:,i)
            ats%el(assignment(i)) = tmpEl(i)
        end do

    end subroutine as_reassignInv

    ! determines if two structures are identical up to the atomic coordinates
    function as_compare(atsA, atsB) result(ret)
        type(atStruct) :: atsA, atsB
        logical :: ret
        integer :: nElA(maxElemNum), nElB(maxElemNum)

        ret = .true.

        if (atsA%nat /= atsB%nat) then
            ret = .false.
            return
        end if

        call as_elemCount(atsA, nElA)
        call as_elemCount(atsB, nElB)

        if (.not. arrayCompare(nElA, nElB)) then
            ret = .false.
            return
        end if

    end function as_compare



end module atomicStructure