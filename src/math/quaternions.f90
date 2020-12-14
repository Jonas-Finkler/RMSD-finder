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



module quaternions
    use atomicStructure
    use precision
    use linalg
    use constants
    implicit none

    real(dp), parameter :: unitQuaternion(4) = [1._dp,0._dp,0._dp,0._dp]

contains
    function buildRotationMatrix(q) result(Rot)
        implicit none
        real(dp), dimension(4), intent(in) :: q
        real(dp) :: a, b, c, d
        real(dp), dimension(3, 3) :: Rot

        a = q(1)
        b = q(2)
        c = q(3)
        d = q(4)

        Rot(1, 1) = a**2 + b**2 - c**2 - d**2
        Rot(2, 1) = 2 * b * c + 2 * a * d
        Rot(3, 1) = 2 * b * d - 2 * a * c

        Rot(1, 2) = 2 * b * c - 2 * a * d
        Rot(2, 2) = a**2 - b**2 + c**2 - d**2
        Rot(3, 2) = 2 * c * d + 2 * a * b

        Rot(1, 3) = 2 * b * d + 2 * a * c
        Rot(2, 3) = 2 * c * d - 2 * a * b
        Rot(3, 3) = a**2 - b**2 - c**2 + d**2

    end function

    ! rotation with quatMult(q,r) is same as if r was applied first and then q
    function quatMult(q, r) result(t) ! computes multiplication of two quaternions
        implicit none
        real(dp), dimension(4), intent(in) :: q, r
        real(dp), dimension(4) :: t

        t(1) = r(1) * q(1) - r(2) * q(2) - r(3) * q(3) - r(4) * q(4)
        t(2) = r(1) * q(2) + r(2) * q(1) - r(3) * q(4) + r(4) * q(3)
        t(3) = r(1) * q(3) + r(2) * q(4) + r(3) * q(1) - r(4) * q(2)
        t(4) = r(1) * q(4) - r(2) * q(3) + r(3) * q(2) + r(4) * q(1)
    end function

    subroutine quatToVectorAndAngle(q, r, a)
        implicit none
        real(dp), dimension(4), intent(in) :: q
        real(dp), dimension(3), intent(out) :: r
        real(dp), intent(out) :: a

        a = 2._dp * datan2(dsqrt(sum(q(2:4)**2._dp)), q(1))
        r = q(2:4) / sum(dsqrt(q(2:4)**2._dp))
    end subroutine

    subroutine quatFromVectorAndAngle(r, a, q)
        real(dp), intent(in) :: r(3)
        real(dp), intent(in) :: a
        real(dp), intent(out) :: q(4)

        q(2:4) = r / sqrt(sum(r**2)) * sin(a / 2._dp)
        q(1) = cos(a / 2._dp)

    end subroutine

    subroutine rotatePointSet(nat, ats, q)
        implicit none
        integer, intent(in) :: nat
        real(dp), intent(inout) :: ats(3,nat)
        real(dp), dimension(4), intent(in) :: q

        real(dp) :: R(3, 3)
        integer :: i

        R = buildRotationMatrix(q)

        do i = 1, nat
            ats(:, i) = matMulVec3(R, ats(:, i))
        end do
    end subroutine

    subroutine randomQuaternion(r)
        use random
        implicit none
        real(dp), dimension(4), intent(out) :: r

        call random_normal_1D(r)

        r = r / dsqrt(sum(r**2))

    end subroutine

    subroutine readQuats(filename, n, quats)
        implicit none
        character(len = *), intent(in) :: filename
        integer, intent(out) :: n
        real(dp), allocatable, intent(out) :: quats(:,:)
        integer :: filestat, unit, i

        open(newUnit = unit, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat /= 0) stop "could not open quaternion file"
        read(unit, *) n
        allocate(quats(4,n))
        do i=1,n
            read(unit, *) quats(:, i)
        end do
        close(unit)
    end subroutine

    subroutine as_quaternionFromAssignment(atsA, atsB, q)
        type(atStruct), intent(in) :: atsA, atsB
        real(dp), intent(out) :: q(4)

        call quaternionFromAssignment(atsA%nat, atsA%ats, atsB%ats, q)

    end subroutine

    subroutine quaternionFromAssignment(nat, atsA, atsB, q)
        integer, intent(in) :: nat
        real(dp), intent(in) :: atsA(3,nat), atsB(3,nat)
        real(dp), intent(out) :: q(4)
        real(dp) :: corr(3, 3), evals(4), evecs(4,4), F(4,4)
        integer :: i, j

        corr = 0._dp
        do i = 1, 3
            do j = 1, 3
                corr(i, j) = sum(atsA(i,:) * atsB(j,:))
            end do
        end do

        F(1, 1) =  corr(1, 1) + corr(2, 2) + corr(3, 3)
        F(2, 1) =  corr(2, 3) - corr(3, 2)
        F(3, 1) =  corr(3, 1) - corr(1, 3)
        F(4, 1) =  corr(1, 2) - corr(2, 1)
        F(2, 2) =  corr(1, 1) - corr(2, 2) - corr(3, 3)
        F(3, 2) =  corr(1, 2) + corr(2, 1)
        F(4, 2) =  corr(1, 3) + corr(3, 1)
        F(3, 3) = -corr(1, 1) + corr(2, 2) - corr(3, 3)
        F(4, 3) =  corr(2, 3) + corr(3, 2)
        F(4, 4) = -corr(1, 1) - corr(2, 2) + corr(3, 3)

        ! dont need upper half
        !F(1,2) = F(2,1)
        !F(1,3) = F(3,1)
        !F(1,4) = F(4,1)
        !F(2,3) = F(3,2)
        !F(2,4) = F(4,2)
        !F(3,4) = F(4,3)

        call eigensystemDiag(4, F, evals, evecs)

        q = evecs(:,4)

    end subroutine


end module quaternions