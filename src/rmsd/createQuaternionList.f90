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

program main
    use precision
    use pointsets
    implicit none
    integer, parameter :: dim=4
    integer :: n
    real(dp), allocatable :: quats(:,:), quatsout(:,:)
    character(len=200) :: fname
    integer :: i
    real(dp), dimension(2,2) :: tt = reshape([1._dp, 1._dp, 1._dp, 1._dp], shape(tt))

    do n=100,2000,100
        allocate(quats(dim, n))
        allocate(quatsout(dim, 2*n))
        call createQuatList(dim, n, quats)

        call sortQuats(dim, n, quats)
        if (dim==3) then
            do i=1,n
                write(fname,'(A,I5.5,A)') 'out/', i, '.xyz'
                quatsout(:,:i) = quats(:,:i)
                quatsout(:,i+1:2*i) = -1._dp * quats(:,:i)
                call writeXYZ(fname, 2*i, quatsout(:,:2*i))
            end do
        end if

        write(fname, '(A,I4.4,A)') 'quats/', n, '.quats'
        open(35, file=fname, action='write')
        write(35, *) n
        do i=1,n
            write(35, '(4F22.18)') quats(:,i)
        end do
        close(35)
        deallocate(quats)
        deallocate(quatsout)
    end do

end program main

subroutine sortQuats(dim, n, quats)
    use precision
    integer, intent(in) :: dim
    integer, intent(in) :: n
    real(dp), intent(inout) :: quats(dim, n)

    integer :: i, j, k
    integer :: order(n)
    real(dp) :: tmpquats(dim, n)
    integer :: best, tmp
    real(dp) :: d, dd, dbest


    do i=1,n
        order(i) = i
    end do

    i = 2
    do i=2,n
        ! pick the one furthest away from all the ones already picked.
        dbest = huge(dbest) !0._dp
        do j=i,n ! over candidates
            d = 0._dp!huge(d)
            do k=1,i-1 ! the ones already selected
                !dd = min(sum((quats(:,order(j)) - quats(:,order(k)))**2),sum((quats(:,order(j)) + quats(:,order(k)))**2))
                !if (dd < d) then
                !    d = dd
                !end if
                d = d + 1._dp / sqrt(sum((quats(:,order(j))-quats(:,order(k)))**2)) &
                        + 1.d0 / sqrt(sum((quats(:,order(j))+quats(:,order(k)))**2))
            end do
            if (d < dbest) then
                dbest = d
                best = j
            end if
        end do
        tmp = order(i)
        order(i) = order(best)
        order(best) = tmp
    end do
    print*, order

    tmpquats = quats
    do i=1,n
        quats(:,i) = tmpquats(:,order(i))
    end do

end subroutine


subroutine createQuatList(dim, n, quats)
    use precision
    use pointsets
    use random
    implicit none
    integer, intent(in) :: dim ! can be used to visually check 3D results
    integer, intent(in) :: n
    real(dp), intent(out) :: quats(dim, n)

    real(dp) :: quatsout(dim, 2*n)
    integer :: i, j
    real(dp) :: e, e2, f(dim, n), fproj(dim, n), dx(dim, n)
    real(dp) :: step
    character(len=200) :: fname

    call random_normal_2D(quats)
    call projectOnSphere(dim, n, quats)

    !call random_normal_2D(dx)
    !dx = dx * 1.e-8_dp
    !call getEnergyAndForce(dim, n, quats, e, f)
    !quats = quats + dx
    !call getEnergyAndForce(dim, n, quats, e2, f)
    !print*, e2 - e, -1._dp * sum(f*dx)

    step = 1.e-5_dp
    call getEnergyAndForce(dim, n, quats, e, f)
    do j=1,n
        fproj(:,j) = f(:,j) - sum(f(:,j) * quats(:,j)) * quats(:,j)
    end do
    e2 = e + 100._dp
    i = 0
    do while(abs(e2-e) > 1.e-5)

        quats = quats + step * fproj
        call projectOnSphere(dim, n, quats)
        e2 = e
        call getEnergyAndForce(dim, n, quats, e, f)
        do j=1,n
            fproj(:,j) = f(:,j) - sum(f(:,j) * quats(:,j)) * quats(:,j)
        end do
        print*, n, i, step, e2-e, sqrt(sum(fproj**2)), e

        if (e < e2) then
            step = step * 1.05
        else
            step = step * 0.5_dp
        end if

        !if (dim==3) then
        !    write(fname,'(A,I5.5,A)') 'out/', i, '.xyz'
        !    quatsout(:,:n) = quats
        !    quatsout(:,n+1:) = -1._dp * quats
        !    call writeXYZ(fname, 2*n, quatsout)
        !end if
        i = i + 1

    end do
    call projectOnSphere(dim, n, quats)


end subroutine createQuatList

subroutine projectOnSphere(dim, n, quats)
    use precision
    implicit none
    integer, intent(in) :: dim
    integer, intent(in) :: n
    real(dp), intent(inout) :: quats(dim, n)
    integer :: i

    do i=1,n
        quats(:,i) = quats(:,i) / sqrt(sum(quats(:,i)**2))
    end do

end subroutine

subroutine getEnergyAndForce(dim, n, quats, e, f)
    use precision
    implicit none
    integer, intent(in) :: dim
    integer, intent(in) :: n
    real(dp), intent(in) :: quats(dim, n)
    real(dp), intent(out) :: e
    real(dp), intent(out) :: f(dim, n)

    real(dp) :: d, dd, ft(dim), ftp(dim)
    integer :: i, j

    f = 0._dp
    e = 0._dp

    do i=1,n
        do j=i+1,n
            d = sqrt(sum((quats(:,i)-quats(:,j))**2))
            dd = sqrt(sum((quats(:,i)+quats(:,j))**2))
            e = e + 1._dp / d + 1.d0 / dd
            ft(:) = (quats(:,i)-quats(:,j)) / d**3
            ftp(:) = (quats(:,i)+quats(:,j)) / dd**3
            f(:,i) = f(:,i) + ft(:) + ftp(:)
            f(:,j) = f(:,j) - ft(:) + ftp(:)
        end do
    end do


end subroutine getEnergyAndForce