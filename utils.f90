module utils

    ! this module saves the subroutines which are less fault dependent

    use const_module
    implicit none

contains

    subroutine locate(ibin, edge,  x)
        ! note edge is the lower bound of a bin
        integer :: ibin, n_edge, i
        real(8) :: x
        real(8), allocatable :: edge(:)
        n_edge = size(edge)
        ibin = 1
        if (x >= edge(n_edge)) then
            ibin = n_edge
            return
        end if

        do i = 1, n_edge - 1
            if (x < edge (i+1)) then
                ibin = i
                exit
                return

            end if
        end do


    end subroutine

    subroutine deg2km_simple(vn, ve,alat_sta, alon_sta, &
            alat_ref, alon_ref)

        real(8) :: ve , vn   !m2f: check dim(:)!m
        real(8) :: alat_sta, alon_sta, alat_ref, alon_ref

        vn = (alat_sta - alat_ref) * DEG2RAD * EARTH_R
        ve = (alon_sta - alon_ref) * DEG2RAD * &
            cos(0.5 * (alat_sta + alat_ref) * DEG2RAD) * EARTH_R

    end subroutine  deg2km_simple

    subroutine delaz2_km (y1, x1, y2, x2, delta, az)

        real(8), intent(out) :: az
        real(8), intent(out) :: delta
        real(8), intent(in) :: x1
        real(8), intent(in) :: x2
        real(8), intent(in) :: y1
        real(8), intent(in) :: y2
        ! Variable declarations
        real(8) :: dx
        real(8) :: dy
        !
        if ((abs(x1 - x2) < 1e-12) .and. (abs(y1 - y2) < 1e-12)) then
            stop 'the two points are too close to calc azimuth'
        end if

        dx = x2 - x1
        dy = y2 - y1

        ! az = mod(2*pi+atan2(dx,dy),(2*pi)); this one gives true north based
        ! azimuth
        az = atan2(dx,dy) ! note the az is not .true. north based azimuth
        delta = sqrt((dx)**2+(dy)**2)



    end subroutine  delaz2_km

    elemental function normcdf(x)

        real(8), intent(in) :: x
        real(8) :: normcdf

        normcdf = 0.5 * erfc(-x * SQRT2_INV)

    end function

    elemental function deltacdf(x)
        real(8), intent(in) :: x
        real(8) :: deltacdf
        if (x < 0.0) then
            deltacdf = 0.0

        else
            deltacdf = 1.0
        end if

    end function

    subroutine truncnormcdf(x, a, b, z)

        ! Arguments declarations
        real(8), intent(in) :: a
        real(8), intent(in) :: b
        real(8), intent(in) :: x
        real(8), intent(out) :: z    !m2f: check dim(:)!m
        ! Variable declarations
        real(8) :: Z0
        !

        Z0 = normcdf(b) - normcdf(a)
        if (x < a) then
            z = 0.0
        else if (x <= b) then
            z = (normcdf(x) - normcdf(a)) / Z0
        else
            z = 1.0
        end if


    end subroutine  truncnormcdf

    function M22DET (A) result (DET)

        double precision, dimension(2,2), intent(in)  :: A
        double precision :: DET


        DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

        return

    end function M22DET

    ! http://web.hku.hk/~gdli/UsefulFiles/matrix/m33det_f90.txt
    ! http://www.davidgsimpson.com/software/m33det_f90.txt

    !***********************************************************************************************************************************
    !  M33DET  -  Compute the determinant of a 3x3 matrix.
    !***********************************************************************************************************************************

    function M33DET (A) result (DET)

        double precision, dimension(3,3), intent(in)  :: A

        double precision :: DET


        DET = A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)



    end function M33DET

    !>
    subroutine pointLineSegDistance(a, b, x, dist)

        real(8), intent(out) :: dist

        real(8), intent(in) :: x(2)
        ! Variable declarations
        real(8) :: a(2)
        real(8), dimension(3,3) :: A0
        real(8) :: b(2)
        real(8) :: d_ab
        real(8) :: d_ax
        real(8) :: d_bx
        !

        d_ab = norm2(a-b)
        d_ax = norm2(a-x)
        d_bx = norm2(b-x)

        if (dot_product(a-b,x-b)*dot_product(b-a,x-a)>=0) then
            A0 = reshape([ a,1.0d0,b,1.0d0,x,1.0d0 ] , [ 3 , 3] )
            dist = abs(M33DET(A0))/d_ab
        else
            dist = min(d_ax, d_bx)
        end if

    end subroutine  pointLineSegDistance

    !>
    subroutine pointTriangleDistance(TRI1, TRI2, TRI3, P, dist)

    ! This subroutine is translated from matlab version in the link below
    ! https://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d
    !
    ! the copyright of the original version is copied hereby

    ! Copyright (c) 2009, Gwendolyn Fischer
    ! All rights reserved.

    ! Redistribution and use in source and binary forms, with or without
    ! modification, are permitted provided that the following conditions are
    ! met:

    !  * Redistributions of source code must retain the above copyright
    !    notice, this list of conditions and the following disclaimer.
    !  * Redistributions in binary form must reproduce the above copyright
    !    notice, this list of conditions and the following disclaimer in
    !    the documentation and/or other materials provided with the distribution

    !    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    !    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    !    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    !    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    !    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    !    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    !    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    !    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    !    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    !    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    !    POSSIBILITY OF SUCH DAMAGE.

        real(8), intent(out) :: dist
        real(8), intent(in) :: P(3)
        real(8), intent(in) :: TRI1(3), TRI2(3), TRI3(3)
        real(8) :: B0(3)
        real(8) :: D0(3)
        real(8) :: E0(3)
        real(8) :: E1(3)

        real(8) :: a
        real(8) :: b
        real(8) :: c
        real(8) :: d
        real(8) :: denom
        real(8) :: det
        real(8) :: e
        real(8) :: f
        real(8) :: invDet
        real(8) :: numer
        real(8) :: s
        real(8) :: sqrDistance
        real(8) :: t
        real(8) :: tmp0
        real(8) :: tmp1

        ! real(8), external :: dot3
        !
        ! calculate distance between a point and a triangle in 3D
        ! SYNTAX
        ! dist = pointTriangleDistance(TRI,P)
        ! [dist,PP0] = pointTriangleDistance(TRI,P)
        !
        ! DESCRIPTION
        ! Calculate the distance of a given point P from a triangle TRI.
        ! Point P is a row vector of the form 1x3. The triangle is a matrix
        ! formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
        ! dist = pointTriangleDistance(TRI,P) returns the distance of the point P
        ! to the triangle TRI.
        ! [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
        ! closest point PP0 to P on the triangle TRI.
        !
        ! Author: Gwendolyn Fischer
        ! Release: 1.0
        ! Release date: 09/02/02
        ! Release: 1.1 Fixed Bug because of normalization
        ! Release: 1.2 Fixed Bug because of typo in region 5 20101013
        ! Release: 1.3 Fixed Bug because of typo in region 2 20101014

        ! Possible extention could be a version tailored not to return the distance
        ! and additionally the closest point, but instead return only the closest
        ! point. Could lead to a small speed gain.

        ! Example:
        ! !! The Problem
        ! P0 = [0.5 -0.3 0.5];
        !
        ! P1 = [0 -1 0];
        ! P2 = [1 0 0];
        ! P3 = [0 0 0];
        !
        ! vertices = [P1; P2; P3];
        ! faces = [1 2 3];
        !
        ! !! The Engine
        ! [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0);
        !
        ! !! Visualization
        ! [x,y,z] = sphere(20);
        ! x = dist*x+P0(1);
        ! y = dist*y+P0(2);
        ! z = dist*z+P0(3);
        !
        ! figure
        ! hold all
        ! patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8);
        ! plot3(P0(1),P0(2),P0(3),'b*');
        ! plot3(PP0(1),PP0(2),PP0(3),'*g')
        ! surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
        ! view(3)

        ! The algorithm is based on
        ! "David Eberly, 'Distance Between Point and Triangle in 3D',
        ! Geometric Tools, LLC, (1999)"
        ! http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
        !
        !
        !        ^t
        !  \     |
        !   \reg2|
        !    \   |
        !     \  |
        !      \ |
        !       \|
        !        *P2
        !        |\
        !        | \
        !  reg3  |  \ reg1
        !        |   \
        !        |reg0\
        !        |     \
        !        |      \ P1
        ! -------*-------*------->s
        !        |P0      \
        !  reg4  | reg5    \ reg6



        ! Do some error checking


        ! ToDo: check for colinearity and/or too small triangles.


        ! rewrite triangle in normal form
        B0 = TRI1
        E0 = TRI2-B0
        !E0 = E0/sqrt(sum(E0.^2)); !normalize vector
        E1 = TRI3-B0
        !E1 = E1/sqrt(sum(E1.^2)); !normalize vector


        D0 = B0 - P
        a = dot3(E0,E0)
        b = dot3(E0,E1)
        c = dot3(E1,E1)
        d = dot3(E0,D0)
        e = dot3(E1,D0)
        f = dot3(D0,D0)

        det = a*c - b*b ! do we have to use abs here?
        s = b*e - c*d
        t = b*d - a*e

        ! Terible tree of conditionals to determine in which region of the diagram
        ! shown above the projection of the point into the triangle-plane lies.
        if ((s+t) <= det) then
            if (s < 0) then
                if (t < 0) then
                    !region4
                    if (d < 0) then
                        t = 0
                        if (-d >= a) then
                            s = 1
                            sqrDistance = a + 2*d + f
                        else
                            s = -d/a
                            sqrDistance = d*s + f
                        end if
                    else
                        s = 0
                        if (e >= 0) then
                            t = 0
                            sqrDistance = f
                        else
                            if (-e >= c) then
                                t = 1
                                sqrDistance = c + 2*e + f
                            else
                                t = -e/c
                                sqrDistance = e*t + f
                            end if
                        end if
                    end if !of region 4
                else
                    ! region 3
                    s = 0
                    if (e >= 0) then
                        t = 0
                        sqrDistance = f
                    else
                        if (-e >= c) then
                            t = 1
                            sqrDistance = c + 2*e +f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        end if
                    end if
                end if !of region 3
            else
                if (t < 0) then
                    ! region 5
                    t = 0
                    if (d >= 0) then
                        s = 0
                        sqrDistance = f
                    else
                        if (-d >= a) then
                            s = 1
                            sqrDistance = a + 2*d + f ! GF 20101013 fixed typo d*s ->2*d
                        else
                            s = -d/a
                            sqrDistance = d*s + f
                        end if
                    end if
                else
                    ! region 0
                    invDet = 1/det
                    s = s*invDet
                    t = t*invDet
                    sqrDistance = s*(a*s + b*t + 2*d) &
                        + t*(b*s + c*t + 2*e) + f
                end if
            end if
        else
            if (s < 0) then
                ! region 2
                tmp0 = b + d
                tmp1 = c + e
                if (tmp1 > tmp0) then ! minimum on edge s+t=1
                    numer = tmp1 - tmp0
                    denom = a - 2*b + c
                    if (numer >= denom) then
                        s = 1
                        t = 0
                        sqrDistance = a + 2*d + f ! GF 20101014 fixed typo 2*b -> 2*d
                    else
                        s = numer/denom
                        t = 1-s
                        sqrDistance = s*(a*s + b*t + 2*d) &
                            + t*(b*s + c*t + 2*e) + f
                    end if
                else ! minimum on edge s=0
                    s = 0
                    if (tmp1 <= 0) then
                        t = 1
                        sqrDistance = c + 2*e + f
                    else
                        if (e >= 0) then
                            t = 0
                            sqrDistance = f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        end if
                    end if
                end if !of region 2
            else
                if (t < 0) then
                    !region6
                    tmp0 = b + e
                    tmp1 = a + d
                    if (tmp1 > tmp0) then
                        numer = tmp1 - tmp0
                        denom = a-2*b+c
                        if (numer >= denom) then
                            t = 1
                            s = 0
                            sqrDistance = c + 2*e + f
                        else
                            t = numer/denom
                            s = 1 - t
                            sqrDistance = s*(a*s + b*t + 2*d) &
                                + t*(b*s + c*t + 2*e) + f
                        end if
                    else
                        t = 0
                        if (tmp1 <= 0) then
                            s = 1
                            sqrDistance = a + 2*d + f
                        else
                            if (d >= 0) then
                                s = 0
                                sqrDistance = f
                            else
                                s = -d/a
                                sqrDistance = d*s + f
                            end if
                        end if
                    end if
                    !end region 6
                else
                    ! region 1
                    numer = c + e - b - d
                    if (numer <= 0) then
                        s = 0
                        t = 1
                        sqrDistance = c + 2*e + f
                    else
                        denom = a - 2*b + c
                        if (numer >= denom) then
                            s = 1
                            t = 0
                            sqrDistance = a + 2*d + f
                        else
                            s = numer/denom
                            t = 1-s
                            sqrDistance = s*(a*s + b*t + 2*d) &
                                + t*(b*s + c*t + 2*e) + f
                        end if
                    end if !of region 1
                end if
            end if
        end if

        ! account for numerical round-off error
        if (sqrDistance < 0) then
            sqrDistance = 0
        end if

        dist = sqrt(sqrDistance)

        !    if (nargout>1) then
        !        PP0 = B + s*E0 + t*E1
        !    end if
    end subroutine  pointTriangleDistance

    real(8) function dot3(x,y)

        implicit none

        real(8) :: x(3), y(3)

        dot3 = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)

    end function



    subroutine dist_rup_seg (Rrup, Rjb, Rx, coor, Ztor, strike, dip, rup_wid)

        real(8) :: Rrup, Rjb, Rx, Rrup1, Rrup2, Rjb1, Rjb2
        real(8) :: coor(2,2), Ztor, strike, dip, rup_wid
        real(8) :: Y1, X1, Y2, X2, Z1, Z2, Z3, Z4
        real(8) :: d_D, Y3, Y4, X3, X4
        real(8) :: site(3), P1(3), P2(3), P3(3), P4(3)
        real(8) :: P1_(3), P2_(3), P3_(3), P4_(3)

        Y1 = coor(1,1); X1 = coor(1,2)
        Y2 = coor(2,1); X2 = coor(2,2)
        Z1 = -Ztor;     Z2 = -Ztor

        Z3 = Z1 - rup_wid * sin(dip * DEG2RAD)
        Z4 = Z3

        d_D = rup_wid * cos(dip * DEG2RAD)

        Y3 = Y2 + d_D * cos(strike * DEG2RAD + PI / 2.0d0)
        Y4 = Y1 + d_D * cos(strike * DEG2RAD + PI / 2.0d0)

        X3 = X2 + d_D * sin(strike * DEG2RAD + PI / 2.0d0)
        X4 = X1 + d_D * sin(strike * DEG2RAD + PI / 2.0d0)

        site = [0.0d0, 0.0d0, 0.0d0]

        P1 = [X1,Y1,Z1]; P1_ = [X1, Y1, 0.0d0]
        P2 = [X2,Y2,Z2]; P2_ = [X2, Y2, 0.0d0]
        P3 = [X3,Y3,Z3]; P3_ = [X3, Y3, 0.0d0]
        P4 = [X4,Y4,Z4]; P4_ = [X4, Y4, 0.0d0]

        call pointTriangleDistance(P1, P2, P3, site, Rrup1)
        call pointTriangleDistance(P3, P4, P1, site, Rrup2)
        Rrup = min(Rrup1, Rrup2)

        if (abs(dip - 90) < 0.0001) then
            call pointLineSegDistance([X1, Y1], [X2, Y2], [0.0d0,0.0d0], Rjb)
        else
            call pointTriangleDistance(P1_, P2_, P3_, site, Rjb1)
            call pointTriangleDistance(P3_, P4_, P1_, site, Rjb2)
            Rjb = min(Rjb1, Rjb2)
        end if
        Rx = cal_Rx(coor)

    end subroutine

    function cal_Rx(coor)

        ! Arguments declarations
        real(8), intent(in) :: coor(2,2)
        ! Variable declarations
        real(8) :: az
        real(8) :: cal_Rx
        real(8) :: x1
        real(8) :: x2
        real(8) :: y1
        real(8) :: y2
        real(8) :: length

        ! the coor is in km
        ! coor = [Y1, X1; Y2, X2]
        ! site coor is (0,0) here

        y1 = coor(1,1)
        x1 = coor(1,2)
        y2 = coor(2,1)
        x2 = coor(2,2)

        call delaz2_km(y1,x1,y2,x2, length, az)
        cal_Rx = -cos(-az)*x1 - sin(-az)*y1

    end function cal_Rx

    subroutine dist_rup_set(Rrup, Rjb, Rx, coor, Ztor, strike, dip, rup_wid)
        real(8) :: Rrup, Rjb, Rx, Ztor, strike, dip, rup_wid
        real(8), allocatable :: coor(:,:)
        real(8) :: coor1(2,2)
        real(8), allocatable :: Rrup1(:), Rjb1(:), Rx1(:)
        integer :: i_rup, n_rup, indx

        n_rup = size(coor) / 2 - 1

        allocate(Rrup1(n_rup))
        allocate(Rjb1(n_rup))
        allocate(Rx1(n_rup))

        Rrup1 = 0.0d0
        Rjb1  = 0.0d0
        Rx1   = 0.0d0


        do i_rup = 1, n_rup

            coor1 = coor(i_rup:(i_rup+1),:)
            call dist_rup_seg (Rrup1(i_rup), Rjb1(i_rup), Rx1(i_rup), &
                coor1, Ztor, strike, dip, rup_wid)

        end do

        Rrup = minval(Rrup1)
        Rjb  = minval(Rjb1)
        indx = minloc(abs(Rx1), dim = 1)
        Rx = Rx1(indx)

    end subroutine

    subroutine interp_coeff (x1,x2,y1,y2,x,y,iflag)

        ! this is copied from Abrahamson haz code

        ! This subroutine will perform the Log-linear interpolation
        ! of the given input values. This routine is used to interpolate
        ! the regression cofficients of the attenuation models for
        ! spectral periods other than those defined in the model.

        integer :: iflag
        real(8) ::  x1, x2, y1, y2, x, y

        ! Check to see if the interpolation period is at an end point.
        ! Return the 'iflag' for output purposes with
        !             iflag = 0  No interpolation
        !                   = 1  Interpolation need.

        if (x .eq. x1 .or. x .eq. x2) then
            iflag = 0
        else
            iflag = 1
        end if

        ! Set the PGA period to 100 Hz (i.e., 0.01 Sec).
        if (x1 .eq. 0.0) then
            x1 = 0.01
        end if

        ! Take the Log of the Period values.
        x1 = log(x1)
        x2 = log(x2)
        x  = log(x)
        ! Perform the log-linear interpolation.
        y = y1 + (y2-y1)*((x-x1)/(x2-x1))

        ! Convert the Log Periods back to period.
        x1 = exp(x1)
        x2 = exp(x2)
        x  = exp(x)

        return
    end

    subroutine prob_exceed (p_exceed, m_eps, m_aleatory_distribution, trunclevel)

        real(8) :: p_exceed, m_eps, trunclevel
        integer :: m_aleatory_distribution

        if (m_aleatory_distribution .eq. NORMAL) then
            p_exceed = 1.0 - normcdf (m_eps)

        else if (m_aleatory_distribution .eq. TRUNC_NORMAL) then
            call truncnormcdf( m_eps,-trunclevel, trunclevel, p_exceed)
            p_exceed = 1.0 - p_exceed

        else if (m_aleatory_distribution .eq. HEAVISIDE) then
            p_exceed = 1.0 - deltacdf(m_eps)
        else
            stop 'error in subroutine prob_exceed '
        end if

    end subroutine


end module utils
