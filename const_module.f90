module const_module

    implicit none

    real(8), parameter :: PI = 3.14159265358979
    real(8), parameter :: DEG2RAD = 0.0174532925199433
    real(8), parameter :: RAD2DEG = 57.2957795130823
    real(8), parameter :: SQRT2_INV = 0.707106781186547
    real(8), parameter :: EARTH_R = 6371.0 ! in km

    integer, parameter :: SS = 1, RV = 2, NM = 3, NA = 4 ! style of faulting
    integer, parameter :: WC94 = 1, PEER = 2, CEUS = 3, POINT = 4 ! magnitude area scaling
    integer, parameter :: EXPONENTIAL = 1, CHARACTERISTIC = 2, DELTA = 3 ! recurrence model
    integer, parameter :: DEG = 1, KM = 2 ! unit
    integer, parameter :: SADIGH97 = 1, CY14 = 2 ! gmpe
    integer, parameter :: UNIFORM = 1, TRIANGULAR = 2 ! depth distribution
    integer, parameter :: NORMAL = 1, TRUNC_NORMAL = 2, HEAVISIDE = 3 ! aleatory distribution

end module const_module
