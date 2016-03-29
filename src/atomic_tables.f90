module atomic_tables

use H5LT
use HDF5
use hdf5_extra

IMPLICIT NONE

interface bt_maxwellian
    module procedure bt_maxwellian_n, bt_maxwellian_n_m
    module procedure bt_maxwellian_q_n, bt_maxwellian_q_n_m
end interface

integer, parameter, private   :: Int32   = 4 !bytes = 32 bits (-2,147,483,648 to 2,147,483,647)
integer, parameter, private   :: Int64   = 8 !bytes = 64 bits (-9,223,372,036,854,775,808 to 9,223,372,036,854,775,807)
integer, parameter, private   :: Float32 = 4 !bytes = 32 bits (1.2E-38 to 3.4E+38) at 6 decimal places
integer, parameter, private   :: Float64 = 8 !bytes = 64 bits (2.3E-308 to 1.7E+308) at 15 decimal places

real(Float64), parameter :: PI = 3.14159265d0

real(Float64), parameter :: e_amu = 5.485799093287202d-4
real(Float64), parameter :: H1_amu = 1.00782504d0
real(Float64), parameter :: H2_amu = 2.0141017778d0
real(Float64), parameter :: B_amu = 10.81d0
real(Float64), parameter :: C_amu = 12.011d0

integer, parameter :: B_q = 5
integer, parameter :: C_q = 6

!Spontaneous de-exitation Table n=1:15,m=1:15 for Hydrogen
!EINSTEIN(initial,final) / EINSTEIN(higher,lower) / EINSTEIN(n,m)
!Wiese, W. L., Smith, M. W., & Glennon, B. M. (1966). ATOMIC TRANSITION PROBABILITIES. VOLUME 1.
!HYDROGEN THROUGH NEON. NATIONAL BUREAU OF STANDARDS WASHINGTON DC INST FOR BASIC STANDARDS.
real(Float64), dimension(15,15), parameter :: EINSTEIN = reshape([ &                                                  !(n,m)
0.d0,4.699d8,5.575d7,1.278d7,4.125d6,1.644d6,7.568d5,3.869d5,2.143d5,1.263d5,7.834d4,5.066d4,3.393d4,2.341d4,1.657d4,&!(:,1) 
0.d0,0.d0   ,4.410d7,8.419d6,2.530d6,9.732d5,4.389d5,2.215d5,1.216d5,7.122d4,4.397d4,2.834d4,1.893d4,1.303d4,9.210d3,&!(:,2)
0.d0,0.d0   ,0.d0   ,8.986d6,2.201d6,7.783d5,3.358d5,1.651d5,8.905d4,5.156d4,3.156d4,2.021d4,1.343d4,9.211d3,6.490d3,&!(:,3)
0.d0,0.d0   ,0.d0   ,0.d0   ,2.699d6,7.711d5,3.041d5,1.424d5,7.459d4,4.235d4,2.556d4,1.620d4,1.069d4,7.288d3,5.110d3,&!(:,4)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.025d6,3.253d5,1.388d5,6.908d4,3.800d4,2.246d4,1.402d4,9.148d3,6.185d3,4.308d3,&!(:,5)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,4.561d5,1.561d5,7.065d4,3.688d4,2.110d4,1.288d4,8.271d3,5.526d3,3.815d3,&!(:,6)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,2.272d5,8.237d4,3.905d4,2.117d4,1.250d4,7.845d3,5.156d3,3.516d3,&!(:,7)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.233d5,4.676d4,2.301d4,1.287d4,7.804d3,5.010d3,3.359d3,&!(:,8)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,7.141d4,2.812d4,1.427d4,8.192d3,5.080d3,3.325d3,&!(:,9)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,4.377d4,1.774d4,9.231d3,5.417d3,3.324d3,&!(:,10)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,2.799d4,1.163d4,6.186d3,3.699d3,&!(:,11)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.857d4,7.884d3,4.271d3,&!(:,12)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,1.271d4,5.496d3,&!(:,13)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,8.933d3,&!(:,14)
0.d0,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ,0.d0   ]&!(:,15)
, [15,15])

contains

!!Collisions with protons
!!Charge Exchange
!!H(+) + H(nl) -> H(m'l') + H(+)
!!Eq. 44 "Collision processes in low-temperature hydrogen plasma" Janev et al
function p_cx_1_janev(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: a = [3.2345d0,  2.3588d2,  2.3713d0, &
                                                   3.8371d-2, 3.8068d-6, 1.1832d-10 ]
    real(Float64), parameter :: n = 1.d0
    real(Float64) :: Ehat
  
    Ehat = Erel * n**2.0
  
    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4)) 
  
end function p_cx_1_janev

function p_cx_2_janev(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: a = [9.2750d-1, 6.5040d3, 2.0699d1, &
                                                   1.3405d-2, 3.0842d-6, 1.1832d-10 ]
    real(Float64), parameter :: n = 2.d0
    real(Float64) :: Ehat
  
    Ehat = Erel * n**2.0
  
    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4)) 
  
end function p_cx_2_janev

function p_cx_3_janev(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: a = [3.7271d-1, 2.7645d6,  1.4857d3, &
                                                   1.5720d-3, 3.0842d-6, 1.1832d-10 ]
  
   
    real(Float64), parameter :: n = 3.d0
    real(Float64) :: Ehat
  
    Ehat = Erel * n**2.0
  
    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4)) 
  
end function p_cx_3_janev
  
function p_cx_4_janev(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: a = [2.1336d-1, 1.0000d10, 1.3426d6, &
                                                   1.8184d-3, 3.0842d-6, 1.1832d-10 ]
   
    real(Float64), parameter :: n = 4.d0
    real(Float64) :: Ehat
  
    Ehat = Erel * n**2.0
    sigma = (1.d-16*a(1)*(n**4))*log(a(2)/Ehat + a(3)) / &
            (1.d0+a(4)*Ehat + a(5)*Ehat**(3.5) + a(6)*Ehat**(5.4)) 
  
end function p_cx_4_janev
  
function p_cx_janev(Erel,n) result(sigma)
    real(Float64), intent(in) :: Erel
    integer, intent(in)       :: n
    real(Float64)             :: sigma

    integer :: i
  
    i = min(n,4)
  
    select case (i)
        case (0) 
            stop
        case (1)
            sigma = p_cx_1_janev(Erel)
        case (2)
            sigma = p_cx_2_janev(Erel)
        case (3) 
            sigma = p_cx_3_janev(Erel)
        case DEFAULT 
            sigma = p_cx_4_janev(Erel)
    end select
  
end function p_cx_janev
  
function p_cx_1_1_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(7), parameter :: a = [-3.496092687d2, 4.724931484d2, &
                                                   -2.720493064d2, 8.158564625d1, &
                                                   -1.339790721d1, 1.138706949d0, &
                                                   -3.914774156d-2 ]
    real(Float64) :: e, ee, fac, l, p
  
    e = Erel*1.d3
    if(e.ge.1.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 1.0d3
        fac = Erel**(-0.2)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0
  
    sigma = fac*(10.d0**p)

end function p_cx_1_1_adas
  
function p_cx_1_2_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(9), parameter :: a = [-4.036239511d3, 6.941235312d3, &
                                                   -5.186974866d3, 2.194885201d3, &
                                                   -5.765960509d2, 9.653534186d1, &
                                                   -1.008066138d1, 6.010731909d-1,&
                                                   -1.567417031d-2 ]
    real(Float64) :: e, ee, fac, l, p
  
    e = Erel*1.d3
    if(e.ge.1.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 1.0d3
        fac = Erel**(0.5)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0
  
    sigma = fac*(10.d0**p)

end function p_cx_1_2_adas
  
function p_cx_1_3_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(10), parameter :: a = [7.037287586d4, -1.479161477d5, &
                                                    1.370120708d5, -7.343180122d4, &
                                                    2.509832081d4, -5.674317075d3, &
                                                    8.487767749d2, -8.102284612d1, &
                                                    4.480007503d0, -1.093512342d-1 ]
    real(Float64) :: e, ee, fac, l, p
  
    e = Erel*1.d3
    if(e.ge.2.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 2.0d3
        fac = (Erel**(1.4))/2.8
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0
  
    sigma = fac*(10.d0**p)

end function p_cx_1_3_adas
  
function p_cx_1_4_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(10), parameter :: a = [6.826447557d4, -1.431980004d5, &
                                                    1.323968679d5, -7.083995050d4, &
                                                    2.417608863d4, -5.458418789d3, &
                                                    8.154875237d2, -7.776012846d1, &
                                                    4.295431731d0, -1.047567211d-1 ] 
    real(Float64) :: e, ee, fac, l, p
  
    e = Erel*1.d3
    if(e.ge.2.d3) then
        ee = max(e,1.0)
        fac = 1.d0
    else
        ee = 2.0d3
        fac = (Erel**(2.0))/4.0
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0
  
    sigma = fac*(10.d0**p)

end function p_cx_1_4_adas
  
function p_cx_1(Erel,m_max) result(sigma)
    real(Float64), intent(in)       :: Erel
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    real(Float64) :: norm_fac
  
    sigma = 0.d0
    sigma(1) = p_cx_1_1_adas(Erel)
    sigma(2) = p_cx_1_2_adas(Erel)
    sigma(3) = p_cx_1_3_adas(Erel)
    sigma(4) = p_cx_1_4_adas(Erel)
  
    !Normalize to Janev
    norm_fac = p_cx_janev(Erel, 1)/sum(sigma)
    sigma = norm_fac*sigma
  
end function p_cx_1
  
function p_cx_2_2_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
   
    real(Float64), dimension(11), parameter :: a2s = [-1.896015167d6, 4.431727330d6, &
                                                      -4.627815357d6, 2.843068107d6, &
                                                      -1.137952956d6, 3.100801094d5, &
                                                      -5.825744660d4, 7.452319142d3, &
                                                      -6.212350647d2, 3.047712749d1, &
                                                      -6.682658463d-1 ]
  
    real(Float64), dimension(11), parameter :: a2p = [-1.614213508d5, 3.772469288d5, &
                                                      -3.924736424d5, 2.393127027d5, &
                                                      -9.470300966d4, 2.541276100d4, &
                                                      -4.682860453d3, 5.851219013d2, &
                                                      -4.744504549d1, 2.254460913d0, &
                                                      -4.767235839d-2 ]
    real(Float64), parameter :: n = 2.d0
  
    real(Float64) :: e, ee, fac, l, sigma2s, sigma2p
  
    e = Erel * 1.d3 * n**2.0
    if(Erel.le.1.5d2) then
        ee = max(e, 1.d3)
        fac = 1.d0
    else
        ee = 1.5e5 * n**2.d0
        fac = 2.d15 * ((e*1.d-3)**(-5.5))
    endif
  
    l = log10(ee)
  
    sigma2s = a2s(1) + a2s(2)*l + a2s(3)*l**2.0 + a2s(4)*l**3.0 + &
              a2s(5)*l**4.0     + a2s(6)*l**5.0 + a2s(7)*l**6.0 + &
              a2s(8)*l**7.0     + a2s(9)*l**8.0 + a2s(10)*l**9.0 + a2s(11)*l**10.0
    sigma2s = 10.d0**(sigma2s)
  
    sigma2p = a2p(1) + a2p(2)*l + a2p(3)*l**2.0 + a2p(4)*l**3.0 + &
              a2p(5)*l**4.0     + a2p(6)*l**5.0 + a2p(7)*l**6.0 + &
              a2p(8)*l**7.0     + a2p(9)*l**8.0 + a2p(10)*l**9.0 + a2p(11)*l**10.0
    sigma2p = 10.d0**(sigma2p)
  
    sigma = fac*(0.25*sigma2s + 0.75*sigma2p)

end function p_cx_2_2_adas
  
function p_cx_2_3_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
   
    real(Float64), dimension(11), parameter :: a2s = [-3.513030327d5, 9.281116596d5, &
                                                      -1.086843398d6, 7.437325055d5, &
                                                      -3.296609685d5, 9.897503768d4, &
                                                      -2.039707143d4, 2.850670244d3, &
                                                      -2.587092857d2, 1.377382945d1, &
                                                      -3.268306303d-1 ]
  
    real(Float64), dimension(11), parameter :: a2p = [-1.901264631d5, 5.124716103d5, &
                                                      -6.101921504d5, 4.234717934d5, &
                                                      -1.899866398d5, 5.764464326d4, &
                                                      -1.199087959d4, 1.689900512d3, &
                                                      -1.545334374d2, 8.285001228d0, &
                                                      -1.978656474d-1 ]
  
    real(Float64), parameter :: n = 2.d0
  
    real(Float64) :: ee, l, sigma2s, sigma2p
  
    ee = max(Erel * 1.d3 * n**2.d0, 1.d3)
  
    l = log10(ee)
  
    sigma2s = a2s(1) + a2s(2)*l + a2s(3)*l**2.0 + a2s(4)*l**3.0 + &
              a2s(5)*l**4.0     + a2s(6)*l**5.0 + a2s(7)*l**6.0 + &
              a2s(8)*l**7.0     + a2s(9)*l**8.0 + a2s(10)*l**9.0 + a2s(11)*l**10.0
    sigma2s = 10.d0**(sigma2s)
  
    sigma2p = a2p(1) + a2p(2)*l + a2p(3)*l**2.0 + a2p(4)*l**3.0 + &
              a2p(5)*l**4.0     + a2p(6)*l**5.0 + a2p(7)*l**6.0 + &
              a2p(8)*l**7.0     + a2p(9)*l**8.0 + a2p(10)*l**9.0 + a2p(11)*l**10.0
    sigma2p = 10.d0**(sigma2p)
  
    sigma = (0.25*sigma2s + 0.75*sigma2p)

end function p_cx_2_3_adas
  
subroutine m_spread(n, m_max, sigma_tot, sigma)
    integer, intent(in)                            :: n
    integer, intent(in)                            :: m_max
    real(Float64), intent(in)                      :: sigma_tot
    real(Float64), dimension(m_max), intent(inout) :: sigma
  
    real(Float64) :: En, Em
    real(Float64) :: norm_fac
    real(Float64), dimension(m_max) :: sigma_m
    integer :: m
  
    sigma_m = 0.d0
    En = 13.6/(real(n)**2.0)
    do m=1,m_max
        Em = 13.6/(real(m)**2.0)
        if(sigma(m).eq.0.d0) then 
            sigma_m(m) = (sigma_tot/sqrt(2.0*PI))*exp(-0.5*(En-Em)**2.0)
        endif
    enddo
    
    norm_fac = sigma_tot/sum(sigma_m)
    do m=1,m_max
        if(sigma(m).eq.0.d0) sigma(m) = sigma_m(m)*norm_fac
        if(sigma(m).ne.sigma(m)) sigma(m) = 0.d0
    enddo

end subroutine m_spread
  
function p_cx_2(Erel,m_max) result(sigma)
    real(Float64), intent(in)       :: Erel
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    real(Float64), parameter :: n2 = 4.d0
    integer :: i
    real(Float64) :: En, Em, sigma_n, norm_fac
  
    sigma = 0.d0
    sigma(1) = p_cx_1_2_adas(Erel*n2)/n2
  
    sigma(2) = p_cx_2_2_adas(Erel)
    sigma(3) = p_cx_2_3_adas(Erel)
  
    sigma_n = max(p_cx_janev(Erel, 2) - sum(sigma), 0.d0)
  
    call m_spread(2,m_max,sigma_n,sigma)
  
    norm_fac = p_cx_janev(Erel, 2)/sum(sigma)
    sigma = sigma*norm_fac
  
end function p_cx_2
  
function p_cx_3_2_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: a = [-1.149224555d6, 2.750368877d6, &
                                                    -2.942222842d6, 1.852584954d6, &
                                                    -7.603284323d5, 2.125284465d5, &
                                                    -4.097580431d4, 5.380901722d3, &
                                                    -4.606297192d2, 2.321345254d1, &
                                                    -5.230186707d-1 ]
  
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p
  
    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 1.d16 * (Erel*n**2.0)**(-5.5)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.d0
  
    sigma = fac*(10.d0**p)

end function p_cx_3_2_adas
  
function p_cx_3_3_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(10), parameter :: a = [-4.302808608d4, 9.499298161d4, &
                                                    -9.264698488d4, 5.236947172d4, &
                                                    -1.890479538d4, 4.519068626d3, &
                                                    -7.152485009d2, 7.227063167d1, &
                                                    -4.230036444d0, 1.092702525d-1 ]
  
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p
  
    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2
        fac = 0.85d16 *(Erel*n**2.0)**(-5.5)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0
  
    sigma = fac*(10.d0**p)
  
end function p_cx_3_3_adas
  
function p_cx_3_4_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(9), parameter :: a = [ 1.705303425d4,-3.316878090d4, &
                                                    2.792556433d4,-1.330264490d4, &
                                                    3.921666688d3,-7.327555138d2, &
                                                    8.476342861d1,-5.551987930d0, &
                                                    1.577120745d-1 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p
  
    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 0.82d16 *(Erel*n**2.0)**(-5.5)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0
  
    sigma = fac*(10.d0**p)
    
end function p_cx_3_4_adas
  
function p_cx_3_5_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: a = [-2.786268232d2, 4.269683825d4, &
                                                    -8.973561028d4, 8.365732310d4, &
                                                    -4.524587937d4, 1.563630402d4, &
                                                    -3.580391824d3, 5.432527332d2, &
                                                    -5.267599631d1, 2.962329657d0, &
                                                    -7.362649692d-2 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p
  
    ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
    fac = 1.d0
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.0
  
    sigma = fac*(10.d0**p)
  
end function p_cx_3_5_adas
  
function p_cx_3_6inf_adas(Erel) result(sigma)
    real(Float64), intent(in) :: Erel
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: a = [ 7.146969470d5,-1.665413326d6, &
                                                     1.735840441d6,-1.065792786d6, &
                                                     4.269334710d5,-1.165954977d5, &
                                                     2.198700496d4,-2.827160468d3, &
                                                     2.372409350d2,-1.173264972d1, &
                                                     2.596865877d-1 ]
    real(Float64), parameter :: n = 3.0
    real(Float64) :: ee, fac, l, p
  
    if(Erel.lt.90.0) then
        ee = max(Erel * 1.d3 * n**2.0, 1.d3) !keV to eV
        fac = 1.d0
    else
        ee = 90.0 * 1.d3 * n**2.0
        fac = 2.d20 *(Erel*n**2.0)**(-7.0)
    endif
  
    l = log10(ee)
  
    p = a(1) + a(2)*l + a(3)*l**2.0 + a(4)*l**3.0 + &
        a(5)*l**4.0   + a(6)*l**5.0 + a(7)*l**6.0 + &
        a(8)*l**7.0   + a(9)*l**8.0 + a(10)*l**9.0 + a(11)*l**10.0
  
    sigma = fac*(10.d0**p)
  
end function p_cx_3_6inf_adas
  
function p_cx_3(Erel,m_max) result(sigma)
    real(Float64), intent(in)       :: Erel
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    real(Float64), parameter :: n2 = 9.d0
    real(Float64) :: eb, En, Em, sigma_m6, norm_fac
    real(Float64), dimension(m_max) :: sigma1
  
    sigma = 0.d0
    sigma1 = 0.d0
    sigma1 = p_cx_1(Erel*n2,m_max)
    sigma(1) = p_cx_1_3_adas(Erel*n2)/n2
  
    sigma(2) = p_cx_3_2_adas(Erel)
    sigma(3) = p_cx_3_3_adas(Erel)
    sigma(4) = p_cx_3_4_adas(Erel)
  
    if(m_max.ge.5) then
        sigma(5) = p_cx_3_5_adas(Erel)
    endif
  
    if(m_max.ge.6) then
        sigma_m6 = p_cx_3_6inf_adas(Erel)
        call m_spread(3, m_max, sigma_m6, sigma)
    endif
  
    norm_fac = p_cx_janev(Erel, 3)/sum(sigma)
    sigma = sigma*norm_fac
  
end function p_cx_3
  
function p_cx_n(Erel, n, m_max) result(sigma)
    real(Float64), intent(in)       :: Erel
    integer, intent(in)             :: n
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    real(Float64), dimension(m_max) :: sigma2,sigma3
    real(Float64) :: sigma_n,e,norm_fac
  
    sigma = 0.d0
    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_cx_1(Erel,m_max)
            return
        case (2)
            sigma = p_cx_2(Erel,m_max)
            return
        case (3)
            sigma = p_cx_3(Erel,m_max)
            return
        case (4)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(1) = p_cx_1_4_adas(e/(1.0**2.0))*(1.d0/n)**2.0
            sigma(2) = sigma2(4)*(2.d0/n)**2.0
            sigma(3) = p_cx_3_4_adas(e/(3.0**2.0))*(3.d0/n)**2.0
        case (5)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(2) = sigma2(5)*(2.d0/n)**2.0
            sigma(3) = p_cx_3_5_adas(e/(3.0**2.0))*(3.d0/n)**2.0
        case (6)
            e = Erel*n**2.0
            sigma2 = p_cx_2(e/(2.0**2.0),m_max)
            sigma(2) = sigma2(6)*(2.d0/n)**2.0
            sigma3 = p_cx_3(e/(3.0**2.0),m_max)*(3.d0/n)**2.0
            sigma(3) = sigma3(6)
        case DEFAULT
    end select 
  
    sigma_n = max(p_cx_janev(Erel,n) - sum(sigma),0.0)
    call m_spread(n, m_max, sigma_n, sigma)
  
    norm_fac = p_cx_janev(Erel, n)/sum(sigma)
    sigma = norm_fac*sigma
  
end function p_cx_n
  
function p_cx_n_m(Erel, n, m) result(sigma)
    real(Float64), intent(in) :: Erel
    integer, intent(in)       :: n
    integer, intent(in)       :: m
    real(Float64)             :: sigma
  
    integer :: m_max = 12
    real(Float64), dimension(12) :: sigma_m
  
    sigma_m = p_cx_n(Erel, n, m_max)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif
  
end function p_cx_n_m
  
function p_cx(Erel, n_max, m_max) result(sigma)
    real(Float64), intent(in)             :: Erel
    integer, intent(in)                   :: n_max
    integer, intent(in)                   :: m_max
    real(Float64), dimension(m_max,n_max) :: sigma
  
    real(Float64), dimension(12,12) :: sigma_full
  
    integer :: n, m
  
    do n=1,12
        sigma_full(n,:) = p_cx_n(Erel, n, 12)
    enddo
  
    sigma = sigma_full(1:n_max,1:m_max)

end function p_cx
  
!Proton impact ionization
function p_ioniz_1_omullane(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [2.0160d-3, 3.7154d0,  &
                                                   3.9890d-2, 3.1413d-1, &
                                                   2.1254d0,  6.3990d3,  &
                                                   6.1897d1,  9.2731d3 ]
  
    real(Float64), parameter :: n2 = 1.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3
  
    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)

end function p_ioniz_1_omullane
  
function p_ioniz_2_omullane(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [3.9330d-3, 1.8188d0,  &
                                                   1.8870d-2, 6.7489d-3, &
                                                   1.3768d0, 6.8852d2,   &
                                                   9.6435d1, 5.6515d23 ]
  
    real(Float64), parameter :: n2 = 4.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3
  
    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)
  
end function p_ioniz_2_omullane
  
function p_ioniz_3_omullane(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [1.1076d-2, 1.6197d0,  &
                                                   6.7154d-3, 5.1188d-3, &
                                                   1.8549d0,  2.3696d2,  &
                                                   7.8286d1,  1.0926d23 ]
  
    real(Float64), parameter :: n2 = 9.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3
  
    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)
  
end function p_ioniz_3_omullane
  
function p_ioniz_4_omullane(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [1.1033d-2, 1.6281d0,  &
                                                   5.5955d-3, 7.2023d-3, &
                                                   1.7358d0,  2.2755d2,  &
                                                   8.6339d1,  3.9151d29 ]
  
    real(Float64), parameter :: n2 = 16.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3
  
    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)
  
end function p_ioniz_4_omullane
  
function p_ioniz_5_omullane(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [1.1297d-2, 1.8685d0,  &
                                                   1.5038d-2, 1.1195d-1, &
                                                   1.0538d0,  8.6096d2,  &
                                                   8.9939d1,  1.9249d4 ]
  
    real(Float64), parameter :: n2 = 25.d0
    real(Float64) :: Ehat
    real(Float64) :: p1, p2, p3
  
    Ehat = eb*n2
    p1 = b(1)*(n2)**2.0
    p2 = Ehat**b(2) * exp(-b(3)*Ehat) / (1.d0 + b(4)*Ehat**b(5))
    p3 = (b(6)* exp(-b(7)/Ehat) *log(1.d0  +b(8)*Ehat) ) /Ehat
    sigma = 1.0d-16 * p1 * (p2 + p3)
  
end function p_ioniz_5_omullane
  
function p_ioniz_n(eb,n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_ioniz_1_omullane(eb)
        case (2)
            sigma = p_ioniz_2_omullane(eb)
        case (3)
            sigma = p_ioniz_3_omullane(eb)
        case (4)
            sigma = p_ioniz_4_omullane(eb)
        case DEFAULT
            sigma = p_ioniz_5_omullane(eb)*(real(n)/5.0)**4.0
    end select

end function p_ioniz_n
  
function p_ioniz(eb,n_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: n_max
    real(Float64), dimension(n_max) :: sigma
  
    integer :: i
  
    do i=1,n_max
        sigma(i) = p_ioniz_n(eb,i)
    enddo

end function p_ioniz
  
  !! Proton impact excitation
function p_excit_1_2_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(10), parameter :: a = [34.433d0, 8.5476d0,  &
                                                    7.8501d0, -9.2217d0, &
                                                    1.8020d-2, 1.6931d0, &
                                                    1.9422d-3, 2.9068d0, &
                                                    44.507d0, 0.56870d0 ]
      
    sigma = 1.d-16 * a(1) * ( a(2)*exp(-a(3)*eb)/(eb**a(4)) + & 
                     a(5)*exp(-a(6)/eb)/(1.+a(7)*eb**a(8))  + &
                     exp(-a(9)/eb)*log(1.+a(10)*eb)/eb )

end function p_excit_1_2_janev
  
function p_excit_1_3_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [ 6.1950d0, 5.5162d-3,  &
                                                    0.29114d0, -4.5264d0, &
                                                    6.0311d0, -2.0679d0,  &
                                                    35.773d0, 0.54818d0 ]
    
    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)
  
end function p_excit_1_3_janev
  
function p_excit_1_4_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [2.0661d0, 5.1335d-4,  &
                                                   0.28953d0, -2.2849d0, &
                                                   0.11528d0, -4.8970d0, &
                                                   34.975d0, 0.91213d0 ]
  
    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)
  
end function p_excit_1_4_janev
  
function p_excit_1_5_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [1.2449d0, 3.0826d-4,   &
                                                   0.31063d0, -2.4161d0,  &
                                                   0.024664d0, -6.3726d0, &
                                                   32.291d0, 0.21176d0 ]
  
    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)
  
end function p_excit_1_5_janev
  
function p_excit_1_6_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: b = [0.63771d0, 3.2949d-4,  &
                                                   0.25757d0, -2.2950d0,  &
                                                   0.050796d0, -5.5986d0, &
                                                   37.174d0, 0.39265d0 ]
  
    sigma = 1.d-16 * b(1) * (b(2)*exp(-b(3)*eb)/ &
                     (eb**b(4)+b(5)*eb**b(6)) + &
                     exp(-b(7)/eb)*log(1.+b(8)*eb)/eb)
  
end function p_excit_1_6_janev
  
function p_excit_1_janev(eb, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    sigma = 0.d0
  
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = p_excit_1_2_janev(eb)
            case (3)
                sigma(3) = p_excit_1_3_janev(eb)
            case (4)
                sigma(4) = p_excit_1_4_janev(eb)
            case (5)
                sigma(5) = p_excit_1_5_janev(eb)
            case (6)
                sigma(6) = p_excit_1_6_janev(eb)
            case DEFAULT
                sigma(m) = p_excit_1_6_janev(eb)*(6.0/real(m))**3.0
        end select
    enddo
  
end function p_excit_1_janev
  
function p_excit_2_3_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [394.51d0, 0.013597d0, &
                                                   0.16565d0, -0.8949d0, &
                                                   21.606d0,  0.62426d0  ]
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_2_3_janev
  
function p_excit_2_4_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [50.744d0, 0.014398d0, &
                                                   0.31584d0, -1.4799d0, &
                                                   19.416d0,   4.0262d0  ]
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_2_4_janev
  
function p_excit_2_5_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [18.264d0, 0.013701d0, &
                                                   0.31711d0, -1.4775d0, &
                                                   18.973d0,   2.9056d0  ]
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_2_5_janev
  
function p_excit_2_6_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 4.61d-1
  
    sigma = A*p_excit_2_5_janev(eb)
  
end function p_excit_2_6_janev
  
function p_excit_2_7_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 2.475d-1
  
    sigma = A*p_excit_2_5_janev(eb)
  
end function p_excit_2_7_janev
  
function p_excit_2_8_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 1.465d-1
  
    sigma = A*p_excit_2_5_janev(eb)
  
end function p_excit_2_8_janev
  
function p_excit_2_9_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 9.2d-2
  
    sigma = A*p_excit_2_5_janev(eb)
  
end function p_excit_2_9_janev
  
function p_excit_2_10_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 6.05d-2
  
    sigma = A*p_excit_2_5_janev(eb)
  
end function p_excit_2_10_janev
  
function p_excit_2_janev(eb, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    sigma = 0.d0
  
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = p_excit_2_3_janev(eb)
            case (4)
                sigma(4) = p_excit_2_4_janev(eb)
            case (5)
                sigma(5) = p_excit_2_5_janev(eb)
            case (6)
                sigma(6) = p_excit_2_6_janev(eb)
            case (7)
                sigma(7) = p_excit_2_7_janev(eb)
            case (8)
                sigma(8) = p_excit_2_8_janev(eb)
            case (9)
                sigma(9) = p_excit_2_9_janev(eb)
            case (10)
                sigma(10) = p_excit_2_10_janev(eb)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/real(m))**3.0
        end select
    enddo
  
end function p_excit_2_janev
  
function p_excit_3_4_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [1247.5d0,  0.068781d0, &
                                                   0.521176d0, -1.2722d0, &
                                                   11.319d0,    2.6235d0  ]  
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_3_4_janev
  
function p_excit_3_5_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [190.59d0, 0.073307d0, &
                                                   0.54177d0, -1.2894d0, &  
                                                   11.096d0,   2.9098d0  ]
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_3_5_janev
  
function p_excit_3_6_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: c = [63.494d0, 0.077953d0, &
                                                   0.53461d0, -1.2881d0, &
                                                   11.507d0,   4.3417d0  ]
  
    sigma = 1.d-16 * c(1)*(c(2)*exp(-c(3)*eb)/(eb**c(4)) + &
                     exp(-c(5)/eb)*log(1.+c(6)*eb)/eb)
  
end function p_excit_3_6_janev
  
function p_excit_3_7_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 4.67d-1
  
    sigma = A*p_excit_3_6_janev(eb)
  
end function p_excit_3_7_janev
  
function p_excit_3_8_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 2.545d-1
  
    sigma = A*p_excit_3_6_janev(eb)
  
end function p_excit_3_8_janev
  
function p_excit_3_9_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 1.54d-1
  
    sigma = A*p_excit_3_6_janev(eb)
  
end function p_excit_3_9_janev
  
function p_excit_3_10_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 1.0d-1
  
    sigma = A*p_excit_3_6_janev(eb)
  
end function p_excit_3_10_janev
  
function p_excit_3_janev(eb, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    sigma = 0.d0
  
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = 0.d0
            case (4)
                sigma(4) = p_excit_3_4_janev(eb)
            case (5)
                sigma(5) = p_excit_3_5_janev(eb)
            case (6)
                sigma(6) = p_excit_3_6_janev(eb)
            case (7)
                sigma(7) = p_excit_3_7_janev(eb)
            case (8)
                sigma(8) = p_excit_3_8_janev(eb)
            case (9)
                sigma(9) = p_excit_3_9_janev(eb)
            case (10)
                sigma(10) = p_excit_3_10_janev(eb)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/real(m))**3.0
        end select
    enddo
  
end function p_excit_3_janev
  
function p_excit_n(eb, n, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: n
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    real(Float64) :: nf, mf, Etil, s, D, A, G, L, F 
    real(Float64) :: y, zpl, zmi, C2pl, C2mi, H
    sigma = 0.d0
  
    select case (n)
        case (0)
            stop
        case (1)
            sigma = p_excit_1_janev(eb,m_max)
        case (2)
            sigma = p_excit_2_janev(eb,m_max)
        case (3)
            sigma = p_excit_3_janev(eb,m_max)
        case DEFAULT
            nf = real(n)
            m_loop: do m=1,m_max
                if(n.ge.m) then
                    sigma(m) = 0.d0
                    cycle m_loop
                endif
                mf = real(m)
                Etil = Eb/25.0
                s = (mf-nf)

                D = exp(-1.0/(nf*mf*Etil**2.0))
                A = 8.0/(3.0*s)*(mf/(s*nf))**3*(0.184-0.04/s**(2.0/3.0)) * &
                    (1.0 - 0.2*s/(nf*mf))**(1.0 + 2.0*s)
                G = 0.5*( Etil*nf**2.0/(mf - 1.0/mf) )**3.
                L = log(1.0 + 0.53*Etil**2.0 * nf*(mf - 2.0/mf)/(1.0 + 0.4*Etil))
                F = ( 1.0 - 0.3*s*D/(nf*mf) )**(1.0 + 2.0*s)
         
                y = 1.0/( 1.0 - D*log(18*s)/(4.0*s) )
                zpl = 2.0/(Etil * nf**2 * ( (2.0 - (nf/mf)**2)**0.5 + 1.0))
                zmi = 2.0/(Etil * nf**2 * ( (2.0 - (nf/mf)**2)**0.5 - 1.0))
                C2pl = zpl**2 * log(1.0 + 2.0*zpl/3.0)/(2.0*y + 3.0*zpl/2.0)
                C2mi = zmi**2 * log(1.0 + 2.0*zmi/3.0)/(2.0*y + 3.0*zmi/2.0)

                H = C2mi - C2pl

                sigma(m) = ((8.8d-17*n**4)/Etil)*(A*L*D + F*G*H)
            enddo m_loop
    end select
  
end function p_excit_n
  
function p_excit_n_m(eb, n, m) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: n
    integer, intent(in)       :: m
    real(Float64)             :: sigma
  
    real(Float64), dimension(12) :: sigma_m
  
    sigma_m = p_excit_n(eb, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif
  
end function p_excit_n_m
  
function p_excit(eb, n_max, m_max) result(sigma)
    real(Float64), intent(in)             :: eb
    integer, intent(in)                   :: m_max
    integer, intent(in)                   :: n_max
    real(Float64), dimension(m_max,n_max) :: sigma
  
    real(Float64), dimension(12,12) :: sigma_full
  
    integer :: n, m
  
    do n=1,12
        sigma_full(n,:) = p_excit_n(eb, n, 12)
    enddo
  
    sigma = sigma_full(1:n_max,1:m_max)

end function p_excit
  
function e_ioniz_1_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    integer, parameter :: n = 1
    real(Float64), dimension(6), parameter :: A = [ 0.18450d0, -0.032226d0, &
                                                   -0.034539d0, 1.4003d0,   &
                                                   -2.8115d0, 2.2986d0 ]
    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s
  
    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV
  
    x = (1.0 - Edn2/e)
    s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
    sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    sigma = max(sigma,0.d0)
  
end function e_ioniz_1_janev
  
function e_ioniz_2_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    integer, parameter :: n = 2
    real(Float64), dimension(6), parameter :: A = [ 0.14784d0, 0.0080871d0, &
                                                   -0.062270d0, 1.9414d0,   & 
                                                   -2.1980d0, 0.95894d0 ]
    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s
  
    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV
  
    x = (1.0 - Edn2/e)
    s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
    sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    sigma = max(sigma, 0.d0)

end function e_ioniz_2_janev
  
function e_ioniz_3_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    integer, parameter :: n = 3
    real(Float64), dimension(6), parameter :: A = [0.058463d0, -0.051272d0, &
                                                   0.85310d0, -0.57014d0,   &
                                                   0.76684d0, 0.00d0 ]
  
    real(Float64) :: Edn2
    real(Float64) :: e, x
    real(Float64) :: s
  
    Edn2 = 13.6/real(n)**2
    e = eb * 1.d3 !keV to eV
  
    if(e.ge.1.5) then
        x = (1.0 - Edn2/e)
        s = A(2)*x + A(3)*(x**2.0) + A(4)*(x**3.0) + A(5)*(x**4.0) + A(6)*(x**5.0)
        sigma = ((1.d-13)/(Edn2*e))*(A(1)*log(e/Edn2) + s)
    else
        sigma = 0.d0
    endif
  
end function e_ioniz_3_janev
  
function e_ioniz_n(eb, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: n
    real(Float64) :: sigma
  
    real(Float64) :: rn, xn, Edn2
    real(Float64) :: g0, g1, g2, An, b, Bn

    select case (n)
        case (0)
            stop
        case (1)
            sigma = e_ioniz_1_janev(eb)
        case (2)
            sigma = e_ioniz_2_janev(eb)
        case (3)
            sigma = e_ioniz_3_janev(eb)
        case DEFAULT
            rn = 1.94/n**1.57
            Edn2 = 13.6/n**2.0
            xn = (eb*1.d3)/Edn2

            g0 = 0.9935 + 0.2328/n - 0.1296/n**2.0
            g1 = -(1.0/n)*(0.6282 - 0.5598/n + 0.5299/n**2.0)
            g2 =  (1.0/n**2.0)*(0.3887 - 1.181/n + 1.47/n**2.0)

            An = 32.0*n/(3.0*sqrt(3.0)*PI)*(g0/3.0 + g1/4.0 + g2/5.0)
            b  = (1.0/n)*(4.0 - 18.63/n + 36.24/n**2.0 - 28.09/n**3.0)
            Bn = (2.0/3.0)*(n**2.0)*(5.0 + b)

            if(xn.gt.1) then
                sigma = 1.76*n**2/xn*(1.0 - exp(-rn*xn)) * &
                        (An*log(xn) + (Bn - An*log(2.0*n**2)) * &
                        (1.0 - 1.0/xn)**2)*1.e-16
            else
                sigma = 0.d0
            endif
    end select
  
end function e_ioniz_n
  
function e_ioniz(eb, n_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: n_max
    real(Float64), dimension(n_max) :: sigma
  
    integer :: i
  
    do i=1,n_max
        sigma(i) = e_ioniz_n(eb,i)
    enddo

end function e_ioniz
  
function e_excit_1_2_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: sigma0=5.984d0
    real(Float64), parameter :: deltaE=10.2d0
    real(Float64), parameter :: a=0.228d0
    real(Float64), parameter :: b=0.1865d0
    real(Float64), parameter :: c=0.5025d0
    real(Float64), dimension(6), parameter :: An = [ 4.4979d0, 1.4182d0, &
                                                    -20.877d0, 49.735d0, &
                                                    -46.249d0, 17.442d0 ]
  
    real(Float64) :: ecoll, x, s
  
    ecoll = eb*1.d3
    x = (ecoll)/deltaE
  
    if((ecoll.gt.10.2).and.(ecoll.le.11.56)) then
        sigma = 1.d-16 * (a + b*(ecoll - deltaE))
        return
    endif
  
    if((ecoll.ge.11.56).and.(ecoll.le.12.23)) then
        sigma = 1.d-16 * c
        return
    endif
  
    if(ecoll.ge.12.23) then
        s = An(2) + An(3)/x + An(4)/x**2.0 + An(5)/x**3.0 + An(6)/x**4.0
        sigma = 1.d-16 * sigma0/(deltaE*x) * (An(1)*log(x) + s)
        return
    endif
  
    if(x.le.1.0) then
        sigma = 0.0
        return
    endif

end function e_excit_1_2_janev
  
function e_excit_1_3_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 12.09d0
    real(Float64), parameter :: alpha = 0.38277d0
    real(Float64), dimension(5), parameter :: A = [ 0.75448d0, 0.42956d0, &
                                                   -0.58288d0, 1.0693d0,  &  
                                                    0.d0 ]     
    real(Float64) :: ecoll, x, s
  
    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0
    
    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)
  
    if(x.le.1.0) sigma = 0.d0
  
end function e_excit_1_3_janev
  
function e_excit_1_4_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 12.75d0
    real(Float64), parameter :: alpha = 0.41844d0
    real(Float64), dimension(5), parameter :: A = [ 0.24300d0, 0.24846d0, &
                                                    0.19701d0, 0.d0,      &
                                                    0.d0 ]
  
    real(Float64) :: ecoll, x, s
  
    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0
    
    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)
  
    if(x.le.1.0) sigma = 0.d0
  
end function e_excit_1_4_janev
  
function e_excit_1_5_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 13.06d0
    real(Float64), parameter :: alpha = 0.45929d0
    real(Float64), dimension(5), parameter :: A = [ 0.11508d0, 0.13092d0, &
                                                    0.23581d0, 0.d0,      &
                                                    0.d0 ]   
  
    real(Float64) :: ecoll, x, s
  
    ecoll = eb*1.d3
    x=ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0
    
    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)
  
    if(x.le.1.0) sigma = 0.d0
  
end function e_excit_1_5_janev
  
function e_excit_f(n, m) result(fnm)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(Float64)       :: fnm
  
    real(Float64), dimension(3) :: g
    real(Float64) :: x, nf, mf, gs
  
    nf = real(n)
    mf = real(m)
    x = 1.0 - (nf/mf)**2.0
  
    select case (n)
        case (1)
            g = [1.133,-0.4059,0.0714]
        case (2)
            g = [1.0785,-0.2319,0.02947]
        case DEFAULT
            g(1) = 0.9935 + 0.2328/nf - 0.1296/nf**2
            g(2) =-1.0/nf * (0.6282 - 0.5598/nf + 0.5299/nf**2)
            g(3) = 1.0/nf**2.0 * (0.3887 - 1.1810/nf + 1.4700/nf**2)
    end select 

    gs = g(1) + g(2)/x + g(3)/x**2
    fnm = 32.0/(3.0*sqrt(3.0)*PI) * nf/mf**3 * 1/x**3 * gs

end function e_excit_f
  
function e_excit_1_janev(eb, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    real(Float64) :: x, y, A, B, deltaE
  
    do m=1,m_max 
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = e_excit_1_2_janev(eb)
            case (3)
                sigma(3) = e_excit_1_3_janev(eb)
            case (4)
                sigma(4) = e_excit_1_4_janev(eb)
            case (5)
                sigma(5) = e_excit_1_5_janev(eb)
            case DEFAULT
                y = 1.0 - (1.d0/m)**2.0
                deltaE = 13.6*y
                x = (eb*1.d3)/deltaE

                A = 2.0 * e_excit_f(1,m)/y
                B = 4.0/(m**3.0 * y)*(1.0 + 4.0/(3.0*y) - 0.603/y**2.0)

                sigma(m) = 1.76e-16/(y*x)*(1.0 - exp(-0.45*y*x))* &
                           (A*(log(x) + 1.0/(2.0*x)) + (B - A*log(2.0/y))* &
                           (1.0 - 1.0/x))

                if(x.le.1.0) sigma(m) = 0.d0
        end select
    enddo
  
end function e_excit_1_janev
  
function e_excit_2_3_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), parameter :: sigma0 = 5.984d0
    real(Float64), parameter :: deltaE = 1.8888d0
    real(Float64), parameter :: alpha = 1.3196d0
    real(Float64), dimension(5), parameter :: A = [ 38.906d0, 5.2373d0, 119.25d0, &
                                                   -595.39d0, 816.71d0]
  
    real(Float64) :: ecoll, x, s
  
    ecoll = eb*1.d3
    x = ecoll/deltaE
    s = A(2) + A(3)/x + A(4)/x**2.0 + A(5)/x**3.0
    
    sigma = 1.d-16 * sigma0/(deltaE*x) * (1.0 - 1.0/x)**alpha * (A(1)*log(x) + s)
  
    if(x.le.1.0) sigma = 0.d0
  
end function e_excit_2_3_janev
  
function e_excit_n(eb, n, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: n
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    real(Float64) :: nf, mf
    real(Float64) :: x, y, A, B, bn, r, deltaE
  
    nf = real(n)
  
    if(n.eq.1) then
        sigma = e_excit_1_janev(eb, m_max)
    else
        m_loop: do m=1,m_max
            mf = real(m)
            if(n.ge.m) then
                sigma(m) = 0.d0
                cycle m_loop
            endif
            if((n.eq.2).and.(m.eq.3)) then
                sigma(m) = e_excit_2_3_janev(eb)
            else
                deltaE=13.6*(1.0/nf**2 - 1.0/mf**2)
                x = (eb*1.d3)/deltaE
                y = 1.0 - (nf/mf)**2
                r = 1.94/nf**1.57

                A = 2.0 * nf**2 * e_excit_f(n,m)/y
                bn = 1.0/nf*(4.0 - 18.63/nf + 36.24/nf**2 - 28.09/nf**3)
                B = 4.0 * nf**4/(mf**3*y**2)*(1.0 + 4.0/(3.0*y) + bn/y**2.0)     

                sigma(m) = 1.76e-16*nf**2/(y*x)*(1.0 - exp(-r*y*x))* &
                           (A*(log(x) + 1.0/(2.0*x)) + (B - A*log(2.0*n**2.0/y))* &
                           (1.0 - 1.0/x))

                if(x.le.1.0) sigma(m) = 0.d0
            endif
        enddo m_loop
    endif
        
end function e_excit_n
  
function e_excit_n_m(eb, n, m) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: n
    integer, intent(in)       :: m
    real(Float64)             :: sigma
  
    real(Float64), dimension(12) :: sigma_m
  
    sigma_m = e_excit_n(eb, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif
  
end function e_excit_n_m
  
function e_excit(eb, n_max, m_max) result(sigma)
    real(Float64), intent(in)             :: eb
    integer, intent(in)                   :: n_max
    integer, intent(in)                   :: m_max
    real(Float64), dimension(n_max,m_max) :: sigma
  
    real(Float64), dimension(12,12) :: sigma_full
  
    integer :: n
  
    do n=1,12
        sigma_full(n,:) = e_excit_n(eb, n, 12)
    enddo
  
    sigma = sigma_full(1:n_max,1:m_max)
  
end function e_excit

!Impurities
!A[q]_cx_[n]_[source]
function B5_cx_1_adas(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(7), parameter :: A = [ 1.174052518d3, -1.793561728d3, &
                                                    1.117522436d3, -3.679435571d2, &
                                                    6.750816878d1, -6.542029074d0, &
                                                    2.614113716d-1 ]
    real(Float64) :: e, l, p
  
    e = max(eb*1.d3,10.0)
    l = log10(e)
  
    if(e.le.4.d5) then
        p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
            A(5)*l**4 + A(6)*l**5 + A(7)*l**6
        sigma = 10.d0**p
    else
        sigma = 0.d0
    endif
  
end function B5_cx_1_adas
  
function B5_cx_2_adas(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(10), parameter :: A = [6.603246818d1, -3.072575676d2, &
                                                    5.030801019d2, -4.585636345d2, &
                                                    2.568666393d2, -9.185150382d1, &
                                                    2.100012584d1, -2.964174788d0, &
                                                    2.346396110d-1, -7.943766873d-3]
    real(Float64) :: e, l, p
  
    e = max(eb*1.d3,10.0)
    l = log10(e)
  
    if(e.le.1.d5) then
        p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
            A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
            A(8)*l**7 + A(9)*l**8 + A(10)*l**9
        sigma = 10.d0**p
    else
        sigma = 0.d0
    endif
  
end function B5_cx_2_adas
  
function C6_cx_1_adas(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(7), parameter :: A = [2.007882674d2, -3.546893286d2, &
                                                   2.381542403d2, -8.355431742d1, &
                                                   1.617519888d1, -1.638152470d0, &
                                                   6.768953863d-2 ]
  
    real(Float64) :: e, l, p
  
    e = max(eb*1.d3,1.5d3)
    l = log10(e)
  
    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6
    sigma = 10.d0**p
  
end function C6_cx_1_adas
  
function C6_cx_2_adas(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: A = [9.151879441d5, -2.134573133d6, &
                                                    2.223792624d6, -1.362648703d6, &
                                                    5.438401343d5, -1.477110500d5, &
                                                    2.764972254d4, -3.522105245d3, &
                                                    2.921934171d2, -1.425552507d1, &
                                                    3.106007048d-1 ]
    real(Float64) :: e, l, p
  
    e = max(eb*1.d3,1.5d3)*2.0**2
    l = log10(e)
  
    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
        A(8)*l**7 + A(9)*l**8 + A(10)*l**9 + A(11)*l**10
    sigma = 10.d0**p
  
end function C6_cx_2_adas
  
function C6_cx_3_adas(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: A = [9.208877916d5, -2.147294379d6, &
                                                    2.236451628d6, -1.370042347d6, &
                                                    5.466461899d5, -1.484338816d5, &
                                                    2.777765778d4, -3.537459450d3, &
                                                    2.933884362d2, -1.430994136d1, &
                                                    3.117002878d-1 ] 
    real(Float64) :: e, l, p
  
    e = max(eb*1.d3,1.5d3)*3.0**2
    l = log10(e)
  
    p = A(1) + A(2)*l + A(3)*l**2 + A(4)*l**3 + &
        A(5)*l**4 + A(6)*l**5 + A(7)*l**6 +     &
        A(8)*l**7 + A(9)*l**8 + A(10)*l**9 + A(11)*l**10
    sigma = 10.d0**p
  
end function C6_cx_3_adas
  
function Aq_cx_n_adas(eb, q, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    sigma = 0.d0
    select case (q)
        case (5)
            if(n.eq.1) sigma = B5_cx_1_adas(eb)
            if(n.eq.2) sigma = B5_cx_2_adas(eb)
        case (6)
            if(n.eq.1) sigma = C6_cx_1_adas(eb)
            if(n.eq.2) sigma = C6_cx_2_adas(eb)
            if(n.eq.3) sigma = C6_cx_3_adas(eb)
        case DEFAULT
            sigma = 0.d0
    end select  
      
end function Aq_cx_n_adas
  
function B5_cx_1_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: A = [31.226d0, 1.1442d0,    &
                                                    4.8372d-8, 3.0961d-10, &
                                                    4.7205d0, 6.2844d-7,   &
                                                    3.1297d0, 0.12556d0,   &
                                                    0.30098d0, 5.9607d-2,  &
                                                   -0.57923d0 ]
  
    sigma = 1.d-16*A(1)*(exp(-A(2)/eb**A(8)) /  &
           (1.0 + A(3)*eb**2 + A(4)*eb**A(5) +  &
            A(6)*eb**A(7)) + A(9)*exp(-A(10)*eb) /eb**A(11))
  
end function B5_cx_1_janev 
  
function C6_cx_1_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(11), parameter :: A = [418.18d0, 2.1585d0,   &
                                                    3.4808d-4, 5.3333d-9, &
                                                    4.6556d0, 0.33755d0,  &
                                                    0.81736d0, 0.27874d0, &
                                                    1.8003d-6, 7.1033d-2, &
                                                    0.53261d0 ]
  
    sigma = 1.d-16*A(1)*(exp(-A(2)/eb**A(8)) /  &
           (1.0 + A(3)*eb**2 + A(4)*eb**A(5) +  &
            A(6)*eb**A(7)) + A(9)*exp(-A(10)*eb) /eb**A(11))
  
end function C6_cx_1_janev 
  
function Aq_cx_n_janev(eb, q, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 1.507d5
    real(Float64), parameter :: B = 1.974d-5
  
    real(Float64) :: etil, nf, qf
  
    nf = real(n)
    qf = real(q)
  
    if((n.eq.1).and.(q.eq.5)) then
        sigma = B5_cx_1_janev(eb)
        return
    endif

    if((n.eq.1).and.(q.eq.6)) then
        sigma = C6_cx_1_janev(eb)
        return
    endif

    if(n.le.1) then
        sigma = 0.d0
        return
    endif
  
    etil = eb*(nf**2.0)/(qf**0.5)
  
    sigma = qf*nf**4 * 7.04d-16 * A/(etil**3.5 * (1.0 + B*etil**2)) * &
            (1.0 - exp(-2.0*etil**3.5 * (1.0 + B*etil**2)/(3.0*A)))
  
end function Aq_cx_n_janev
  
function Aq_cx_n(eb, q, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    sigma = Aq_cx_n_adas(eb, q, n)
    if(sigma.eq.0.d0) then
        sigma = Aq_cx_n_janev(eb, q, n)
    endif
  
end function Aq_cx_n
  
function Aq_cx(eb, q, n_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: n_max
    real(Float64), dimension(n_max) :: sigma
  
    integer :: n
  
    do n=1,n_max
        sigma(n) = Aq_cx_n(eb, q, n)
    enddo
  
end function Aq_cx
  
!Impurity impact ionization
function B5_ioniz_1_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: A = [351.52d0, 233.63d0,   &
                                                   3.2952d3, 5.3787d-6,  &
                                                   1.8834d-2, -2.2064d0, &
                                                   7.2074d0, -3.78664d0 ]
  
    sigma = 1.d-16*A(1)*(exp(-A(2)/eb)*log(1 + A(3)*eb)/eb &
          + A(4)*exp(-A(5)*eb)/((eb**A(6)) + A(7)*(eb**A(8))))
  
end function B5_ioniz_1_janev
  
function C6_ioniz_1_janev(eb) result(sigma)
    real(Float64), intent(in) :: eb
    real(Float64)             :: sigma
  
    real(Float64), dimension(8), parameter :: A = [ 438.36d0, 327.10d0,    &
                                                    1.4444d5, 3.5212d-3,   &
                                                    8.3031d-3, -0.63731d0, &
                                                    1.9116d4, -3.1003d0 ]
  
    sigma = 1.d-16*A(1)*(exp(-A(2)/eb)*log(1 + A(3)*eb)/eb &
          + A(4)*exp(-A(5)*eb)/((eb**A(6)) + A(7)*(eb**A(8))))
  
end function C6_ioniz_1_janev
  
function Aq_ioniz_n_janev(eb, q, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    real(Float64), parameter :: M = 0.283d0
    real(Float64), parameter :: B = 4.04d0
    real(Float64), parameter :: c = 137.d0
    real(Float64), parameter :: g = 0.662d0
    real(Float64), parameter :: lambda = 0.76d0
  
    real(Float64) :: nf, qf, u, v, sigma_b
  
    nf = real(n)
    qf = real(q)
    v = sqrt(eb/25.)
    u = nf*v
  
    sigma_b = 3.52d-16 * (nf**4) * (qf**2)/(u**2) * &
              (M * (log((u**2)/(c**2 - u**2)) - (u**2)/(c**2)) + B - g/u**2)
    sigma_b = max(sigma_b,0.d0)
  
    sigma = exp(-lambda*qf/u**2)*sigma_b
  
end function Aq_ioniz_n_janev
  
function Aq_ioniz_n(eb, q, n) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    real(Float64)             :: sigma
  
    if((q.eq.5).and.(n.eq.1)) then
        sigma = B5_ioniz_1_janev(eb)
        return
    endif
  
    if((q.eq.6).and.(n.eq.1)) then
        sigma = C6_ioniz_1_janev(eb)
        return
    endif
  
    sigma = Aq_ioniz_n_janev(eb, q, n)
  
end function Aq_ioniz_n
  
function Aq_ioniz(eb, q, n_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: n_max
    real(Float64), dimension(n_max) :: sigma
  
    integer :: n
  
    do n=1,n_max
        sigma(n) = Aq_ioniz_n(eb, q, n)
    enddo
  
end function Aq_ioniz
  
!Impurity impact excitation
function Aq_excit_1_2_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [38.738d0, 37.033d0,   &
                                                   0.39862d0, 7.7582d-5, &
                                                   0.25402d0, -2.7418d0 ]
    real(Float64) :: Etil, xsi, qf
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_1_2_janev
  
function Aq_excit_1_3_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [4.3619d0, 57.451d0,  &
                                                   21.001d0, 2.3292d-4, &
                                                   0.083130d0, -2.2364d0 ]
    real(Float64) :: Etil, xsi, qf
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_1_3_janev
  
function Aq_excit_1_4_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [1.3730d0, 60.710d0,  &
                                                   31.797d0, 2.0207d-4, &
                                                   0.082513d0, -2.3055d0 ]
    real(Float64) :: Etil, xsi, qf
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_1_4_janev
  
function Aq_excit_1_5_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [0.56565d0, 67.333d0, &
                                                   55.290d0, 2.1595d-4, &
                                                   0.081624d0, -2.1971d0 ]
    real(Float64) :: Etil, xsi, qf
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.d-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_1_5_janev
  
function Aq_excit_1_janev(eb, q, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
  
    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = Aq_excit_1_2_janev(eb, q)
            case (3)
                sigma(3) = Aq_excit_1_3_janev(eb, q)
            case (4)
                sigma(4) = Aq_excit_1_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_1_5_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(5)*(5.0/m)**3.0
        end select
    enddo
  
end function Aq_excit_1_janev
  
function Aq_excit_2_3_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [358.03d0, 25.283d0,   &
                                                   1.4726d0, 0.014398d0, &
                                                   0.12207d0, -0.86210d0 ]
    real(Float64) :: etil, qf, xsi
  
    qf = real(q) 
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_2_3_janev
  
function Aq_excit_2_4_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [50.744d0, 19.416d0,   &
                                                   4.0262d0, 0.014398d0, &
                                                   0.31584d0, -1.4799d0 ]
    real(Float64) :: etil, qf, xsi
  
    qf = real(q) 
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_2_4_janev
  
function Aq_excit_2_5_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [18.264d0, 18.973d0,   &
                                                   2.9056d0, 0.013701d0, &
                                                   0.31711d0, -1.4775d0 ]
    real(Float64) :: etil, qf, xsi
  
    qf = real(q) 
    etil = eb/qf
    xsi = 2.0**(0.5238*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1 + A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_2_5_janev
  
function Aq_excit_2_6_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.4610d0
  
    real(Float64) :: hi
  
    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))
  
    sigma = A*hi*Aq_excit_2_5_janev(eb, q)

end function Aq_excit_2_6_janev
  
function Aq_excit_2_7_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.2475d0
  
    real(Float64) :: hi
  
    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))
  
    sigma = A*hi*Aq_excit_2_5_janev(eb, q)
  
end function Aq_excit_2_7_janev
  
function Aq_excit_2_8_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.1465d0
  
    real(Float64) :: hi
  
    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))
  
    sigma = A*hi*Aq_excit_2_5_janev(eb, q)
  
end function Aq_excit_2_8_janev
  
function Aq_excit_2_9_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.092d0
  
    real(Float64) :: hi
  
    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))
  
    sigma = A*hi*Aq_excit_2_5_janev(eb, q)
  
end function Aq_excit_2_9_janev
  
function Aq_excit_2_10_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.0605d0
  
    real(Float64) :: hi
  
    hi = 2.0**(0.397*(1.0 - sqrt(2.0/q)))
  
    sigma = A*hi*Aq_excit_2_5_janev(eb, q)
  
end function Aq_excit_2_10_janev
  
function Aq_excit_2_janev(eb, q, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
  
    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = Aq_excit_2_3_janev(eb, q)
            case (4)
                sigma(4) = Aq_excit_2_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_2_5_janev(eb, q)
            case (6)
                sigma(6) = Aq_excit_2_6_janev(eb, q)
            case (7)
                sigma(7) = Aq_excit_2_7_janev(eb, q)
            case (8)
                sigma(8) = Aq_excit_2_8_janev(eb, q)
            case (9)
                sigma(9) = Aq_excit_2_9_janev(eb, q)
            case (10)
                sigma(10) = Aq_excit_2_10_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/m)**3.0
        end select
    enddo
  
end function Aq_excit_2_janev
  
function Aq_excit_3_4_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [1247.5d0, 11.319d0,   &
                                                   2.6235d0, 0.068781d0, &
                                                   0.521176d0, -1.2722d0 ]
  
    real(Float64) :: Etil, qf, xsi
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_3_4_janev
  
function Aq_excit_3_5_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [190.59d0, 11.096d0,   &
                                                   2.9098d0, 0.073307d0, &
                                                   0.54177d0, -1.2894d0 ]
  
    real(Float64) :: Etil, qf, xsi
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_3_5_janev
  
function Aq_excit_3_6_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), dimension(6), parameter :: A = [63.494d0, 11.507d0,   &
                                                   4.3417d0, 0.077953d0, &
                                                   0.53461d0, -1.2881d0 ]
  
    real(Float64) :: Etil, qf, xsi
  
    qf = real(q)
    etil = eb/qf
    xsi = 2.0**(0.397*(1 - sqrt(2.0/qf)))
    sigma = qf*1.e-16*xsi*A(1)*(exp(-A(2)/etil)*log(1+A(3)*etil)/etil &
          + A(4)*exp(-A(5)*etil)/etil**A(6))
  
end function Aq_excit_3_6_janev
  
function Aq_excit_3_7_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.4670d0
  
    real(Float64) :: hi
  
    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)
  
end function Aq_excit_3_7_janev
  
function Aq_excit_3_8_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
  
    real(Float64), parameter :: A = 0.2545d0
  
    real(Float64) :: hi
  
    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)
  
end function Aq_excit_3_8_janev
  
function Aq_excit_3_9_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 0.1540d0
  
    real(Float64) :: hi
  
    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)
  
end function Aq_excit_3_9_janev
  
function Aq_excit_3_10_janev(eb, q) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    real(Float64)             :: sigma
    
    real(Float64), parameter :: A = 0.1d0
  
    real(Float64) :: hi
  
    hi=2.0**(0.397*(1.0 - sqrt(2.0/q)))
    sigma=hi*A*Aq_excit_3_6_janev(eb, q)
  
end function Aq_excit_3_10_janev
  
function Aq_excit_3_janev(eb, q, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
  
    sigma = 0.d0
    do m=1,m_max
        select case (m)
            case (1)
                sigma(1) = 0.d0
            case (2)
                sigma(2) = 0.d0
            case (3)
                sigma(3) = 0.d0
            case (4)
                sigma(4) = Aq_excit_3_4_janev(eb, q)
            case (5)
                sigma(5) = Aq_excit_3_5_janev(eb, q)
            case (6)
                sigma(6) = Aq_excit_3_6_janev(eb, q)
            case (7)
                sigma(7) = Aq_excit_3_7_janev(eb, q)
            case (8)
                sigma(8) = Aq_excit_3_8_janev(eb, q)
            case (9)
                sigma(9) = Aq_excit_3_9_janev(eb, q)
            case (10)
                sigma(10) = Aq_excit_3_10_janev(eb, q)
            case DEFAULT
                sigma(m) = sigma(10)*(10.0/m)**3.0
        end select
    enddo
  
end function Aq_excit_3_janev
  
function Aq_excit_n_janev(eb, q, n, m_max) result(sigma)
    real(Float64), intent(in)       :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: n
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    integer :: m
    real(Float64) :: nf, mf, qf, etil, hi, s
    real(Float64) :: D, A, G, L, F, H, y, zpl, zmi, C2pl, C2mi
  
    nf = real(n)
    qf = real(q)
    sigma = 0.d0
    m_loop: do m=1,m_max
        mf = real(m)
        if(n.ge.m) then
            sigma(m) = 0.d0
            cycle m_loop
        endif
        etil = eb/(25.0*qf)
        hi = 2.0**(0.322*(1.0 - sqrt(2.0/qf)))
        s = (mf - nf)

        D = exp(-1.0/(nf*mf*etil**2))
        A = 8.0/(3.0*s) * (mf/(s*nf))**3 * (0.184 - 0.04/s**(2.0/3.0)) * &
           (1.0 - 0.2*s/(nf*mf))**(1.0 + 2.0*s)
        G = 0.5*(etil*nf**2.0 / (mf - 1.0/mf))**3.0
        L = log(1.0 + 0.53*etil**2.0 * nf*(mf - 2.0/mf )/(1.0 + 0.4*etil))
        F = (1.0 - 0.3*s*D/(nf*mf))**(1.0 + 2.0*s)
        
        y = 1.0/(1.0 - D*log(18*s)/(4.0*s))
        zpl = 2.0/(etil*nf**2*(sqrt(2.0 - nf**2/mf**2) + 1.0))
        zmi = 2.0/(etil*nf**2*(sqrt(2.0 - nf**2/mf**2) - 1.0))
        C2pl = zpl**2*log(1.0 + 2.0*zpl/3.0)/(2.0*y + 3.0*zpl/2.0)
        C2mi = zmi**2*log(1.0 + 2.0*zmi/3.0)/(2.0*y + 3.0*zmi/2.0)
        H = C2mi - C2pl

        sigma(m) = q*hi*8.86e-17*nf**4/etil*(A*D*L+F*G*H)
    enddo m_loop
  
end function Aq_excit_n_janev
  
function Aq_excit_n(eb, q, n, m_max) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)             :: q
    integer, intent(in)             :: n
    integer, intent(in)             :: m_max
    real(Float64), dimension(m_max) :: sigma
  
    select case (n)
        case (0)
            stop
        case (1)
            sigma = Aq_excit_1_janev(eb, q, m_max)
        case (2)
            sigma = Aq_excit_2_janev(eb, q, m_max)
        case (3)
            sigma = Aq_excit_3_janev(eb, q, m_max)
        case DEFAULT
            sigma = Aq_excit_n_janev(eb, q, n, m_max)
    end select
  
end function Aq_excit_n
  
function Aq_excit_n_m(eb ,q, n, m) result(sigma)
    real(Float64), intent(in) :: eb
    integer, intent(in)       :: q
    integer, intent(in)       :: n
    integer, intent(in)       :: m
    real(Float64)             :: sigma
  
    real(Float64), dimension(12) :: sigma_m
  
    sigma_m = Aq_excit_n(eb, q, n, 12)
    if(m.le.0) then
        sigma = sum(sigma_m)
    else
        sigma = sigma_m(m)
    endif
  
end function Aq_excit_n_m
  
function Aq_excit(eb, q, n_max, m_max) result(sigma)
    real(Float64), intent(in)              :: eb
    integer, intent(in)                    :: q
    integer, intent(in)                    :: n_max
    integer, intent(in)                    :: m_max
    real(Float64), dimension(n_max, m_max) :: sigma
  
    real(Float64), dimension(12,12) :: sigma_full
    integer :: n, m
  
    do n=1,12
        sigma_full(n,:) = Aq_excit_n(eb, q, n, 12)
    enddo
  
    sigma = sigma_full(1:n_max,1:m_max)
  
end function Aq_excit
  
function simpsons_rule(f, dx) result(I)
    real(Float64), dimension(:), intent(in) :: f
    real(Float64), intent(in)               :: dx
    real(Float64)                           :: I
  
    integer :: s, ii
  
    s = size(f)
    I = 0.d0
    if(mod(s,2).eq.1) then
        write(*,'(a)') "Length of array must be even"
        return
    endif
  
    I = f(1)
    do ii=2,s-1
        if(mod(ii,2).eq.1) then
            I = I + 4.0*f(ii)
        else
            I = I + 2.0*f(ii)
        endif
    enddo
    I = I + f(s)
    I = (dx/3.0)*I
    
end function simpsons_rule
  
subroutine bt_maxwellian_n(fn, T, eb, am, ab, n, rate)
    interface
        function fn(a, b)
            real(8)              :: fn !sigma
            real(8), intent(in)  :: a !eb
            integer, intent(in)  :: b !n
        end function fn
    end interface
    real(Float64), intent(in)  :: T
    real(Float64), intent(in)  :: eb
    real(Float64), intent(in)  :: am
    real(Float64), intent(in)  :: ab
    integer, intent(in)        :: n
    real(Float64), intent(out) :: rate
  
    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max,dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz 
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, v_therm, dE
   
    integer :: i, j
  
    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo
  
    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo
  
    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
    dE = (13.6d-3)/(n**2.0)
  
    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0
  
    fz = 0.d0
    fr = 0.d0
  
    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(erel.ge.dE) then
                sig = fn(erel/ared,n)
            else
                sig = 0.d0
            endif
            fr(j) = sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo
  
    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
    
end subroutine bt_maxwellian_n
  
subroutine bt_maxwellian_q_n(fqn, q, T, eb, am, ab, n, rate)
    interface
        function fqn(a, b, c)
            real(8)             :: fqn !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !q
            integer, intent(in) :: c !n
        end function fqn
    end interface
    integer, intent(in)        :: q
    real(Float64), intent(in)  :: T
    real(Float64), intent(in)  :: eb
    real(Float64), intent(in)  :: am
    real(Float64), intent(in)  :: ab
    integer, intent(in)        :: n
    real(Float64), intent(out) :: rate
  
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz 
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, v_therm, dE
   
    integer :: i, j
  
    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo
  
    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo
  
    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
    dE = (13.6d-3)/(n**2.0)
  
    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0
  
    fz = 0.d0
    fr = 0.d0
  
    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(erel.ge.dE) then
                sig = fqn(erel/ared, q, n)
            else
                sig = 0.d0
            endif
            fr(j) = sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo
  
    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
  
end subroutine bt_maxwellian_q_n
  
subroutine bt_maxwellian_n_m(fnm, T, eb, am, ab, n, m, rate, deexcit)
    interface
        function fnm(a, b, c)
            real(8)             :: fnm !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !n
            integer, intent(in) :: c !m
        end function fnm
    end interface
    real(Float64), intent(in)     :: T
    real(Float64), intent(in)     :: eb
    real(Float64), intent(in)     :: am
    real(Float64), intent(in)     :: ab
    integer, intent(in)           :: n
    integer, intent(in)           :: m
    real(Float64), intent(out)    :: rate
    logical, intent(in), optional :: deexcit
  
    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz 
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, dE, factor, En, Em, v_therm
   
    integer :: i, j
  
    if(present(deexcit)) then
        dxc = deexcit
    else
        dxc = .False.
    endif
  
    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo
  
    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo
  
    En = (13.6d-3)*(1.0 - (1.d0/n)**2.0)
    Em = (13.6d-3)*(1.0 - (1.d0/m)**2.0)
    dE = Em - En
  
    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
  
    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0
  
    fz = 0.d0
    fr = 0.d0
  
    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(dxc) then
                factor = (erel + dE)/erel
                erel = erel + dE
            else
                factor = 1.0
            endif
            if(erel.ge.dE) then
                sig = fnm(erel/ared, n, m)
            else
                sig = 0.d0
            endif
            fr(j) = factor*sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo
  
    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
    if(dxc) rate = rate*(real(n)/real(m))**2.0
  
end subroutine bt_maxwellian_n_m
  
subroutine bt_maxwellian_q_n_m(fqnm, q, T, eb, am, ab, n, m, rate, deexcit)
    interface
        function fqnm(a, b, c, d)
            real(8)             :: fqnm !sigma
            real(8), intent(in) :: a !eb
            integer, intent(in) :: b !q
            integer, intent(in) :: c !n
            integer, intent(in) :: d !m
        end function fqnm
    end interface
    integer, intent(in) :: q
    real(Float64), intent(in)     :: T
    real(Float64), intent(in)     :: eb
    real(Float64), intent(in)     :: am
    real(Float64), intent(in)     :: ab
    integer, intent(in)           :: n
    integer, intent(in)           :: m
    real(Float64), intent(out)    :: rate
    logical, intent(in), optional :: deexcit
  
    logical :: dxc
    integer :: n_vr
    real(Float64) :: vr_max, dvr
    real(Float64), dimension(32) :: vr
    real(Float64), dimension(32) :: fr
    integer :: n_vz
    real(Float64) :: vz_max, dvz
    real(Float64), dimension(62) :: vz
    real(Float64), dimension(62) :: fz 
    real(Float64) :: T_per_amu, eb_per_amu, ared, sig, sig_eff
    real(Float64) :: zb, u2_to_erel, u2, erel, dE, factor, En, Em, v_therm
   
    integer :: i, j
  
    if(present(deexcit)) then
        dxc = deexcit
    else
        dxc = .False.
    endif
  
    n_vr = 32
    vr_max = 4.d0
    dvr = vr_max/(n_vr - 1.d0)
    do i=1,n_vr
        vr(i) = (i-1)*dvr
    enddo
  
    n_vz = 62
    vz_max = 4.d0
    dvz = 2.0*vz_max/(n_vz - 1.d0)
    do i=1,n_vz
        vz(i) = (i-1)*dvz - vz_max
    enddo
  
    En = (13.6d-3)*(1.0 - (1.d0/n)**2.0)
    Em = (13.6d-3)*(1.0 - (1.d0/m)**2.0)
    dE = Em - En
  
    T_per_amu = max(T, 1.d-6)/am
    eb_per_amu = eb/ab
    ared = am*ab/(am + ab)
  
    v_therm = 1.384d6 * sqrt(T_per_amu*1.d3)
    zb = sqrt(eb_per_amu/T_per_amu)
    u2_to_erel = ared*T_per_amu
    if(ared.lt.0.5) ared = 1.0
  
    fz = 0.d0
    fr = 0.d0
  
    do i=1,n_vz
        do j=1,n_vr
            u2 = (zb - vz(i))**2.0 + vr(j)**2.0
            erel = u2_to_erel*u2
            if(dxc) then
                factor = (erel + dE)/erel
                erel = erel + dE
            else
                factor = 1.0
            endif
            if(erel.ge.dE) then
                sig = fqnm(erel/ared, q, n, m)
            else
                sig = 0.d0
            endif
            fr(j) = factor*sig*sqrt(u2)*exp(-(vz(i)**2.0 + vr(j)**2.0))*vr(j)
        enddo
        fz(i) = simpsons_rule(fr, dvr)
    enddo
  
    sig_eff = (2.0/sqrt(PI))*simpsons_rule(fz, dvz)
    rate = sig_eff*v_therm
    if(dxc) rate = rate*(real(n)/real(m))**2.0

end subroutine bt_maxwellian_q_n_m
  
subroutine write_einstein(id, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max
  
    real(Float64), dimension(n_max,m_max) :: ein
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer :: error
  
    ein(:,:) = EINSTEIN(1:n_max,1:m_max)
   
    call h5gcreate_f(id, "spontaneous", gid, error)
  
    dim1 = [1]
    dim2 = [n_max, m_max]
  
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
  
    call h5ltmake_compressed_dataset_double_f(gid, "einstein", 2, dim2, ein, error)
  
    call h5ltset_attribute_string_f(id, "spontaneous", "description", &
         "Atomic rates for spontaneous emission/deexcitation", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
  
    call h5ltset_attribute_string_f(gid, "einstein", "description", &
         "n/m resolved einstein coefficients: einstein(n,m)", error)
    call h5ltset_attribute_string_f(gid, "einstein", "units", "1/s", error)
    call h5ltset_attribute_string_f(gid, "einstein", "reaction", &
         "H(n) -> H(m) + ph, n > m", error)
  
    call h5gclose_f(gid, error)
  
end subroutine write_einstein
  
subroutine write_bb_H_H(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy
  
    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64), dimension(:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:), allocatable :: cx, excit
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
  
    integer :: i, cnt, error
    logical :: exis

    NAMELIST /H_H_cross/ nenergy, emin, emax

    nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_H: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_H_cross)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(ioniz(n_max,nenergy))
    allocate(cx(n_max,m_max,nenergy))
    allocate(excit(n_max,m_max,nenergy))

    ebarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0
  
    write(*,'(a)') "---- H-H cross sections settings ----"
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,*) ''

    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=1, nenergy
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb
  
        cx(:,:,i) = p_cx(eb, n_max, m_max)
        excit(:,:,i) = p_excit(eb, n_max, m_max)
        ioniz(:,i) = p_ioniz(eb, n_max)
        cnt = cnt + 1
        WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO
  
    call h5gcreate_f(id, "H_H", gid, error)
  
    dim1 = [1]
    dim2 = [n_max, nenergy]
    dim3 = [n_max, m_max, nenergy]
  
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
  
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "cx", 3, dim3, cx, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)
  
    call h5ltset_attribute_string_f(id, "H_H", "description", &
         "Cross sections for Hydrogen-Hydrogen interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)  
  
    call h5ltset_attribute_string_f(gid, "cx", "description", &
         "n/m resolved charge exchange cross sections: cx(n,m,energy)", error)
    call h5ltset_attribute_string_f(gid, "cx", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "cx", "reaction", &
         "H(+) + H(n) -> H(m) + H(+)", error)
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "H(+) + H(n) -> H(+) + H(m), m > n", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization cross sections: ionization(n,energy)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "H(+) + H(n) -> H(+) + H(+) + e-", error)
  
    call h5gclose_f(gid, error)
  
    deallocate(ebarr, cx, excit, ioniz)

end subroutine write_bb_H_H
  
subroutine write_bb_H_e(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max
  
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64), dimension(:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:), allocatable :: excit
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
  
    integer :: i, cnt, error
    logical :: exis
  
    NAMELIST /H_e_cross/ nenergy, emin, emax

    nenergy = 200; emin = 1.d-3 ; emax = 8.d2

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_E: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_e_cross)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(ioniz(n_max,nenergy))
    allocate(excit(n_max,m_max,nenergy))

    ebarr = 0.d0
    ioniz = 0.d0
    excit = 0.d0
  
    write(*,'(a)') "---- H-e cross sections settings ----"
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,*) ''
  
    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=1, nenergy
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb
  
        excit(:,:,i) = e_excit(eb, n_max, m_max)
        ioniz(:,i) = e_ioniz(eb, n_max)
  
        cnt = cnt + 1
        WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy),char(13)
    enddo
    !$OMP END PARALLEL DO
  
    call h5gcreate_f(id, "H_e", gid, error)
  
    dim1 = [1]
    dim2 = [n_max, nenergy]
    dim3 = [n_max, m_max, nenergy]
  
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
  
  
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)
  
    call h5ltset_attribute_string_f(id, "H_e", "description", &
         "Cross sections for Hydrogen-Electron interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum Energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum Energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)  
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "e- + H(n) -> e- + H(m), m > n", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization cross sections: ionization(n,energy)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "e- + H(n) -> e- + H(+) + e-", error)
  
    call h5gclose_f(gid, error)

    deallocate(ebarr, ioniz, excit)

end subroutine write_bb_H_e
  
subroutine write_bb_H_Aq(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max

    integer :: q
    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy
  
    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64), dimension(:,:), allocatable :: cx, ioniz
    real(Float64), dimension(:,:,:), allocatable :: excit
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
  
    character(len=10) :: aname
    character(len=5) :: asym
    integer :: i, cnt, error
    logical :: exis
  
    NAMELIST /H_Aq_cross/ q, nenergy, emin, emax

    nenergy = 200; emin = 1.d-3 ; emax = 8.d2
    q = 6

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BB_H_Aq: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_Aq_cross)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(ioniz(n_max,nenergy))
    allocate(cx(n_max,nenergy))
    allocate(excit(n_max,m_max,nenergy))

    ebarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0

    select case (q)
        case (5)
            aname = "Boron"
            asym = "H_B5"
        case (6)
            aname = "Carbon"
            asym = "H_C6"
        case DEFAULT
            aname = "Impurity"
            asym = "H_Aq"
    end select
  
    write(*,'(a)') "---- H-"//trim(adjustl(aname))//" cross sections settings ----"
    write(*,'(T2,"q = ", i2)'), q
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,*) ''
  
    cnt = 0
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    !$OMP PARALLEL DO private(i, eb)
    do i=1, nenergy
        eb = 10.d0**(log10(emin) + (i-1)*dlogE)
        ebarr(i) = eb
  
        cx(:,i) = Aq_cx(eb, q, n_max)
        ioniz(:,i) = Aq_ioniz(eb, q, n_max)
        excit(:,:,i) = Aq_excit(eb, q, n_max, m_max)
        cnt = cnt + 1
        WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy),char(13)
    enddo
  
    call h5gcreate_f(id, trim(adjustl(asym)), gid, error)
  
    dim1 = [1]
    dim2 = [n_max, nenergy]
    dim3 = [n_max, m_max, nenergy]
  
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
  
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "cx", 2, dim2, cx, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 2, dim2, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 3, dim3, excit, error)
  
    call h5ltset_attribute_string_f(id, trim(adjustl(asym)), "description", &
         "Cross sections for Hydrogen-"//trim(adjustl(aname))//" interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Nucleon energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV/amu", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV/amu)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV/amu", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV/amu", error)  
  
    call h5ltset_attribute_string_f(gid, "cx", "description", &
         "n resolved charge exchange / electron capture cross sections: cx(n,energy)", error)
    call h5ltset_attribute_string_f(gid, "cx", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "cx", "reaction", &
         "A(q+) + H(n) -> A((q-1)+) + H(+)", error)
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved excitation cross sections: excitation(n,m,energy)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "A(q+) + H(n) -> A(q+) + H(m), m > n", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization cross sections: ionization(n,energy)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^2", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "A(q+) + H(n) -> A(q+) + H(+) + e-", error)
  
    call h5gclose_f(gid, error)
  
    deallocate(ebarr, ioniz, cx, excit)

end subroutine write_bb_H_Aq
  
subroutine write_bt_H_H(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp
  
    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr
  
    real(Float64), dimension(:,:,:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit, cx
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
  
    integer :: ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 4
    real(Float64), dimension(2,n_bt_amu) :: a
    logical :: exis

    NAMELIST /H_H_rates/ nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_H: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_H_rates)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
    allocate(cx(n_max,m_max,nenergy,ntemp,n_bt_amu))
    allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))

    ebarr = 0.d0
    tarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0
    a(:,1) = [H1_amu, H1_amu]
    a(:,2) = [H1_amu, H2_amu]
    a(:,3) = [H2_amu, H1_amu]
    a(:,4) = [H2_amu, H2_amu]
     
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo
  
    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo
  
    write(*,'(a)') "---- H-H reaction rates settings ----"
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
    write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
    write(*,'(T2,"Ntemp = ", i4)') ntemp
    write(*,*) ''
  
    cnt = 0
    !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
    do ie=1, nenergy
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            do ia=1, n_bt_amu
                do n=1, n_max
                    do m=1, m_max
                        call bt_maxwellian(p_cx_n_m, ti, eb, &
                                           a(2,ia), a(1,ia), n, m, rate)
                        cx(n,m,ie,it,ia) = rate
                        if(m.gt.n) then
                            call bt_maxwellian(p_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, rate)
                            excit(n,m,ie,it,ia) = rate

                            call bt_maxwellian(p_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, &
                                               rate, deexcit=.True.)
                            excit(m,n,ie,it,ia) = rate
                        endif
                    enddo
                    call bt_maxwellian(p_ioniz_n, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    ioniz(n,ie,it,ia) = rate
                enddo
            enddo
            cnt = cnt + 1
            WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO
  
    call h5gcreate_f(id, "H_H", gid, error)
  
    dim1 = [1]
    dim4 = [n_max, nenergy, ntemp, n_bt_amu]
    dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
    call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
    call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
    call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
    call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)
  
    dim2 = [2,n_bt_amu]
    call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error) 
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    dim1 = [ntemp]
    call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "cx", 5, dim5, cx, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)
  
    call h5ltset_attribute_string_f(id, "H_H", "description", &
         "Reaction rates for Hydrogen(beam)-Hydrogen(target) interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(gid, "ntemp", "description", &
         "Number of target temperature values", error)
    call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
         "Number of beam-target amu combinations", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
  
    call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
         "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "temperature", "description", &
         "Target temperature values", error)
    call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "description", &
         "Temperature spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "tmin","description", &
         "Minimum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "tmax","description", &
         "Maximum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "cx", "description", &
         "n/m resolved charge exchange reaction rates: cx(n,m,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "cx", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "cx", "reaction", &
         "H(+) + H(n) -> H(m) + H(+)", error)
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "H(+) + H(n) -> H(+) + H(m); m > n excitation, m < n de-excitation", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "H(+) + H(n) -> H(+) + H(+) + e-", error)
  
    call h5gclose_f(gid, error)

    deallocate(ebarr, tarr, cx, excit, ioniz)

end subroutine write_bt_H_H
  
subroutine write_bt_H_e(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp
  
    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr
  
    real(Float64), dimension(:,:,:,:), allocatable :: ioniz
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
  
    integer :: ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 2
    real(Float64), dimension(2,n_bt_amu) :: a
    logical :: exis

    NAMELIST /H_e_rates/ nenergy, emin, emax, ntemp, tmin, tmax

    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_E: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_e_rates)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
    allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))
  
    ebarr = 0.d0
    ioniz = 0.d0
    excit = 0.d0
    a(:,1) = [H1_amu, e_amu]
    a(:,2) = [H2_amu, e_amu]
     
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo
  
    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo
  
    write(*,'(a)') "---- H-e reaction rates settings ----"
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
    write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
    write(*,'(T2,"Ntemp = ", i4)') ntemp
    write(*,*) ''
  
    cnt = 0
    !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
    do ie=1, nenergy
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            do ia=1, n_bt_amu
                do n=1, n_max
                    do m=1, m_max
                        if(m.gt.n) then
                            call bt_maxwellian(e_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, rate)
                            excit(n,m,ie,it,ia) = rate

                            call bt_maxwellian(e_excit_n_m, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, &
                                               rate, deexcit=.True.)
                            excit(m,n,ie,it,ia) = rate
                        endif
                    enddo
                    call bt_maxwellian(e_ioniz_n, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    ioniz(n,ie,it,ia) = rate
                enddo
            enddo
            cnt = cnt + 1
            WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO
  
    call h5gcreate_f(id, "H_e", gid, error)
  
    dim1 = [1]
    dim4 = [n_max, nenergy, ntemp, n_bt_amu]
    dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
    call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
    call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
    call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
    call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)
  
    dim2 = [2,n_bt_amu]
    call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error) 
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    dim1 = [ntemp]
    call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)
  
    call h5ltset_attribute_string_f(id, "H_e", "description", &
         "Reaction rates for Hydrogen(beam)-Electron(target) interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(gid, "ntemp", "description", &
         "Number of target temperature values", error)
    call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
         "Number of beam-target amu combinations", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
  
    call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
         "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "temperature", "description", &
         "Target temperature values", error)
    call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "description", &
         "Temperature spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "tmin","description", &
         "Minimum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "tmax","description", &
         "Maximum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "e- + H(n) -> e- + H(m); m > n excitation, m < n de-excitation", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "e- + H(n) -> e- + H(+) + e-", error)
  
    call h5gclose_f(gid, error)
  
    deallocate(ebarr, tarr, excit, ioniz)

end subroutine write_bt_H_e
  
subroutine write_bt_H_Aq(id, namelist_file, n_max, m_max)
    integer(HID_T), intent(inout) :: id
    character(len=*), intent(in)  :: namelist_file
    integer, intent(in)           :: n_max
    integer, intent(in)           :: m_max

    integer :: q
    real(Float64) :: mass

    real(Float64) :: emin
    real(Float64) :: emax
    integer :: nenergy

    real(Float64) :: tmin
    real(Float64) :: tmax
    integer :: ntemp
  
    real(Float64) :: eb
    real(Float64) :: dlogE
    real(Float64), dimension(:), allocatable :: ebarr
  
    real(Float64) :: ti
    real(Float64) :: dlogT
    real(Float64), dimension(:), allocatable :: tarr
  
    real(Float64), dimension(:,:,:,:), allocatable :: ioniz, cx
    real(Float64), dimension(:,:,:,:,:), allocatable :: excit
  
    integer(HID_T) :: gid
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(2) :: dim2
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
  
    integer :: ie, it, ia, n, m, error, cnt
    real(Float64) :: rate
    integer, parameter :: n_bt_amu = 2
    real(Float64), dimension(2,n_bt_amu) :: a
  
    character(len=10) :: aname
    character(len=5) :: asym
    logical :: exis

    NAMELIST /H_Aq_rates/ q, mass, nenergy, emin, emax, ntemp, tmin, tmax

    q = 6 ; mass = C_amu
    nenergy = 100; emin = 1.d-3 ; emax = 4.d2
    ntemp = 100; tmin = 1.d-3 ; tmax = 2.d1
    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'WRITE_BT_H_Aq: Input file does not exist: ',trim(namelist_file)
        write(*,'(a)') 'Continuing with default settings...'
    else
        open(13,file=namelist_file)
        read(13,NML=H_Aq_rates)
        close(13)
    endif

    allocate(ebarr(nenergy))
    allocate(tarr(ntemp))
    allocate(ioniz(n_max,nenergy,ntemp,n_bt_amu))
    allocate(cx(n_max,nenergy,ntemp,n_bt_amu))
    allocate(excit(n_max,m_max,nenergy,ntemp,n_bt_amu))
  
  
    select case (q)
      case (5)
          aname = "Boron"
          asym = "H_B5"
      case (6)
          aname = "Carbon"
          asym = "H_C6"
      case DEFAULT
          aname = "Impurity"
          asym = "H_Aq"
    end select
  
    ebarr = 0.d0
    ioniz = 0.d0
    cx = 0.d0
    excit = 0.d0
    a(:,1) = [H1_amu, mass]
    a(:,2) = [H2_amu, mass]
     
    dlogE = (log10(emax) - log10(emin))/(nenergy - 1)
    do ie=1, nenergy
        ebarr(ie) = 10.d0**(log10(emin) + (ie-1)*dlogE)
    enddo
  
    dlogT = (log10(tmax) - log10(tmin))/(ntemp - 1)
    do it=1, ntemp
        tarr(it) = 10.d0**(log10(tmin) + (it-1)*dlogT)
    enddo
    write(*,'(a)') "---- H-"//trim(adjustl(aname))//" reaction rates settings ----"
    write(*,'(T2,"q = ", i2)'), q
    write(*,'(T2,"mass = ",f7.2, " amu")') mass
    write(*,'(T2,"Emin = ",e9.2, " keV")') emin
    write(*,'(T2,"Emax = ",e9.2, " keV")') emax
    write(*,'(T2,"Nenergy = ", i4)') nenergy
    write(*,'(T2,"Tmin = ",e9.2, " keV")') tmin
    write(*,'(T2,"Tmax = ",e9.2, " keV")') tmax
    write(*,'(T2,"Ntemp = ", i4)') ntemp
    write(*,*) ''
  
    cnt = 0
    !$OMP PARALLEL DO private(ie, it, ia, n, m, eb, ti, rate)
    do ie=1, nenergy
        eb = ebarr(ie)
        do it=1, ntemp
            ti = tarr(it)
            do ia=1, n_bt_amu
                do n=1, n_max
                    do m=1, m_max
                        if(m.gt.n) then
                            call bt_maxwellian(Aq_excit_n_m, q, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, rate)
                            excit(n,m,ie,it,ia) = rate

                            call bt_maxwellian(Aq_excit_n_m, q, ti, eb, &
                                               a(2,ia), a(1,ia), n, m, &
                                               rate, deexcit=.True.)
                            excit(m,n,ie,it,ia) = rate
                        endif
                    enddo
                    call bt_maxwellian(Aq_cx_n, q, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    cx(n,ie,it,ia) = rate

                    call bt_maxwellian(Aq_ioniz_n, q, ti, eb, &
                                       a(2,ia), a(1,ia), n, rate)
                    ioniz(n,ie,it,ia) = rate
                enddo
            enddo
            cnt = cnt + 1
            WRITE(*,'(f7.2,"%",a,$)') 100*cnt/real(nenergy*ntemp),char(13)
        enddo
    enddo
    !$OMP END PARALLEL DO
  
    call h5gcreate_f(id, trim(adjustl(asym)), gid, error)
  
    dim1 = [1]
    dim4 = [n_max, nenergy, ntemp, n_bt_amu]
    dim5 = [n_max, m_max, nenergy, ntemp, n_bt_amu]
    call h5ltmake_dataset_int_f(gid, "nenergy", 0, dim1, [nenergy], error)
    call h5ltmake_dataset_int_f(gid, "ntemp", 0, dim1, [ntemp], error)
    call h5ltmake_dataset_int_f(gid, "n_bt_amu", 0, dim1, [n_bt_amu], error)
    call h5ltmake_dataset_int_f(gid, "n_max", 0, dim1, [n_max], error)
    call h5ltmake_dataset_int_f(gid, "m_max", 0, dim1, [m_max], error)
    call h5ltmake_dataset_double_f(gid, "dlogE", 0, dim1, [dlogE], error)
    call h5ltmake_dataset_double_f(gid, "emin", 0, dim1, [emin], error)
    call h5ltmake_dataset_double_f(gid, "emax", 0, dim1, [emax], error)
    call h5ltmake_dataset_double_f(gid, "dlogT", 0, dim1, [dlogT], error)
    call h5ltmake_dataset_double_f(gid, "tmin", 0, dim1, [tmin], error)
    call h5ltmake_dataset_double_f(gid, "tmax", 0, dim1, [tmax], error)
  
    dim2 = [2,n_bt_amu]
    call h5ltmake_compressed_dataset_double_f(gid, "bt_amu", 2, dim2, a, error) 
    dim1 = [nenergy]
    call h5ltmake_compressed_dataset_double_f(gid, "energy", 1, dim1, ebarr, error)
    dim1 = [ntemp]
    call h5ltmake_compressed_dataset_double_f(gid, "temperature", 1, dim1, tarr, error)
    call h5ltmake_compressed_dataset_double_f(gid, "cx", 4, dim4, cx, error)
    call h5ltmake_compressed_dataset_double_f(gid, "ionization", 4, dim4, ioniz, error)
    call h5ltmake_compressed_dataset_double_f(gid, "excitation", 5, dim5, excit, error)
  
    call h5ltset_attribute_string_f(id, trim(adjustl(asym)), "description", &
         "Reaction rates for Hydrogen(beam)-"//trim(adjustl(aname))// &
         "(target) interactions", error)
    call h5ltset_attribute_string_f(gid, "nenergy", "description", &
         "Number of energy values", error)
    call h5ltset_attribute_string_f(gid, "ntemp", "description", &
         "Number of target temperature values", error)
    call h5ltset_attribute_string_f(gid, "n_bt_amu", "description", &
         "Number of beam-target amu combinations", error)
    call h5ltset_attribute_string_f(gid, "n_max", "description", &
         "Number of initial energy levels", error)
    call h5ltset_attribute_string_f(gid, "m_max", "description", &
         "Number of final energy levels", error)
  
    call h5ltset_attribute_string_f(gid, "bt_amu", "description", &
         "Combinations of beam-target amu's e.g. b_amu, t_amu = bt_amu[:,i]", error)
    call h5ltset_attribute_string_f(gid, "energy", "description", &
         "Energy values", error)
    call h5ltset_attribute_string_f(gid, "energy", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "description", &
         "Energy spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogE", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "emin","description", &
         "Minimum energy", error)
    call h5ltset_attribute_string_f(gid, "emin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "emax","description", &
         "Maximum energy", error)
    call h5ltset_attribute_string_f(gid, "emax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "temperature", "description", &
         "Target temperature values", error)
    call h5ltset_attribute_string_f(gid, "temperature", "units", "keV", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "description", &
         "Temperature spacing in log-10", error)
    call h5ltset_attribute_string_f(gid, "dlogT", "units", "log10(keV)", error)
    call h5ltset_attribute_string_f(gid, "tmin","description", &
         "Minimum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmin", "units", "keV", error)  
    call h5ltset_attribute_string_f(gid, "tmax","description", &
         "Maximum temperature", error)
    call h5ltset_attribute_string_f(gid, "tmax", "units", "keV", error)  
  
    call h5ltset_attribute_string_f(gid, "cx", "description", &
         "n-resolved charge exchange reaction rates: cx(n,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "cx", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "cx", "reaction", &
         "A(q+) + H(n) -> A((q-1)+) + H(+)", error)
  
    call h5ltset_attribute_string_f(gid, "excitation", "description", &
         "n/m resolved (de-)excitation reaction rates: excitation(n,m,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "excitation", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "excitation", "reaction", &
         "A(q+) + H(n) -> A(q+) + H(m); m > n excitation, m < n de-excitation", error)
  
    call h5ltset_attribute_string_f(gid, "ionization", "description", &
         "n resolved ionization reaction rates: ionization(n,energy,temp,bt_amu)", error)
    call h5ltset_attribute_string_f(gid, "ionization", "units", "cm^3/s", error)
    call h5ltset_attribute_string_f(gid, "ionization", "reaction", &
         "A(q+) + H(n) -> A(q+) + H(+) + e-", error)
  
    call h5gclose_f(gid, error)
  
    deallocate(ebarr, tarr, excit, ioniz)

end subroutine write_bt_H_Aq

end module atomic_tables

program generate_tables
    use atomic_tables
    use H5LT
    use HDF5
    use hdf5_extra

#ifdef _OMP
  use omp_lib
#endif

    character(len=200) :: namelist_file
    character(len=3) :: arg

    character(len=200) :: tables_file = ''
    integer :: n_max, m_max

    integer, dimension(8) :: time_arr, time_start, time_end
    integer :: hour, minu, sec
    integer :: argc, max_threads, nthreads
    integer(HID_T) :: fid, gid
    integer :: error
    logical :: exis

    NAMELIST /general_settings/ n_max, m_max, tables_file

    argc = command_argument_count()
    if(argc.ge.1) then
        call get_command_argument(1, namelist_file)
    endif

    inquire(file=namelist_file,exist=exis)
    if(.not.exis) then
        write(*,'(a,a)') 'Input file does not exist: ',trim(namelist_file)
        stop
    else
        open(13,file=namelist_file)
        read(13,NML=general_settings)
        close(13)
        n_max = min(n_max,15)
        m_max = min(m_max,15)
    endif

    write(*,'(a)') "---- General settings ----"
    write(*,'(T2,"n_max = ",i2)') n_max
    write(*,'(T2,"m_max = ",i2)') m_max
    write(*,'(T2,"Tables File: ",a)') trim(tables_file)
    write(*,*) ''

#ifdef _OMP
    max_threads = OMP_get_num_procs()
    if(argc.ge.2) then
        call get_command_argument(2,arg)
        read(arg,'(i3)') nthreads
    else
        nthreads = max_threads
    endif
    max_threads = min(nthreads,max_threads)
    write(*,'(a)') "---- OpenMP settings ----"
    write(*,'(T2,"Number of threads: ",i2)') max_threads
    write(*,*) ''
    call OMP_set_num_threads(max_threads)
#endif

    ! Check if compression is possible
    call check_compression_availability()

    ! measure time
    call date_and_time (values=time_start)

    call h5open_f(error)
  
    call h5fcreate_f(tables_file, H5F_ACC_TRUNC_F, fid, error)
  
    call h5gcreate_f(fid, "cross", gid, error)
  
    call date_and_time(values=time_arr)
    write(*,"(A,I2,A,I2.2,A,I2.2)") 'Cross Sections:   ',time_arr(5),':',time_arr(6),':',time_arr(7)
    call write_bb_H_H(gid, namelist_file, n_max, m_max)
    call write_bb_H_e(gid, namelist_file, n_max, m_max)
    call write_bb_H_Aq(gid, namelist_file, n_max, m_max)
  
    call h5gclose_f(gid, error)
  
    call h5gcreate_f(fid, "rates", gid, error)
    
    call date_and_time(values=time_arr)
    write(*,"(A,I2,A,I2.2,A,I2.2)") 'Reaction Rates:   ',time_arr(5),':',time_arr(6),':',time_arr(7)
    call write_bt_H_H(gid, namelist_file, n_max, m_max)
    call write_bt_H_e(gid, namelist_file, n_max, m_max)
    call write_bt_H_Aq(gid, namelist_file, n_max, m_max)
    call write_einstein(gid, n_max, m_max)
  
    call h5gclose_f(gid, error)
  
    call h5ltset_attribute_string_f(fid, "/", "description", &
         "Atomic Cross Sections and Rates", error)
    call h5ltset_attribute_string_f(fid, "cross", "description", &
         "Atomic Cross Sections", error)
    call h5ltset_attribute_string_f(fid, "rates", "description", &
         "Atomic Reaction Rates", error)
  
    call h5fclose_f(fid, error)
    call h5close_f(error)
  
    write(*,'(a)') "Atomic tables written to "//trim(tables_file)
    write(*,*) ''

    call date_and_time (values=time_arr)
    write(*,'(A,I2,":",I2.2,":",I2.2)') 'END: hour, minute, second: ',time_arr(5), time_arr(6),time_arr(7)

    call date_and_time (values=time_end)
    hour = time_end(5) - time_start(5)
    minu = time_end(6) - time_start(6)
    sec  = time_end(7) - time_start(7)
    if (minu.lt.0.) then
        minu = minu +60
        hour = hour -1
    endif
    if (sec.lt.0.) then
        sec  = sec +60
        minu = minu -1
    endif

    write(*,'(A,18X,I2,":",I2.2,":",I2.2)') 'duration:',hour,minu,sec
end program
