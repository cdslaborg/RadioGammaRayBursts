module NicoleRadio_mod

    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@NicoleRadio_mod"

    integer(IK), parameter  :: NUM_NICOLE_RADIO_DARK = 41, NUM_NICOLE_RADIO_BRIGHT = 78

#ifdef kfacOneThird
    real(RK)   , parameter :: TIME_DILATION_EXPO = 0.666666666666667_RK
#else
#error "kfac must be oneThird in NicholeRadio_mod.f90. Otherwise, remove this error command."
#endif

    type ThetaJ_type
        real(RK), allocatable :: Maria(:), JetBreak(:)
    end type ThetaJ_type

    type NicoleRadio_type
        real(RK)                        :: kfac
        integer(IK)                     :: count, fileUnit
        character(len=:), allocatable   :: filePath
        character(len=7), allocatable   :: GrbId(:)
        real(RK), allocatable           :: LogEiso(:), LogT90z(:), LogZone(:), Zone(:), LogLisoLogPbolDiff(:), LogEisoLogSbolDiff(:)
        type(ThetaJ_type)               :: ThetaJ
    end type NicoleRadio_type
    !type(NicoleRadio_type) :: RadioDark, RadioBright

    interface NicoleRadio_type
        module procedure :: constructNicoleRadio
    end interface NicoleRadio_type

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function constructNicoleRadio(nsample,filePath) result(NicoleRadio)

        use Cosmology_mod, only: LOGMPC2CMSQ4PI, getLogLumDisWicMpc
        use Constants_mod, only: IK, RK, LN10
        implicit none
        integer(IK), intent(in)     :: nsample
        character(*), intent(in)    :: filePath
        type(NicoleRadio_type)      :: NicoleRadio
        integer                     :: i
        real(RK)                    :: twiceLogLumDisMpc
        NicoleRadio%filePath = filePath
        NicoleRadio%count = nsample
        allocate( NicoleRadio%Zone(nsample) &
                , NicoleRadio%GrbId(nsample) &
                , NicoleRadio%LogEiso(nsample) &
                , NicoleRadio%LogT90z(nsample) &
                , NicoleRadio%LogZone(nsample) &
                , NicoleRadio%LogLisoLogPbolDiff(nsample) &
                , NicoleRadio%LogEisoLogSbolDiff(nsample) &
                , NicoleRadio%ThetaJ%JetBreak(nsample) &
                , NicoleRadio%ThetaJ%Maria(nsample) &
                )
        open( newunit   = NicoleRadio%fileUnit &
            , file      = NicoleRadio%filePath &
            , status    = "old" &
            )
        read( NicoleRadio%fileUnit , * )
        do i = 1, nsample
            read(NicoleRadio%fileUnit,* ) NicoleRadio%GrbId(i) &
                                        , NicoleRadio%LogEiso(i) &
                                        , NicoleRadio%LogT90z(i) &
                                        , NicoleRadio%Zone(i) &
                                        , NicoleRadio%ThetaJ%Maria(i) &
                                        , NicoleRadio%ThetaJ%JetBreak(i)

            ! convert to proper units and frames of reference

            NicoleRadio%LogEiso(i) = log( NicoleRadio%LogEiso(i) ) + 52._RK * log(1.e1_RK)
            NicoleRadio%Zone(i) = 1._RK + NicoleRadio%Zone(i)
            NicoleRadio%LogZone(i) = log( 1._RK + NicoleRadio%Zone(i) )

            twiceLogLumDisMpc = 2 * getLogLumDisWicMpc(NicoleRadio%Zone(i))
            NicoleRadio%LogLisoLogPbolDiff(i) = LOGMPC2CMSQ4PI + twiceLogLumDisMpc
            NicoleRadio%LogEisoLogSbolDiff(i) = LOGMPC2CMSQ4PI + twiceLogLumDisMpc - NicoleRadio%LogZone(i)

#ifdef kfacOneThird
            NicoleRadio%kfac = TIME_DILATION_EXPO
            NicoleRadio%LogT90z(i) = log(NicoleRadio%LogT90z(i)) - NicoleRadio%LogZone(i) * NicoleRadio%kfac;
#else
            NicoleRadio%LogT90z(i) = log(NicoleRadio%LogT90z(i)) - NicoleRadio%LogZone(i);
#endif
        end do

!write(*,*) sum(NicoleRadio%LogEiso) / NicoleRadio%count
!read(*,*)

    end function constructNicoleRadio

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module NicoleRadio_mod