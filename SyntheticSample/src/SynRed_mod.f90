module SynRed_mod

    use Err_mod, only: Err_type
    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@SynRed_mod"

    type :: RedshiftSample_type
        real(RK) :: z, logZone, logLisoLogPbolDiff, logEisoLogSbolDiff
    end type RedshiftSample_type

    type :: SynRed_type
        integer(IK) :: count
        character(:), allocatable :: filePath
        type(Err_type) :: Err
        type(RedshiftSample_type), allocatable :: Sample(:)
    end type SynRed_type

    interface SynRed_type
        module procedure :: constructSynRed
    end interface

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! All parameters are assumed to be in log Neper (not log10) wherever needed.
    function constructSynRed(redshiftChainFilePath) result(SynRed)

        !use, intrinsic :: iso_fortran_env, only: output_unit
        use Cosmology_mod, only: LOGMPC2CMSQ4PI, getLogLumDisWicMpc
        use Constants_mod, only: IK, RK
        implicit none

        character(*), intent(in)        :: redshiftChainFilePath
        type(SynRed_type)               :: SynRed

        character(*), parameter         :: PROCEDURE_NAME = MODULE_NAME//"@getModelIntegral()"
        integer(IK)                     :: isample, fileUnit
        real(RK)                        :: zone, twiceLogLumDisMpc, Dummy(6)

        SynRed%filePath = trim(adjustl(redshiftChainFilePath))

        SynRed%count = 277000
        if (allocated(SynRed%Sample)) deallocate(SynRed%Sample); allocate(SynRed%Sample(SynRed%count))

        open(newunit = fileUnit, file = SynRed%filePath, status = "old")
        read(fileUnit,*)
        do isample = 1,SynRed%count
            read(fileUnit,*) Dummy(1:6), SynRed%Sample(isample)%z
            zone = SynRed%Sample(isample)%z + 1._RK
            twiceLogLumDisMpc = 2 * getLogLumDisWicMpc(zone)
            SynRed%Sample(isample)%logZone = log(zone)
            SynRed%Sample(isample)%logLisoLogPbolDiff  = LOGMPC2CMSQ4PI + twiceLogLumDisMpc
            SynRed%Sample(isample)%logEisoLogSbolDiff  = LOGMPC2CMSQ4PI + twiceLogLumDisMpc - SynRed%Sample(isample)%logZone
        end do

    end function constructSynRed

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SynRed_mod