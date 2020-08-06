program SynSam_prog

use, intrinsic :: iso_fortran_env, only: output_unit
use System_mod, only: CmdArg_type
use Constants_mod, only: IK, RK
use ParaPost_mod, only: ParaPost_type, NSAMPLE_DEFAULT
use NicoleRadio_mod, only: NicoleRadio_type, NUM_NICOLE_RADIO_DARK, NUM_NICOLE_RADIO_BRIGHT
use SynRed_mod, only: SynRed_type
use SynSam_mod, only: genSynSam

implicit none

integer(IK) , parameter         :: N_SFR_MODEL = 2
character(:), allocatable       :: StarFormationModel(:)
character(:), allocatable       :: SynRedFilePath(:), ParaPostFilePath(:)
character(:), allocatable       :: nicoleRadioFilePathDark, nicoleRadioFilePathBright
character(:), allocatable       :: outPathBase, synSamOutPath, synSamOutPathDark, synSamOutPathBright   !, detector
type(SynRed_type)               :: SynRed
type(NicoleRadio_type)          :: NicoleRadioDark, NicoleRadioBright
type(ParaPost_type)             :: ParaPost
type(CmdArg_type)               :: CmdArg
integer(IK)                     :: inFileUnit, imodel, numSynSam
logical                         :: nicoleEisoEnabled, nicoleZoneEnabled
real(RK)                        :: Log10LisoRange(2)

namelist /InputData/ StarFormationModel, SynRedFilePath, Log10LisoRange !, detector
namelist /InputData/ ParaPostFilePath, outPathBase, numSynSam, nicoleRadioFilePathDark, nicoleRadioFilePathBright, nicoleEisoEnabled, nicoleZoneEnabled

nicoleEisoEnabled = .true.
nicoleZoneEnabled = .true.
Log10LisoRange = [52._RK, 60._RK]

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%% query input data file name from the command line
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call CmdArg%query()
if (CmdArg%count/=1) then
    write(output_unit,"(*(g0))")
    write(output_unit,"(*(g0))") "FATAL: Invalid number of command-line arguments: ", CmdArg%count
    write(output_unit,"(*(g0))") "       Use the following example syntax to invoke the program: "
    write(output_unit,"(*(g0))") "       a.exe <input file path: ../in/SynSam.nml>"
end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%% read simulation input data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

open( newunit = inFileUnit, file = CmdArg%Arg(1)%record, status="old" )
    !allocate(character(10)      :: detector)
    allocate(character(1000)    :: outPathBase)
    allocate(character(1000)    :: nicoleRadioFilePathDark, nicoleRadioFilePathBright)
    allocate(character(3)       :: StarFormationModel(N_SFR_MODEL))
    allocate( character(2047)   :: SynRedFilePath(N_SFR_MODEL), ParaPostFilePath(N_SFR_MODEL) )
    StarFormationModel = "NUL"
    read(inFileUnit,nml=InputData)
    outPathBase = trim(adjustl(outPathBase))
    !detector = trim(adjustl(detector))
close(inFileUnit)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%% Read the Nicole Lloyd's Radio LGRB data
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nicoleEisoEnabled) then

    NicoleRadioDark     = NicoleRadio_type  ( nsample   = NUM_NICOLE_RADIO_DARK &
                                            , filePath  = trim(adjustl(nicoleRadioFilePathDark)) )

    NicoleRadioBright   = NicoleRadio_type  ( nsample   = NUM_NICOLE_RADIO_Bright &
                                            , filePath  = trim(adjustl(nicoleRadioFilePathBright)) )

end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%% simulate for different GRB formation models
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do imodel = 1, N_SFR_MODEL

    write(output_unit,"(*(g0))"); write(output_unit,"(*(g0))")
    write(output_unit,"(*(g0))") "Generating GRB redshift sample based on the SFR model of ", StarFormationModel(imodel)
    write(output_unit,"(*(g0))"); write(output_unit,"(*(g0))")

    ! estimate redshifts

    call execute_command_line("mkdir "//trim(adjustl(outPathBase)))
    synSamOutPath         = outPathBase // "SynSam" // trim(adjustl(StarFormationModel(imodel))) // ".csv"
    synSamOutPathDark     = outPathBase // "SynSam" // trim(adjustl(StarFormationModel(imodel))) // "Dark.csv"
    synSamOutPathBright   = outPathBase // "SynSam" // trim(adjustl(StarFormationModel(imodel))) // "Bright.csv"

    if ( &
        trim(adjustl(StarFormationModel(imodel))) == "H06" &
        .or. &
        !trim(adjustl(StarFormationModel(imodel))) == "L08" &
        .or. &
        trim(adjustl(StarFormationModel(imodel))) == "B10" &
        ) then

        ! Read the LGRB world model parameters

        ParaPost = ParaPost_type( nsample = NSAMPLE_DEFAULT &
                                , paraPostFilePath = trim(adjustl(ParaPostFilePath(imodel))) &
                                )

        if (nicoleEisoEnabled .and. nicoleZoneEnabled) then

            ! generate Radio Dark sample

            call genSynSam  ( ParaPost          = ParaPost              &
                            , outFilePath       = synSamOutPathDark     &
                            , numSynSam         = numSynSam             &
                            , NicoleRadio       = NicoleRadioDark       &
                            )

            ! generate Radio Bright sample

            call genSynSam  ( ParaPost          = ParaPost              &
                            , outFilePath       = synSamOutPathBright   &
                            , numSynSam         = numSynSam             &
                            , NicoleRadio       = NicoleRadioBright     &
                            )

        elseif (nicoleEisoEnabled) then

            ! Read the synthetic redshifts

            SynRed = SynRed_type ( redshiftChainFilePath = trim(adjustl(SynRedFilePath(imodel))) )

            ! generate Radio Dark sample

            call genSynSam  ( ParaPost          = ParaPost              &
                            , outFilePath       = synSamOutPathDark     &
                            , numSynSam         = numSynSam             &
                            , NicoleRadio       = NicoleRadioDark       &
                            , SynRed            = SynRed                &
                            )

            ! generate Radio Bright sample

            call genSynSam  ( ParaPost          = ParaPost              &
                            , outFilePath       = synSamOutPathBright   &
                            , numSynSam         = numSynSam             &
                            , NicoleRadio       = NicoleRadioBright     &
                            , SynRed            = SynRed                &
                            )

        else

            ! generate Radio sample irrespective of the Nicole Radio samples

            call genSynSam  ( ParaPost          = ParaPost              &
                            , outFilePath       = synSamOutPath         &
                            , numSynSam         = numSynSam             &
                            , Log10LisoRange    = Log10LisoRange        &
                           !, NicoleRadio       = NicoleRadioBright     &
                            , SynRed            = SynRed                &
                            )

        end if

    elseif (StarFormationModel(imodel) /= "NUL") then
    
        write(output_unit,"(*(g0))")
        write(output_unit,"(*(g0))") "FATAL: Invalid star formation model on input: ", trim(adjustl(StarFormationModel(imodel)))
        write(output_unit,"(*(g0))") "       Use the one of following supported models: "
       !write(output_unit,"(*(g0))") "       H06   L08   B10"
        write(output_unit,"(*(g0))") "       H06   B10"
        write(output_unit,"(*(g0))")

    end if

end do

end program SynSam_prog
