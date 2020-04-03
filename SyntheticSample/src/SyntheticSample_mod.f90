module SyntheticSample_mod

    use Constants_mod, only: IK, RK

    implicit none

    character(*), parameter :: MODULE_NAME = "@SyntheticSample_mod"

    ! the exponent of zplus1 in time-dilation translation of T90 to T90z
#ifdef KFAC_ONETHIRD_ENABLED
    real(RK)   , parameter :: TIME_DILATION_EXPO = 0.666666666666667_RK
#endif

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine generateSyntheticSample  ( SynRed &
                                        , ParaPost &
                                        , outFilePath &
                                        , nsim &
                                        , NicoleRadio &
                                        , Log10LisoRange &
                                       !, detector &
                                       !, nicoleEisoUsed &
                                        )

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK, RK, SPR, CARRIAGE_RETURN, CLOCK_TICK, LN10, SQRT2, LOG_ERG2KEV
        use SyntheticRedshift_mod, only: SyntheticRedshift_type
        use NicoleRadio_mod, only: NicoleRadio_type !, NUM_NICOLE_RADIO_DARK, NUM_NICOLE_RADIO_BRIGHT
        use ParaPost_mod, only: ParaPost_type, Thresh_type
        use CorrCoef_mod, only: CorrCoefSpearman_type
        use Statistics_mod, only: getRandInt, getRandMVN, getMean, getVariance
        use BandModel_mod, only: getPhotonFluenceFromEnergyFluence
        use TranGaus_mod, only: getTranGaus
        use Timer_mod, only: Timer_type
        use Batse_mod, only: getLogPF53
        use Err_mod, only: Err_type

        implicit none

        type(SyntheticRedshift_type), intent(in)        :: SynRed
        type(ParaPost_type), intent(inout)              :: ParaPost
       !logical, intent(in)                             :: nicoleEisoUsed
        character(*), intent(in)                        :: outFilePath  !, detector
        integer(IK) , intent(in)                        :: nsim
        type(NicoleRadio_type), intent(in), optional    :: NicoleRadio
        real(RK), intent(in), optional                  :: Log10LisoRange(2)

        type(Err_type)                                  :: Err
        type(CorrCoefSpearman_type)                     :: CorZoneEiso, CorZoneDurz, CorEisoDurz
        logical                                         :: nicoleEisoNotUsed, nicoleEisoUsed
        real(RK)                                        :: Log10LisoLimit(2)

        real(RK), parameter     :: LOGEISO_LOWER_LIM = 52._RK * LN10, LOGEISO_UPPER_LIM = 60._RK * LN10
        real(RK)                :: logPF53, normedLogPF53, unifrnd
        real(RK), allocatable   :: GrbIntDep(:), NicoleSynSam(:,:)
        integer(IK)             :: isim, isample, iParaPost, iRedshift, fileUnit, itick, inicole, ivar, sampleSize!, fileUnitSample
        integer(IK)             :: numSumStat   ! number of event variables whose statistics will be computed
        real(RK), allocatable   :: Vector(:), SampleMean(:), SampleVariance(:)
        type(Timer_type)        :: Timer

#ifdef SWIFT_ENABLED

        ! Swift detection parameters
        real(RK), parameter :: SWIFT_INV_STD_LOG_THRESH_SQRT2 = 1._RK / (0.1_RK * LN10 * SQRT2)
        real(RK), parameter :: SWIFT_AVG_LOG_THRESH_B07 = log(3._RK)
        real(RK), parameter :: SWIFT_AVG_LOG_THRESH_B10 = log(0.15_RK)    ! = log( butlerCeffmin / sqrt(avgSwiftPartialCoding*medSwiftT90OverTr45) )
        real(RK)            :: invSqrtT90Obs
        real(RK)            :: probDetectionSwiftB07, normedEffectivePhotonFluenceB07, photonFluence
        real(RK)            :: probDetectionSwiftB10, normedEffectivePhotonFluenceB10, photonFluence15150
        !real(RK), parameter     :: LOG_T90_OVER_TR45 = 1.728056032899616_RK
        !real(RK), parameter     :: LOG_T90_OVER_TR45 = 1.728056032899616_RK
        !SwiftLogThresh%invStdSqrt2 = 1._RK / (0.1_RK*LN10*SQRT2)  ! taken from Butler 2010 Appendix: sigmaLogThresh = 0.1
        !SwiftLogThresh%avg = log(0.257039578276886_RK)  ! taken from Butler 2010 Appendix: avgLogThresh = log(0.257039578276886)
#endif

        if (present(Log10LisoRange)) then
            Log10LisoLimit = Log10LisoRange * LN10
        else
            Log10LisoLimit(1) = LOGEISO_LOWER_LIM
            Log10LisoLimit(2) = LOGEISO_UPPER_LIM
        end if
        write(*,"(*(g0,:,' '))")
        write(*,"(*(g0,:,' '))") "Log10LisoLimit:", Log10LisoLimit
        write(*,"(*(g0,:,' '))")

        nicoleEisoUsed = present(NicoleRadio)
        nicoleEisoNotUsed = .not. nicoleEisoUsed
        sampleSize = 60_IK  ! average of Dark and Bright sample sizes
        if (nicoleEisoUsed) sampleSize = NicoleRadio%count

        itick = 0
        call Timer%tic()

        open(newunit=fileUnit,file=outFilePath,status="replace")
!        write(fileUnit,"(*(g0,:,','))"  ) &!"meanLogLiso",          &
!                                        !, "meanLogEpkz"             &
!                                        !, "meanLogDurz"             &
!                                           "meanLogDurz"             &
!                                        !, "meanLogEiso"             &
!                                        !, "meanLogPbol"             &
!                                        !, "meanLogEpko"             &
!                                        !, "meanLogDuro"             &
!                                        !, "meanLogSbol"             &
!                                         , "meanLogZone"             &
!                                         , "meanLogProbDetection"    &
!                                        !, "varLogLiso"              &
!                                        !, "varLogEpkz"              &
!                                         , "varLogDurz"              &
!                                        !, "varLogEiso"              &
!                                        !, "varLogPbol"              &
!                                        !, "varLogEpko"              &
!                                        !, "varLogDuro"              &
!                                        !, "varLogSbol"              &
!                                         , "varLogZone"              &
!                                         , "varLogProbDetection"     &
!                                         , "corZoneEiso"             &
!                                         , "corZoneDurz"             &
!                                         , "corEisoDurz"

        write(fileUnit,"(*(g0,:,','))"  )  "meanLogLiso"            &
                                         , "meanLogEpkz"            &
                                         , "meanLogDurz"            &
                                         , "meanLogEiso"            &
                                         , "meanLogPbol"            &
                                         , "meanLogEpko"            &
                                         , "meanLogDuro"            &
                                         , "meanLogSbol"            &
                                         , "meanLogZone"            &
                                         , "meanProbDetection"      &
                                         , "varLogLiso"             &
                                         , "varLogEpkz"             &
                                         , "varLogDurz"             &
                                         , "varLogEiso"             &
                                         , "varLogPbol"             &
                                         , "varLogEpko"             &
                                         , "varLogDuro"             &
                                         , "varLogSbol"             &
                                         , "varLogZone"             &
                                         , "varProbDetection"       &
                                         , "corZoneEiso"            &
                                         , "corZoneDurz"            &
                                         , "corEisoDurz"

        !open(newunit=fileUnitSample,file="SynSam.csv",status="replace")
        !write(fileUnitSample,"(*(g0))") "LogLiso,&
        !                          &LogEpkz,&
        !                          &LogDurz,&
        !                          &LogEiso,&
        !                          &LogPbol,&
        !                          &LogEpko,&
        !                          &LogDuro,&
        !                          &LogSbol,&
        !                          &LogZone,&
        !                          &LogProbDetection"

        ! assign the Eiso values to the synthetic Nicole sample size
        numSumStat = ParaPost%ndim * 2 + 2

        allocate(Vector(numSumStat))
        allocate(SampleMean(numSumStat))
        allocate(SampleVariance(numSumStat))

        ! first determine the conditional mean of [LogLiso,LogEpkz,LogDurz] given LogEiso
        if (allocated(NicoleSynSam)) deallocate(NicoleSynSam)
        allocate(NicoleSynSam(numSumStat,sampleSize))
        do iParaPost = 1, ParaPost%count
            if (allocated(ParaPost%Sample(iParaPost)%ConMean)) deallocate(ParaPost%Sample(iParaPost)%ConMean)
            allocate(ParaPost%Sample(iParaPost)%ConMean(ParaPost%numDepVar,sampleSize))
            if (nicoleEisoUsed) then
                NicoleSynSam(4,:) = NicoleRadio%LogEiso
                do isample = 1, sampleSize
                    ParaPost%Sample(iParaPost)%ConMean(:,isample)   = ParaPost%Sample(iParaPost)%Avg(1:ParaPost%numDepVar) &
                                                                    + ParaPost%Sample(iParaPost)%RegresCoefMat(1:ParaPost%numDepVar,1) &
                                                                    * ( NicoleSynSam(4,isample) - ParaPost%Sample(iParaPost)%Avg(4) )
                end do
            end if
        end do

        ! simulate samples
        allocate(GrbIntDep(ParaPost%numDepVar))
        do isim = 1, nsim

            iParaPost = getRandInt(1,ParaPost%count)

!write(*,"(*(g0))")
!write(*,"(*(g0))") "Avg: ", ParaPost%Sample(iParaPost)%Avg
!write(*,"(*(g0))") "Low: ", ParaPost%Sample(iParaPost)%CholFacLower
!write(*,"(*(g0))") "Dia: ", ParaPost%Sample(iParaPost)%CholFacDiag
!write(*,"(*(g0))")

            loopNicole: do inicole = 1, sampleSize

                if (nicoleEisoNotUsed) then
                    do isample = 1, sampleSize
                        NicoleSynSam(4,isample) = getTranGaus   ( Log10LisoLimit(1) &
                                                                , Log10LisoLimit(2) &
                                                                , ParaPost%Sample(iParaPost)%Avg(4) &
                                                                , ParaPost%Sample(iParaPost)%Std(4) &
                                                                )
                        ParaPost%Sample(iParaPost)%ConMean(:,isample)   = ParaPost%Sample(iParaPost)%Avg(1:ParaPost%numDepVar) &
                                                                        + ParaPost%Sample(iParaPost)%RegresCoefMat(1:ParaPost%numDepVar,1) &
                                                                        * ( NicoleSynSam(4,isample) - ParaPost%Sample(iParaPost)%Avg(4) )
                    end do
                end if

                ! detect an event
                loopDetection: do
                    GrbIntDep = getRandMVN  ( nd = ParaPost%numDepVar &
                                            , MeanVec = ParaPost%Sample(iParaPost)%ConMean(:,inicole) &
                                            , CholeskyLower = ParaPost%Sample(iParaPost)%CholFacLower &
                                            , Diagonal = ParaPost%Sample(iParaPost)%CholFacDiag &
                                            )
                    iRedshift = getRandInt(1,SynRed%count)

                    ! rest-frame properties
                    NicoleSynSam(1:3,inicole) = GrbIntDep

                    ! observer-frame properties
                    NicoleSynSam(5,inicole) = NicoleSynSam(1,inicole) - SynRed%Sample(iRedshift)%logLisoLogPbolDiff ! LogLiso -> LogPbol
                    NicoleSynSam(6,inicole) = NicoleSynSam(2,inicole) - SynRed%Sample(iRedshift)%logzplus1          ! LogEpkz -> LogEpk

                    NicoleSynSam(8,inicole) = NicoleSynSam(4,inicole) - SynRed%Sample(iRedshift)%logEisoLogSbolDiff ! LogEiso -> LogSbol
#ifdef KFAC_ONETHIRD_ENABLED
                    NicoleSynSam(7,inicole) = NicoleSynSam(3,inicole) + SynRed%Sample(iRedshift)%logzplus1 * TIME_DILATION_EXPO ! LogDurz -> LogDur
#else
                    NicoleSynSam(7,inicole) = NicoleSynSam(3,inicole) + SynRed%Sample(iRedshift)%logzplus1  ! LogDurz -> LogDur
#endif


                    ! compute the probability of detection
#ifdef SWIFT_ENABLED


                    invSqrtT90Obs = 1._RK / sqrt(exp(NicoleSynSam(7,inicole)))
                    ! get the bolometric photon fluence


                    ! According to Butler et al 2007/2010, detection quantity for Swift is Nbol / sqrt(T90)
                    call getPhotonFluenceFromEnergyFluence  ( energyFluence = exp(NicoleSynSam(8,inicole)+LOG_ERG2KEV)  &
                                                            , lowerLim      = 1.e-1_RK                                  &
                                                            , upperLim      = 2.e+4_RK                                  &
                                                            , epk           = exp(NicoleSynSam(6,inicole))              &
                                                            , alpha         = -1.1_RK                                   &
                                                            , beta          = -2.3_RK                                   &
                                                            , tolerance     = 1.e-5_RK                                  &
                                                            , photonFluence = photonFluence                             &
                                                            , Err           = Err                                       &
                                                            )
                    if (Err%occurred) then
                        write(*,"(*(g0,:,' '))") Err%msg
                        write(*,"(*(g0,:,' '))") "Err%stat: ", Err%stat
                        error stop
                    end if
                    ! get probability of detection by Swift according to Butler 2007
                    normedEffectivePhotonFluenceB07 = SWIFT_INV_STD_LOG_THRESH_SQRT2 * ( log(photonFluence*invSqrtT90Obs) - SWIFT_AVG_LOG_THRESH_B07 )     ! Nbol / sqrt(T90)
                    probDetectionSwiftB07 = 0.5_RK + 0.5_RK * erf( real( normedEffectivePhotonFluenceB07 , kind=SPR ) )
                    NicoleSynSam(10,inicole) = 0.5_RK + 0.5_RK * erf( real( probDetectionSwiftB07 , kind=SPR ) ) ! probDetection


                    !! get the bolometric photon fluence
                    !call getPhotonFluenceFromEnergyFluence  ( energyFluence = exp(NicoleSynSam(8,inicole)+LOG_ERG2KEV)  &
                    !                                        , lowerLim      = 1.e-1_RK                                  &
                    !                                        , upperLim      = 2.e+4_RK                                  &
                    !                                        , epk           = exp(NicoleSynSam(6,inicole))              &
                    !                                        , alpha         = -1.1_RK                                   &
                    !                                        , beta          = -2.3_RK                                   &
                    !                                        , tolerance     = 1.e-5_RK                                  &
                    !                                        , photonFluence = photonFluence15150                        &
                    !                                        , Err           = Err                                       &
                    !                                        , lowerLimNew   = 15._RK                                    &
                    !                                        , upperLimNew   = 150._RK                                   &
                    !                                        )
                    !if (Err%occurred) then
                    !    write(*,"(*(g0,:,' '))") Err%msg
                    !    write(*,"(*(g0,:,' '))") "Err%stat: ", Err%stat
                    !    error stop
                    !end if
                    !! get probability of detection by Swift according to Butler 2010
                    !normedEffectivePhotonFluenceB10 = SWIFT_INV_STD_LOG_THRESH_SQRT2 * ( log(photonFluence15150*invSqrtT90Obs) - SWIFT_AVG_LOG_THRESH_B10 )     ! Nbol / sqrt(T90)
                    !probDetectionSwiftB10 = 0.5_RK + 0.5_RK * erf( real( normedEffectivePhotonFluenceB10 , kind=SPR ) )
                    !NicoleSynSam(10,inicole) = 0.5_RK + 0.5_RK * erf( real( probDetectionSwiftB10 , kind=SPR ) ) ! probDetection

#else

                    ! get the photon flux in 50-300 BATSE detection energy range
                    logPF53 = getLogPF53(logEpk=NicoleSynSam(6,inicole),logPbol=NicoleSynSam(5,inicole))
                    normedLogPF53 = ParaPost%Sample(iParaPost)%Thresh%invStdSqrt2 * ( logPF53 - ParaPost%Sample(iParaPost)%Thresh%avg )
                    NicoleSynSam(10,inicole) = 0.5_RK + 0.5_RK * erf( real( normedLogPF53 , kind=SPR ) ) ! probDetection

#endif


                    ! check if it is detectable
                    call random_number(unifrnd)
                    if (NicoleSynSam(10,inicole)>unifrnd) exit loopDetection
                end do loopDetection

                ! convert detection prob to log scale
                !NicoleSynSam(10,inicole) = log( NicoleSynSam(10,inicole) )
                NicoleSynSam(9,inicole) = SynRed%Sample(iRedshift)%logzplus1    ! log(redshift + 1) or LogZone

                !write(fileUnitSample,"(*(g0.8,:,','))") NicoleSynSam(:,inicole)

            end do loopNicole

            ! compute the correlations of the detected sample
            call CorZoneEiso%get( ndata             = sampleSize                    &
                                , Data1             = NicoleSynSam(9,1:sampleSize)  &
                                , Data2             = NicoleSynSam(4,1:sampleSize)  &
                                , rho               = CorZoneEiso%rho               &
                                , rhoProb           = CorZoneEiso%rhoProb           &
                                , dStarStar         = CorZoneEiso%dStarStar         &
                                , dStarStarSignif   = CorZoneEiso%dStarStarSignif   &
                                , dStarStarProb     = CorZoneEiso%dStarStarProb     &
                                )

            ! compute the correlations of the detected sample
            call CorZoneDurz%get( ndata             = sampleSize                    &
                                , Data1             = NicoleSynSam(9,1:sampleSize)  &
                                , Data2             = NicoleSynSam(3,1:sampleSize)  &
                                , rho               = CorZoneDurz%rho               &
                                , rhoProb           = CorZoneDurz%rhoProb           &
                                , dStarStar         = CorZoneDurz%dStarStar         &
                                , dStarStarSignif   = CorZoneDurz%dStarStarSignif   &
                                , dStarStarProb     = CorZoneDurz%dStarStarProb     &
                                )

            ! compute the correlations of the detected sample
            call CorEisoDurz%get( ndata             = sampleSize                    &
                                , Data1             = NicoleSynSam(4,1:sampleSize)  &
                                , Data2             = NicoleSynSam(3,1:sampleSize)  &
                                , rho               = CorEisoDurz%rho               &
                                , rhoProb           = CorEisoDurz%rhoProb           &
                                , dStarStar         = CorEisoDurz%dStarStar         &
                                , dStarStarSignif   = CorEisoDurz%dStarStarSignif   &
                                , dStarStarProb     = CorEisoDurz%dStarStarProb     &
                                )


            ! compute the statistics of the detected sample
            SampleMean = getMean( nd = numSumStat &
                                , np = sampleSize &
                                , Point = NicoleSynSam &
                                )
            do ivar = 1, numSumStat
                Vector = NicoleSynSam(ivar,1:sampleSize)
                SampleVariance(ivar) = getVariance( np = sampleSize, mean = SampleMean(ivar), Point = Vector )
            end do
            !write(fileUnit,"(*(g0.8,:,','))") SampleMean(3),SampleMean(9),SampleMean(10) &
            !                                , sqrt(SampleVariance(3)),sqrt(SampleVariance(9)),sqrt(SampleVariance(10)) &
            !                                , CorZoneEiso%rho, CorZoneDurz%rho, CorEisoDurz%rho
            write(fileUnit,"(*(g0.8,:,','))") (SampleMean(ivar),ivar=1,numSumStat) &
                                            , (SampleVariance(ivar),ivar=1,numSumStat) &
                                            , CorZoneEiso%rho, CorZoneDurz%rho, CorEisoDurz%rho

            if (mod(nsim,10000)==0) then
                call Timer%toc()
                itick = itick + 1
                write(output_unit,"(*(' ',g0))", advance="no"   ) CARRIAGE_RETURN, CLOCK_TICK(itick) &
                                                                , isim, " out of  ", nsim &
                                                                , "GRB samples generated in" &
                                                                , Timer%Time%total, "seconds."
                flush(output_unit)
                if (itick==4) itick = 0
            end if

        end do
        close(fileUnit)
        !close(fileUnitSample)

        write(output_unit,"(*(g0))")
        write(output_unit,"(*(g0))")

        !call Timer%toc()
        !write(output_unit,"(*(g0))")
        !write(output_unit,"(*(g0))")
        !write(output_unit,"(*(g0))") "Total Time: ", Timer%Time%total, " seconds."
        !write(output_unit,"(*(g0))")

    end subroutine generateSyntheticSample

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SyntheticSample_mod