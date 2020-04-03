module ParaPost_mod

    use Constants_mod, only: IK, RK, PI

    implicit none

    character(*), parameter :: MODULE_NAME = "@ParaPost_mod"

    integer(IK) , parameter :: NSAMPLE_DEFAULT = 2500_IK
    integer(IK) , parameter :: NPAR = 16    ! number of world model's parameters
    integer(IK) , parameter :: NVAR = 4     ! number of GRB attributes used in the world model

    ! BATSE threshold parameters
    type :: Thresh_type
        real(RK) :: avg, invStdSqrt2
    end type Thresh_type

    ! Posterior sample
    type :: ParaPostSample_type
        ! order of variables: logLiso, logEpkz, logEiso, logT90z
        type(Thresh_type)       :: Thresh
        real(RK)                :: logFunc
        real(RK), allocatable   :: Avg(:), Std(:), CholFacLower(:,:), CholFacDiag(:)
        real(RK), allocatable   :: RegresCoefMat(:,:), SchurComplement(:,:)
        real(RK), allocatable   :: ConMean(:,:)   ! conditional mean (1:NVAR-1,1:ndata)
    end type ParaPostSample_type

    type :: ParaPost_type
        integer(IK) :: count, ndim, fileUnit
        integer(IK) :: numDepVar    ! number of dependent variables
        character(:), allocatable :: filePath
        type(ParaPostSample_type), allocatable :: Sample(:)
    end type ParaPost_type

    interface ParaPost_type
        module procedure :: constructParaPost
    end interface ParaPost_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! All parameters are assumed to be read in log Neper (not log10) from the input file, wherever needed.
    function constructParaPost(nsample,paraPostFilePath) result(ParaPost)

        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK, RK, SQRT2, SQRT2PI, CARRIAGE_RETURN, CLOCK_TICK
        use Matrix_mod, only: getCholeskyFactor, sortPosDefMat, getRegresCoef
        use Timer_mod, only: Timer_type
        implicit none

        integer(IK), intent(in)                 :: nsample
        character(*), intent(in)                :: paraPostFilePath
        type(ParaPostSample_type), allocatable  :: ParaPostSample(:)    ! This is the original data read from input file
        type(ParaPost_type)                     :: ParaPost             ! This is the ordered data: Eiso goes in place of Durz

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@getModelIntegral()"

        integer(IK)                             :: i, j, isample, itick, Indx(1), IndxMap(1)   !, ierr
        real(RK), parameter                     :: INV_SQRT2 = 1._RK / SQRT2
        type(Timer_type)                        :: Timer

        ! Note ParaPost object will contain the shuffled sample. ParaPostSample object will contain the original sample from the input file.
        ParaPost%ndim = NVAR
        ParaPost%count = nsample
        ParaPost%numDepVar = ParaPost%ndim - 1_IK
        ParaPost%filePath = trim(adjustl(paraPostFilePath))
        if (allocated(ParaPost%Sample)) deallocate(ParaPost%Sample)
        allocate(ParaPost%Sample(ParaPost%count))
        if (allocated(ParaPostSample)) deallocate(ParaPostSample)
        allocate(ParaPostSample(ParaPost%count))

        write(output_unit,"(*(g0))")

        itick = 0
        call Timer%tic()

        open(newunit=ParaPost%fileUnit,file=ParaPost%filePath,status="old")
        read(ParaPost%fileUnit,*)

        ! This is the map from Eiso column to Durz in the covariance matrix
        Indx = [3]
        IndxMap = [4]

        isample = 0_IK
        do isample = 1, ParaPost%count

            !*******************************************************************************************************************
            !****                                         begin reading sample data                                         ****
            !*******************************************************************************************************************

            ! now read ParaPost sample properties
            allocate(ParaPostSample(isample)%Avg(ParaPost%ndim))
            allocate(ParaPostSample(isample)%Std(ParaPost%ndim))
            allocate(ParaPostSample(isample)%CholFacLower(ParaPost%ndim,ParaPost%ndim))
           !allocate(ParaPostSample(isample)%CholFacDiag(ParaPost%ndim))

           !read(ParaPost%fileUnit,"(*(g0.8,:,','))")  ParaPostSample(isample)%logFunc &
            read(ParaPost%fileUnit,*)  ParaPostSample(isample)%logFunc &
                                    ,  ParaPostSample(isample)%Avg(1:ParaPost%ndim) &
                                    , (ParaPostSample(isample)%CholFacLower(i,i),i=1,ParaPost%ndim) &  ! log of standard deviations
                                    ,  ParaPostSample(isample)%CholFacLower(2:4,1) &                   ! Rho: LisoEpkz, LisoEiso, LisoT90z
                                    ,  ParaPostSample(isample)%CholFacLower(3:4,2) &                   ! Rho: EpkzEiso, EpkzT90z
                                    ,  ParaPostSample(isample)%CholFacLower(4:4,3) &                   ! Rho: EisoT90z
                                    ,  ParaPost%Sample(isample)%Thresh%avg &                           ! log of threshold mean
                                    ,  ParaPost%Sample(isample)%Thresh%invStdSqrt2                     ! log of threshold standard deviation

            ParaPostSample(isample)%CholFacLower = transpose(ParaPostSample(isample)%CholFacLower)

            !*******************************************************************************************************************
            !****                             begin transform of the sample to the original scales                          ****
            !*******************************************************************************************************************

            ! threshold parameter: convert standard deviation to invStdSqrt2
            ParaPost%Sample(isample)%Thresh%invStdSqrt2 = INV_SQRT2 / exp(ParaPostSample(isample)%Thresh%invStdSqrt2)

            ! compute the standard deviations and Covariance elements
            do i = 1, ParaPost%ndim
                ParaPostSample(isample)%Std(i) = exp(ParaPostSample(isample)%CholFacLower(i,i))
                ParaPostSample(isample)%CholFacLower(i,i) = ParaPostSample(isample)%Std(i)**2  ! This is now squared standard deviation
            end do

            ! compute the covariance matrix's upper-triangle elements (excluding CholFacDiag variance elements)
            do j = 2, ParaPost%ndim
                do i = 1,j-1

                    ! convert the lower-triangle elements to correlation coefficients
                    ParaPostSample(isample)%CholFacLower(i,j) = tanh(ParaPostSample(isample)%CholFacLower(i,j))

                    ! compute the upper-triangle elements of the covariance matrix
                    ParaPostSample(isample)%CholFacLower(i,j)   = ParaPostSample(isample)%CholFacLower(i,j) &
                                                                * ParaPostSample(isample)%Std(i) &
                                                                * ParaPostSample(isample)%Std(j)

                end do
            end do

            !*******************************************************************************************************************
            !****                                 begin swapping the position of Eiso and Durz                              ****
            !*******************************************************************************************************************

            ! allocate ParaPost sample properties
            allocate(ParaPost%Sample(isample)%Avg(ParaPost%ndim))
            allocate(ParaPost%Sample(isample)%Std(ParaPost%ndim))
            allocate(ParaPost%Sample(isample)%CholFacLower(ParaPost%ndim,ParaPost%ndim))

            ! now reshuffle the variables so that Eiso in column 3 goes in place of Durz in column 4
            ParaPost%Sample(isample)%Avg(1:2)   = ParaPostSample(isample)%Avg(1:2)
            ParaPost%Sample(isample)%Avg(3)     = ParaPostSample(isample)%Avg(4)    ! Durz goes in place of Eiso
            ParaPost%Sample(isample)%Avg(4)     = ParaPostSample(isample)%Avg(3)    ! Eiso goes in place of Durz
            ParaPost%Sample(isample)%Std(1:2)   = ParaPostSample(isample)%Std(1:2)
            ParaPost%Sample(isample)%Std(3)     = ParaPostSample(isample)%Std(4)    ! Durz goes in place of Eiso
            ParaPost%Sample(isample)%Std(4)     = ParaPostSample(isample)%Std(3)    ! Eiso goes in place of Durz

            ! now reshuffle the variables so that Eiso in column 3 goes in place of Durz in column 4 in the covariance matrix
            ParaPost%Sample(isample)%CholFacLower = sortPosDefMat(ParaPost%ndim,ParaPostSample(isample)%CholFacLower,1_IK,Indx,IndxMap)
            do j = 2,ParaPost%ndim
                do i = 1, j-1
                    ParaPost%Sample(isample)%CholFacLower(j,i) = ParaPost%Sample(isample)%CholFacLower(i,j)
                end do
            end do

!write(*,"(*(g0.8,:,','))")
!write(*,"(*(g0.8,:,','))") "Original"
!do i = 1,ParaPost%ndim
!    block
!        real(RK) :: Vec(NVAR)
!        Vec = ParaPostSample(isample)%CholFacLower(i,:)
!        write(*,"(*(F15.8,:))") Vec
!    end block
!end do
!
!write(*,"(*(g0.8,:,','))") "Modified"
!do i = 1,ParaPost%ndim
!    block
!        real(RK) :: Vec(NVAR)
!        Vec = ParaPost%Sample(isample)%CholFacLower(i,:)
!        write(*,"(*(F15.8,:))") Vec
!    end block
!end do
!write(*,"(*(g0.8,:,','))")
!read(*,*)

            ! get conditional covariance matrix
            if (allocated(ParaPost%Sample(isample)%RegresCoefMat)) deallocate(ParaPost%Sample(isample)%RegresCoefMat)
            allocate(ParaPost%Sample(isample)%RegresCoefMat(3_IK,1_IK))
            allocate(ParaPost%Sample(isample)%SchurComplement(3_IK,3_IK))
            call getRegresCoef  ( rankPDM           = NVAR &
                                , rankS11           = ParaPost%numDepVar &
                                , rankS22           = ParaPost%ndim - ParaPost%numDepVar &
                                , PosDefMat         = ParaPost%Sample(isample)%CholFacLower     &
                                , RegresCoefMat     = ParaPost%Sample(isample)%RegresCoefMat    &
                                , SchurComplement   = ParaPost%Sample(isample)%SchurComplement  &
                                )
            if (ParaPost%Sample(isample)%RegresCoefMat(1,1)<0._RK) then
                write(*,"(*(g0,:,' '))") MODULE_NAME//": RegresCoefMat(1,1)<0._RK"
                error stop
            end if

!write(*,"(*(g0,:,' '))")
!write(*,"(*(g0,:,' '))") "CovMat:", shape(ParaPost%Sample(isample)%CholFacLower)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(1,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(2,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(3,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(4,:)
!read(*,*)

!write(*,"(*(g0,:,' '))")
!write(*,"(*(g0,:,' '))") "Schur:", shape(ParaPost%Sample(isample)%SchurComplement)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%SchurComplement(1,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%SchurComplement(2,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%SchurComplement(3,:)
!read(*,*)

            ! get the Cholesky factor
            allocate(ParaPost%Sample(isample)%CholFacDiag(ParaPost%numDepVar))
            deallocate(ParaPost%Sample(isample)%CholFacLower)
            allocate( ParaPost%Sample(isample)%CholFacLower, source = ParaPost%Sample(isample)%SchurComplement )

!write(*,"(*(g0,:,' '))")
!write(*,"(*(g0,:,' '))") "Chol=Schur:", shape(ParaPost%Sample(isample)%CholFacLower)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(1,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(2,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(3,:)
!write(*,"(*(g0,:,' '))") "Diag:", ParaPost%Sample(isample)%CholFacDiag
!read(*,*)

            call getCholeskyFactor  ( nd = ParaPost%numDepVar &
                                    , PosDefMat = ParaPost%Sample(isample)%CholFacLower &
                                    , Diagonal  = ParaPost%Sample(isample)%CholFacDiag &
                                    )

!write(*,"(*(g0,:,' '))")
!write(*,"(*(g0,:,' '))") "Chol:", shape(ParaPost%Sample(isample)%CholFacLower)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(1,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(2,:)
!write(*,"(*(g0,:,' '))") ParaPost%Sample(isample)%CholFacLower(3,:)
!write(*,"(*(g0,:,' '))") "Diag:", ParaPost%Sample(isample)%CholFacDiag
!write(*,"(*(g0,:,' '))")
!read(*,*)

            if (ParaPost%Sample(isample)%CholFacDiag(1)<0._RK) then
                write(*,*) "FATAL: Covariance Matrix corresponding to the parameter sample ", isample, "not positive-definite."
                write(*,*) "ParaPost%Sample(isample)%CholFacDiag(1): ", ParaPost%Sample(isample)%CholFacDiag(1)
                error stop
            end if

            !*******************************************************************************************************************
            !****                                           report time and progress                                        ****
            !*******************************************************************************************************************

            call Timer%toc()
            itick = itick + 1
            write(output_unit,"(*(' ',g0))", advance="no"   ) CARRIAGE_RETURN, CLOCK_TICK(itick) &
                                                            , isample, " out of  ", ParaPost%count &
                                                            , "parameter posterior samples processed in" &
                                                            , Timer%Time%total, "seconds."
            flush(output_unit)
            if (itick==4) itick = 0

        end do
        write(output_unit,"(*(g0))")
        close(ParaPost%fileUnit)

    end function constructParaPost


!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module ParaPost_mod