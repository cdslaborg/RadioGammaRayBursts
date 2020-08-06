classdef BandSpectrum

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        KEV2ERG = 1.60217662080000e-9;
    end

    properties(Access = public)
        alpha = -1.1;
        beta = -2.3;
        epk
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = public)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = BandSpectrum(epk,alpha,beta)
            if nargin<4
                self.beta = beta;
                if nargin<3
                    self.alpha = alpha;
                    if nargin<2
                        self.epk = epk;
                    end
                end
            else
                error("The correct calling syntax is BandSpectrum(epk,alpha,beta).");
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function BandFunc = getBandFunc(energy,epk,alpha,beta)
            % returns the value of the Band differential spectrum for the given energy and input parameters.
            % It is expected that the input energy is in units of KeV, although it does not affect the computations here.
            % NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.
            %
            % Parameters
            %
            %   energy
            %       a scalar or vector of real number at which we want to calculate the Band spectrum
            %
            %   epk
            %       scalar real number representing the spectral peak energy of the Band function
            %
            % Returns
            %
            %   BandFunc
            %       a scalar or vector of the same length as energy
            %

            % check if the photon indices are consistent with the mathematical rules

            BandSpectrum.verifyPhotonIndices(alpha,beta);

            alphaPlus2 = alpha + 2.;
            alphaMinusBeta = alpha - beta;
            ebrk = epk * alphaMinusBeta / alphaPlus2;
            coef = ebrk^alphaMinusBeta * exp(-alphaMinusBeta);

            % compute the spectrum

            if energy <= ebrk
                BandFunc = energy^alpha * exp(-energy*alphaPlus2/epk);
            else
                BandFunc = coef * energy^beta;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function PhotonFlux = getPhotonFlux(energy,epk,alpha,beta,ebrk,coef,alphaPlusTwo)
            % returns the value the Band differential spectrum for the given energy and input parameters.
            % It is expected that the input energy is in units of KeV, although it does not affect the computations here.
            % Under normal conditions, output is always positive.

            ! check if the photon indices are consistent with the mathematical rules

            BandSpectrum.verifyPhotonIndices(alpha,beta);

            if (alpha<beta .or. alpha<-2._RK) then
                PhotonFlux = -HUGE_RK
                return
            end if

            ! compute the spectrum
            if (energy<=ebrk) then
                PhotonFlux = energy**alpha * exp(-energy*alphaPlusTwo/epk)
            else
                PhotonFlux = coef * energy**beta
            end if
        end function getPhotonFlux

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function verifyPhotonIndices(alpha,beta)
            if alpha < beta || alpha < -2.; error("alpha<beta .or. alpha<-2: " + string(alpha) + " " + string(beta) ); end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end