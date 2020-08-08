function BandFunc = getBandFunc(energy,epk,alpha,beta)
% returns the value of the Band differential spectrum for the given energy and input parameters.
% It is expected that the input energy is in units of KeV, although it does not affect the computations here.
% NOTE: A negative huge output value is used to signal error has occurred. Under normal conditions, output is always positive.

    KEV2ERG = 1.60217662080000e-9;

    if nargin==2
        alpha = -1.1;
        beta = -2.3;
    elseif nargin==3
        beta = -2.3;
    elseif nargin~=4
        error('Incorrect input arguments...');
    end

    % check if the photon indices are consistent with the mathematical rules

    if alpha<beta || alpha<-2.
        warning('alpha<beta .or. alpha<-2.');
    end

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

