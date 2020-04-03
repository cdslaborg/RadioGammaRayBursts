function Func = integrand(Energy)
    global epk alpha beta
    Func = zeros(1,length(Energy));
    for i = 1:length(Energy)
        Func(i) = Energy(i) * getBandFunc(Energy(i),epk,alpha,beta);
    end
    %disp([epk, alpha, beta])
end