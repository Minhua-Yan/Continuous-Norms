%   This program investigates how the norm in a focal population evolves 
%   under natural selection and social selection, when the cognitive process has attractors

%   The function is used in "EQ_to_D_inv"
%   This function takes D_inv, Opt, S_inv, E_sd, f, R and h as inputs (see below for their meanings) and gives the equilibeium 



%   Input parameters:
%           D_inv: the strength of direct effects (natural selection) on the normative traits
%           Opt: the optimal normative trait in direct effects (natural selection)
%           f: initial norm (mean of the initial trait distribution) of the focal population
%           S_inv: the strength of social selection as coordination
%                  S_inv is the inverse of S in the base model in the main text
%           E_sd: standard deviation of the error in social learning
%           R: distance between neighboring cognitive attractors, 
%              operated as the absolute value difference between two neighboring cognitive attractors' centers
%              b/binWid needs to be even
%           h: cognitive fitness difference between cognitive attractors and non-cognitive attractors 
%                (highest cognitive fitness (reached at cognitive attractors) is set to be 1, lowest cognitive fitness is set to be 1-h       
              



function eq = Function_CA_Selection_D_inv(D_inv, Opt, f, S, E_sd, p)
    if S > 0
        V_sd= sqrt(E_sd^2 + sqrt(E_sd^4+E_sd^2*S));                       % equilibrium variance with only transission error and norm selection, see "eq V no mig calculation.pdf"                                                       
    else
        error('S must be positive')
    end                    

    if D_inv<0 || E_sd <=0 || p <=0 
        error('error: one or more parameter values are not biologically sensible')
    end


%% initialize vector of norm trait values

    % since we cannot maintain normality, we solve for the distribution analytically by approxiamating the distribution of the norm trait and of the error, with bins

    binRange= 20*V_sd;                                                                % the trait range on one side of f we keep track of
    binWid = 0.01;                                                              % the width of the bins we cut the distributions into
    nbins = 2*round( binRange /binWid ) +1;                                      % make the total number of bins an odd number, so the center of the center bin is exactly "f"

    range_min = -nbins/2*binWid + f;                                              % updated lower boundary of the trait range

    xU = (1:nbins)*binWid + range_min;                                           % xU is a vector that gives the upper bound of each bin
    xL = xU - binWid;                                                            % xL is a vector that gives the lower bound of each bin
    xM = xU - binWid/2;                                                          % xM is a vector that gives the center of each bin

    %% natural selection fitness function
    NSFit = exp((-D_inv*(xM-Opt).^2)/2);


    %% cognitive fitness function
    nCog = round(binRange/p);
    
    maxCog = f + (nCog+5)*p;                                             % make sure the calculated range of cognitive attractors is big enough so that the bin range is well covered
    minCog = f - (nCog+5)*p;
    CogFit = zeros(1,nbins);
    for cogAtr = minCog:p:maxCog
        newAtr = exp((-1*(xM-cogAtr).^2)/2);
        CogFit = CogFit + newAtr;
    end
    maxCogFit = max(CogFit);
    minCogFit = min(CogFit);
    CogFit = (CogFit-minCogFit)/(maxCogFit-minCogFit);



    %% initialize copying error
    
    % the copying error follows a normal distribution with mean == 0 and standard deviation E_sd
    % error is also approxiamited with bins of width "binWid" 
    
    errRange = 6*E_sd;                                                            % the error range on one side of 0 we keep track of
    nEbins = 2* round( errRange/binWid ) + 1;                                     % number of error bins, odd so the center of the center bin is exactly 0
    minErrRange = -nEbins/2*binWid;                                               % updated lower boundary of the error range
    xU_Err = ((1:nEbins)*binWid + minErrRange);                                   % upper bounds of the error bins
    xL_Err = xU_Err - binWid;                                                     % lower bounds of the error bins
    
    %errFreq(i) is a vector of the probability of the error falling in the ith bin's range: [minErr +(i-1)*binWid, minErr +i*binWid]
    errFreq = normcdf(xU_Err,0,E_sd) - normcdf(xL_Err,0,E_sd);                    % a row, density distribution of the error
    errFreq = errFreq/sum(errFreq);      


    %% transition matrix for errorous social learning process
    % If a cultural parent will produce cultural children out of the trait range we keep track of, 
    %    we assume this cultural parent does not reproduce.
    %    therefore, the probability of the first and last (nEbins-1)/2 bins transitioning to any other bin is 0
    % To make sure the results are accurate, binRange has to be big enough
    
    % The element on the ith row and jth column gives the probability 
    %    that a cultural parent in the ith bin produces a cultural child in the jth bin in xM
    
    Trans = zeros(nbins, nbins);
    n0 = 0;
    for row = ((nEbins+1)/2) : (nbins- (nEbins-1)/2)
        Trans(row,:) = [zeros(1,n0),errFreq,zeros(1,nbins-nEbins-n0)];
        n0 = n0+1;
    end

    %% matrix of distance from each bin to each bin-- to be used in norm selection
    dis = zeros(nbins, nbins);
    for row=1:nbins
        dis(row,:) = abs((0:(nbins-1))*binWid-(row-1)*binWid);
    end


    %% initialize focal population N(f,V)
        % because the initial distribution is not normal under biased cognitive process, we approximate it with the corresponding normal distribution
        % and then let the population evolve under SS and SL error to get the real initial distribution
    fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                              % a row, density of the island tra
    fFreq = fFreq/sum(fFreq);
    sd_current = getsd(fFreq, xM);
    
    
    for pre_gen = 1:3000
        sd_bf =  sd_current;                                                           % the sd of the population distribution in the generation before
        
        
        %% social selection 
        fit_Soc = zeros(1,nbins);
        for i = 1: nbins
            fit_Soc(i) = sum(exp(-1 *dis(i,:).*dis(i,:)/(2*S)).*fFreq);  
        end       
        
        fFreq = fFreq.*fit_Soc;
        fFreq = fFreq/sum(fFreq);
    
    
        %% unbiased social learning with error
        fFreq = fFreq*Trans;                                                    
        fFreq = fFreq/sum(fFreq);
        
        %% biased cultural transmission
        fFreq = fFreq.*CogFit;
        fFreq = fFreq/sum(fFreq);                                               
        sd_current = getsd(fFreq, xM);
    
    
        % jumping out of the loop after the population has already stablized  
        if (abs(sd_bf - sd_current) < 0.00000001)
            break
        end
    
        if pre_gen==3000
            error('initial population distribution has not stablized')
        end
    end

    %% prepare matrices to contain information through the evolutionary process (for plotting)

    gen = 50000000;                                                        % the number of generations to run the evolution for
    
    fit_Soc = zeros(1, nbins);                                             % initiate Social fitness vector for the bins
                        
    Mean_Gen = zeros(gen, 1);    Mean_Gen(1,1) = sum(fFreq.*xM);
    %% focal population distribution's evolutionary process
    for g = 2:gen
        g;
      % natural selection
        fFreq = fFreq.*NSFit;
        fFreq = fFreq/sum(fFreq);
        
      % social selection 
        for i = 1: nbins 
            fit_Soc(i) = sum(exp(-1*dis(i,:).*dis(i,:)/(2*S)).*fFreq);  
        end    
        fFreq = fFreq.*fit_Soc;
        fFreq = fFreq/sum(fFreq);
    
      % unbiased social learning
        fFreq = fFreq*Trans;                                                    
        fFreq = fFreq/sum(fFreq);
  
    
      % biased cultural transmission
        fFreq = fFreq.*CogFit;
        fFreq = fFreq/sum(fFreq);   

        Mean_Gen(g,1) = sum(fFreq.*xM);
        Mean_Gen(g,1);
    
        %% jumping out of the loop after the population has already stablized 
        if g>101
            if (abs(Mean_Gen(g,1) - Opt) < 0.000001 && abs(Mean_Gen(g-100,1) - Opt) < 0.000001) || (abs(Mean_Gen(g,1) - Mean_Gen(g-1,1)) < 0.000000001 && abs(Mean_Gen(g,1) - Mean_Gen(g-100,1)) < 0.0000001)
               break
            end
        end    
    end
    % the evolutionary process ends

    if g==gen
        error('Error occurred: generations not enough')
    end
    
    
    eq = Mean_Gen(g,1);
end

function sd=getsd(mypdf, xvec)
    mypdf = max(mypdf, zeros(1,length(xvec)));
    mycdf = cumsum(mypdf);
    mycdf = min(mycdf,ones(1,length(xvec)));   
    mycdf(end) = 1; mycdf(1) = 0;
    pd = makedist('PieceWiselinear','X',xvec,'FX',mycdf);
    sd = sqrt(var(pd));
end




