%   This program investigates how the norm in a focal population evolves under migration and social selection on the norm's category
%   It is the leaner version of "discretization" as a function

%   The function is used in "EQ_to_m_cat_constant", "EQ_to_c_cat_constant"
%   This function takes m, S, E_sd, step, f, and z as inputs (see below for their meanings) and gives 
%   1) the effect of social selection on changing norm mean in the first generation with migration
%   2) the equilibrium in four forms of representation
%             1. mean of the equilibrium distribution of true trait values
%             2. mean of the equilibrium distribution of categories
%             3. the category with the highest fitness in coordination at equilibrium
%             4. the most common type after social selection as coordination at equilibrium


% Input arameters:
%     m: migration rate
%     S: strength of social selection as pairwise coordination in the focal population
%     E_sd: standard deviation of the error in social learning
%     step: how many base bins there are in each category
%     f: initial norm (mean of the initial trait distribution) of the focal population
%     z: norm (mean of the trait distribution) of the source population




function result=Function_discretization_migration_fast(m,S_inv,E_sd,cbn,f,z,binWid)
    if S_inv > 0
        V_sd= sqrt(E_sd^2 + sqrt(E_sd^4+E_sd^2/S_inv));                       % equilibrium variance with only transission error and norm selection, see "eq V no mig calculation.pdf"
        U_sd= V_sd;                                                           % for simplicity and empirical realism, we assume the source population has the same S and E as the focal
    elseif S_inv == 0
        error('S_inv must be positive')
    end

    %% initialize vector of norm trait values
    binRange= 20*V_sd;                                                           % the width of the bins we cut the trait distribution into
    nbins = 2*round( binRange /binWid ) +1;                                    % make the total number of bins an odd number, so the center of the center bin is exactly "f"
        
    range_min = -nbins/2*binWid + f;
    
    xM = (1:nbins)*binWid + range_min- binWid/2;                                    % xM is a vector that gives the center of each bin
        
        
    
    %% add in bins at left and right of the range we consider so that the bins can be turned into complete categories
        
    if rem(nbins, cbn) ~= 0
        for j = 1:cbn-1
            if rem((nbins-cbn)/2, cbn) == j
                nbins = nbins + 2*(cbn-j);
                xM = [xM(1)-binWid*((cbn-j):-1:1),xM,xM(end)+binWid*(1:(cbn-j))];
            end
        end
    end
    xU = xM+ binWid/2;                                                          %update xU
    xL = xM- binWid/2;                                                          %update xL

    nCat = nbins/cbn;
    catM = (1:nCat)*cbn*binWid + (xM(1)-binWid/2) - cbn*binWid/2;      %the center values of the categories
    
        
    
     
    

    %% initialize copying error, a normal distribution with mean == 0 and standard deviation E_sd
    % error is also approxiamited with bins of width "binWid"       
    errRange = 10*E_sd;                                                           % the error range on one side of 0 we keep track of
    nEbins = 2* round( errRange/binWid ) + 1;                                     % number of error bins, odd so the center of the center bin is exactly 0
    minErrRange = -nEbins/2*binWid;                                               % updated lower boundary of the error range
    xU_Err = ((1:nEbins)*binWid + minErrRange);                                   % upper bounds of the error bins
    xL_Err = xU_Err - binWid;                                                     % lower bounds of the error bins
    
    %errFreq(i) is a vector of the probability of the error falling in the ith bin's range: [minErr +(i-1)*binWid, minErr +i*binWid]
    errFreq = normcdf(xU_Err,0,E_sd) - normcdf(xL_Err,0,E_sd);                    % a row, density distribution of the error
    errFreq = errFreq/sum(errFreq);      


    
    %% transition matrix for errorous social learning process if error is continuous
        
    binErrTrans = zeros(nbins, nbins);                                               %(again, if a bin will produce cultural children that go out of the range of our tracked bins, those are treated as not reproducing)
    n0 = 0;
    for i = ((nEbins+1)/2) : (nbins- (nEbins-1)/2)
        binErrTrans(i,:) = [zeros(1,n0),errFreq,zeros(1,nbins-nEbins-n0)];
        n0 = n0+1;
    end
    %(note: if a bin will produce cultural children that go out of the range of our tracked bins, 
    % those are treated as not reproducing
    % accuracy is ensured by making the range of trait bins sufficiently large, especially compared to the range of error)
    

    
    %% matrix of distance from each category to each category-- to be used in social selection
    dis_cat = zeros(nCat, nCat);
    for i=1:nCat
        dis_cat(i,:) = abs((0:(nCat-1))*cbn*binWid-(i-1)*cbn*binWid);
    end
            

    %% initialize focal population N(f,V)
        % because the initial distribution is not normal under discretization, we approximate it with the corresponding normal distribution
        % and then let the population evolve under SS and SL error to get the real initial distribution
    fBinFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                              % a row, density of the island tra
    fBinFreq = fBinFreq/sum(fBinFreq);
    sd_current = getsd(fBinFreq, xM);

    pre_gen = 3000;

    for pre_g = 1:pre_gen
        sd_bf =  sd_current;                                                           % the sd of the population distribution in the generation before
        
        %% SS
        catFit=zeros(1,nCat);
        binFit=zeros(1,nbins);

        fCatFreq = get_catPdf(fBinFreq, cbn, nCat);

        for i=1:nCat
            catFit(1,i) = sum(exp(-S_inv*dis_cat(i,:).*dis_cat(i,:)/2).*fCatFreq);      % this is the fitness value of each category when coordination is by category
            for j=0:cbn-1
                binFit(1,i*cbn-j) = catFit(1,i);
            end
        end   
        fBinFreq = fBinFreq.*binFit;
        fBinFreq = fBinFreq/sum(fBinFreq);

        %% SL no error
    
        % version 1: trying to learn the true norm value of a single model
        fBinFreq = fBinFreq;
            
    %     % version 2: trying to learn the category value of a single model
    %         for i= (0:nCat-1)*cbn + (cbn+1)/2
    %             for j = 1:((cbn-1)/2)
    %                 fBinFreq(i) = fBinFreq(i) + fBinFreq(i-j) + fBinFreq(i+j);
    %                 fBinFreq(i-j) = 0;
    %                 fBinFreq(i+j) = 0;
    %             end
    %         end
    %         fBinFreq = fBinFreq/sum(fBinFreq);
    %         binPdf_SL_Gen(g,:) = fBinFreq;
    %         binMean_Gen(g,2) = sum(fBinFreq.*xM);


        %% SL error
        fBinFreq = fBinFreq*binErrTrans;
        fBinFreq = fBinFreq/sum(fBinFreq);

        sd_current = getsd(fBinFreq,xM);
    
    
        % jumping out of the loop after the population has already stablized  
        if (abs(sd_bf - sd_current) < 0.00000001)
            break
        end
    
        if pre_g==pre_gen
            error('initial population distribution has not stablized')
        end
    end


 

     %% initialize source population N(c,U)
    zFreq = normcdf(xU,z,U_sd) - normcdf(xL,z,U_sd);                                  % a row, density of the continent trait values between [minX +(i-1)*binWid, minX +i*binWid]
    zFreq = zFreq/sum(zFreq);
    z = sum(zFreq.*xM);      
        

    %% prepare matrices to contain information through the evolutionary process (for plotting)
    
    gen = 200000;                                                                 % the number of generations to run the evolution for
    
              
    binFit = zeros(1,nbins);                                           % gth row: fitness values for each base bin in social selection in generation g
    catFit = zeros(1, nCat);                                           % gth row: fitness values for each category in social selection in generation g
    
    
        
    binMean_Gen = zeros(gen, 3);                                             % mean values of the base bin distribution-- 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with err

    binMean_Gen(1,:) = ones(1,3)*sum(fBinFreq.*xM);


    %% focal population distribution's evolutionary process
    for g = 2:gen
        g;    
       %% evolutionary section: migration
        % distribution of true trait values
        fBinFreq = (1-m).*fBinFreq + m.*zFreq;                                       
        fBinFreq = fBinFreq/sum(fBinFreq);
        binMean_Gen(g,1) = sum(fBinFreq.*xM);
        
        %% analytical section: distribution of categories after migration
        if cbn==1
            fCatFreq = fBinFreq;
        else
            fCatFreq = sum(reshape(fBinFreq,cbn,nCat));  
        end
      
        

        %% evolutionary section: social selection
         % distribution of the true continuous trait values
        for i=1:nCat
            catFit(1,i) = sum(exp(-S_inv*dis_cat(i,:).*dis_cat(i,:)/2).*fCatFreq);      % this is the fitness value of each category when coordination is by category
            for j=0:cbn-1
                binFit(1,i*cbn-j) = catFit(1,i);
            end
        end   
        fBinFreq = fBinFreq.*binFit(1,:);
        fBinFreq = fBinFreq/sum(fBinFreq);
        binMean_Gen(g,2) = sum(fBinFreq.*xM);

            

    
        %% evolutionary section: social learning with no error
        %% version 1: trying to learn the true norm value of a single model
        fBinFreq = fBinFreq;
        
%         % version 2: trying to learn the category value of a single model
%         for i= (0:nCat-1)*cbn + (cbn+1)/2
%             for j = 1:((cbn-1)/2)
%                 fBinFreq(i) = fBinFreq(i) + fBinFreq(i-j) + fBinFreq(i+j);
%                 fBinFreq(i-j) = 0;
%                 fBinFreq(i+j) = 0;
%             end
%         end
%         fBinFreq = fBinFreq/sum(fBinFreq);
%         binPdf_SL_Gen(g,:) = fBinFreq;
%         binMean_Gen(g,2) = sum(fBinFreq.*xM);

                   
    
        %% unbiased continuous social learning error

        fBinFreq = fBinFreq*binErrTrans;
        fBinFreq = fBinFreq/sum(fBinFreq);
        binMean_Gen(g,3) = sum(fBinFreq.*xM);
        binMean_Gen(g,3);
            

        %% jumping out of the loop after the population has already stablized 
        if g>101
            if (abs(binMean_Gen(g,3) - z) < 0.000001 && abs(binMean_Gen(g-100,3) - z) < 0.000001) || (abs(binMean_Gen(g,3) - binMean_Gen(g-1,3)) < 0.000000001 && abs(binMean_Gen(g,3) - binMean_Gen(g-100,3)) < 0.0000001)
               break
            end
        end    

    end    
    % the evolutionary process ends


    if g==gen
        error('Generations not enough, focal population norm not stablized')
    end
    


    %% calculate the change of mean norm from social selection
    SS_effect_Gen = binMean_Gen(1:g,2) - binMean_Gen(1:g,1);
    SS_effect_gen2 = SS_effect_Gen(2,1);

    if z < f
        SS_effect_max = max(SS_effect_Gen);
    else
        SS_effect_max = min(SS_effect_Gen);
    end


    
    %% calculate the equilibrium true mean norm
    binMean_Eq = binMean_Gen(g,3);


    %% calculate the category with the highest fitness at equilibrium
    maxfit = max(catFit);
    cat_mostFit = catM(catFit==maxfit);

    result=[binMean_Eq, cat_mostFit, SS_effect_gen2, SS_effect_max];


end

function sd=getsd(mypdf, xvec)
    mypdf = max(mypdf, zeros(1,length(xvec)));
    mycdf = cumsum(mypdf);
    mycdf = min(mycdf,ones(1,length(xvec)));   
    mycdf(end) = 1; mycdf(1) = 0;
    pd = makedist('PieceWiselinear','X',xvec,'FX',mycdf);
    sd = sqrt(var(pd));
end

function catPdf = get_catPdf(binFreq, cbn, nCat)
    if cbn==1   
        catPdf = binFreq;
    else
        catPdf = sum(reshape(binFreq,cbn,nCat));  
    end
end
