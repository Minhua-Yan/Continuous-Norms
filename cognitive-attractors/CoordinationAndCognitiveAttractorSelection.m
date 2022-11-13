% Compared to the base model, here we add a cognitive process that issimilar to selection, 
% where the norm variants' frequency change because they have different cognitive fitness

% Parameters:
%     D: variance-like coefficient of individual natural selection fitness function
%     Opt: The optimal norm trait for individual natural selection
%     f: initial norm (mean of the initial trait distribution) of the focal population
%     S_inv: strength of social selection as pairwise coordination in the focal population
%        S_inv is the inverse of S in the base model in the main text
%     E_sd: standard deviation of the error in social learning
%     R: The distance between neighbor congnitive attractors; we assume that the original focal population norm is a cognitive attractor

%% assign values to parameters

D = 0.5;    Opt = 2.39;      
f = 0;                                                                
S = 4;    E_sd = 1;
binWid = 0.01; 

p = 1;                                                                     % The distance between neighbor congnitive attractors; we assume that the original island norm is a cognitive attractor


if S > 0
    V_sd= sqrt(E_sd^2 + sqrt(E_sd^4+E_sd^2*S));                        % equilibrium variance with only transission error and norm selection, see "eq V no mig calculation.pdf"                                                          
else
    error('S must be positive')
end

%% initialize vector of norm trait values

binRange= 20*V_sd;                                                           % the trait range on one side of f we keep track of
    % reminder: this may need to be changed depending on Opt and D_inv
D_sd = sqrt(D);
if f+binRange < Opt+6*D_sd || f-binRange > Opt-6*D_sd
    error('the range of bins we keep track of needs to be increased to comfortably cover direct effects')
end
                                                              % the width of the bins we cut the distributions into
nbins = 2*round( binRange /binWid ) +1;                                      % make the total number of bins an odd number, so the center of the center bin is exactly "f"


range_min = -nbins/2*binWid + f;                                             % updated lower boundary of the trait range

xU = (1:nbins)*binWid + range_min;                                           % xU is a vector that gives the upper bound of each bin
xL = xU - binWid;                                                            % xL is a vector that gives the lower bound of each bin
xM = xU - binWid/2;                                                          % xM is a vector that gives the center of each bin


%% natural selection fitness function
NSFit = exp((-(xM-Opt).^2)/(2*D));

%% cognitive fitness function
nCog = round(binRange/p);

maxCog = f + (nCog+5)*p;                                             % make sure the calculated range of cognitive attractors is big enough so that the bin range is well covered
minCog = f - (nCog+5)*p;
CogFit = zeros(1,nbins);
for cogAtr = minCog:p:maxCog
    newAtr = exp((-4*(xM-cogAtr).^2)/2);
    CogFit = CogFit + newAtr;
end

maxCogFit = max(CogFit);
minCogFit = min(CogFit);
if maxCogFit == minCogFit
    error('the width of each cognitive attractor needs to be adjusted')
end
CogFit = (CogFit-minCogFit)/(maxCogFit-minCogFit);


figure;
plot(xM, CogFit,'LineWidth',2)
xlim([-15 15])
ylim([0 1])
xlabel('trait value')
ylabel('cognitive fitness')
ax = gca;
ax.FontSize = 12;

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
        fit_Soc(i) = sum(exp(-1*dis(i,:).*dis(i,:)/(2*S)).*fFreq);  
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

gen = 5000;                                                                  % the number of generations to run the evolution for

pdf_NS_Gen = zeros(gen, nbins);        pdf_NS_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after migration
pdf_SS_Gen = zeros(gen, nbins);        pdf_SS_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after social selection (SS)
pdf_Err_Gen = zeros(gen, nbins);       pdf_Err_Gen(1,:) = fFreq;             % gth row: trait distribution in generation g after SL error
pdf_CA_Gen = zeros(gen, nbins);        pdf_CA_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after cognitive fitness process
fit_Soc_Gen = zeros(gen,nbins);                                              % gth row: fitness values for xM in social selection in generation g
for i = 1: nbins
        fit_Soc_Gen(1,i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
end
Mean_Gen = zeros(gen, 4);              Mean_Gen(1,:) = ones(1,4)*f;          % 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with err





        
%-------------------------------------------------------------------------- % evolutionary process of the island under the influence of migration from the continent, norm selection against deviants from island mean and transmission
for g = 2:gen
    g
  % natural selection
    fFreq = fFreq.*NSFit;
    fFreq = fFreq/sum(fFreq);
    pdf_NS_Gen(g,:) = fFreq; 
    Mean_Gen(g,1) = sum(fFreq.*xM);

  % social selection 
    fit_Soc = zeros(1,nbins);
    for i = 1: nbins
        fit_Soc(i) = sum(exp(-1 *dis(i,:).*dis(i,:)/(2*S)).*fFreq);  
    end
    fit_Soc_Gen(g,:) = fit_Soc;         
    
    fFreq = fFreq.*fit_Soc;
    fFreq = fFreq/sum(fFreq);
    pdf_SS_Gen(g,:) = fFreq;
    Mean_Gen(g,2) = sum(fFreq.*xM);
    
  % unbiased social learning with error
    fFreq = fFreq*Trans;                                                    
    fFreq = fFreq/sum(fFreq);
    pdf_Err_Gen(g,:) = fFreq;
    Mean_Gen(g,3) = sum(fFreq.*xM);
    
  % biased cultural transmission
    fFreq = fFreq.*CogFit;
    fFreq = fFreq/sum(fFreq);                                               
    pdf_CA_Gen(g,:) = fFreq;
    Mean_Gen(g,4) = sum(fFreq.*xM);  
    Mean_Gen(g,4)
    
end
figure
hold on;
gen_plot = 5;
plot(xM, pdf_NS_Gen(gen_plot,:),xM, pdf_SS_Gen(gen_plot,:),xM, pdf_CA_Gen(gen_plot,:),'LineWidth',2);
title(['population distribution after different stages in generation', num2str(gen_plot)])
legend('after natural selection','after social selection','after cultural transmission')
xlabel('trait value')
ylabel('density')
xlim([-15 15])
ax = gca;
ax.FontSize = 14;

figure;
hold on;
tplot = 0:niter-1;
plot(tplot,Mean_Gen(:,3),'LineWidth',2)
xlabel('generation')
ylabel('mean of norm trait distribution')
line([0,100],[3,3],'color','black','LineWidth',1)
xlim([0 100])
ylim([0 3])
ax = gca;
ax.FontSize = 14;

function sd=getsd(mypdf, xvec)
    mypdf = max(mypdf, zeros(1,length(xvec)));
    mycdf = cumsum(mypdf);
    mycdf = min(mycdf,ones(1,length(xvec)));   
    mycdf(end) = 1; mycdf(1) = 0;
    pd = makedist('PieceWiselinear','X',xvec,'FX',mycdf);
    sd = sqrt(var(pd));
end


