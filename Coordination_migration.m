%   Author: Robert Boyd, Minhua Yan
%   This program investigates how the norm in a focal population evolves under migration and social selection
%      The initial norm trait distribution in the focal population follows N(f,V), i.e., the initial norm is f and the initial variance is V
%         where V is the equilibrium distribution variance with no migration, given the strength of social selection and magnitude of social learning error
%      In each generation, the population takes in immigrants from a source population with a constant distribution of N(c, U)
%         such that m of the post-migration population is immigrants
%      The population then go through social selection, where each member's fitness is the expected payoff from pairwise coordination with other members
%      The next generation is composed of all naive individuals who have socially learned the norm trait 
%         with probability proportional to the social learning model's fitness, and with a error that follows N(0, E)
%   Here we have migration happen before social selection
%   Reversing their order gives the same results

% Parameters:
%     m: migration rate
%     z: norm (mean of the trait distribution) of the source population
%     U_sd: standard deviation of the source population's trait distribution
%     f: initial norm (mean of the initial trait distribution) of the focal population
%     S: strength of social selection as pairwise coordination in the focal population
%     E_sd: standard deviation of the error in social learning

%% ----------------------------------------------------------------- assign values to parameters

m = 0.0000000;               
z = 3;      f = 0;                                                                
U_sd = 2;   S = 2;    E_sd = 1;



%% ----------------------------------------------------------------- initialize vector of norm trait values

% since we cannot maintain normality, we solve for the distribution analytically by approxiamating the distribution of the norm trait and of the error, with bins

binRange= 30;                                                                % the trait range on one side of f we keep track of
binWid = 0.01;                                                               % the width of the bins we cut the distributions into
nbins = 2*round( binRange /binWid ) +1;                                      % make the total number of bins an odd number, so the center of the center bin is exactly "f"

minRange = -nbins/2*binWid + f;                                              % updated lower boundary of the trait range

xU = (1:nbins)*binWid + minRange;                                            % xU is a vector that gives the upper bound of each bin
xL = xU - binWid;                                                            % xL is a vector that gives the lower bound of each bin
xM = xU - binWid/2;                                                          % xM is a vector that gives the center of each bin



%% initialize focal population N(f,V)
V = E_sd^2 + sqrt(E_sd^4+E_sd^2/S);                                          % equilibrium variance with only transission error and norm selection, see "eq V no mig calculation.pdf"
V_sd = sqrt(V);                                                              % equilibrium standard deviation
iFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             % a row, density of the initial island norm traits; the ith element gives the frequency of the ith bin
iFreq = iFreq/sum(iFreq);


%% initialize source population N(c,U)
zFreq = normcdf(xU,z,U_sd) - normcdf(xL,z,U_sd);                             % a row, density of the continent trait values between [minX +(i-1)*binWid, minX +i*binWid]
zFreq = zFreq/sum(zFreq);


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


%% transition matrix for errorous social learning process
% If a cultural parent will produce cultural children out of the trait range we keep track of, 
%    we assume this cultural parent does not reproduce.
%    therefore, the probability of the first and last (nEbins-1)/2 bins transitioning to any other bin is 0

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

   
%% prepare matrices to contain information through the evolutionary process (for plotting)

gen = 30000;                                                                     % the number of generations to run the evolution for

pdf_Mig_Gen = zeros(gen, nbins);                                             % gth row: trait distribution in generation g after migration
pdf_SS_Gen = zeros(gen, nbins);                                              % gth row: trait distribution in generation g after social selection (SS)
pdf_SL_Gen = zeros(gen, nbins);                                              % gth row: trait distribution in generation g after social learning (SL) w/o error
pdf_Err_Gen = zeros(gen, nbins);                                             % gth row: trait distribution in generation g after SL error

Fit_Gen = zeros(gen,nbins);                                                  % gth row: fitness values for xM in social selection in generation g
Mean_Gen = zeros(gen, 3);                                                    % 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with err


%% focal population distribution's evolutionary process

for g = 2:gen

  % migration
    iFreq = (1-m).*iFreq + m.*zFreq;    
    iFreq = iFreq/sum(iFreq);
    pdf_Mig_Gen(g,:) = iFreq; 
    Mean_Gen(g,1) = sum(iFreq.*xM);
    

  % social selection as coordination
    fit = zeros(1,nbins);
    for i = 1: nbins
        fit(i) = sum(exp(-S*dis(i,:).*dis(i,:)/2).*iFreq);  
    end
    Fit_Gen(g,:) = fit;         
    
    iFreq = iFreq.*fit;
    iFreq = iFreq/sum(iFreq);
    pdf_SS_Gen (g,:) = iFreq;
    Mean_Gen(g,2) = sum(iFreq.*xM);
 
    
%   % trying to learn the mean of two cultural parents
%     iCombM = transpose(iFreq)*iFreq;                                        %a matrix containing probabilities of different combinations                                                 %the number of averages of combinations
%     iCombMean = zeros(1,2*nbins-1);                                         %a rwo containing the probabilities of different averages of combinations
%     for i = 1:nbins
%         for j = 1:i
%            iCombMean(1,i) = iCombMean(1,i) +iCombM(j,i+1-j); 
%         end
%     end
%     for i = nbins+1 : 2*nbins-1
%         for j = i+1-nbins :nbins
%             iCombMean(1,i) = iCombMean(1,i) +iCombM(j,i+1-j);
%         end
%     end
%     iFreq = iCombMean(1:2:end) + [1/2*iCombMean(2:2:end),0]+[0,1/2*iCombMean(2:2:end)];   
%     Dis_SL_Gen(t,:)=iFreq; 
%     Mean_Gen(g,3) = sum(iFreq.*xM);    
      
    
    
    % social learning error
    iFreq = iFreq*Trans;
    iFreq = iFreq/sum(iFreq);
    pdf_Err_Gen(g,:) = iFreq;
    Mean_Gen(g,3) = sum(iFreq.*xM);                                          
    
end

%% figure 1: the norm trait distribution after different stages (migration, social selection, social learning w/ error) in one generation 
figure; hold on;
g_plot=2;                 % the generation we plot the distributions for
plot(xM, pdf_Mig_Gen(g_plot,:), xM, pdf_SS_Gen(g_plot,:),  xM, pdf_Err_Gen(g_plot,:),  'LineWidth',1);
legend('after migration','after social selection','after social learning w/ error')
xlabel('trait value');
ylabel('density');
xlim([-10 10])
title(['norm trait distribution after different stages in generation ', num2str(g_plot)]);
ax = gca;
ax.FontSize = 16;

%% figure 2: the norm evolutionary trajectory
figure; hold on;
plot(1:gen, Mean_Gen(:,3),'LineWidth',2)
xlabel('generation')
ylabel('mean of norm trait distribution')
line([0,300],[z,z],'color','black','LineWidth',1)
title('norm evolution trajectory in the focal population')
ax = gca;
ax.FontSize = 16;



% Mean_each_step = zeros(1,3*g);
% for g = 1:gen
%     Mean_each_step(1,3*g-2)=Mean_Gen(g,1);     %after migration
%     Mean_each_step(1,3*g-1)=Mean_Gen(g,2);     %after social selection
%     Mean_each_step(1,3*g)=Mean_Gen(g,3);       %after social learning with error
% end
% 
% figure;
% hold on;
% plot(1:3*gen,  Mean_each_step,   'lineWidth',0.01)
% xlabel('generation')
% ylabel('mean of norm trait distribution')
% title('norm evolution trajectory in the focal population')


% figure;
% hold on;
% [ax, h1, h2] = plotyy(xM(500:2502), iNSpostMig(2,500:2502),xM(500:2502),iNSfit(2,500:2502));
% h1.LineStyle = '--';
% h1.LineWidth = 2;
% h2.LineWidth = 2;
% xlim([-10 10])
% ylim(ax(1),[0 0.007]);
% ylim(ax(2),[0 0.8]);
% line([iNSMean(2,2),iNSMean(2,2)],[0,0.2],'color','black','LineWidth',1.5)
% legend('population distribution','population mean','social selection fitness function')
% x3=xlabel('trait value','Position',[0 -0.00045 -1]);
% h= title('social selection fitness function for a population distribution');
% ax = gca;
% ax.FontSize = 16;
