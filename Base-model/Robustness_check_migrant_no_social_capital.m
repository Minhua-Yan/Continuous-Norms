%   This program investigates how the norm in a focal population evolves under migration and social selection
%      The initial norm trait distribution in the focal population follows N(f,V), i.e., the initial norm is f and the initial variance is V
%         where V is the equilibrium distribution variance with no migration, given the strength of social selection and magnitude of social learning error
%      In each generation, the population takes in immigrants from a source population with a constant distribution of N(c, U)
%         such that m of the post-migration population is immigrants
%      The population then go through social selection
%      Compared to the base model, the migrants have no social capital,
%         so a norm trait's expected payoff in similarity-based norm effects is determined only by the residents' norm trait distribution
%      The next generation is composed of all naive individuals who have socially learned the norm trait 
%         with probability proportional to the social learning model's fitness, and with a error that follows N(0, E)
%   Here we have migration happen before social selection
%   Reversing their order gives the same results

% Parameters:
%     m: migration rate
%     z: norm (mean of the trait distribution) of the source population
%     U_sd: standard deviation of the source population's trait distribution
%     f: initial norm (mean of the initial trait distribution) of the focal population
%     S_inv: strength of social selection as pairwise coordination in the focal population
%        S_inv is the inverse of S in the base model in the main text
%     E_sd: standard deviation of the error in social learning

%% ----------------------------------------------------------------- assign values to parameters

m = 0.01;               
z = 3;      f = 0;                                                                
S_inv = 0.5;    E_sd = 1;

binWid = 0.01;                                                            % the width of the bins we cut the distributions into
    % if the migration rate is low or if V_sd is small, the binWid has to be set smaller
    % Otherwise we get the same effect as discretization and the norm does not converge to the source population norm

if S_inv > 0
    V_sd= sqrt(E_sd^2 + sqrt(E_sd^4+E_sd^2/S_inv));                        % equilibrium variance with only transission error and norm selection, see "eq V no mig calculation.pdf"                                                          
elseif S_inv == 0
    error('S_inv must be positive')
end
U_sd = 0.1;                                                 


%% ----------------------------------------------------------------- initialize vector of norm trait values

% since we cannot maintain normality, we solve for the distribution analytically by approxiamating the distribution of the norm trait and of the error, with bins

binRange= 20*V_sd;                                                          % the trait range on one side of f we keep track of
if f+binRange < z+6*U_sd || f-binRange > z-6*U_sd
    error('the range of bins we keep track of needs to be increased to comfortably cover the migrants')
end
nbins = 2*round( binRange /binWid ) +1;                                    % make the total number of bins an odd number, so the center of the center bin is exactly "f"

range_min = -nbins/2*binWid + f;                                           % updated lower boundary of the trait range

xU = (1:nbins)*binWid + range_min;                                            % xU is a vector that gives the upper bound of each bin
xL = xU - binWid;                                                            % xL is a vector that gives the lower bound of each bin
xM = xU - binWid/2;                                                          % xM is a vector that gives the center of each bin



%% initialize focal population N(f,V)
fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             % a row, density of the initial island norm traits; the ith element gives the frequency of the ith bin
fFreq = fFreq/sum(fFreq);
f = sum(fFreq.*xM);                                                          % make sure this is close enough to the original value assigned to f



%% initialize source population N(z,U)
zFreq = normcdf(xU,z,U_sd) - normcdf(xL,z,U_sd);                             % a row, density of the continent trait values between [minX +(i-1)*binWid, minX +i*binWid]
zFreq = zFreq/sum(zFreq);


%% initialize copying error, a normal distribution with mean == 0 and standard deviation E_sd
% error is also approxiamited with bins of width "binWid"            
errRange = 15*E_sd;                                                            % the error range on one side of 0 we keep track of
if f+errRange < z+6*U_sd || f-errRange > z-6*U_sd
    error('the range of SL error we keep track of needs to be increased to comfortably cover the migrants')
end
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

   
%% prepare matrices to contain information through the evolutionary process (for plotting)

gen = 30000;                                                                   % the number of generations to run the evolution for

pdf_Mig_Gen = zeros(gen, nbins);       pdf_Mig_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after migration
pdf_SS_Gen = zeros(gen, nbins);        pdf_SS_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after social selection (SS)
pdf_Err_Gen = zeros(gen, nbins);       pdf_Err_Gen(1,:) = fFreq;             % gth row: trait distribution in generation g after SL error
fit_Soc_Gen = zeros(gen,nbins);                                              % gth row: fitness values for xM in social selection in generation g
for i = 1: nbins
        fit_Soc_Gen(1,i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
end
Mean_Gen = zeros(gen, 3);              Mean_Gen(1,:) = ones(1,3)*f;          % 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with err


%% focal population distribution's evolutionary process

for g = 2:gen

    g

  %% calculate the norm payoff function before migration happens, since the migrants have no social capital
    fit_Soc = zeros(1,nbins);
    for i = 1: nbins
        fit_Soc(i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
    end
    fit_Soc_Gen(g,:) = fit_Soc;   

  %% migration
    fFreq = (1-m).*fFreq + m.*zFreq;    
    fFreq = fFreq/sum(fFreq);
    pdf_Mig_Gen(g,:) = fFreq; 
    Mean_Gen(g,1) = sum(fFreq.*xM);
    

  %% social selection as coordination
    fFreq = fFreq.*fit_Soc;
    fFreq = fFreq/sum(fFreq);
    pdf_SS_Gen (g,:) = fFreq;
    Mean_Gen(g,2) = sum(fFreq.*xM);
  
  %% social learning without error
  % trying to learn the norm of one cultural parent, fFreq remains the same

    
  %% social learning error
    fFreq = fFreq*Trans;
    fFreq = fFreq/sum(fFreq);
    pdf_Err_Gen(g,:) = fFreq;
    Mean_Gen(g,3) = sum(fFreq.*xM);
    Mean_Gen(g,3) 

    
end



%% figure 1: the norm trait distribution after different stages (migration, norm effects, social learning w/ error) in one generation 
figure; hold on;
g_plot=2;                 % the generation we plot the distributions for
plot(xM, pdf_Mig_Gen(g_plot,:), 'color', [0.9290 0.6940 0.1250], 'LineWidth',1);
plot(xM, pdf_SS_Gen(g_plot,:), 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
plot(xM, pdf_Err_Gen(g_plot,:), 'color', [0 0.4470 0.7410], 'LineStyle','--', 'LineWidth',1);
legend('After migration','After similarity-based norm effects','After social learning with error','FontSize',12)
xlabel('Norm trait value');
ylabel('Distribution density');
xlim([-10 10])
ylim([0 0.0032])
%title(['norm trait distribution after different stages in generation ', num2str(g_plot)]);
ax = gca;
ax.FontSize = 12;

%% figure 2: the norm evolutionary trajectory
figure; hold on;
plot(1:100:401, z*ones(1,5), 'LineStyle','--', 'Marker','*', 'color','black', 'LineWidth',0.2)
plot(1:401, Mean_Gen(1:401,3),'color',[0.4940 0.1840 0.5560],'LineWidth',2)
xlabel('Generation')
ylabel('Mean norm value')
legend('Source population','Focal population','FontSize',12)
xlim([1 400])

%title('norm evolution trajectory in the focal population')
ax = gca;
ax.FontSize = 12;


%% figure 3: the norm evolution in each stage
Mean_each_step = zeros(1,3*g);
for g = 1:gen
    Mean_each_step(1,3*g-2)=Mean_Gen(g,1);     %after migration
    Mean_each_step(1,3*g-1)=Mean_Gen(g,2);     %after social selection
    Mean_each_step(1,3*g)=Mean_Gen(g,3);       %after social learning with error
end

figure;
hold on;
plot(1:3*gen,  Mean_each_step,   'lineWidth',0.01)
xlabel('generation')
ylabel('mean of norm trait distribution')
title('by stage norm evolution trajectory in the focal population')


%% figure 4: population distribution after migration and social selection fitness function
figure;
hold on;
g_plot=150;                 % the generation we plot the distribution and fitness function for

yyaxis left
plot(xM,fit_Soc_Gen(g_plot,:),'Color',[0.8500, 0.3250, 0.0980]);
ylabel('fitness','FontSize',20);

yyaxis right
plot(xM, pdf_Mig_Gen(g_plot,:),'Color','black');
ylim([0 0.005])
ylabel('density','FontSize',20);

line([Mean_Gen(g_plot,1),Mean_Gen(g_plot,1)],[0,0.2],'color','black','LineWidth',1.5)    % after migration distribution mean

%xlim([-10 10])

leg = legend('social selection fitness function','norm distribution density function','FontSize',12);
title('social selection fitness function for a population distribution');
ax = gca;   ax.FontSize = 16;
ax.YAxis(1).Color = [0.8500, 0.3250, 0.0980];
ax.YAxis(2).Color = 'black';
xlabel('trait value','FontSize',20);


