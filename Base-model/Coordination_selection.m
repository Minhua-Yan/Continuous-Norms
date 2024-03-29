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
%     D_inv: strength of individual natural selection 
%        D_inv is the inverse of D in the base model in the main text
%     Opt: The optimal norm trait for individual natural selection
%     f: initial norm (mean of the initial trait distribution) of the focal population
%     S_inv: strength of social selection as pairwise coordination in the focal population
%        S_inv is the inverse of S in the base model in the main text
%     E_sd: standard deviation of the error in social learning

%% ----------------------------------------------------------------- assign values to parameters

D = 0.2;    Opt = 6.12;      
f = 0;                                                                
S_inv = 1;    E_sd = 1;



%% ----------------------------------------------------------------- initialize vector of norm trait values

% since we cannot maintain normality, we solve for the distribution analytically by approxiamating the distribution of the norm trait and of the error, with bins

binRange= 30;                                                                % the trait range on one side of f we keep track of
% reminder: this may need to be changed depending on E, S and z
binWid = 0.01;                                                               % the width of the bins we cut the distributions into
nbins = 2*round( binRange /binWid ) +1;                                      % make the total number of bins an odd number, so the center of the center bin is exactly "f"

range_min = -nbins/2*binWid + f;                                              % updated lower boundary of the trait range

xU = (1:nbins)*binWid + range_min;                                            % xU is a vector that gives the upper bound of each bin
xL = xU - binWid;                                                            % xL is a vector that gives the lower bound of each bin
xM = xU - binWid/2;                                                          % xM is a vector that gives the center of each bin



%% initialize focal population N(f,V)
V = E_sd^2 + sqrt(E_sd^4+E_sd^2/S_inv);                                      % equilibrium variance with only 1-parent social learning, transission error N~(0,E) and norm selection with strength S_inv, see "eq V with social selection and SL error.pdf"
V_sd = sqrt(V);                                                              % equilibrium standard deviation
fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             % a row, density of the initial island norm traits; the ith element gives the frequency of the ith bin
fFreq = fFreq/sum(fFreq);
f = sum(fFreq.*xM);                                                          % make sure this is close enough to the original value assigned to f



%% initialize copying error, a normal distribution with mean == 0 and standard deviation E_sd
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

   
%% prepare matrices to contain information through the evolutionary process (for plotting)

gen = 500;                                                                   % the number of generations to run the evolution for

pdf_NS_Gen = zeros(gen, nbins);        pdf_NS_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after migration
pdf_SS_Gen = zeros(gen, nbins);        pdf_SS_Gen(1,:) = fFreq;              % gth row: trait distribution in generation g after social selection (SS)
pdf_Err_Gen = zeros(gen, nbins);       pdf_Err_Gen(1,:) = fFreq;             % gth row: trait distribution in generation g after SL error
fit_Soc_Gen = zeros(gen,nbins);                                              % gth row: fitness values for xM in social selection in generation g
for i = 1: nbins
        fit_Soc_Gen(1,i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
end
Mean_Gen = zeros(gen, 3);              Mean_Gen(1,:) = ones(1,3)*f;          % 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with err


%% focal population distribution's evolutionary process

for g = 2:gen

  %% migration
    fFreq = (1-m).*fFreq + m.*zFreq;    
    fFreq = fFreq/sum(fFreq);
    pdf_NS_Gen(g,:) = fFreq; 
    Mean_Gen(g,1) = sum(fFreq.*xM);
    

  %% social selection as coordination
    fit_Soc = zeros(1,nbins);
    for i = 1: nbins
        fit_Soc(i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
    end
    fit_Soc_Gen(g,:) = fit_Soc;         
    
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

    
end

eq_cdf = zeros(1,nbins+1);
for loc = 1:nbins
    eq_cdf(loc+1) = sum(fFreq(1:loc));
     if eq_cdf(loc+1) >1
        eq_cdf(loc+1) =1;
     end
end
    eq_cdf=[eq_cdf,1];
    dist_x = [xM(1)-binWid,xM,xM(end)+binWid];
    
    for pdf_pos = 1:length(eq_cdf)-1
        if eq_cdf(pdf_pos) > eq_cdf(pdf_pos+1)
           msg2 = 'Error occurred: eq_cdf is not a cdf';
           error(msg2)
        end
    end
    pd = makedist('PieceWiselinear','X',dist_x,'FX',eq_cdf);
    sd = sqrt(var(pd));
    

%% figure 1: the norm trait distribution after different stages (migration, social selection, social learning w/ error) in one generation 
figure; hold on;
g_plot=9;                 % the generation we plot the distributions for
plot(xM, pdf_Mig_Gen(g_plot,:), xM, pdf_SS_Gen(g_plot,:),  xM, pdf_Err_Gen(g_plot,:),  'LineWidth',1);
legend('After migration','After social selection','After social learning with error','FontSize',12)
xlabel('Norm trait value');
ylabel('Distribution density');
xlim([-10 10])
%title(['norm trait distribution after different stages in generation ', num2str(g_plot)]);
ax = gca;
ax.FontSize = 16;

%% figure 2: the norm evolutionary trajectory
figure; hold on;
plot(1:gen, Mean_Gen(:,3),'LineWidth',2)
xlabel('Generation')
ylabel('Mean norm value')
line([0,gen],[z,z],'color','black','LineWidth',1)
%title('norm evolution trajectory in the focal population')
ax = gca;
ax.FontSize = 16;


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


