%   See "Coordination_migration.m" for a description of the file and the parameters

%   The only difference this file has compared to "Coordination_migration.m" 
%     is that there are 2 social learning parents
%     and the child tries to learn the average of their norm traits. 
%     If the average is not at the center of a base bin,
%     half of the children acquire the trait of the bin on the left of the average, 
%     and the other half of the children acquire the trait of the bin on the right of the average.

%% ----------------------------------------------------------------- assign values to parameters

m = 0.1;               
z = 6.12;      f = 0;                                                                
U_sd = 0.001;   S_inv = 1;    E_sd = 1;



%% ----------------------------------------------------------------- initialize vector of norm trait values

binRange= 30;                                                                
% reminder: this may need to be changed depending on E, S and z
binWid = 0.01;                                                               
nbins = 2*round( binRange /binWid ) +1;       

range_min = -nbins/2*binWid + f;         

xU = (1:nbins)*binWid + range_min;                                           
xL = xU - binWid;                                                           
xM = xU - binWid/2;                                                          




%% initialize source population N(z,U)
zFreq = normcdf(xU,z,U_sd) - normcdf(xL,z,U_sd);                             % a row, density of the continent trait values between [minX +(i-1)*binWid, minX +i*binWid]
zFreq = zFreq/sum(zFreq);


%% initialize copying error, a normal distribution with mean == 0 and standard deviation E_sd
% error is also approxiamited with bins of width "binWid"            
errRange = 6*E_sd;                                                           % the error range on one side of 0 we keep track of
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
V = E_sd^2 + sqrt(E_sd^4+E_sd^2/S_inv);                                      % equilibrium variance with only 1-parent social learning, transission error N~(0,E) and norm selection with strength S_inv, see "eq V with social selection and SL error.pdf"
V_sd = sqrt(V);                                                              % fake equilibrium standard deviation
fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             % fake population density function
fFreq = fFreq/sum(fFreq);
f = sum(fFreq.*xM);                                                          % make sure this is close enough to the original value assigned to f


% to get the equilibrium variance under 2-parent social learning, I run the
% population with only social selection and SL error until the variance stabalizes

gen_V = 1000;

V_vec = zeros(1,gen_V);
V_vec(1) = V;

for g_V = 2:gen_V
    g_V    

  %% social selection as coordination
    fit_Soc_V = zeros(1,nbins);
    for i = 1: nbins
        fit_Soc_V(i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
    end    
    
    fFreq = fFreq.*fit_Soc_V;
    fFreq = fFreq/sum(fFreq);
  
  %% social learning without error
    
  % trying to learn the mean of two cultural parents
    fCombM_V = transpose(fFreq)*fFreq;                                        %a matrix containing probabilities of different combinations                                                 %the number of averages of combinations
    fCombMean_V = zeros(1,2*nbins-1);                                         %a rwo containing the probabilities of different averages of combinations
    for i = 1:nbins
        for j = 1:i
           fCombMean_V(1,i) = fCombMean_V(1,i) +fCombM_V(j,i+1-j); 
        end
    end
    for i = nbins+1 : 2*nbins-1
        for j = i+1-nbins :nbins
            fCombMean_V(1,i) = fCombMean_V(1,i) +fCombM_V(j,i+1-j);
        end
    end
    fFreq = fCombMean_V(1:2:end) + [1/2*fCombMean_V(2:2:end),0]+[0,1/2*fCombMean_V(2:2:end)];   
      
    
   %% social learning error
    fFreq = fFreq*Trans;
    fFreq = fFreq/sum(fFreq);
    
   %% calculate the distribution variance and end the loop if it has stablized
    fFreq_Cdf = zeros(1,nbins+1);
    for loc = 1:nbins
         fFreq_Cdf(loc+1) = sum(fFreq(1:loc));
         if fFreq_Cdf(loc+1) >1
            fFreq_Cdf(loc+1) =1;
         end
    end
    fFreq_Cdf=[fFreq_Cdf,1];
    dist_x = [xM(1)-binWid,xM,xM(end)+binWid];
    
    for pdf_pos = 1:length(fFreq_Cdf)-1
        if fFreq_Cdf(pdf_pos) > fFreq_Cdf(pdf_pos+1)
           msg2 = 'Error occurred: eq_cdf is not a cdf';
           error(msg2)
        end
    end
    pd = makedist('PieceWiselinear','X',dist_x,'FX',fFreq_Cdf);
    V_vec(1,g_V) = var(pd);
    
    if g_V > 10 && abs(V_vec(1,g_V)-V_vec(1,g_V-1)) < 0.000001 && abs(V_vec(1,g_V)-V_vec(1,g_V-10)) < 0.00001
        V = V_vec(1,g_V);
        break;
    end
end

% calculate the actual focal population distribution before migration starts
V_sd = sqrt(V);                                                              % actual equilibrium standard deviation
fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             % actual focal population distribution density function
fFreq = fFreq/sum(fFreq);
f = sum(fFreq.*xM);                                                          % make sure this is close enough to the original value assigned to f




   
%% prepare matrices to contain information through the evolutionary process (for plotting)

gen = 500;                                                                     % the number of generations to run the evolution for

pdf_Mig_Gen = zeros(gen, nbins);        pdf_Mig_Gen(1,:) = fFreq;            % gth row: trait distribution in generation g after migration
pdf_SS_Gen = zeros(gen, nbins);         pdf_SS_Gen(1,:) = fFreq;             % gth row: trait distribution in generation g after social selection (SS)
pdf_SL_Gen = zeros(gen, nbins);         pdf_SL_Gen(1,:) = fFreq;             % gth row: trait distribution in generation g after social learning (SL) w/o error
pdf_Err_Gen = zeros(gen, nbins);        pdf_Err_Gen(1,:) = fFreq;            % gth row: trait distribution in generation g after SL error

fit_Soc_Gen = zeros(gen,nbins);                                                  % gth row: fitness values for xM in social selection in generation g
for i = 1: nbins
        fit_Soc_Gen(1,i) = sum(exp(-S_inv *dis(i,:).*dis(i,:)/2).*fFreq);  
end
Mean_Gen = zeros(gen, 4);                                                    % 1 column: mean after Mig, 2nd column: mean after SS, 3rd column: mean after SL with 2 parents, 4th column: mean after SL error


%% focal population distribution's evolutionary process

for g = 2:gen

  %% migration
    fFreq = (1-m).*fFreq + m.*zFreq;    
    fFreq = fFreq/sum(fFreq);
    pdf_Mig_Gen(g,:) = fFreq; 
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
    
  % trying to learn the mean of two cultural parents
    fCombM = transpose(fFreq)*fFreq;                                        %a matrix containing probabilities of different combinations                                                 %the number of averages of combinations
    fCombMean = zeros(1,2*nbins-1);                                         %a rwo containing the probabilities of different averages of combinations
    for i = 1:nbins
        for j = 1:i
           fCombMean(1,i) = fCombMean(1,i) +fCombM(j,i+1-j); 
        end
    end
    for i = nbins+1 : 2*nbins-1
        for j = i+1-nbins :nbins
            fCombMean(1,i) = fCombMean(1,i) +fCombM(j,i+1-j);
        end
    end
    fFreq = fCombMean(1:2:end) + [1/2*fCombMean(2:2:end),0]+[0,1/2*fCombMean(2:2:end)];   
    pdf_SL_Gen(g,:)=fFreq; 
    Mean_Gen(g,3) = sum(fFreq.*xM);    
      
    
    
    % social learning error
    fFreq = fFreq*Trans;
    fFreq = fFreq/sum(fFreq);
    pdf_Err_Gen(g,:) = fFreq;
    Mean_Gen(g,4) = sum(fFreq.*xM);
    
end

    

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
Mean_each_step = zeros(1,4*g);
for g = 1:gen
    Mean_each_step(1,4*g-3)=Mean_Gen(g,1);     %after migration
    Mean_each_step(1,4*g-2)=Mean_Gen(g,2);     %after social selection
    Mean_each_step(1,4*g-1)=Mean_Gen(g,3);     %after social learning w/o error
    Mean_each_step(1,4*g)=Mean_Gen(g,4);         %after social learning error
end

figure;
hold on;
plot(1:4*gen,  Mean_each_step,   'lineWidth',0.01)
xlabel('generation')
ylabel('mean of norm trait distribution')
title('norm evolution trajectory in the focal population given the stages')


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


