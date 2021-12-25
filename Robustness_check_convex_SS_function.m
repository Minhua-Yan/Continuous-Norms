%   See "Coordination_migration.m" for a description of the file and the parameters

%   The only difference this file has compared to "Coordination_migration.m" 
%      is this file uses a different social selection fitness function: W = exp(-S_inv *sqrt(x_i-x_j)/2)


%% ----------------------------------------------------------------- assign values to parameters

m = 0.1;               
z = 6.12;      f = 0;                                                                
U_sd = 1;   S_inv = 1;    E_sd = 1;



%% ----------------------------------------------------------------- initialize vector of norm trait values

binRange= 30;                                                                
% reminder: this may need to be changed depending on E, S and z
binWid = 0.01;                                                               
nbins = 2*round( binRange /binWid ) +1;      

range_min = -nbins/2*binWid + f;                                              

xU = (1:nbins)*binWid + range_min;                                            
xL = xU - binWid;                                                            
xM = xU - binWid/2;                                                          



%% initialize focal population N(f,V)
V = E_sd^2 + sqrt(E_sd^4+E_sd^2/S_inv);                                     
V_sd = sqrt(V);                                                              
fFreq = normcdf(xU,f,V_sd) - normcdf(xL,f,V_sd);                             
fFreq = fFreq/sum(fFreq);


%% initialize source population N(z,U)
zFreq = normcdf(xU,z,U_sd) - normcdf(xL,z,U_sd);                             
zFreq = zFreq/sum(zFreq);


%% initialize copying error, a normal distribution with mean == 0 and standard deviation E_sd
            
errRange = 6*E_sd;                                                           
nEbins = 2* round( errRange/binWid ) + 1;                                     
minErrRange = -nEbins/2*binWid;                                              
xU_Err = ((1:nEbins)*binWid + minErrRange);                                   
xL_Err = xU_Err - binWid;                                                     

errFreq = normcdf(xU_Err,0,E_sd) - normcdf(xL_Err,0,E_sd);                    
errFreq = errFreq/sum(errFreq);      


%% transition matrix for errorous social learning process

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

gen = 500;                                                                     

pdf_Mig_Gen = zeros(gen, nbins);                                            
pdf_SS_Gen = zeros(gen, nbins);                                              
pdf_Err_Gen = zeros(gen, nbins);                                             

fit_Soc_Gen = zeros(gen,nbins);                                                 
Mean_Gen = zeros(gen, 3);                                                    


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
        fit_Soc(i) = sum(exp(-S_inv *sqrt(dis(i,:))/2).*fFreq);  
        % this is the only difference compared to "Coordination_migration"
    end
    fit_Soc_Gen(g,:) = fit_Soc;         
    
    fFreq = fFreq.*fit_Soc;
    fFreq = fFreq/sum(fFreq);
    pdf_SS_Gen (g,:) = fFreq;
    Mean_Gen(g,2) = sum(fFreq.*xM);
  
  %% social learning without error

    
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


