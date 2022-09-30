% This program plots the equilibrium norm of two focal population, one
% starting at norm f1, another starting at norm f2, under migration from the
% same continent population with a norm z, against the migration rate

% the factor of interest here is category width

% Here, the categories that people are grouped into in social
% judgment (norm coordination) stay constant through out the evolutionary
% history

S_inv = 0.1;
E_sd = 1;
binWid = 0.01;
cbn_base = 141;
cbn_A = cbn_base;
cbn_B = 3*cbn_base;                     % this makes it easier to select a migrant norm that is at the center of some category for all category width that we investigate
z = cbn_B*binWid;                     % make sure the migrant norm is at the center of a category for all category width that we investigate

m_grid_base = 0.01;                     % this is how finely we plot the equilibria against m
mA_vec_base = 0:m_grid_base:10*m_grid_base;
mB_vec_base = 0:m_grid_base:10*m_grid_base;
% A has cbn_A as the cbn; B has cbn_B as the cbn; each has two eq-m plots
%     for two focal population starting with different initial norms, f1 and f2
mA1_vec = mA_vec_base;
mA2_vec = mA_vec_base;
mB1_vec = mB_vec_base;
mB2_vec = mB_vec_base;

f1=0;
f2=2*z;                                   % the two populations are symmetrical around the source population

%% cbn_A, f1
    resultA1_vec= zeros(length(mA1_vec),4);
    k=1;       % loop counter
    for mA1 = mA1_vec
        mA1
        resultA1_vec(k,:) = Function_discretization_migration_fast(mA1,S_inv,E_sd,cbn_A,f1,z,binWid);
        k=k+1;    
    end
    
    % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
    finer = 0;
    m_grid_A1 = m_grid_base;
    while finer <  2                  % this is how many times we refine the m vector before giving up and plotting the eq-m graph as discontinuous
                                      % note that when m is small, binWid may need to be adjusted smaller to produce an accurate equilibrium, 
                                      % so binWid and m_grid_base restricts how big finer can be
                               
                                     
        finer = finer +1;
        m_grid_A1 = m_grid_A1/10;
                
        % check whether there is (still) a jump
        eqA1_i=1;
        while eqA1_i < length(mA1_vec)
           if abs(resultA1_vec(eqA1_i,1) - resultA1_vec((eqA1_i+1),1)) > abs(z-f1)/2
               mA1_middle_vec = mA1_vec(eqA1_i):m_grid_A1:mA1_vec(eqA1_i+1); 
               mA1_middle_vec = mA1_middle_vec(2:end-1);

               % calculate the equilibria for the newly added m values in the finer range
               k=1;
               resultA1_middle_vec= zeros(length(mA1_middle_vec),4);

               for mA1 = mA1_middle_vec
                   mA1
                   resultA1_middle_vec(k,:)=Function_discretization_migration_fast(mA1,S_inv,E_sd,cbn_A,f1,z,binWid);
                   k=k+1;
               end
               % update the m vector and the equilibrium vector
               mA1_vec = [mA1_vec(1: eqA1_i),mA1_middle_vec,mA1_vec(eqA1_i+1:end)];% this is still the original m1_vec if there is no jump
               resultA1_vec = [resultA1_vec(1:eqA1_i,:); resultA1_middle_vec; resultA1_vec(eqA1_i+1:end,:)]; 

           end
           eqA1_i = eqA1_i+1;
        end
    end
    
    


    %% cbn_A, f2
    resultA2_vec= zeros(length(mA2_vec),4); 
    k=1;       % loop counter
    for mA2 = mA2_vec
        mA2
        resultA2_vec(k,:)=Function_discretization_migration_fast(mA2,S_inv,E_sd,cbn_A,f2,z,binWid);
        k=k+1;
    end
    
    % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
    finer = 0;
    m_grid_A2 = m_grid_base;
    while finer <  2                  % this is how many times we refine the m vector before giving up and plotting the eq-m graph as discontinuous
        finer = finer +1;
        m_grid_A2 = m_grid_A2/10;
               
        % check whether there is still a jump
        eqA2_i=1;
        while eqA2_i < length(mA2_vec)
           if abs(resultA2_vec(eqA2_i,1) - resultA2_vec((eqA2_i+1),1)) > abs(z-f2)/2
               mA2_middle_vec = mA2_vec(eqA2_i):m_grid_A2:mA2_vec(eqA2_i+1); 
               mA2_middle_vec = mA2_middle_vec(2:end-1);

               % calculate the equilibria for the newly added m values during the finer range
               k=1;
               resultA2_middle_vec= zeros(length(mA2_middle_vec),4);

               for mA2 = mA2_middle_vec
                   mA2
                   resultA2_middle_vec(k,:)=Function_discretization_migration_fast(mA2,S_inv,E_sd,cbn_A,f2,z,binWid);
                   k=k+1;
               end
               % update the m vector and the equilibrium vector
               mA2_vec = [mA2_vec(1: eqA2_i),mA2_middle_vec,mA2_vec(eqA2_i+1:end)];% this is still the original m1_vec if there is no jump
               resultA2_vec = [resultA2_vec(1:eqA2_i,:); resultA2_middle_vec; resultA2_vec(eqA2_i+1:end,:)]; 
           end
           eqA2_i = eqA2_i+1;
        end
    end
    


%% cbn_B, f1
    resultB1_vec= zeros(length(mB1_vec),4);
    k=1;       % loop counter
    for mB1 = mB1_vec
        mB1
        resultB1_vec(k,:)=Function_discretization_migration_fast(mB1,S_inv,E_sd,cbn_B,f1,z,binWid);
        k=k+1;    
    end
    
    % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
    finer = 0;
    m_grid_B1 = m_grid_base;
    while finer <  2                  % this is how many times we refine the m vector before giving up and plotting the eq-m graph as discontinuous
                                      % note that when m is small, binWid may need to be adjusted smaller to produce an accurate equilibrium, 
                                      % so binWid and m_grid_base restricts how big finer can be
                               
                                     
        finer = finer +1;
        m_grid_B1 = m_grid_B1/10;
                
        % check whether there is (still) a jump
        eqB1_i=1;
        while eqB1_i < length(mB1_vec)
           if abs(resultB1_vec(eqB1_i,1) - resultB1_vec((eqB1_i+1),1)) > abs(z-f1)/2
               mB1_middle_vec = mB1_vec(eqB1_i):m_grid_B1:mB1_vec(eqB1_i+1); 
               mB1_middle_vec = mB1_middle_vec(2:end-1);

               % calculate the equilibria for the newly added m values in the finer range
               k=1;
               resultB1_middle_vec= zeros(length(mB1_middle_vec),4);

               for mB1 = mB1_middle_vec
                   mB1
                   resultB1_middle_vec(k,:)=Function_discretization_migration_fast(mB1,S_inv,E_sd,cbn_B,f1,z,binWid);
                   k=k+1;
               end
               % update the m vector and the equilibrium vector
               mB1_vec = [mB1_vec(1: eqB1_i),mB1_middle_vec,mB1_vec(eqB1_i+1:end)];% this is still the original m1_vec if there is no jump
               resultB1_vec = [resultB1_vec(1:eqB1_i,:); resultB1_middle_vec; resultB1_vec(eqB1_i+1:end,:)];

           end
           eqB1_i = eqB1_i+1;
        end
    end
    
    


    %% cbn_B, f2
    resultB2_vec= zeros(length(mB2_vec),4);
    k=1;       % loop counter
    for mB2 = mB2_vec
        mB2
        resultB2_vec(k,:)=Function_discretization_migration_fast(mB2,S_inv,E_sd,cbn_B,f2,z,binWid);
        k=k+1;
    end
    
    % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
    finer = 0;
    m_grid_B2 = m_grid_base;
    while finer <  2                  % this is how many times we refine the m vector before giving up and plotting the eq-m graph as discontinuous
        finer = finer +1;
        m_grid_B2 = m_grid_B2/10;
               
        % check whether there is still a jump
        eqB2_i=1;
        while eqB2_i < length(mB2_vec)
           if abs(resultB2_vec(eqB2_i,1) - resultB2_vec((eqB2_i+1),1)) > abs(z-f2)/2
               mB2_middle_vec = mB2_vec(eqB2_i):m_grid_B2:mB2_vec(eqB2_i+1); 
               mB2_middle_vec = mB2_middle_vec(2:end-1);

               % calculate the equilibria for the newly added m values during the finer range
               k=1;
               resultB2_middle_vec= zeros(length(mB2_middle_vec),4);

               for mB2 = mB2_middle_vec
                   mB2
                   resultB2_middle_vec(k,:)=Function_discretization_migration_fast(mB2,S_inv,E_sd,cbn_B,f2,z,binWid);
                   k=k+1;
               end
               % update the m vector and the equilibrium vector
               mB2_vec = [mB2_vec(1: eqB2_i),mB2_middle_vec,mB2_vec(eqB2_i+1:end)];% this is still the original m1_vec if there is no jump
               resultB2_vec = [resultB2_vec(1:eqB2_i,:); resultB2_middle_vec; resultB2_vec(eqB2_i+1:end,:)];
           end
           eqB2_i = eqB2_i+1;
        end
    end

save eq-m.mat mA1_vec resultA1_vec mA2_vec resultA2_vec mB1_vec resultB1_vec mB2_vec resultB2_vec 
    
eqA1_vec = resultA1_vec(:,1);
eqA2_vec = resultA2_vec(:,1);
eqB1_vec = resultB1_vec(:,1);
eqB2_vec = resultB2_vec(:,1);



%% figure 1: equilibrium norm as mean of the equilibrium distribution of true trait values - migration rate
figure;
hold on;
subplot(1, 2, 1)
hold on;


plot_z = line([mA1_vec(1),mA1_vec(end)],[z,z],'LineStyle','--','color','black','LineWidth',0.2);

    %% plot the equilibrium-m for cbn_A, f1
    k=1;
    while k < length(mA1_vec)
        if abs(eqA1_vec(k,1) - eqA1_vec(k+1,1)) > abs(z-f1)/2               % check whether the equilibrium vector still has a jump in it
            mA1_vec_1 = mA1_vec(1:k);      eqA1_vec_1 = eqA1_vec(1:k,:);
            mA1_vec_2 = mA1_vec(k+1:end);  eqA1_vec_2 = eqA1_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(mA1_vec)                                                 % if there is no jump, we plot the equilibrium vector as continuous by connecting all the dots
        plotA1 = plot(mA1_vec, eqA1_vec(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    elseif length(mA1_vec_1) == 1                                           % if there is a jump and the initial equilibrium is only maintained under m=0, plot the rest as continuous and m=0 as a separate dot
        scatter(mA1_vec_1(1), eqA1_vec_1(1,1), 50, [0, 0.4470, 0.7410], 'filled');   plotA1= plot(mA1_vec_2, eqA1_vec_2(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    else                                                                    % if there is a jump in the middle, plot the two parts separately as continuous but do not connect them
        plot(mA1_vec_1, eqA1_vec_1(:,1),'color',[0, 0.4470, 0.7410], 'LineWidth',3); plotA1= plot(mA1_vec_2, eqA1_vec_2(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    end 

    %% plot the equilibrium-m for step_A, f2
    k=1;
    while k < length(mA2_vec)
        if abs(eqA2_vec(k,1) - eqA2_vec(k+1,1)) > abs(z-f2)/2
            mA2_vec_1 = mA2_vec(1:k);     eqA2_vec_1 = eqA2_vec(1:k,:);
            mA2_vec_2 = mA2_vec(k+1:end); eqA2_vec_2 = eqA2_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(mA2_vec)
        plotA2= plot(mA2_vec, eqA2_vec(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    elseif length(mA2_vec_1) == 1
        scatter(mA2_vec_1(1), eqA2_vec_1(1,1), 50, [0.8500, 0.3250, 0.0980], 'filled');    plotA2= plot(mA2_vec_2, eqA2_vec_2(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    else
        plot(mA2_vec_1, eqA2_vec_1(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1); plotA2= plot(mA2_vec_2, eqA2_vec_2(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    end 


xlim([mA1_vec(1) mA1_vec(end)])    
ylim([f1 f2])
lgd = legend([plot_z plotA1 plotA2], 'migrant norm', ['initial norm =',num2str(f1)],  ['initial norm =',num2str(f2)]);
lgd.FontSize = 11;
xlabel('{\itm}: migration rate')
ylabel('equilibrium  norm')
title(['{\itc}=',num2str(0.01*cbn_A), ', {\itS}=',num2str(1/S_inv), ', {\itE}=',num2str(E_sd^2), ', {\itz}=',num2str(z)]);
ax = gca;   ax.FontSize = 12;




subplot(1, 2, 2)
hold on;

plot_z = line([mB1_vec(1),mB1_vec(end)],[z,z],'LineStyle','--','color','black','LineWidth',0.2);


    %% plot the equilibrium-m for step_B, f1
    k=1;
    while k < length(mB1_vec)
        if abs(eqB1_vec(k,1) - eqB1_vec(k+1,1)) > abs(z-f1)/2
            mB1_vec_1 = mB1_vec(1:k);     eqB1_vec_1 = eqB1_vec(1:k,:);
            mB1_vec_2 = mB1_vec(k+1:end); eqB1_vec_2 = eqB1_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(mB1_vec)
        plotB1 = plot(mB1_vec, eqB1_vec(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    elseif length(mB1_vec_1) == 1
        scatter(mB1_vec_1(1), eqB1_vec_1(1,1), 50, [0, 0.4470, 0.7410], 'filled');    plotB1= plot(mB1_vec_2, eqB1_vec_2(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    else
        plot(mB1_vec_1, eqB1_vec_1(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3); plotB1= plot(mB1_vec_2, eqB1_vec_2(:,1), 'color',[0, 0.4470, 0.7410], 'LineWidth',3);
    end 
    

    %% plot the equilibrium-m for cbn_B, f2
    k=1;
    while k < length(mB2_vec)
        if abs(eqB2_vec(k,1) - eqB2_vec(k+1,1)) > abs(z-f2)/2
            mB2_vec_1 = mB2_vec(1:k);     eqB2_vec_1 = eqB2_vec(1:k,:);
            mB2_vec_2 = mB2_vec(k+1:end); eqB2_vec_2 = eqB2_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(mB2_vec)
        plotB2= plot(mB2_vec, eqB2_vec(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    elseif length(mB2_vec_1) == 1
        scatter(mB2_vec_1(1), eqB2_vec_1(1,1), 50, [0.8500, 0.3250, 0.0980], 'filled');       plotB2= plot(mB2_vec_2, eqB2_vec_2(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    else
        plot(mB2_vec_1, eqB2_vec_1(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);  plotB2= plot(mB2_vec_2, eqB2_vec_2(:,1), 'color',[0.8500, 0.3250, 0.0980], 'LineWidth',1);
    end 


xlim([mB1_vec(1) mB1_vec(end)])    
ylim([f1 f2])
%legend([plot_z plotB1 plotB2], 'migrant norm', ['initial norm =',num2str(f1)],  ['initial norm =',num2str(f2)]);
xlabel('{\itm}: migration rate')
ylabel('equilibrium  norm')
title(['{\itc}=',num2str(0.01*cbn_B), ', {\itS}=',num2str(1/S_inv), ', {\itE}=',num2str(E_sd^2), ', {\itz}=',num2str(z)]);
%xlim(0.01*[start_width,end_width])
ax = gca;     ax.FontSize = 12;     
