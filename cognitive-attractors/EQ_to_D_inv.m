% % This program plots the equilibrium norm of two focal population, one
% % starting at norm f1, another starting at norm f2, under natural selection
% % of different strengths, when the cognitive process involves attractors
% 
% % the x axis here is natural selection strength, D_inv
% S = 1;
% E_sd = 1;
% Opt = 4;
% p_A = 1.5;
% p_B = 6;
% 
% 
% D_inv_grid_base = 0.1;
% D_inv_vec_base = 0:D_inv_grid_base:5*D_inv_grid_base;
% 
% % D_inv_A1_vec = D_inv_vec_base;
% % D_inv_A2_vec = D_inv_vec_base;
% D_inv_B1_vec = D_inv_vec_base;
% D_inv_B2_vec = D_inv_vec_base;
% 
% 
% % the norm equilibria for different natural selection strength when the population starts with the norm 0
% f1 = 0;
% f2 = 6;
% 
% %% p_A, f1
% 
%     % calculate the equilibrium for each value of I and store them in a vector eq_vec
%     eqA1_vec= zeros(length(D_inv_A1_vec),1);
%     k=1;       % loop counter
%     for D_inv_A1 = D_inv_A1_vec
%         D_inv_A1
%         eqA1_vec(k,1)=Function_CA_Selection_D_inv(D_inv_A1, Opt, f1, S, E_sd, p_A);
%         eqA1_vec(k,1)
%         k=k+1;
%     end
% 
%     
%     % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
%     finer = 0;
%     D_inv_grid_A1 = D_inv_grid_base;
%     while finer < 2 
%         finer = finer +1;
%         D_inv_grid_A1 = D_inv_grid_A1/10;
% 
%         % check whether there is (still) a jump
%         eqA1_i=1;
%         while eqA1_i < length(D_inv_A1_vec)
%             if abs(eqA1_vec(eqA1_i,1) - eqA1_vec(eqA1_i+1,1)) > abs(Opt-f1)/2
%                D_inv_A1_middle_vec = D_inv_A1_vec(eqA1_i):D_inv_grid_A1:D_inv_A1_vec(eqA1_i+1);     
%                D_inv_A1_middle_vec = D_inv_A1_middle_vec(2:end-1);   
%                % calculate the equilibria for the newly added D_inv values in the finer range
%                k=1;
%                eqA1_middle_vec = zeros(length(D_inv_A1_middle_vec),1);
%                for D_inv_A1 = D_inv_A1_middle_vec
%                    D_inv_A1
%                    eqA1_middle_vec(k,1)=Function_CA_Selection_D_inv(D_inv_A1, Opt, f1, S, E_sd, p_A);
%                    eqA1_middle_vec(k,1)
%                    k=k+1;
%                end
%                %n update the D_inv vector and the equilibria vector
%                D_inv_A1_vec = [D_inv_A1_vec(1:eqA1_i), D_inv_A1_middle_vec, D_inv_A1_vec(eqA1_i+1 :end)];
%                eqA1_vec = [eqA1_vec(1:eqA1_i,1); eqA1_middle_vec; eqA1_vec(eqA1_i+1 :end,1)];
%             end
%             eqA1_i = eqA1_i+1;
%         end
%     end
%     
%  
% 
% 
% %% p_A, f2
% 
%     % calculate the equilibrium for each value of I and store them in a vector eq_vec
%     eqA2_vec= zeros(length(D_inv_A2_vec),1);
%     k=1;       % loop counter
%     for D_inv_A2 = D_inv_A2_vec
%         D_inv_A2
%         eqA2_vec(k,1)=Function_CA_Selection_D_inv(D_inv_A2, Opt, f2, S, E_sd, p_A);
%         eqA2_vec(k,1)
%         k=k+1;
%     end
% 
%     
%     % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
%     finer = 0;
%     D_inv_grid_A2 = D_inv_grid_base;
%     while finer < 2 
%         finer = finer +1;
%         D_inv_grid_A2 = D_inv_grid_A2/10;
% 
%         % check whether there is (still) a jump
%         eqA2_i=1;
%         while eqA2_i < length(D_inv_A2_vec)
%             if abs(eqA2_vec(eqA2_i,1) - eqA2_vec(eqA2_i+1,1)) > abs(Opt-f2)/2
%                D_inv_A2_middle_vec = D_inv_A2_vec(eqA2_i):D_inv_grid_A2:D_inv_A2_vec(eqA2_i+1);     
%                D_inv_A2_middle_vec = D_inv_A2_middle_vec(2:end-1);   
%                % calculate the equilibria for the newly added D_inv values in the finer range
%                k=1;
%                eqA2_middle_vec = zeros(length(D_inv_A2_middle_vec),1);
%                for D_inv_A2 = D_inv_A2_middle_vec
%                    D_inv_A2
%                    eqA2_middle_vec(k,1)=Function_CA_Selection_D_inv(D_inv_A2, Opt, f2, S, E_sd, p_A);
%                    eqA2_middle_vec(k,1)
%                    k=k+1;
%                end
%                %n update the D_inv vector and the equilibria vector
%                D_inv_A2_vec = [D_inv_A2_vec(1:eqA2_i), D_inv_A2_middle_vec, D_inv_A2_vec(eqA2_i+1 :end)];
%                eqA2_vec = [eqA2_vec(1:eqA2_i, 1); eqA2_middle_vec; eqA2_vec(eqA2_i+1 :end, 1)];
%             end
%             eqA2_i = eqA2_i+1;
%         end
%     end
% 
% 
% %% p_B, f1
% 
%     % calculate the equilibrium for each value of I and store them in a vector eq_vec
%     eqB1_vec= zeros(length(D_inv_B1_vec),1);
%     k=1;       % loop counter
%     for D_inv_B1 = D_inv_B1_vec
%         D_inv_B1
%         eqB1_vec(k,1)=Function_CA_Selection_D_inv(D_inv_B1, Opt, f1, S, E_sd, p_B);
%         eqB1_vec(k,1)
%         k=k+1;
%     end
% 
%     
%     % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
%     finer = 0;
%     D_inv_grid_B1 = D_inv_grid_base;
%     while finer < 2 
%         finer = finer +1;
%         D_inv_grid_B1 = D_inv_grid_B1/10;
% 
%         % check whether there is (still) a jump
%         eqB1_i=1;
%         while eqB1_i < length(D_inv_B1_vec)
%             if abs(eqB1_vec(eqB1_i,1) - eqB1_vec(eqB1_i+1,1)) > abs(Opt-f1)/2
%                D_inv_B1_middle_vec = D_inv_B1_vec(eqB1_i):D_inv_grid_B1:D_inv_B1_vec(eqB1_i+1);     
%                D_inv_B1_middle_vec = D_inv_B1_middle_vec(2:end-1);   
%                % calculate the equilibria for the newly added D_inv values in the finer range
%                k=1;
%                eqB1_middle_vec = zeros(length(D_inv_B1_middle_vec),1);
%                for D_inv_B1 = D_inv_B1_middle_vec
%                    D_inv_B1
%                    eqB1_middle_vec(k,1)=Function_CA_Selection_D_inv(D_inv_B1, Opt, f1, S, E_sd, p_B);
%                    eqB1_middle_vec(k,1)
%                    k=k+1;
%                end
%                %n update the D_inv vector and the equilibria vector
%                D_inv_B1_vec = [D_inv_B1_vec(1:eqB1_i), D_inv_B1_middle_vec, D_inv_B1_vec(eqB1_i+1 :end)];
%                eqB1_vec = [eqB1_vec(1:eqB1_i,1); eqB1_middle_vec; eqB1_vec(eqB1_i+1 :end,1)];
%             end
%             eqB1_i = eqB1_i+1;
%         end
%     end
%     
%  
% 
% 
% %% p_A, f2
% 
%     % calculate the equilibrium for each value of I and store them in a vector eq_vec
%     eqB2_vec= zeros(length(D_inv_B2_vec),1);
%     k=1;       % loop counter
%     for D_inv_B2 = D_inv_B2_vec
%         D_inv_B2
%         eqB2_vec(k,1)=Function_CA_Selection_D_inv(D_inv_B2, Opt, f2, S, E_sd, p_B);
%         eqB2_vec(k,1)
%         k=k+1;
%     end
% 
%     
%     % find the region of jump in the equilibria vector to segment more finely for a more accurate plot
%     finer = 0;
%     D_inv_grid_B2 = D_inv_grid_base;
%     while finer < 2 
%         finer = finer +1;
%         D_inv_grid_B2 = D_inv_grid_B2/10;
% 
%         % check whether there is (still) a jump
%         eqB2_i=1;
%         while eqB2_i < length(D_inv_B2_vec)
%             if abs(eqB2_vec(eqB2_i,1) - eqB2_vec(eqB2_i+1,1)) > abs(Opt-f2)/2
%                D_inv_B2_middle_vec = D_inv_B2_vec(eqB2_i):D_inv_grid_B2:D_inv_B2_vec(eqB2_i+1);     
%                D_inv_B2_middle_vec = D_inv_B2_middle_vec(2:end-1);   
%                % calculate the equilibria for the newly added D_inv values in the finer range
%                k=1;
%                eqB2_middle_vec = zeros(length(D_inv_B2_middle_vec),1);
%                for D_inv_B2 = D_inv_B2_middle_vec
%                    D_inv_B2
%                    eqB2_middle_vec(k,1)=Function_CA_Selection_D_inv(D_inv_B2, Opt, f2, S, E_sd, p_B);
%                    eqB2_middle_vec(k,1)
%                    k=k+1;
%                end
%                %n update the D_inv vector and the equilibria vector
%                D_inv_B2_vec = [D_inv_B2_vec(1:eqB2_i), D_inv_B2_middle_vec, D_inv_B2_vec(eqB2_i+1 :end)];
%                eqB2_vec = [eqB2_vec(1:eqB2_i, 1); eqB2_middle_vec; eqB2_vec(eqB2_i+1 :end, 1)];
%             end
%             eqB2_i = eqB2_i+1;
%         end
%     end
% 
% 
% 
% save eq-D-CA.mat D_inv_A1_vec eqA1_vec D_inv_A2_vec eqA2_vec D_inv_B1_vec eqB1_vec D_inv_B2_vec eqB2_vec 


%% figure 1: equilibrium norm as mean of the equilibrium distribution of true trait values - migration rate
figure;
hold on;


subplot(1, 2, 1)
hold on;

    %% subplot 1
plot_Opt = plot(D_inv_vec_base, Opt*ones(1,length(D_inv_vec_base)),'LineStyle','--','Marker','*','color','black','LineWidth',0.2);

    %% plot the equilibrium-m for cbn_A, f1

    k=1;
    while k < length(D_inv_A1_vec)
        if abs(eqA1_vec(k,1) - eqA1_vec(k+1,1)) > abs(Opt-f1)/2               % check whether the equilibrium vector still has a jump in it
            D_A1_vec_1 = D_inv_A1_vec(1:k);      eqA1_vec_1 = eqA1_vec(1:k,:);
            D_A1_vec_2 = D_inv_A1_vec(k+1:end);  eqA1_vec_2 = eqA1_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(D_inv_A1_vec)                                                 % if there is no jump, we plot the equilibrium vector as continuous by connecting all the dots
        plotA1 = plot(D_inv_A1_vec, eqA1_vec(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    elseif length(D_A1_vec_1) == 1                                           % if there is a jump and the initial equilibrium is only maintained under m=0, plot the rest as continuous and m=0 as a separate dot
        scatter(D_A1_vec_1(1), eqA1_vec_1(1,1), 50, [0.3010 0.7450 0.9330], 'filled');   plotA1= plot(D_A1_vec_2, eqA1_vec_2(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    else                                                                    % if there is a jump in the middle, plot the two parts separately as continuous but do not connect them
        plot(D_A1_vec_1, eqA1_vec_1(:,1),'color',[0.3010 0.7450 0.9330], 'LineWidth',3); plotA1= plot(D_A1_vec_2, eqA1_vec_2(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    end 

    %% plot the equilibrium-m for step_A, f2
    k=1;
    while k < length(D_inv_A2_vec)
        if abs(eqA2_vec(k,1) - eqA2_vec(k+1,1)) > abs(Opt-f2)/2
            D_A2_vec_1 = D_inv_A2_vec(1:k);     eqA2_vec_1 = eqA2_vec(1:k,:);
            D_A2_vec_2 = D_inv_A2_vec(k+1:end); eqA2_vec_2 = eqA2_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(D_inv_A2_vec)
        plotA2= plot(D_inv_A2_vec, eqA2_vec(:,1), 'color',[1,0,0], 'LineWidth',1);
    elseif length(D_A2_vec_1) == 1
        scatter(D_A2_vec_1(1), eqA2_vec_1(1,1), 50, [1,0,0], 'filled');    plotA2= plot(D_A2_vec_2, eqA2_vec_2(:,1), 'color',[1,0,0], 'LineWidth',1);
    else
        plot(D_A2_vec_1, eqA2_vec_1(:,1), 'color',[1,0,0], 'LineWidth',1); plotA2= plot(D_A2_vec_2, eqA2_vec_2(:,1), 'color',[1,0,0], 'LineWidth',1);
    end 


xlim([D_inv_A1_vec(1) D_inv_A1_vec(end)])    
ylim([f1 f2])
lgd = legend([plot_Opt plotA1 plotA2], 'Direct effect optimum', ['initial norm = ',num2str(f1)],  ['initial norm = ',num2str(f2)]);
lgd.FontSize = 11;
xlabel('Direct effect strength (^{1}/_{D})')
ylabel('Equilibrium mean norm')
title(['{\itp}=',num2str(p_A), ', {\itS}=',num2str(S), ', {\itE}=',num2str(E_sd^2), ', {\theta}=',num2str(Opt)]);
ax = gca;   ax.FontSize = 12;




subplot(1, 2, 2)
hold on;

plot_Opt = plot(D_inv_vec_base, Opt*ones(1,length(D_inv_vec_base)),'LineStyle','--','Marker','*','color','black','LineWidth',0.2);


    %% plot the equilibrium-m for step_B, f1
    k=1;
    while k < length(D_inv_B1_vec)
        if abs(eqB1_vec(k,1) - eqB1_vec(k+1,1)) > abs(Opt-f1)/2
            D_inv_B1_vec_1 = D_inv_B1_vec(1:k);     eqB1_vec_1 = eqB1_vec(1:k,:);
            D_inv_B1_vec_2 = D_inv_B1_vec(k+1:end); eqB1_vec_2 = eqB1_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(D_inv_B1_vec)
        plotB1 = plot(D_inv_B1_vec, eqB1_vec(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    elseif length(D_inv_B1_vec_1) == 1
        scatter(D_inv_B1_vec_1(1), eqB1_vec_1(1,1), 50, [0.3010 0.7450 0.9330], 'filled');    plotB1= plot(D_inv_B1_vec_2, eqB1_vec_2(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    else
        plot(D_inv_B1_vec_1, eqB1_vec_1(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3); plotB1= plot(D_inv_B1_vec_2, eqB1_vec_2(:,1), 'color',[0.3010 0.7450 0.9330], 'LineWidth',3);
    end 
    

    %% plot the equilibrium-m for cbn_B, f2
    k=1;
    while k < length(D_inv_B2_vec)
        if abs(eqB2_vec(k,1) - eqB2_vec(k+1,1)) > abs(Opt-f2)/2
            D_inv_B2_vec_1 = D_inv_B2_vec(1:k);     eqB2_vec_1 = eqB2_vec(1:k,:);
            D_inv_B2_vec_2 = D_inv_B2_vec(k+1:end); eqB2_vec_2 = eqB2_vec(k+1:end,:);
            break;
        else
            k = k+1;
        end
    end
    
    if k == length(D_inv_B2_vec)
        plotB2= plot(D_inv_B2_vec, eqB2_vec(:,1), 'color',[1,0,0], 'LineWidth',1);
    elseif length(D_inv_B2_vec_1) == 1
        scatter(D_inv_B2_vec_1(1), eqB2_vec_1(1,1), 50, [1,0,0], 'filled');       plotB2= plot(D_inv_B2_vec_2, eqB2_vec_2(:,1), 'color',[1,0,0], 'LineWidth',1);
    else
        plot(D_inv_B2_vec_1, eqB2_vec_1(:,1), 'color',[1,0,0], 'LineWidth',1);    plotB2= plot(D_inv_B2_vec_2, eqB2_vec_2(:,1), 'color',[1,0,0], 'LineWidth',1);
    end 


xlim([D_inv_B1_vec(1) D_inv_B1_vec(end)])    
ylim([f1 f2])
%legend([plot_z plotB1 plotB2], 'migrant norm', ['initial norm = ',num2str(f1)],  ['initial norm = ',num2str(f2)]);
xlabel('Direct effect strength (^{1}/_{D})')
ylabel('Equilibrium mean norm')
title(['{\itp}=',num2str(p_B), ', {\itS}=',num2str(S), ', {\itE}=',num2str(E_sd^2), ', {\theta}=',num2str(Opt)]);
%xlim(0.01*[start_width,end_width])
ax = gca;     ax.FontSize = 12;     
