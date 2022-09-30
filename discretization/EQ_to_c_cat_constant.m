% This program plots the equilibrium norm of two focal population, one
% starting at norm f1, another starting at norm f2, under migration with migration rate m from the
% same continent population with a norm z, as a function of different category width values

% Here the categories are symmetrical around the focal population's original norm and do not evolve



S_inv=0.2;
E_sd=1;
z=4.23;
m=0.05;
binWid = 0.01;

cbn_grid_base = 10;                      % this is how finely we plot the equilibria against cbn
cbn_vec_base = 1:cbn_grid_base:80*cbn_grid_base+1;

cbn_vec = cbn_vec_base;

f=0;


%% f1
    result_vec= zeros(length(cbn_vec),4);
    k=1;       % loop counter
    for cbn = cbn_vec
        cbn
        result_vec(k,:) =Function_discretization_migration_fast(m,S_inv,E_sd,cbn,f,z,binWid);
        k=k+1;    
    end

    % find the region of jump in the equilibria vector to segment more finely for a more accurate plot

    cbn_grid_finer = max(cbn_grid_base/10 , 2);

                
    % check whether there is a jump
    eq_i=1;
    while eq_i < length(cbn_vec)
       if abs(result_vec(eq_i,1) - result_vec((eq_i+1),1)) > abs(z-f)/2
           cbn_middle_vec = cbn_vec(eq_i):cbn_grid_finer:cbn_vec(eq_i+1); 
           cbn_middle_vec = cbn_middle_vec(2:end-1);

           % calculate the equilibria for the newly added cbn values during the finer range
           k=1; 
        
           result_middle_vec= zeros(length(cbn_middle_vec),4);

           for cbn = cbn_middle_vec
               cbn
               result_middle_vec(k,:)=Function_discretization_migration_fast(m,S_inv,E_sd,cbn,f,z,binWid);
               k=k+1;
           end
           % update the m vector and the equilibrium vector
           cbn_vec = [cbn_vec(1:eq_i), cbn_middle_vec, cbn_vec(eq_i+1:end)];% this is still the original m1_vec if there is no jump
           result_vec = [result_vec(1:eq_i,:); result_middle_vec; result_vec(eq_i+1:end,:)]; 
       end
       eq_i = eq_i+1;
    end


save eq-c.mat cbn_vec result_vec

%% plots

eq_vec = result_vec(:,1);
mostFitCat_vec = result_vec(:,2);
SS_effect_gen2_vec = result_vec(:,3);
SS_effect_max_vec = result_vec(:,4);

sz=10;
plot_num = length(cbn_vec);



c=jet(3*plot_num);
cmap = colormap(c(1:3:end,:));


%% figure 1: equilibrium mean of true norm values - category width
figure;
plot(cbn_vec*binWid,eq_vec(:,1));
xlabel('category width'); ylabel('equilibrium  mean  norm')
%title(['E=', num2str(E_sd^2), ', S=', num2str(1/S_inv), ', m=', num2str(m), ', z=', num2str(z), ', f=', num2str(f)])

% figure 2: equilibrium mean of true norm values - effect of SS on cat on changing mean norm in gen 2
figure; 
scatter(SS_effect_gen2_vec,eq_vec,c);
xlabel('effect of SS on cat on changing mean norm in gen 2'); ylabel('equilibrium mean of true norm values')
%title(['E=', num2str(E_sd^2), ', S=', num2str(1/S_inv), ', m=', num2str(m), ', z=', num2str(z), ', f=', num2str(f)])




%% figure 3: equilibrium mean of true norm values - max effect of SS on cat on changing mean norm
figure; 
scatter(SS_effect_max_vec,eq_vec,sz,cmap);
xlabel('max mean norm change due to discretized social effects'); ylabel('equilibrium mean norm')
%title(['E=', num2str(E_sd^2), ', S=', num2str(1/S_inv), ', m=', num2str(m), ', z=', num2str(z), ', f=', num2str(f)])

%% figure 4: max effect of SS on cat on changing mean norm - category width
figure; 
scatter(cbn_vec*binWid,SS_effect_max_vec,sz,cmap);
xlabel('{\itc}: category width'); ylabel('max  mean  norm  change  due  to  discretized  social  effects')
%title(['E=', num2str(E_sd^2), ', S=', num2str(1/S_inv), ', m=', num2str(m), ', z=', num2str(z), ', f=', num2str(f)])




