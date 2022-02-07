function error = fitness(para)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for objective functions of parameter optimizations intumorigenesis  
% regulatory network model
% This script needs 'model.m' for running.
% Input para: 27 parameters need to be optimized
% Output error: one objective function value
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

global EGF ERBB2 EGF_ERBB2 MAPK1 FOXO3 TGFB2 TGFBR1 TGFB2_TGFBR1 SMAD23 pSMAD23 MYC ERBB2_DRUG DRUG

% y index for each species
EGF = 1;
ERBB2 = 2;
EGF_ERBB2 = 3;
MAPK1 = 4;
FOXO3 = 5;
TGFB2 = 6;
TGFBR1 = 7;
TGFB2_TGFBR1 = 8;
SMAD23 = 9;
pSMAD23 = 10;
MYC = 11;
ERBB2_DRUG = 12;
DRUG = 13;


t_data = [0 1 3 48];

% empirical data
FOXO3_data = [1.5403 1.7199 1.7289];
TGFBR1_data = [1.3959 1.6867 1.6087];
FOXO3_change = [0.1166 0.1224];
TGFBR1_change = [0.2083 0.1524];

TGFB2_change_48 = [0.12];
pSMAD23_change_48 = [0.18];

% CMYC_change = [0.3 0.4];
% MYC_change = [0.3 0.2];
MAPK1_change = [-0.3 -0.4];

tstable = [0:1:1000];
tspan = [0:1:50];
t_index = t_data + 1;

y0 = ones(1,13);
y0(ERBB2_DRUG) = 0;
y0(DRUG) = 0;

odefun = @(t, y) model(t, y, para);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[~, YOUT1] = ode23tb(odefun, tstable, y0, options);

if size(YOUT1,1) ~= length(tstable)
    error = 999;  return
elseif sum(abs(YOUT1(end-10,:) - YOUT1(end,:))./YOUT1(end,:))>0.01 || sum(any(YOUT1<0))>=1 %check stable and negative cases
    error = 999;  return
end

y0 = YOUT1(end,:);
y0(end) = 10;
[~, YOUT2] = ode23tb(odefun, tspan, y0, options);

%check stable and negative cases
if size(YOUT2,1) ~= length(tspan) || sum(any(YOUT2<0))>=1 
    error = 999; 
else
    Y_SIM = [YOUT2(t_index,:) ...
             YOUT2(t_index,TGFBR1)+YOUT2(t_index,TGFB2_TGFBR1) ... % total TGFBR1
             YOUT2(t_index,TGFB2)+YOUT2(t_index,TGFB2_TGFBR1)];  % total TGFB2
    
    Y_rel = [Y_SIM(2,:)./Y_SIM(1,:)-1; Y_SIM(3,:)./Y_SIM(1,:)-1; ...
             Y_SIM(4,:)./Y_SIM(1,:)-1]; % relative change
    
    error = sumabs(Y_rel(3,TGFB2)-TGFB2_change_48)...
            + sumabs(Y_rel(1:2,TGFBR1)-TGFBR1_change) + sumabs(Y_rel(1:2,FOXO3)-FOXO3_change)...
            + sumabs(Y_rel(3,pSMAD23) - pSMAD23_change_48) + sumabs(Y_rel(1:2,MAPK1)-MAPK1_change); %+ sumabs(Y_rel(1:2,MYC)-MYC_change);
        
    error = error * 100;
     
%     if error < 20            % + sumabs(Y_rel(:,end)-TGFBR_change) ...
%         disp(error)
%         disp(para)
%     end

end

end