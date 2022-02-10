%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for plotting simulated levels of signature genes pSMAD23, FOXO3,
% MYC, EGF-ERBB2, MAPK1, and TGFβ-TGFβR1 under one, two, and three 
% doses of DEX treatment.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clc; clear; close all;
load('cancer_ga_outputs53.mat', "value", "error");
load('gene_percentile.mat');

Inits = value;

global EGF ERBB2 EGF_ERBB2 MAPK1 FOXO3 TGFB2 TGFBR1 TGFB2_TGFBR1 SMAD23 pSMAD23 MYC DRUG ERBB2_DRUG

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

t_data = [0 1 3];
FOXO3_data = [1.5403 1.7199 1.7289];
TGFBR1_data = [1.3959 1.6867 1.6087];

FOXO3_scaled = FOXO3_data./FOXO3_data(1);
TGFB2_scaled = [1 1.12];
TGFBR1_scaled = TGFBR1_data./TGFBR1_data(1);
pSMAD23_scaled = [1 1.18];

tstable = [0 1000];
tspan = [0 50];
t_index = t_data + 2 + tstable(end);

y0 = ones(1,13);
y0(ERBB2_DRUG) = 0;
y0(DRUG) = 0;

odefun = @(t, y) model(t, y, Inits);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[TOUT1, YOUT1] = ode23tb(odefun, tstable, y0, options);
%% one dose applied 
y0 = YOUT1(end,:);
y0(end) = 10;
[TOUT_d, YOUT_d] = ode23tb(odefun, tspan, y0, options);

TOUT = [TOUT1; TOUT_d+1000];
YOUT = [YOUT1; YOUT_d];
YOUT_d = [YOUT_d YOUT_d(:,TGFBR1)+YOUT_d(:,TGFB2_TGFBR1) YOUT_d(:,TGFB2)+YOUT_d(:,TGFB2_TGFBR1)];
YOUTd_scaled = YOUT_d./YOUT_d(1,:);

%% two doses applied
y0 = YOUT1(end,:);
y0(end) = 20;
[TOUT_2d, YOUT_2d] = ode23tb(odefun, tspan, y0, options);

TOUT2d = [TOUT1; TOUT_2d+1000];
YOUT2d = [YOUT1; YOUT_2d];
YOUT_2d = [YOUT_2d YOUT_2d(:,TGFBR1)+YOUT_2d(:,TGFB2_TGFBR1) YOUT_2d(:,TGFB2)+YOUT_2d(:,TGFB2_TGFBR1)];
YOUT2d_scaled = YOUT_2d./YOUT_2d(1,:);

%% three doses applied 
y0 = YOUT1(end,:);
y0(end) = 30;
[TOUT_3d, YOUT_3d] = ode23tb(odefun, tspan, y0, options);

TOUT3 = [TOUT1; TOUT_3d+1000];
YOUT3 = [YOUT1; YOUT_3d];
YOUT_3d = [YOUT_3d YOUT_3d(:,TGFBR1)+YOUT_3d(:,TGFB2_TGFBR1) YOUT_3d(:,TGFB2)+YOUT_3d(:,TGFB2_TGFBR1)];
YOUT3d_scaled = YOUT_3d./YOUT_3d(1,:);

% settings for figures 
extraInputs = {'fontsize',20};
path ='./diffDose';

%% pSMAD23
ytime = 1:1:length(TOUT_d);
ytime02 = 1:1:length(TOUT_2d);
ytime04 = 1:1:length(TOUT_3d);
figure(1)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, pSMAD23),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, pSMAD23),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, pSMAD23),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0.8 4])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
ax = gca;ax.FontSize = 20;
r = rectangle('Position',[0 0 50 4]');
title('pSMAD2/3');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['pSMAD23','.png'] ));

%% FOXO3
figure(2)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, FOXO3),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, FOXO3),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, FOXO3),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0.9 1.5])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.5]');
ax = gca;ax.FontSize = 20;
title('FOXO3');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['FOXO3','.png'] ));

%% MYC
figure(3)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, MYC),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, MYC),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, MYC),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0.4 1.2])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.2]');
ax = gca;ax.FontSize = 20;
title('MYC');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'location','Southeast');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['MYC','.png'] ));
%% EGF-ERBB2
figure(4)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, EGF_ERBB2),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, EGF_ERBB2),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, EGF_ERBB2),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0 1.2])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.2]');
ax = gca;ax.FontSize = 20;
title('EGF-ERBB2');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'location','Southeast');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['EGF_ERBB2','.png'] ));
%% MAPK1
figure(5)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, MAPK1),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, MAPK1),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, MAPK1),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0.3 1.2])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
r = rectangle('Position',[0 0 50 1.2]');
ax = gca;ax.FontSize = 20;
title('MAPK1');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'location','Southeast');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['MAPK1','.png'] ));
%% TGFβ-TGFβR1
figure(6)
hold on;
plot(TOUT_d(ytime),YOUTd_scaled(ytime, TGFB2_TGFBR1),'k','LineWidth',3,'color','#DE7878')
plot(TOUT_2d(ytime02),YOUT2d_scaled(ytime02, TGFB2_TGFBR1),'k','LineWidth',3,'color','#EEDF96')
plot(TOUT_3d(ytime04),YOUT3d_scaled(ytime04, TGFB2_TGFBR1),'k','LineWidth',3,'color','#9CB4E0')

axis([0 50,0 8])
xlabel('Time (hour)',extraInputs{:});
ylabel('Level',extraInputs{:});
r = rectangle('Position',[0 0 50 8]');
ax = gca;ax.FontSize = 20;
title('TGFβ-TGFβR1');
s_leg =legend('1x DEX','2x DEX','3x DEX');
set(s_leg,'Units','Normalized','FontUnits','Normalized')
saveas(gcf, fullfile(path, ['TGFB_TGFBR1','.png'] ));