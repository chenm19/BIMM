%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for ordinary differential equations of tumorigenesis regulatory 
% network model
% Input para: 27 optimized parameters 
% In the article, pSMAD23 stands for pSMAD2/3; SMAD23 stands for SMAD2/3 
% TGFB2 stands for TGFβ2; TGFBR1 stands for TGFβR1;
% TGFB2_TGFBR1 stands for TGFβ2_TGFβR1
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function dy = model(t, y, parameters)
global EGF ERBB2 EGF_ERBB2 MAPK1 FOXO3 TGFB2 TGFBR1 TGFB2_TGFBR1 SMAD23 pSMAD23 MYC DRUG ERBB2_DRUG

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
dy = zeros(13,1);
ks_ef = parameters(1);
kd_ef = 1;
kb_efeb = parameters(2);
ku_efeb = parameters(3);
ks_eb = parameters(4);
kd_eb = 1;
kd_efeb = 1;
ks_mk = parameters(5);
kd_mk = 1;
alpha_efeb = parameters(6);
beta_mcmk = parameters(7);
km_mcmk = parameters(8);
ks_fo = parameters(9);
kd_fo = 1;
beta_mkfo = parameters(10);
km_mkfo = parameters(11);
ks_tb = parameters(12);
kd_tb = 1;
kb_tbtr = parameters(13);
ku_tbtr = parameters(14);
kb_dgeb = parameters(15);
ku_dgeb = parameters(16);
ks_tr = parameters(17);
kd_tr = 1;
kd_tbtr = 1;
ks_sd = parameters(18);
kd_sd = 1;
alpha_tbtr = parameters(19);
ks_mc = parameters(20);
kd_mc = 1;
beta_sdmc = parameters(21);
km_psdmc = parameters(22);
kd_dg = log(2)/24;
kpho_sd = parameters(23);
kdepho_sd = parameters(24);
alpha_tbdg = parameters(25);
alpha_mkmc = parameters(26);
alpha_trdg = parameters(27);

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
% Differential equations of the model
dy(EGF) = ks_ef - kd_ef * y(EGF) - kb_efeb * y(EGF) * y(ERBB2) + ku_efeb * y(EGF_ERBB2) ;
dy(ERBB2) = ks_eb - kd_eb * y(ERBB2) - kb_efeb * y(EGF) * y(ERBB2) + ku_efeb * y(EGF_ERBB2) ...
          - kb_dgeb * y(DRUG) * y(ERBB2) + ku_dgeb * y(ERBB2_DRUG);
dy(EGF_ERBB2) =  - kd_efeb * y(EGF_ERBB2) + kb_efeb * y(EGF) * y(ERBB2) - ku_efeb * y(EGF_ERBB2);
dy(MAPK1) = ks_mk * (1 + alpha_efeb * y(EGF_ERBB2) + beta_mcmk /(km_mcmk + y(MYC)) ) - kd_mk * y(MAPK1);
dy(FOXO3) = ks_fo * (1 + beta_mkfo/(km_mkfo + y(MAPK1)) ) - kd_fo * y(FOXO3); % + ksd_fo * y(SMAD23)

dy(TGFB2) =  ks_tb * (1 + alpha_tbdg * y(DRUG))- kd_tb * y(TGFB2) - kb_tbtr * y(TGFB2) * y(TGFBR1) + ku_tbtr * y(TGFB2_TGFBR1);
dy(TGFBR1) = ks_tr * (1 + alpha_trdg * y(DRUG)) - kd_tr * y(TGFBR1) - kb_tbtr * y(TGFB2) * y(TGFBR1) + ku_tbtr * y(TGFB2_TGFBR1);
dy(TGFB2_TGFBR1) = - kd_tbtr * y(TGFB2_TGFBR1) + kb_tbtr * y(TGFB2) * y(TGFBR1) - ku_tbtr * y(TGFB2_TGFBR1);
dy(SMAD23) = ks_sd - kd_sd * y(SMAD23) - kpho_sd * (1 + alpha_tbtr * y(TGFB2_TGFBR1)) * y(SMAD23) + kdepho_sd * y(pSMAD23);
dy(pSMAD23) = - kd_sd * y(pSMAD23) + kpho_sd * (1 + alpha_tbtr * y(TGFB2_TGFBR1)) * y(SMAD23) - kdepho_sd * y(pSMAD23);
dy(MYC) = ks_mc * (1 + alpha_mkmc * y(MAPK1) + beta_sdmc/(km_psdmc+y(pSMAD23))) - kd_mc * y(MYC);
dy(ERBB2_DRUG) = - kd_dg * y(ERBB2_DRUG) + kb_dgeb * y(DRUG) * y(ERBB2) - ku_dgeb * y(ERBB2_DRUG);
dy(DRUG) = - kd_dg * y(DRUG) - kb_dgeb * y(DRUG) * y(ERBB2) + ku_dgeb * y(ERBB2_DRUG);
% end of equations
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%

end
