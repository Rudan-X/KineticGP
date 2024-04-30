function [param_names]=load_parameter_name()
reactions=cell(57,1);
enzymes=cell(57,1);
vmax=cell(53,1);

reactions{1}='1. CO2->HCO3 [MC]';
enzymes{1}='4.2.1.1 Carbonic anhydrase (CA)';

reactions{2}='2. HCO3->PEP+OAA [MC]';
enzymes{2}='4.1.1.31 Phosphoenolpyruvate carboxylase (PEPC)';

reactions{3}='3. OAA+NADPH->MAL+NADP [MCchl]';
enzymes{3}='1.1.1.82 Malate dehydrogenase (MDH)';

reactions{4}='4. MAL+NADP->CO2+PYR+NADPH [BCchl]';
enzymes{4}='1.1.1.40 NADP-Malic enzyme (ME)';

reactions{5}='5. PYR+ATP->PEP [MCchl]';
enzymes{5}='2.7.9.1 Pyruvate, phosphate dikinase (PPDK)';

reactions{6}='6. CO2+RuBP->2PGA [BCchl]';
enzymes{6}='4.1.1.39 Rubisco';

reactions{7}='7. PGA+ATP->ADP+DPGA [BCchl]';
enzymes{7}='2.7.2.3 (PGAK)';

reactions{8}='8. DPGA+NADPH->GAP+Pi+NADP [BCchl]';
enzymes{8}='1.2.1.13 ()';

reactions{9}='9. GAP<->DHAP [BCchl]';
enzymes{9}='5.3.1.1 ()';

reactions{10}='10. DHAP+GAP<->FBP [BCchl]';
enzymes{10}='4.1.2.13 ()';

reactions{11}='11. FBP<->F6P+Pi [BCchl]';
enzymes{11}='3.1.3.11 ()';

reactions{12}='12. E4P+DHAP<->SBP [BCchl]';
enzymes{12}='4.1.2.13 ()';

reactions{13}='13. SBP<->S7P+Pi [BCchl]';
enzymes{13}='3.1.3.37 ()';

reactions{14}='14. F6P+GAP<->E4P+Xu5P [BCchl]';
enzymes{14}='2.2.1.1X ()';

reactions{15}='15. S7P+GAP<->Ri5P+Xu5P [BCchl]';
enzymes{15}='2.2.1.1R ()';

reactions{16}='16. Ri5P<->Ru5P [BCchl]';
enzymes{16}='5.3.1.6 ()';

reactions{17}='17. Xu5P<->Ru5P [BCchl]';
enzymes{17}='5.1.3.1 ()';

reactions{18}='18. Ru5P+ATP<->RuBP+ADP [BCchl]';
enzymes{18}='2.7.1.19 ()';

reactions{19}='19. PGA+ATP->ADP+DPGA [MCchl]';
enzymes{19}='2.7.2.3 ()';

reactions{20}='20. DPGA+NADPH->GAP+Pi+NADP [MCchl]';
enzymes{20}='1.2.1.13 ()';

reactions{21}='21. ADPG+ATP->Starch+ADP [MCchl]';
enzymes{21}='2.4.1.21 starch synthase ()';

reactions{22}='22. PGA->sink ';
enzymes{22}='PGA sink ()';

reactions{23}='23. DHAP+GAP<->FBP [MC]';
enzymes{23}='4.1.2.13 ()';

reactions{24}='24. FBP<->F6P+Pi [MC]';
enzymes{24}='3.1.3.11 ()';

reactions{25}='25. F6P<->G6P, G6P<->G1P [MC]'; %Suc6
enzymes{25}='5.3.1.9, 5.4.2.2 ()';

reactions{26}='26. G1P+UTP<->UDPG+Pi [MC]'; %Suc7
enzymes{26}='2.7.7.9 ()';

reactions{27}='27. UDPG+F6P<->SUCP+UDP [MC]'; %Suc8
enzymes{27}='2.4.1.14 ()';

reactions{28}='28. SUCP<->Pi+SUC [MC]';
enzymes{28}='3.1.3.24 ()';

reactions{29}='29. SUC->sink [MC]';
enzymes{29}='SUC sink ()';

reactions{30}='30. F6P+ATP->F26BP+ADP [MC]';
enzymes{30}='2.7.1.105 ()';

reactions{31}='31. F26BP[c]->F6P[c]+Pi[c] [MC]';
enzymes{31}='3.1.3.46 ()';

reactions{32}='32. ADP+Pi<->ATP [MCthy]'; % light reactions
enzymes{32}='3.6.3.14 ()';

reactions{33}='33. NADP<->NADPH [MCthy]'; % light reactions
enzymes{33}='1.18.1.2 ()';

reactions{34}='34. ADP+Pi<->ATP [BCthy]'; % light reactions
enzymes{34}='3.6.3.14 ()';

reactions{35}='35. Metabolite transport through plasmodesmata';
enzymes{35}='Metabolite transport through plasmodesmata';

reactions{36}='36. Pi equilibrium';% 
enzymes{36}='KEPi ()';

reactions{37}='37. NADP<->NADPH [BCthy]'; % light reactions
enzymes{37}='1.18.1.2 ()';


% Photorespiration
reactions{38}='38. O2+RuBP<->PGCA+PGA [BCchl]';
enzymes{38}='4.1.1.39';

reactions{39}='39. PGCA->Pi+GCA [BCchl]';
enzymes{39}='3.1.3.18';

reactions{40}='40. GCA+O2<->H2O2+GOA [BCmit]'; % GLX==GOA
enzymes{40}='1.1.3.15';

reactions{41}='41. GOA+GLU<->GLY+KG [BCmit]';
enzymes{41}='2.6.1.4';

reactions{42}='42. GLY+NAD<->SER+NADH+NH3 [BCmit]';
enzymes{42}='Gly_ser';

reactions{43}='43. SER+GOA<->HPR+GLY [BCmit]';
enzymes{43}='2.6.1.45';

reactions{44}='44. HPR+NAD<->GCEA+NADH [BCmit]';
enzymes{44}='1.1.1.29';

reactions{45}='45. GCEA+ATP<->PGA+ADP [BCmit]';
enzymes{45}='2.7.1.31';

reactions{46}='46. GCA[BCmit]<->GCA [BCchl]';
enzymes{46}='GCA transport';

reactions{47}='47. GCEA[BCmit]<->GCEA [BCchl]';
enzymes{47}='GCEA transport';

reactions{48}='48. PGA<->PEP [MC]';
enzymes{48}='PGA mutase enolase';

reactions{49}='49. PPDK inactivation';
enzymes{49}='';

reactions{50}='50. PPDK activation';
enzymes{50}='';

reactions{51}='51. G1P + ATP <-> ADP + ADPG';
enzymes{51}='2.7.7.27 (ADPG pyrophosphorylase)';

reactions{52}='52. PPi + H2O <-> 2 Pi [Bchl]';
enzymes{52}='3.6.1.1 (inorganic diphosphatase)';

reactions{53}='53. ADPG <-> Starch';
enzymes{53}='2.4.1.21 (starch synthase)';

reactions{54}='54. Hexose phosphate rate';
enzymes{54}='';

reactions{55}='55. PGA, GAP and DHAP transport [BC]';
enzymes{55}='TPT';

reactions{56}='56. PGA, GAP and DHAP transport [MC]';
enzymes{56}='TPT';

reactions{57}='57. Gs dynamics [MC]';
enzymes{57}='';

%% VMAX
vmax{1}='Vmax[1]';
vmax{2}='Vpmax';
vmax{3}='Vmax[3]';
vmax{4}='Vmax[4]';
vmax{5}='Vmax[5]';
vmax{6}='Vmax';
vmax{7}='Vmax[7_8]';
% vmax{8}='Vmax[8]';
vmax{8}='Vmax[10]';
vmax{9}='Vmax[11]';
vmax{10}='Vmax[12]';
vmax{11}='Vmax[13]';
vmax{12}='Vmax[14]';
vmax{13}='Vmax[15]';
vmax{14}='Vmax[18]';
vmax{15}='Vmax[19_20]'; 
% vmax{17}='Vmax[20]';
% vmax{18}='Vmax[21]'; starch synthesis in Mchl
vmax{16}='Vmax[22]';
vmax{17}='Vmax[23]';
vmax{18}='Vmax[24]';
vmax{19}='Vmax[26]';
vmax{20}='Vmax[27]';
vmax{21}='Vmax[28]';
vmax{22}='Vmax[29]';
vmax{23}='Vmax[30]';
vmax{24}='Vmax[31]';
vmax{25}='Jmax';
% ATP and NADPH synthesis in BC and MC
vmax{26}='Vmax[32]';
vmax{27}='Vmax[33]';
vmax{28}='Vmax[34]';
vmax{29}='Vmax[37]';
% photorespiration
vmax{30}='ARratio';
vmax{31}='Vmax[39]';
vmax{32}='Vmax[40]';
vmax{33}='Vmax[41]';
vmax{34}='Vmax[42]';
vmax{35}='Vmax[43]';
vmax{36}='Vmax[44]';
vmax{37}='Vmax[45]';
vmax{38}='Vmax[46]';
vmax{39}='Vmax[47]';
vmax{40}='Vmax[48]';

vmax{41}='Vtp_Bchl';
vmax{42}='Vtp_Mchl';
vmax{43}='Vm_Sta1';
vmax{44}='Vm_Sta2';
vmax{45}='Vm_Sta3';
vmax{46}='Vmax_OAA_MC_MCchl[35]';
vmax{47}='Vmax_PYR_BC[35]';
vmax{48}='Vmax_PYR_MC[35]';
vmax{49}='Vmax_PEP_MC[35]';
vmax{50}='Vmax_PEP_BC[35]';
vmax{51}='Vmax_MAL_BC[35]';
vmax{52}='Vmax_MAL_MC[35]';
vmax{53}='Vmax_Hep';

%% KM

km=cell(50,10);
km{1,1}='KmCO2';  km{1,2}='Ke[1]';
km{2,1}='KmHCO3';  km{2,2}='KmPEP[2]';   km{2,3}='KiMAL[2]'; km{2,4}='KiMALn[2]';
km{3,1}='KmNADPH';  km{3,2}='KmOAA[3]';  km{3,3}='KmNADP[3]';  km{3,4}='KmMAL[3]';  km{3,5}='Ke[3]';
km{4,1}='KmCO2';  km{4,2}='KmNADP[4]';  km{4,3}='KmNADPH[4]';  km{4,4}='KmPyr[4]';  km{4,5}='KmMAL[4]';  km{4,6}='Ke[4]';
km{5,1}='KiPEP';  km{5,2}='KmATP[5]';  km{5,3}='KmPyr[5]';
km{6,1}='KmCO2';  km{6,2}='KmO2[6]';  km{6,3}='KmRuBP[6]';  km{6,4}='KiPGA[6]';  km{6,5}='KiFBP[6]';  km{6,6}='KiSBP[6]';  km{6,7}='KiPi[6]';  km{6,8}='KiNADPH[6]';

% km{7,1}='KmATP[78]';  km{7,2}='KmPGA[78]'; 
% km{8,1}='KmNADPH[78]';
km{7,1}='KmATP[7_8]';  km{7,2}='KmPGA[7_8]'; 
km{8,1}='KmNADPH[7_8]';

km{9,1}='Ke'; 
km{10,1}='KmDHAP';  km{10,2}='KmFBP[10]';  km{10,3}='KmGAP[10]';  km{10,4}='Ke[10]';
km{11,1}='KiF6P';  km{11,2}='KiPi[11]';  km{11,3}='KmFBP[11]';  km{11,4}='Ke[11]';
km{12,1}='KmDHAP';  km{12,2}='KmE4P[12]';  km{12,3}='Ke[12]';
km{13,1}='KiPi';  km{13,2}='KmSBP[13]';  km{13,3}='Ke[13]';
km{14,1}='KmE4P';  km{14,2}='KmF6P[14]';  km{14,3}='KmGAP[14]';  km{14,4}='KmXu5P]';  km{14,5}='Ke[14]';
km{15,1}='KmGAP';  km{15,2}='KmRi5P[15]';  km{15,3}='KmS7P[15]';  km{15,4}='KmXu5P[15]';  km{15,5}='Ke[15]';
km{16,1}='Ke';
km{17,1}='Ke';
km{18,1}='KiADP';  km{18,2}='KiADP[18]';  km{18,3}='KiPGA[18]';  km{18,4}='KiPi[18]';  km{18,5}='KiRuBP[18]';  km{18,6}='KmATP[18]';  km{18,7}='KmRu5P[18]';  km{18,8}='Ke[18]';

% km{19,1}='KmADP';  km{19,2}='KmATP[19]';  km{19,3}='KmPGA[19]';
% km{20,1}='KmDPGA';  km{20,1}='KmNADPH[20]';

km{19,1}='KmATP[19_20]';  km{19,2}='KmPGA[19_20]';
km{20,1}='KmNADPH[19_20]';

% km{21,1}='KiADP';  km{21,2}='KmATP[21]';  km{21,3}='KmG1P[21]';  km{21,4}='KaF6P[21]';  km{21,5}='KaFBP[21]';  km{21,6}='KaPGA[21]';  km{21,7}='Ke[21]';   km{21,8}='Ke[21]';
km{21,1}='Ke_Starch1';   km{21,2}='Ke_Starch2[21]';

km{22,1}='KmPGA';
km{23,1}='KmDHAP';  km{23,2}='KmGAP[23]';  km{23,3}='KmFBP[23]';  km{23,4}='Ke[23]';
km{24,1}='KiF26BP';  km{24,2}='KiF6P[24]';  km{24,3}='KiPi[24]';  km{24,4}='KmFBP[24]';  km{24,5}='Ke[24]';
km{25,1}='Ke';  km{25,2}='Ke[25]';
km{26,1}='KmG1P';  km{26,2}='KmPPi[26]';  km{26,3}='KmUDPG[26]';  km{26,4}='KmUTP[26]';  km{26,5}='Ke[26]';
km{27,1}='KiFBP';  km{27,2}='KiPi[27]';  km{27,3}='KiSuc[27]';  km{27,4}='KiSucP[27]';  km{27,5}='KiUDP[27]';  km{27,6}='KmF6P[27]';  km{27,7}='KmUDPG[27]';  km{27,8}='Ke[27]'; 
km{28,1}='KmSuc';  km{28,2}='KmSucP[28]';  km{28,3}='Ke[28]';
km{29,1}='KmSuc';
km{30,1}='KiADP';  km{30,2}='KIDHAP[30]';  km{30,3}='KmATP[30]';  km{30,4}='KmF26BP[30]';  km{30,5}='KmF6P[30]';  km{30,6}='Ke[30]';
km{31,1}='KiF6P';  km{31,2}='KiPi[31]';  km{31,3}='KmF26BP[31]';  
km{36,1}='KePi';

km{32,1}='KmADP';  km{32,2}='KmATP[32]';  km{32,3}='KmPi[32]';  km{32,4}='X[32]';  km{32,5}='Y[32]';  km{32,6}='F[32]';  km{32,7}='Q[32]';  km{32,8}='D[32]'; km{32,9}='Ke[32]';
km{33,1}='KmNADP';  km{33,2}='KmNADPH[33]'; km{33,3}='Ke[33]';  km{33,4}='E[33]'; 
km{34,1}='KmADP';  km{34,2}='KmPi[34]';   km{34,3}='KmATP[34]';  km{34,4}='Ke[34]';   km{34,5}='G[34]';
km{37,1}='KmNADP';  km{37,2}='KmNADPH[37]'; km{37,3}='Ke[37]'; 

% km{35,1}='Voaa';  km{35,2}='VMAL[35]';  km{35,3}='Vpyr[35]';  km{35,4}='Vpep[35]';  km{35,5}='Vt[35]';  km{35,6}='Vleak[35]'; km{35,7}='Vpga[35]';
km{35,1}='Km_OAA_M';  km{35,2}='Kimal_OAA_M';  km{35,3}='Km_MAL_M';  km{35,4}='KiOAA_MAL_M'; km{35,5}='Km_MAL_B'; km{35,6}='Km_PYR_B';  km{35,7}='Km_PYR_M'; km{35,8}='Km_PEP_M';

km{38,1}='KmCO2'; km{38,2}='KmO2[38]';  km{38,3}='KmRuBP[38]';  km{38,4}='KiPGA[38]';  km{38,5}='KiFBP[38]';  km{38,6}='KiSBP[38]';  km{38,7}='KiPi[38]';  km{38,8}='KiNADPH[38]';
km{39,1}='KmPGCA';  km{39,2}='KiPI[39]';  km{39,3}='KiGCA[39]';
km{40,1}='KmGCA[40]';
km{41,1}='Ke';  km{41,2}='KmGOA[41]';  km{41,3}='KmGLU[41]';  km{41,4}='KiGLY[41]';
km{42,1}='KmGLY';  km{42,2}='KiSER[42]';
km{43,1}='Ke';  km{43,2}='KmGOA[43]';  km{43,3}='KmSER[43]';  km{43,4}='KmGLY[43]';
km{44,1}='Ke';  km{44,2}='KiHPR[44]';  km{44,3}='KmHPR[44]';
km{45,1}='Ke';  km{45,2}='KmATP[45]';  km{45,3}='KmGCEA[45]';  km{45,4}='KiPGA[45]';
km{46,1}='KmGCA';  km{46,2}='KiGCEA[46]';
km{47,1}='KmGCEA';  km{47,2}='KiGCA[47]';
km{48,1}='KmPGA'; km{48,2}='KmPEP[48]'; km{48,3}='Ke[48]';

% new:
km{49,1}='Kcat_EA_PPDKRP_I'; km{49,2}='Km_EA_PPDKRP_I_ADP[49]'; km{49,3}='Ki_EA_PPDKRP_I_Pyr[49]'; km{49,4}='Km_EA_PPDKRP_I_E[49]';
km{50,1}='Kcat_EA_PPDKRP_A'; km{50,2}='Km_EA_PPDKRP_A_Pi[50]'; km{50,3}='Ki_EA_PPDKRP_A_AMP[50]'; km{50,4}='Km_EA_PPDKRP_A_EP[50]'; km{50,5}='Ki_EA_PPDKRP_A_ADP[50]'; km{50,6}='Ki_EA_PPDKRP_A_PPI[50]';

km{51,1}='KaPGA_Sta1'; km{51,2}='KmG1P_Sta1[51]'; km{51,3}='KmATP_Sta1[51]'; km{51,4}='KIAPi_ATP_Sta1[51]'; km{51,5}='KmPPi_Sta1[51]'; km{51,6}='KICPP1_ATP_Sta1[51]'; km{51,7}='KmADPG_Sta1[51]'; km{51,8}='KIAADP_ATP_Sta1[51]'; km{51,9}='Ke_Sta1[51]'; 
km{52,1}='KmPPi_Sta2'; km{52,2}='Ke_Sta2[52]';
km{53,1}='KmADPG_Sta3'; 
km{54,1}='Kmpi_hexp'; km{54,2}='Kmhexp_hexp[54]';

km{55,1}='KmPGA_B'; km{55,2}='KmGAP_B[55]'; km{55,3}='KmDHAP_B[55]'; 
km{56,1}='KmPGA_M'; km{56,2}='KmGAP_M[56]'; km{56,3}='KmDHAP_M[56]';
km{57,1}='Ki[57]'; km{57,2}='Kd[57]';

%%
final_km=[""];
global KVlen
for i=1:size(km,1)
    temp=km(i,1:KVlen(i));
    temp{1}=strcat(reactions{i}," ",temp{1});
    final_km=[final_km;temp'];
end
final_km=final_km(2:length(final_km));

%%
act_rate=cell(10,1);
act_rate{1}='tao_ActPEPC';
act_rate{2}='tao_ActFBPase';
act_rate{3}='tao_ActSBPase';
act_rate{4}='tao_ActATPsynthase';
act_rate{5}='tao_ActGAPDH';
act_rate{6}='tao_ActPRK';
act_rate{7}='tao_ActNADPMDH';
act_rate{8}='KaRac';
act_rate{9}='tao_ActRubisco';
act_rate{10}='tao_ActRca';


%%
perm=cell(10,1);
perm{1}='Perm-MAL';
perm{2}='Perm-PYR';
perm{3}='Perm-CO2';
perm{4}='Perm-PGA';
perm{5}='Perm-CO2[BC]';
perm{6}='Bper_GLU';
perm{7}='Bper_KG';
perm{8}='Bper_NADH';
perm{9}='Bper_NAD';
perm{10}='gm';


% params=["BBslope","BBintercept","factorvp","factorvc","[PDRP]","MRd"];
params=["BBslope","BBintercept","factorvp","factorvc","[PDRP]","MRd"];

% tempvmax=vmax(setdiff(1:length(vmax),[2,6]));
param_names=[params';final_km;vmax;act_rate;perm];