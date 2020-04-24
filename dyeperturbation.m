
%dye perturbation test
% A function to test the perturbation of purifed mCP on the sample pH
%      This function calculate the extra TA contributed by the dye addition
%      and the equilibrium move with the mCP indicator addition. The mCP
%      characteriztion parameters are from :
%%%%   Mueller, J.D. and Rehder, G., 2018. Metrology of pH measurements 
%      in Brackish Waters—Part 2: experimental characterization of 
%      purified meta-Cresol purple for spectrophotometric pHT measurements. 
%      Frontiers in Marine Science, 5, p.177.
%%%%   Liu, X., Patsavas, M.C. and Byrne, R.H., 2011. Purification and characterization of 
%      meta-cresol purple for spectrophotometric seawater pH measurements.
%      Environmental science & technology, 45(11), pp.4862-4868.
%%%%   Lai, C.Z., DeGrandpre, M.D., Wasser, B.D., Brandon, T.A., Clucas, D.S., 
%      Jaqueth, E.J., Benson, Z.D., Beatty, C.M. and Spaulding, R.S., 2016.
%      Spectrophotometric measurement of freshwater pH with purified meta?cresol purple
%      and phenol red. Limnology and Oceanography: Methods, 14(12), pp.864-873.
% This function is modified from CO2SYS MATLAB version, which calculates 
%      and returns the state of the carbonate system of oceanographic water samples
%      if supplied with enough input. For more information about CO2SYS
%      MATLAB version, please look at:
%     Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
%     CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
%     Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
%     Oak Ridge, Tennessee. 
%     http://cdiac.ornl.gov/oceans/co2rprt.html
%**************************************************************************
%
%  **** SYNTAX:
%  [DATA,HEADERS,NICEHEADERS]= dyeperturbation(PAR1,PAR2,SAL,TEMPIN,SI,PO4,pHSCALEIN,
%                K1K2CONSTANTS,KSO4CONSTANTS,DYEC,DYEpH,DYES,SAMPM,DYEM,MCPCONSTANTS)
% 
%  **** SYNTAX EXAMPLES:
%  [Result]                     = dyeperturbation(2500,2000,35,25,0,0,1,9,1,2500,8,35,0.015,2e-5,2)
%  [Result,Headers]             = dyeperturbation(2500,2000,15,25,0,0,1,9,1,2000,8,0,0.015,2e-5,1)
%  [Result,Headers,Niceheaders] = dyeperturbation(2500,2000,0,25,0,0,1,9,1,2000,8,0,0.015,2e-5,3)
%  [A]                          = dyeperturbation(2500,2000,0,25,0,0,1,9,1,2000,8,0,0.015,2e-5,3)
%
%**************************************************************************
%
% INPUT:
%
%   PAR1(TA, umol/kg) : scalar or vector of size n
%   PAR2(DIC, umol/kg): scalar or vector of size n
%   SAL            () : scalar or vector of size n
%   TEMPIN  (degr. C) : scalar or vector of size n 
%   SI    (umol/kgSW) : scalar or vector of size n
%   PO4   (umol/kgSW) : scalar or vector of size n
%   pHSCALEIN         : scalar or vector of size n (*)
%   K1K2              : scalar or vector of size n (**)
%   KSO4              : scalar or vector of size n (***)
%   DYEC(umol/kg)     : scalar or vector of size n (This is the indicator concentration)
%   DYEpH             : scalar or vector of size n (This is the indicator original pH)
%   DYES              : scalar or vector of size n (This is the indicator salinity)
%   SAMPM(kg)         : scalar or vector of size n (This is the sample mass)
%   DYEM (kg)         : scalar or vector of size n (This is the indicator mass)
%   MCPCONSTANTS      : scalar or vector of size n (****)


%  (*) Each element must be an integer, 
%       indicating that the pH-input (PAR1 or PAR2, if any) is at:
%  1 = Total scale
%  2 = Seawater scale
%  3 = Free scale
%  4 = NBS scale
% 
%  (**) Each element must be an integer, 
%        indicating the K1 K2 dissociation  that are to be used:
%   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0. 
%   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
% 
%  (***) Each element must be an integer that 
%         indicates the KSO4 dissociation constants that are to be used,
%         in combination with the formulation of the borate-to-salinity ratio to be used.
%         Having both these choices in a single argument is somewhat awkward, 
%         but it maintains syntax compatibility with the previous version.
%  1 = KSO4 of Dickson & TB of Uppstrom 1979  (PREFERRED) 
%  2 = KSO4 of Khoo    & TB of Uppstrom 1979
%  3 = KSO4 of Dickson & TB of Lee 2010
%  4 = KSO4 of Khoo    & TB of Lee 2010
%
%   (**) Each element must be an integer, 
%        indicating the mCP K1 K2 dissociation and e1 e2 e3 molar absorptivity that are to be used:
%   1 = Müller and Rehder,2018						        T:    5-35  S: 0-40. Total scale. Purified mCP.
%   2 = Liu et al., 2011								    T:    5-35  S: 25-40. Total scale. Purified mCP
%   3 = Lai et al., 2016                                  	T:    8-30  S: 0. Total scale.  Purified mCP
%**************************************************************************%
%
% OUTPUT: * an array containing the following parameter values (one row per sample):
%         *  a cell-array containing crudely formatted headers
%         *  a cell-array containing nicely formatted headers
%
%     POS - PARAMETER         UNIT
%
%     '01 - delta pH                     '; The pH with dye addition - pH without dye addition
%     '02 - TAlk             (umol/kgSW) '; Input TA
%     '03 - TCO2             (umol/kgSW) '; Input DIC
%     '04 - pCO2             (uatm)      '; pCO2 calculated from Input TA and DIC
%     '05 - fCO2             (uatm)      '; fCO2 calculated from Input TA and DIC
%     '06 - pHin(Total)                  ';
%     '07 - pHin(SWS)                    ';
%     '08 - pHin(FREE)                   ';
%     '09 - pHin(NBS)                    ';
%     '10 - HCO3in           (umol/kgSW) ';
%     '11 - CO3in            (umol/kgSW) ';
%     '12 - CO2in            (umol/kgSW) ';
%     '13 - BAlkin           (umol/kgSW) ';
%     '14 - OHin             (umol/kgSW) ';
%     '15 - PAlkin           (umol/kgSW) ';
%     '16 - SiAlkin          (umol/kgSW) ';
%     '17 - Hfreein          (umol/kgSW) ';
%     '18 - pHmix(Total)                 ';
%     '19 - pHmix(SWS)                   ';
%     '20 - pHmix(FREE                   ';
%     '21 - pHmix(NBS)                   ';
%     '22 - HCO3mix          (umol/kgSW) ';
%     '23 - CO3mix           (umol/kgSW) ';
%     '24 - CO2mix           (umol/kgSW) ';
%     '25 - BAlmix           (umol/kgSW) ';
%     '26 - OHmix            (umol/kgSW) ';
%     '27 - PAlkmix          (umol/kgSW) ';
%     '28 - SiAlmix          (umol/kgSW) ';
%     '29 - IAlk             (umol/kgSW) '; %the dye alkalinity in the mixed solution
%     '30 - Hfreemix         (umol/kgSW) ';
%     '31 - Dye concentration(umol/kgSW) '; %the mCP stock concentration (before mixing)
%     '32 - H2I              (umol/kgSW) '; %the mCP species concentration (before mixing)
%     '33 - HI               (umol/kgSW) '; %the mCP species concentration (before mixing)
%     '34 - I2               (umol/kgSW) '; %the mCP species concentration (before mixing)
%     '35 - Dye TA           (umol/kgSW) '; %the TA of the dye (in the stock solution)
%     '36 - Absorbance ratio             '; %the theoretical R with dye addition
%     '37 - Dye pH                       '; %the mCP stock pH
%     '38 - Dye Salinity                 '; %the mCP stock salinity
%     '39 - Sample Mass                  ';
%     '40 - Dye Mass                     ';
%     '41 - Input temperature            '; 
%     '42 - Sample salinity              '; %Input sample salinity
%     '43 - PO4                          ';
%     '44 - Si                           ';
%     '45 - K1K2CONSTANTS    ()          ';
%     '46 - KSO4CONSTANTS    ()          ';
%     '47 - pHSCALEIN        ()          ';
%     '48 - MCPCONSTANTS     ()          ';
%     '49 - Revelle Factors  ()          ';
%     '50 - TB               (umol/kgSW) ';
%     '51 - TF               (umol/kgSW) ';
%     '52 - TS               (umol/kgSW) ';
%     '53 - TI_m             (umol/kgSW) '; %TI_m, the concentration of mCP in the mixed solution
%
%
%  % Questions, bug reports to xinyuli@udel.edu
%**************************************************************************
% NOTHING BELOW THIS SHOULD REQUIRE EDITING BY USER!
%**************************************************************************
function [DATA,HEADERS,NICEHEADERS]=dyeperturbation(PAR1,PAR2,SAL,TEMPIN,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS,DYEC,DYEpH,DYES,SAMPM,DYEM,MCPCONSTANTS)%**************************************************************************

global pHScale WhichKs WhoseKSO4 Whichmcp Pbar
global Sal sqrSal TempK logTempK TempCi;
global FugFac VPFac PengCorrection ntps RGasConstant;
global fH RT;
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KI2 KI1;
global TB TF TS TP TSi F TI pHI SalI samplms dyems;

% Input conditioning
% Determine lengths of input vectors
veclengths=[length(PAR1) length(PAR2) ...
            length(SAL) length(TEMPIN)...
            length(SI) length(PO4) length(pHSCALEIN)...
            length(K1K2CONSTANTS) length(KSO4CONSTANTS) length(DYEC)...
            length(DYEpH) length(DYES) length(SAMPM) length(DYEM) length(MCPCONSTANTS)];

if length(unique(veclengths))>2
	disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
end

% Make column vectors of all input vectors
PAR1         =PAR1         (:);
PAR2         =PAR2         (:);
SAL          =SAL          (:);
TEMPIN       =TEMPIN       (:);
SI           =SI           (:);
PO4          =PO4          (:);
pHSCALEIN    =pHSCALEIN    (:);
K1K2CONSTANTS=K1K2CONSTANTS(:);
KSO4CONSTANTS=KSO4CONSTANTS(:);
DYEC         =DYEC         (:);
DYEpH        =DYEpH        (:);
DYES         =DYES         (:);
SAMPM        =SAMPM        (:);
DYEM         =DYEM         (:);
PRES         =1;
DYEN         =0;       
MCPCONSTANTS =MCPCONSTANTS (:);

% Find the longest column vector:
ntps = max(veclengths);

% Populate column vectors
PAR1(1:ntps,1)          = PAR1(:)          ;
PAR2(1:ntps,1)          = PAR2(:)          ;
SAL(1:ntps,1)           = SAL(:)           ;
TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
SI(1:ntps,1)            = SI(:)            ;
PO4(1:ntps,1)           = PO4(:)           ;
pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
KSO4CONSTANTS(1:ntps,1) = KSO4CONSTANTS(:) ;
DYEC(1:ntps,1)          = DYEC(:)          ;
DYEpH(1:ntps,1)         = DYEpH(:)         ;
DYES(1:ntps,1)          = DYES(:)          ;
SAMPM(1:ntps,1)         = SAMPM(:)         ;
DYEM(1:ntps,1)          = DYEM(:)          ;
PRES(1:ntps,1)          = PRES(:)          ;
DYEN(1:ntps,1)          = DYEN(:)          ;       
MCPCONSTANTS(1:ntps,1)  = MCPCONSTANTS(:)  ;

% Assign input to the variable.
pHScale      = pHSCALEIN;
WhichKs      = K1K2CONSTANTS;
WhoseKSO4    = KSO4CONSTANTS;
TempCi       = TEMPIN;
Sal          = SAL;
sqrSal       = sqrt(SAL);
TP           = PO4;
TSi          = SI;
RGasConstant = 83.1451;  % ml bar-1 K-1 mol-1, DOEv2
%RGasConstant = 83.14472; % ml bar-1 K-1 mol-1, DOEv3
pHI          = DYEpH;
SalI         = DYES;
samplms      = SAMPM;
dyems        = DYEM;
Whichmcp     = MCPCONSTANTS;

% Generate empty vectors...
TA  = nan(ntps,1); % Talk
TC  = nan(ntps,1); % DIC
PH  = nan(ntps,1); % pH
PC  = nan(ntps,1); % pCO2
FC  = nan(ntps,1); % fCO2

% Assign values to empty vectors.
TA=PAR1./1e6; % Convert from micromol/kg to mol/kg
TC=PAR2./1e6; % Convert from micromol/kg to mol/kg

% Generate the columns holding Si, Phos and Sal.
% Pure Water case:
F=(WhichKs==8);
Sal(F) = 0;
% GEOSECS and Pure Water:
F=(WhichKs==8 | WhichKs==6);  
TP(F)  = 0;
TSi(F) = 0;
% All other cases
F=~F;                         
TP(F)  = TP(F)./1e6;
TSi(F) = TSi(F)./1e6;

% PengCorrection is 0 for all cases where WhichKs is not 7
PengCorrection=zeros(ntps,1); F=(WhichKs==7); PengCorrection(F)=TP(F);

% Generate vector for results, and copy the raw input values into them. This
TAc  = TA;
TCc  = TC;
TAnc = TA;
TCnc = TC;
PHic = PH;
PHnc = PH;
PCic = PC;
FCic = FC;

%*****************Calculation process has several steps 1-4*************%

%%%%%%1. The original pH and other values.%%%%%%
% input TA, TC
    Constants(TempCi,PRES);
    TI = DYEN;
    F=true(ntps,1);
if any(F)
    [PHic(F) FCic(F)] = CalculatepHfCO2fromTATC(TAc(F)-PengCorrection(F), TCc(F));
end
    PCic = FCic./FugFac;

% CalculateOtherParamsAtInputConditions:
[HCO3ip CO3ip BAlkip...
    OHip PAlkip...
    SiAlkip Hfreeip ...
    HSO4ip HFip DYEALKip]      = CalculateAlkParts(PHic, TCc);
PAlkip                         = PAlkip+PengCorrection;
CO2ip                          = TCc - CO3ip - HCO3ip;
[Revelleip]                    = RevelleFactor(TAc-PengCorrection, TCc);

% convert pH at input conditions to the other scales 
[pHicT pHicS pHicF pHicN]=FindpHOnAllScales(PHic);

% Merge the Ks at input into an array. Ks at output will be glued to this later.
KIVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KI1 KI2];

TVEC =[TB TF TS];


%%%%%%2. The effect of the dye addition%%%%%%%
%input Temperature, dye salinity, pH and Total alkalinity
   DC = DYEC./1e6;
   [H2I HI I2 TAd] = DyeSalinity(TempCi,SalI,pHI,DC,MCPCONSTANTS);
   TAnc  = (TAc.*(samplms./(samplms+dyems))) + (TAd.*(dyems./(samplms+dyems)));
   TCnc =  TCc.*(samplms./(samplms+dyems));
   Newsal =  (Sal.*(samplms./(samplms+dyems)))+(SalI.*(dyems./(samplms+dyems)));
   TI           =  DC.*(dyems./(samplms+dyems));
   Sal          =  Newsal;
  sqrSal        = sqrt(Sal);
  
  
%%%%%%3. The new pH and other variables
    Constants(TempCi,PRES);
    F=true(ntps,1); 
    if any(F)
    [PHnc(F) FCnc(F)] = CalculatepHfCO2fromTATC(TAnc(F)-PengCorrection(F), TCnc(F));
    end

% CalculateOtherParamsAtInputConditions:
[HCO3np CO3np BAlknp...
    OHnp PAlknp...
    SiAlknp Hfreenp ...
    HSO4np HFnp DYEALKnp]      = CalculateAlkParts(PHnc, TCnc);
    PAlknp                     = PAlknp+PengCorrection;
    CO2np                      = TCnc - CO3np - HCO3np;

% Just for reference, convert pH at input conditions to the other scales, too. 
[pHncT pHncS pHncF pHncN]=FindpHOnAllScales(PHnc);


%%%%%%4.calculate the difference between dye addition and without dye
%%%%%addition
deltapHT = pHncT - pHicT;


%%%%%%5. Calculate R from the pH values
[Ratio] = calculateRfrompH(pHncT,TempCi,SAL,MCPCONSTANTS);


% % %test if Newton's method iterated to a value within the desired accuracy
% % % itenumber= testCalculatepHfromTATC(TAc, TCc);
% % % if abs(itenumber) == 1000
% % %    warning('Reached maximum iteration. The desired accracy may not attained')
% % % end

% Saving data in array, 53 columns, as many rows as samples input
DATA=[deltapHT ...
      TAc*1e6        TCc*1e6        PCic*1e6      FCic*1e6...
      pHicT          pHicS          pHicF         pHicN...     
      HCO3ip*1e6     CO3ip*1e6      CO2ip*1e6     BAlkip*1e6...   
      OHip*1e6       PAlkip*1e6     SiAlkip*1e6   Hfreeip*1e6...  
      pHncT          pHncS          pHncF         pHncN ...
      HCO3np*1e6     CO3np*1e6      CO2np*1e6     BAlknp*1e6...
      OHnp*1e6       PAlknp*1e6     SiAlknp*1e6   DYEALKnp*1e6     Hfreenp*1e6...
      DYEC           H2I*1e6        HI*1e6        I2*1e6...
      TAd*1e6                                     Ratio... 
      DYEpH          DYES           SAMPM         DYEM...
      TEMPIN         SAL            PO4           SI...
      K1K2CONSTANTS  KSO4CONSTANTS  pHSCALEIN     MCPCONSTANTS...
      Revelleip      TVEC*1e6       TI*1e6 ];

HEADERS={'delta_pH';'TAlk';'TCO2';'pCO2in';'fCO2in';...
    'pHinTOTAL';'pHinSWS';'pHinFREE';'pHinNBS';...
    'HCO3in';'CO3in';'CO2in';'BAlkin';...
    'OHin';'PAlkin';'SiAlkin';'Hfreein';...
    'pHmixTOTAL';'pHmixSWS';'pHmixFREE';'pHmixNBS';...
    'HCO3mix';'CO3mix';'CO2mix';'BAlkmix';...
    'OHmix';'PAlkmix';'SiAlkmix';'IAlk'; 'Hfreemix';...
    'Dyeconcentration';'H2I';'HI';'I2';...
    'DyeTA'; 'R';...
    'DyepH';'DyeSal';'SampleMass';'DyeMass';...
    'Tempinput';'Salinput';'PO4';'SI';...
    'K1K2CONSTANTS';'KSO4CONSTANTS';'pHSCALEIN';'mCPCONSTANTS';...
    'Revelle';'TB';'TF';'TS';'TI';};

NICEHEADERS={...
    '01 - delta pH                     '; %The pH with dye addition - pH without dye addition
    '02 - TAlk             (umol/kgSW) '; %Input TA
    '03 - TCO2             (umol/kgSW) '; %Input DIC
    '04 - pCO2             (uatm)      '; %pCO2 calculated from Input TA and DIC
    '05 - fCO2             (uatm)      '; %fCO2 calculated from Input TA and DIC
    '06 - pHin(Total)                  ';
    '07 - pHin(SWS)                    ';
    '08 - pHin(FREE)                   ';
    '09 - pHin(NBS)                    ';
    '10 - HCO3in           (umol/kgSW) ';
    '11 - CO3in            (umol/kgSW) ';
    '12 - CO2in            (umol/kgSW) ';
    '13 - BAlkin           (umol/kgSW) ';
    '14 - OHin             (umol/kgSW) ';
    '15 - PAlkin           (umol/kgSW) ';
    '16 - SiAlkin          (umol/kgSW) ';
    '17 - Hfreein          (umol/kgSW) ';
    '18 - pHmix(Total)                 ';
    '19 - pHmix(SWS)                   ';
    '20 - pHmix(FREE                   ';
    '21 - pHmix(NBS)                   ';
    '22 - HCO3mix          (umol/kgSW) ';
    '23 - CO3mix           (umol/kgSW) ';
    '24 - CO2mix           (umol/kgSW) ';
    '25 - BAlmix           (umol/kgSW) ';
    '26 - OHmix            (umol/kgSW) ';
    '27 - PAlkmix          (umol/kgSW) ';
    '28 - SiAlmix          (umol/kgSW) ';
    '29 - IAlk             (umol/kgSW) '; %the dye alkalinity in the mixed solution
    '30 - Hfreemix         (umol/kgSW) ';
    '31 - Dye concentration(umol/kgSW) '; %the mCP stock concentration (before mixing)
    '32 - H2I              (umol/kgSW) '; %the mCP species concentration (before mixing)
    '33 - HI               (umol/kgSW) '; %the mCP species concentration (before mixing)
    '34 - I2               (umol/kgSW) '; %the mCP species concentration (before mixing)
    '35 - Dye TA           (umol/kgSW) '; %the TA of the dye (in the stock solution)
    '36 - Absorbance ratio             '; %the theoretical R with dye addition
    '37 - Dye pH                       '; %the mCP stock pH
    '38 - Dye Salinity                 '; %the mCP stock salinity
    '39 - Sample Mass                  ';
    '40 - Dye Mass                     ';
    '41 - Input temperature            '; 
    '42 - Sample salinity              '; %Input sample salinity
    '43 - PO4                          ';
    '44 - Si                           ';
    '45 - K1K2CONSTANTS    ()          ';
    '46 - KSO4CONSTANTS    ()          ';
    '47 - pHSCALEIN        ()          ';
    '48 - MCPCONSTANTS     ()          ';
    '49 - Revelle Factors  ()          ';
    '50 - TB               (umol/kgSW) ';
    '51 - TF               (umol/kgSW) ';
    '52 - TS               (umol/kgSW) ';
    '53 - TI_m             (umol/kgSW) ';} %TI_m, the concentration of mCP in the mixed solution

clear global pHScale WhichKs WhoseKSO4 Pbar
clear global Sal sqrSal TempK logTempK TempCi;
clear global FugFac VPFac PengCorrection ntps RGasConstant;
clear global fH RT;
clear global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KI2 KI1;
clear global TB TF TS TP TSi F TI pHI SalI samplms dyems;
	
end % end main function


%**************************************************************************
% Subroutines:
%**************************************************************************


function Constants(TempC,Pdbar)
global pHScale WhichKs WhoseKSO4 Whichmcp sqrSal Pbar RT;
global K0 fH FugFac VPFac ntps TempK logTempK;
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KI1 KI2;
global TB TF TS TP TSi RGasConstant Sal DYE;

% SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
% Outputs: K0, K(), T(), fH, FugFac, VPFac
% This finds the Constants of the CO2 system in seawater or freshwater,
% corrects them for pressure, and reports them on the chosen pH scale.
% The process is as follows: the Constants (except KS, KF which stay on the
% free scale - these are only corrected for pressure) are
%       1) evaluated as they are given in the literature
%       2) converted to the SWS scale in mol/kg-SW or to the NBS scale
%       3) corrected for pressure
%       4) converted to the SWS pH scale in mol/kg-SW
%       5) converted to the chosen pH scale
%
%       PROGRAMMER'S NOTE: all logs are log base e
%       PROGRAMMER'S NOTE: all Constants are converted to the pH scale
%               pHScale% (the chosen one) in units of mol/kg-SW
%               except KS and KF are on the free scale
%               and KW is in units of (mol/kg-SW)^2
TempK    = TempC + 273.15;
RT       = RGasConstant.*TempK;
logTempK = log(TempK);
Pbar     = Pdbar ./ 10;

% Generate empty vectors for holding results(They are all consertative elements)
TB = nan(ntps,1);
TF = nan(ntps,1);
TS = nan(ntps,1);

% CalculateTB - Total Borate:
F=(WhichKs==8); % Pure water case.
if any(F)
    TB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    TB(F) = 0.0004106.*Sal(F)./35; % in mol/kg-SW
    % this is .00001173.*Sali
    % this is about 1% lower than Uppstrom's value
    % Culkin, F., in Chemical Oceanography,
    % ed. Riley and Skirrow, 1965:
    % GEOSECS references this, but this value is not explicitly
    % given here
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8); % All other cases
if any(F)
	FF=F&(WhoseKSO4==1|WhoseKSO4==2); % If user opted for Uppstrom's values:
	if any(FF)
	    % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	    % this is .000416.*Sali./35. = .0000119.*Sali
		% TB(FF) = (0.000232./10.811).*(Sal(FF)./1.80655); % in mol/kg-SW
	    TB(FF) =  0.0004157.*Sal(FF)./35; % in mol/kg-SW
	end
	FF=F&(WhoseKSO4==3|WhoseKSO4==4); % If user opted for the new Lee values:
	if any(FF)
		% Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.	
	 	% Geochimica Et Cosmochimica Acta 74 (6): 1801â€“1811.
		TB(FF) =  0.0004326.*Sal(FF)./35; % in mol/kg-SW
	end
end

% CalculateTF;
% Riley, J. P., Deep-Sea Research 12:219-220, 1965:
% this is .000068.*Sali./35. = .00000195.*Sali
TF = (0.000067./18.998).*(Sal./1.80655); % in mol/kg-SW

% CalculateTS ;
% Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
% this is .02824.*Sali./35. = .0008067.*Sali
TS = (0.14./96.062).*(Sal./1.80655); % in mol/kg-SW

% CalculateK0:
% Weiss, R. F., Marine Chemistry 2:203-215, 1974.
TempK100  = TempK./100;
lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + Sal .*...
    (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
K0   = exp(lnK0);                  % this is in mol/kg-SW/atm

% CalculateIonS:
% This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
IonS         = 19.924 .* Sal ./ (1000 - 1.005   .* Sal);

% CalculateKS: (HSO4-)
lnKS = nan(ntps,1); pKS  = nan(ntps,1); KS   = nan(ntps,1);
F=(WhoseKSO4==1|WhoseKSO4==3);
if any(F)
    % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
    % The goodness of fit is .021.
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % TYPO on p. 121: the constant e9 should be e8.
    % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
  lnKS(F) = -4276.1./TempK(F) + 141.328 - 23.093.*logTempK(F) +...             
      (-13856./TempK(F) + 324.57 - 47.986.*logTempK(F)).*sqrt(IonS(F)) +...     
      (35474./TempK(F) - 771.54 + 114.723.*logTempK(F)).*IonS(F) +...           
      (-2698./TempK(F)).*sqrt(IonS(F)).*IonS(F) + (1776./TempK(F)).*IonS(F).^2; 
	KS(F) = exp(lnKS(F))...            % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005 .* Sal(F));   % convert to mol/kg-SW
end
F=(WhoseKSO4==2|WhoseKSO4==4);
if any(F)
    % Khoo, et al, Analytical Chemistry, 49(1):29-34, 1977
    % KS was found by titrations with a hydrogen electrode
    % of artificial seawater containing sulfate (but without F)
    % at 3 salinities from 20 to 45 and artificial seawater NOT
    % containing sulfate (nor F) at 16 salinities from 15 to 45,
    % both at temperatures from 5 to 40 deg C.
    % KS is on the Free pH scale (inherently so).
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % He finds log(beta) which = my pKS;
    % his beta is an association constant.
    % The rms error is .0021 in pKS, or about .5% in KS.
    % This is equation 20 on p. 33:
    pKS(F) = 647.59 ./ TempK(F) - 6.3451 + 0.019085.*TempK(F) - 0.5208.*sqrt(IonS(F));
    KS(F) = 10.^(-pKS(F))...          % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005.*Sal(F));    % convert to mol/kg-SW
end

% CalculateKF:
% Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
lnKF = 1590.2./TempK - 12.641 + 1.525.*IonS.^0.5;
KF   = exp(lnKF)...                 % this is on the free pH scale in mol/kg-H2O
    .*(1 - 0.001005.*Sal);          % convert to mol/kg-SW
% Another expression exists for KF: Perez and Fraga 1987. Not used here since ill defined for low salinity. (to be used for S: 10-40, T: 9-33)
% Nonetheless, P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
% lnKF = 874./TempK - 9.68 + 0.111.*Sal.^0.5; 
% KF   = exp(lnKF);                   % this is on the free pH scale in mol/kg-SW

% CalculatepHScaleConversionFactors:
%       These are NOT pressure-corrected.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

% CalculatefH
fH = nan(ntps,1);
% Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
F=(WhichKs==8);
if any(F)
    fH(F) = 1; % this shouldn't occur in the program for this case
end
F=(WhichKs==7);
if any(F)
    fH(F) = 1.29 - 0.00204.*  TempK(F) + (0.00046 -...
        0.00000148.*TempK(F)).*Sal(F).*Sal(F);
    % Peng et al, Tellus 39B:439-458, 1987:
    % They reference the GEOSECS report, but round the value
    % given there off so that it is about .008 (1%) lower. It
    % doesn't agree with the check value they give on p. 456.
end
F=(WhichKs~=7 & WhichKs~=8);
if any(F)
    fH(F) = 1.2948 - 0.002036.*TempK(F) + (0.0004607 -...
        0.000001475.*TempK(F)).*Sal(F).^2;
    % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
    % v. 3, 1982 (p. 80);
end

% CalculateKB:
KB      = nan(ntps,1); logKB   = nan(ntps,1);
lnKBtop = nan(ntps,1); lnKB    = nan(ntps,1);
F=(WhichKs==8); % Pure water case
if any(F)
    KB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % This is for GEOSECS and Peng et al.
    % Lyman, John, UCLA Thesis, 1957
    % fit by Li et al, JGR 74:5507-5525, 1969:
    logKB(F) = -9.26 + 0.00886.*Sal(F) + 0.01.*TempC(F);
    KB(F) = 10.^(logKB(F))...  % this is on the NBS scale
        ./fH(F);               % convert to the SWS scale
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
    lnKBtop(F) = -8966.9 - 2890.53.*sqrSal(F) - 77.942.*Sal(F) +...
        1.728.*sqrSal(F).*Sal(F) - 0.0996.*Sal(F).^2;
    lnKB(F) = lnKBtop(F)./TempK(F) + 148.0248 + 137.1942.*sqrSal(F) +...
        1.62142.*Sal(F) + (-24.4344 - 25.085.*sqrSal(F) - 0.2474.*...
        Sal(F)).*logTempK(F) + 0.053105.*sqrSal(F).*TempK(F);
    KB(F) = exp(lnKB(F))...    % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);         % convert to SWS pH scale
end

% CalculateKW:
lnKW = nan(ntps,1); KW = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-79.2447 + 3298.72./TempK(F) + 12.0408.*logTempK(F)).*...
        sqrSal(F) - 0.019813.*Sal(F);
end
F=(WhichKs==8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    % refit data of Harned and Owen, The Physical Chemistry of
    % Electrolyte Solutions, 1958
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    % his check value of 1.6 umol/kg-SW should be 6.2
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-5.977 + 118.67./TempK(F) + 1.0495.*logTempK(F)).*...
        sqrSal(F) - 0.01615.*Sal(F);
end
KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
F=(WhichKs==6);
if any(F)
    KW(F) = 0; % GEOSECS doesn't include OH effects
end

% CalculateKP1KP2KP3KSi:
KP1      = nan(ntps,1); KP2      = nan(ntps,1);
KP3      = nan(ntps,1); KSi      = nan(ntps,1);
lnKP1    = nan(ntps,1); lnKP2    = nan(ntps,1);
lnKP3    = nan(ntps,1); lnKSi    = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    KP1(F) = 0.02;
    % Peng et al don't include the contribution from this term,
    % but it is so small it doesn't contribute. It needs to be
    % kept so that the routines work ok.
    % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
    % Limnology and Oceanography 12:243-252, 1967:
    % these are only for sals 33 to 36 and are on the NBS scale
    KP2(F) = exp(-9.039 - 1450./TempK(F))... % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    KP3(F) = exp(4.466 - 7276./TempK(F))...  % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
    % The Chemical Society (London), Special Publ. 17:751, 1964:
    KSi(F) = 0.0000000004...              % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
end
F=(WhichKs==6 | WhichKs==8);
if any(F)
    KP1(F) = 0; KP2(F) = 0; KP3(F) = 0; KSi(F) = 0;
    % Neither the GEOSECS choice nor the freshwater choice
    % include contributions from phosphate or silicate.
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
    % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
    % KSi was given on the SWS pH scale in molal units.
    lnKP1(F) = -4576.752./TempK(F) + 115.54 - 18.453.*logTempK(F) + (-106.736./TempK(F) +...
        0.69171).*sqrSal(F) + (-0.65643./TempK(F) - 0.01844).*Sal(F);
    KP1(F) = exp(lnKP1(F));
    lnKP2(F) = -8814.715./TempK(F) + 172.1033 - 27.927.*logTempK(F) + (-160.34./TempK(F) +...
        1.3566).*sqrSal(F) + (0.37335./TempK(F) - 0.05778).*Sal(F);
    KP2(F) = exp(lnKP2(F));
    lnKP3(F) = -3070.75./TempK(F) - 18.126 + (17.27039./TempK(F) + 2.81197).*sqrSal(F) +...
        (-44.99486./TempK(F) - 0.09984).*Sal(F);
    KP3(F) = exp(lnKP3(F));
    lnKSi(F) = -8904.2./TempK(F) + 117.4 - 19.334.*logTempK(F) + (-458.79./TempK(F) +...
        3.5913).*sqrt(IonS(F)) + (188.74./TempK(F) - 1.5998).*IonS(F) +...
        (-12.1652./TempK(F) + 0.07871).*IonS(F).^2;
    KSi(F) = exp(lnKSi(F))...                % this is on the SWS pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F));        % convert to mol/kg-SW
end

% CalculateKI1,KI2:
KI1      = nan(ntps,1); KI2      = nan(ntps,1);
PKI1    = nan(ntps,1);  PKI2    = nan(ntps,1);

F=(Whichmcp==1);
if any(F)
   PKI1(F) = -782.62./TempK(F) + 1.1131;
   KI1(F) = 10.^(PKI1(F)); % this is on the total pH scale in mol/kg-H2O
   KI1(F) = KI1(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
   KI1(F) = KI1(F)./SWStoTOT(F); % convert to SWS pH scale
    a0 = 1.08071477e3;
    a1 = -1.35394946e-1;
    a2 = -1.98063716e2;
    a3 = 6.31924397e1;
    a4 = -5.18141866;
    b1 = -2.66457425e4;
    b2 = 5.08796578e3;
    b3 = -1.62454827e3;
    b4 = 1.33276788e2;
    c1 = -1.89671212e2;
    c2 = 3.49038762e1;
    c3 = -1.11336508e1;
    c4 = 9.12761930e-1;
    d1 = 3.27430677e-1;
    d2 = -7.51448528e-4;
    d3 = 3.94838229e-4;
    d4 = -6.00237846e-2;
    d5 = 1.90997693e-2;
    d6 = -1.56396488e-3;
    PKI2e2(F)= a0 + (a1.*Sal(F).^0.5) + (a2.*Sal(F).^1.5) +(a3.*Sal(F).^2)+(a4.*Sal(F).^2.5)...
       + b1.*TempK(F).^-1+b2.*Sal(F).^1.5.*TempK(F).^-1+ b3.*Sal(F).^2.*TempK(F).^-1+b4.*Sal(F).^2.5.*TempK(F).^-1....
       + c1.*log(TempK(F))+c2.*Sal(F).^1.5.*log(TempK(F))+c3.*Sal(F).^2.*log(TempK(F))+c4.*Sal(F).^2.5.*log(TempK(F))...
       + d1.*TempK(F)+d2.*Sal(F).^0.5.*TempK(F) + d3.*Sal(F).*TempK(F) + d4.*Sal(F).^1.5.*TempK(F) + d5.*Sal(F).^2.*TempK(F) + d6.*Sal(F).^2.5.*TempK(F);
   e2(F) = ((2.22-2.306)/35)*Sal(F)+2.306;
   lge2(F) = log10(e2(F));
   %linear relationship between S = 35, e2 =2.22 and S =0,e2= 2.306
   PKI2(F) = PKI2e2(F) + lge2(F);
   KI2(F) = 10.^(-PKI2(F)); %%%%%I temporally assume it's e2 = 2.222
   KI2(F) = KI2(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
   KI2(F) = KI2(F)./SWStoTOT(F); % convert to SWS pH scale
end

F=(Whichmcp==2);
if any(F)
   PKI1(F) = -782.62./TempK(F) + 1.1131;
   KI1(F) = 10.^(PKI1(F)); % this is on the total pH scale in mol/kg-H2O
   KI1(F) = KI1(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
   KI1(F) = KI1(F)./SWStoTOT(F); % convert to SWS pH scale
    a = -246.64209+0.315971.*Sal(F)+2.8855.*10^-4.*Sal(F).^2;
    b = 7229.23864 - 7.098137.*Sal(F) - 0.057034.*Sal(F).^2;
    c = 44.493382 - 0.052711.*Sal(F);
    d = 0.0781344;
    PKI2e2(F) = a+b./TempK(F)+c.*log(TempK(F))-d.*TempK(F);
    e2(F) = ((2.22-2.306)/35)*Sal(F)+2.306;
    lge2(F) = log10(e2(F));
   %linear relationship between S = 35, e2 =2.22 and S =0,e2= 2.306
    PKI2(F) = PKI2e2(F) + lge2(F);
    KI2(F) = 10.^(-PKI2(F)); %%%%%I temporally assume it's e2 = 2.222
    KI2(F) = KI2(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
    KI2(F) = KI2(F)./SWStoTOT(F); % convert to SWS pH scale
end
  
F=(Whichmcp==3);
if any(F)
   PKI1(F) = -782.62./TempK(F) + 1.1131;
   KI1(F) = 10.^(PKI1(F)); % this is on the total pH scale in mol/kg-H2O
   KI1(F) = KI1(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
   KI1(F) = KI1(F)./SWStoTOT(F); % convert to SWS pH scale
    a = 2.0129398*1e3;
    b = 6.1409*1e-1;
    c = -5.024240*1e4;
    d = -3.543347*1e2;
    PKI2(F) = a + b.*TempK(F) + c./TempK(F) + d.*log(TempK(F));
    KI2(F) = 10.^(-PKI2(F));
    KI2(F) = KI2(F).*(1-0.001005.*Sal(F)); % convert to mol/kg-SW
    KI2(F) = KI2(F)./SWStoTOT(F); % convert to SWS pH scale
end

% CalculateK1K2:
logK1    = nan(ntps,1); lnK1     = nan(ntps,1);
pK1      = nan(ntps,1); K1       = nan(ntps,1);
logK2    = nan(ntps,1); lnK2     = nan(ntps,1);
pK2      = nan(ntps,1); K2       = nan(ntps,1);
F=(WhichKs==1);
if any(F)
    % ROY et al, Marine Chemistry, 44:249-267, 1993
    % (see also: Erratum, Marine Chemistry 45:337, 1994
    % and Erratum, Marine Chemistry 52:183, 1996)
    % Typo: in the abstract on p. 249: in the eq. for lnK1* the
    % last term should have S raised to the power 1.5.
    % They claim standard deviations (p. 254) of the fits as
    % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
    % They also claim (p. 258) 2s precisions of .004 in pK1 and
    % .006 in pK2. These are consistent, but Andrew Dickson
    % (personal communication) obtained an rms deviation of about
    % .004 in pK1 and .003 in pK2. This would be a 2s precision
    % of about 2% in K1 and 1.5% in K2.
    % T:  0-45  S:  5-45. Total Scale. Artificial sewater.
    % This is eq. 29 on p. 254 and what they use in their abstract:
    lnK1(F) = 2.83655 - 2307.1266./TempK(F) - 1.5529413.*logTempK(F) +...
        (-0.20760841 - 4.0484./TempK(F)).*sqrSal(F) + 0.08468345.*Sal(F) -...
        0.00654208.*sqrSal(F).*Sal(F);
    K1(F) = exp(lnK1(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
    % This is eq. 30 on p. 254 and what they use in their abstract:
    lnK2(F) = -9.226508 - 3351.6106./TempK(F) - 0.2005743.*logTempK(F) +...
        (-0.106901773 - 23.9722./TempK(F)).*sqrSal(F) + 0.1130822.*Sal(F) -...
        0.00846934.*sqrSal(F).*Sal(F);
    K2(F) = exp(lnK2(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
end
F=(WhichKs==2);
if any(F)
    % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
    % The 2s precision in pK1 is .011, or 2.5% in K1.
    % The 2s precision in pK2 is .02, or 4.5% in K2.
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK1(F) = 812.27./TempK(F) + 3.356 - 0.00171.*Sal(F).*logTempK(F)...
        + 0.000091.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK2(F) = 1450.87./TempK(F) + 4.604 - 0.00385.*Sal(F).*logTempK(F)...
        + 0.000182.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==3);
if any(F)
    % HANSSON refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
    % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
    % on the SWS pH scale in mol/kg-SW.
    % Hansson gave his results on the Total scale (he called it
    % the seawater scale) and in mol/kg-SW.
    % Typo in DM on p. 1739 in Table 4: the equation for pK2*
    % for Hansson should have a .000132 *S^2
    % instead of a .000116 *S^2.
    % The 2s precision in pK1 is .013, or 3% in K1.
    % The 2s precision in pK2 is .017, or 4.1% in K2.
    % This is from Table 4 on p. 1739.
    pK1(F) = 851.4./TempK(F) + 3.237 - 0.0106.*Sal(F) + 0.000105.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is from Table 4 on p. 1739.
    pK2(F) = -3885.4./TempK(F) + 125.844 - 18.141.*logTempK(F)...
        - 0.0192.*Sal(F) + 0.000132.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==4);
if any(F)
    % MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Mehrbach et al gave results on the NBS scale.
    % The 2s precision in pK1 is .011, or 2.6% in K1.
    % The 2s precision in pK2 is .020, or 4.6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 4 on p. 1739.
    pK1(F) = 3670.7./TempK(F) - 62.008 + 9.7944.*logTempK(F)...
             - 0.0118.*Sal(F) + 0.000116.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 4 on p. 1739.
    pK2(F) = 1394.7./TempK(F) + 4.777 - 0.0184.*Sal(F) + 0.000118.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==5);
if any(F)
    % HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
    % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
    % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Typo in DM on p. 1740 in Table 5: the second equation
    % should be pK2* =, not pK1* =.
    % The 2s precision in pK1 is .017, or 4% in K1.
    % The 2s precision in pK2 is .026, or 6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 5 on p. 1740.
    pK1(F) = 845./TempK(F) + 3.248 - 0.0098.*Sal(F) + 0.000087.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1740.
    pK2(F) = 1377.3./TempK(F) + 4.824 - 0.0185.*Sal(F) + 0.000122.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
    % Limnology and Oceanography, 18(6):897-907, 1973.
	% I.e., these are the original Mehrbach dissociation constants.
    % The 2s precision in pK1 is .005, or 1.2% in K1.
    % The 2s precision in pK2 is .008, or 2% in K2.
    pK1(F) = - 13.7201 + 0.031334.*TempK(F) + 3235.76./TempK(F)...
        + 1.3e-5*Sal(F).*TempK(F) - 0.1032.*Sal(F).^0.5;
    K1(F) = 10.^(-pK1(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
    pK2(F) = 5371.9645 + 1.671221.*TempK(F) + 0.22913.*Sal(F) + 18.3802.*log10(Sal(F))...
             - 128375.28./TempK(F) - 2194.3055.*log10(TempK(F)) - 8.0944e-4.*Sal(F).*TempK(F)...
             - 5617.11.*log10(Sal(F))./TempK(F) + 2.136.*Sal(F)./TempK(F); % pK2 is not defined for Sal=0, since log10(0)=-inf
    K2(F) = 10.^(-pK2(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
end
F=(WhichKs==8);
if any(F)	
	% PURE WATER CASE
    % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
    % K1 from refit data from Harned and Davis,
    % J American Chemical Society, 65:2030-2037, 1943.
    % K2 from refit data from Harned and Scholes,
    % J American Chemical Society, 43:1706-1709, 1941.
	% This is only to be used for Sal=0 water (note the absence of S in the below formulations)
    % These are the thermodynamic Constants:
    lnK1(F) = 290.9097 - 14554.21./TempK(F) - 45.0575.*logTempK(F);
    K1(F) = exp(lnK1(F));
    lnK2(F) = 207.6548 - 11843.79./TempK(F) - 33.6485.*logTempK(F);
    K2(F) = exp(lnK2(F));
end
F=(WhichKs==9);
if any(F)
    % From Cai and Wang 1998, for estuarine use.
	% Data used in this work is from:
	% K1: Merhback (1973) for S>15, for S<15: Mook and Keone (1975)
	% K2: Merhback (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	% Sigma of residuals between fits and above data: Â±0.015, +0.040 for K1 and K2, respectively.
	% Sal 0-40, Temp 0.2-30
    % Limnol. Oceanogr. 43(4) (1998) 657-668
	% On the NBS scale
	% Their check values for F1 don't work out, not sure if this was correctly published...
	F1 = 200.1./TempK(F) + 0.3220;
	pK1(F) = 3404.71./TempK(F) + 0.032786.*TempK(F) - 14.8435 - 0.071692.*F1.*Sal(F).^0.5 + 0.0021487.*Sal(F);
    K1(F)  = 10.^-pK1(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential);
	F2 = -129.24./TempK(F) + 1.4381;
	pK2(F) = 2902.39./TempK(F) + 0.02379.*TempK(F) - 6.4980 - 0.3191.*F2.*Sal(F).^0.5 + 0.0198.*Sal(F);
    K2(F)  = 10.^-pK2(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential); 
end
F=(WhichKs==10);
if any(F)
    % From Lueker, Dickson, Keeling, 2000
	% This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
    % Mar. Chem. 70 (2000) 105-119
    % Total scale and kg-sw
    pK1(F) = 3633.86./TempK(F)-61.2172+9.6777.*log(TempK(F))-0.011555.*Sal(F)+0.0001152.*Sal(F).^2;
	K1(F)  = 10.^-pK1(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
    pK2(F) = 471.78./TempK(F)+25.929 -3.16967.*log(TempK(F))-0.01781 .*Sal(F)+0.0001122.*Sal(F).^2;
	K2(F)  = 10.^-pK2(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end
F=(WhichKs==11);
if any(F)
	% Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
	% sigma for pK1 is reported to be 0.0056
	% sigma for pK2 is reported to be 0.010
	% This is from the abstract and pages 2536-2537
    pK1 =  -43.6977 - 0.0129037.*Sal(F) + 1.364e-4.*Sal(F).^2 + 2885.378./TempK(F) +  7.045159.*log(TempK(F));
    pK2 = -452.0940 + 13.142162.*Sal(F) - 8.101e-4.*Sal(F).^2 + 21263.61./TempK(F) + 68.483143.*log(TempK(F))...
				+ (-581.4428.*Sal(F) + 0.259601.*Sal(F).^2)./TempK(F) - 1.967035.*Sal(F).*log(TempK(F));
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==12);
if any(F)
	% Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	% Calculated from overdetermined WOCE-era field measurements 
	% sigma for pK1 is reported to be 0.005
	% sigma for pK2 is reported to be 0.008
	% This is from page 1715
    pK1 =  6.359 - 0.00664.*Sal(F) - 0.01322.*TempC(F) + 4.989e-5.*TempC(F).^2;
    pK2 =  9.867 - 0.01314.*Sal(F) - 0.01904.*TempC(F) + 2.448e-5.*TempC(F).^2;
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==13);
if any(F)
    % From Millero 2006 work on pK1 and pK2 from titrations
	% Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
    % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	pK1_0 = -126.34048 + 6320.813./TempK(F) + 19.568224*log(TempK(F));
	A_1   = 13.4191*Sal(F).^0.5 + 0.0331.*Sal(F) - 5.33e-5.*Sal(F).^2;
	B_1   = -530.123*Sal(F).^0.5 - 6.103.*Sal(F);
	C_1   = -2.06950.*Sal(F).^0.5;
	pK1(F)= A_1 + B_1./TempK(F) + C_1.*log(TempK(F)) + pK1_0; % pK1 sigma = 0.0054
    K1(F) = 10.^-(pK1(F));
	pK2_0= -90.18333 + 5143.692./TempK(F) + 14.613358*log(TempK(F));	
	A_2   = 21.0894*Sal(F).^0.5 + 0.1248.*Sal(F) - 3.687e-4.*Sal(F).^2;
	B_2   = -772.483*Sal(F).^0.5 - 20.051.*Sal(F);
	C_2   = -3.3336.*Sal(F).^0.5;
	pK2(F)= A_2 + B_2./TempK(F) + C_2.*log(TempK(F)) + pK2_0; %pK2 sigma = 0.011
    K2(F) = 10.^-(pK2(F));
end
F=(WhichKs==14);
if any(F)
    % From Millero, 2010, also for estuarine use.
	% Marine and Freshwater Research, v. 61, p. 139â€“142.
	% Fits through compilation of real seawater titration results:
	% Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	% Constants for K's on the SWS;
	% This is from page 141
	pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	% This is from their table 2, page 140.
	A1 = 13.4038.*Sal(F).^0.5 + 0.03206.*Sal(F) - 5.242e-5.*Sal(F).^2;
	B1 = -530.659.*Sal(F).^0.5 - 5.8210.*Sal(F);
	C1 = -2.0664*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1;
	% This is from page 141
	pK20 =  -90.18333 + 5143.692./TempK(F) + 14.613358.*log(TempK(F));
	% This is from their table 3, page 140.
	A2 = 21.3728.*Sal(F).^0.5 + 0.1218.*Sal(F) - 3.688e-4.*Sal(F).^2;
	B2 = -788.289.*Sal(F).^0.5 - 19.189.*Sal(F);
	C2 = -3.374.*Sal(F).^0.5;
	pK2 = pK20 + A2 + B2./TempK(F) + C2.*log(TempK(F));
	K2(F) = 10.^-pK2;
end

%***************************************************************************
%CorrectKsForPressureNow:
% Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
%       the free scale) are on the SWS scale.
%       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
%       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
%       (the pH scales are the same in this case); the other Ks don't matter.
%
%
% No salinity dependence is given for the pressure coefficients here.
% It is assumed that the salinity is at or very near Sali = 35.
% These are valid for the SWS pH scale, but the difference between this and
% the total only yields a difference of .004 pH units at 1000 bars, much
% less than the uncertainties in the values.
%****************************************************************************
% The sources used are:
% Millero, 1995:
%       Millero, F. J., Thermodynamics of the carbon dioxide system in the
%       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
%       See table 9 and eqs. 90-92, p. 675.
%       TYPO: a factor of 10^3 was left out of the definition of Kappa
%       TYPO: the value of R given is incorrect with the wrong units
%       TYPO: the values of the a's for H2S and H2O are from the 1983
%                values for fresh water
%       TYPO: the value of a1 for B(OH)3 should be +.1622
%        Table 9 on p. 675 has no values for Si.
%       There are a variety of other typos in Table 9 on p. 675.
%       There are other typos in the paper, and most of the check values
%       given don't check.
% Millero, 1992:
%       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
%       CRC Press, 1992. See chapter 6.
%       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
%               79, and 96 have typos).
% Millero, 1983:
%       Millero, Frank J., Influence of pressure on chemical processes in
%       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
%       Chester, R., Academic Press, 1983.
%       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
%       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
%       these two are necessary to match the values given in Table 43.24
% Millero, 1979:
%       Millero, F. J., The thermodynamics of the carbon dioxide system
%       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
%       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
% Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
%       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
%       This matches the GEOSECS results and is in Edmond and Gieskes.
% Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
%       boric acid, and the pH of seawater, Limnology and Oceanography
%       13:403-417, 1968.
% Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
%       seawater with respect to calcium carbonate under in situ conditions,
%       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
%****************************************************************************
% These references often disagree and give different fits for the same thing.
% They are not always just an update either; that is, Millero, 1995 may agree
%       with Millero, 1979, but differ from Millero, 1983.
% For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
%       KP3, and KSi as for the other cases. Peng et al didn't consider the
%       case of P different from 0. GEOSECS did consider pressure, but didn't
%       include Phos, Si, or OH, so including the factors here won't matter.
% For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
%       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
%       including the factors won't matter.
%****************************************************************************
%       deltaVs are in cm3/mole
%       Kappas are in cm3/mole/bar
%****************************************************************************

%CorrectK1K2KBForPressure:
deltaV    = nan(ntps,1); Kappa     = nan(ntps,1);
lnK1fac   = nan(ntps,1); lnK2fac   = nan(ntps,1);
lnKBfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    %***PressureEffectsOnK1inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -30.54 + 0.1849 .*TempC(F) - 0.0023366.*TempC(F).^2;
    Kappa(F)   = (-6.22 + 0.1368 .*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
    lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %***PressureEffectsOnK2inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -29.81 + 0.115.*TempC(F) - 0.001816.*TempC(F).^2;
    Kappa(F)   = (-5.74 + 0.093.*TempC(F) - 0.001896.*TempC(F).^2)./1000;
    lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    lnKBfac(F) = 0 ;%; this doesn't matter since TB = 0 for this case
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
    %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
    %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
    %               but the fits are the same as those in
    %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
    %               who in turn quote Li, personal communication
    lnK1fac(F) = (24.2 - 0.085.*TempC(F)).*Pbar(F)./RT(F);
    lnK2fac(F) = (16.4 - 0.04 .*TempC(F)).*Pbar(F)./RT(F);
    %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
    %               and matches the GEOSECS results
    lnKBfac(F) = (27.5 - 0.095.*TempC(F)).*Pbar(F)./RT(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    %***PressureEffectsOnK1:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -25.5 + 0.1271.*TempC(F);
    %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (-3.08 + 0.0877.*TempC(F))./1000;
    %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
 	lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fits given in Millero, 1983 are somewhat different.
    
    %***PressureEffectsOnK2:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -15.82 - 0.0219.*TempC(F);
    %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (1.13 - 0.1475.*TempC(F))./1000;
    %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
	lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fit given in Millero, 1983 is different.
    %               Not by a lot for deltaV, but by much for Kappa. %
    
    %***PressureEffectsOnKB:
    %               This is from Millero, 1979.
    %               It is from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -29.48 + 0.1622.*TempC(F) - 0.002608.*TempC(F).^2;
    %               Millero, 1983 has:
    %                 'deltaV = -28.56 + .1211.*TempCi - .000321.*TempCi.*TempCi
    %               Millero, 1992 has:
    %                 'deltaV = -29.48 + .1622.*TempCi + .295.*(Sali - 34.8)
    %               Millero, 1995 has:
    %                 'deltaV = -29.48 - .1622.*TempCi - .002608.*TempCi.*TempCi
    %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = -2.84./1000; % Millero, 1979
    %               Millero, 1992 and Millero, 1995 also have this.
    %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
    %               Millero, 1983 has:
    %                 'Kappa = (-3 + .0427.*TempCi)./1000
    lnKBfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
end

% CorrectKWForPressure:
lnKWfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    % PressureEffectsOnKWinFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  =  -25.6 + 0.2324.*TempC(F) - 0.0036246.*TempC(F).^2;
    Kappa(F)   = (-7.33 + 0.1368.*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
 	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);

    %               NOTE the temperature dependence of KappaK1 and KappaKW
    %               for fresh water in Millero, 1983 are the same.
end
F=(WhichKs~=8);
if any(F)
    % GEOSECS doesn't include OH term, so this won't matter.
    % Peng et al didn't include pressure, but here I assume that the KW correction
    %       is the same as for the other seawater cases.
    % PressureEffectsOnKW:
    %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
    deltaV(F)  = -20.02 + 0.1119.*TempC(F) - 0.001409.*TempC(F).^2;
    %               Millero, 1992 and Millero, 1995 have:
    Kappa(F)   = (-5.13 + 0.0794.*TempC(F))./1000; % Millero, 1983
    %               Millero, 1995 has this too, but Millero, 1992 is different.
	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               Millero, 1979 does not list values for these.
end

% PressureEffectsOnKF:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KF is on the free pH scale.
deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;
Kappa = (-3.91 + 0.054.*TempC)./1000;
lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKS:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KS is on the free pH scale.
deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;
Kappa = (-4.53 + 0.09.*TempC)./1000;
lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKP1KP2KP3KSiForPressure:
% These corrections don't matter for the GEOSECS choice (WhichKs% = 6) and
%       the freshwater choice (WhichKs% = 8). For the Peng choice I assume
%       that they are the same as for the other choices (WhichKs% = 1 to 5).
% The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
%       same as Millero, 1983.
% PressureEffectsOnKP1:
deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;
Kappa  = (-2.67 + 0.0427.*TempC)./1000;
lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP2:
deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;
Kappa  = (-5.15 + 0.09  .*TempC)./1000;
lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP3:
deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;
Kappa  = (-4.08 + 0.0714.*TempC)./1000;
lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKSi:
%  The only mention of this is Millero, 1995 where it is stated that the
%    values have been estimated from the values of boric acid. HOWEVER,
%    there is no listing of the values in the table.
%    I used the values for boric acid from above.
deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
Kappa  = -2.84./1000;
lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKsForPressureHere:
K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
KSifac = exp(lnKSifac); KSi = KSi.*KSifac;

% CorrectpHScaleConversionsForPressure:
% fH has been assumed to be independent of pressure.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

%  The values KS and KF are already pressure-corrected, so the pH scale
%  conversions are now valid at pressure.

% FindpHScaleConversionFactor:
% this is the scale they will be put on
pHfactor  = nan(ntps,1);
F=(pHScale==1); %Total
pHfactor(F) = SWStoTOT(F);
F=(pHScale==2); %SWS, they are all on this now
pHfactor(F) = 1;
F=(pHScale==3); %pHfree
pHfactor(F) = SWStoTOT(F)./FREEtoTOT(F);
F=(pHScale==4); %pHNBS
pHfactor(F) = fH(F);

% ConvertFromSWSpHScaleToChosenScale:
K1  = K1.* pHfactor; K2  = K2.* pHfactor;
KW  = KW.* pHfactor; KB  = KB.* pHfactor;
KP1 = KP1.*pHfactor; KP2 = KP2.*pHfactor;
KP3 = KP3.*pHfactor; KSi = KSi.*pHfactor;
KI1 = KI1.*pHfactor; KI2 = KI2.*pHfactor;
% CalculateFugacityConstants:
% This assumes that the pressure is at one atmosphere, or close to it.
% Otherwise, the Pres term in the exponent affects the results.
%       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
%       Delta and B in cm3/mol
FugFac=ones(ntps,1);
Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
% For a mixture of CO2 and air at 1 atm (at low CO2 concentrations);
P1atm = 1.01325; % in bar
FugFac = exp((b + 2.*Delta).*P1atm./RT);
F=(WhichKs==6 | WhichKs==7); % GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
FugFac(F) = 1;
% CalculateVPFac:
% Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
%       seawater, Marine Chemistry 8:347-359, 1980.
% They fit the data of Goff and Gratch (1946) with the vapor pressure
%       lowering by sea salt as given by Robinson (1954).
% This fits the more complicated Goff and Gratch, and Robinson equations
%       from 273 to 313 deg K and 0 to 40 Sali with a standard error
%       of .015%, about 5 uatm over this range.
% This may be on IPTS-29 since they didn't mention the temperature scale,
%       and the data of Goff and Gratch came before IPTS-48.
% The references are:
% Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
%       to 212 deg F, Transactions of the American Society of Heating and
%       Ventilating Engineers 52:95-122, 1946.
% Robinson, Journal of the Marine Biological Association of the U. K.
%       33:449-455, 1954.
%       This is eq. 10 on p. 350.
%       This is in atmospheres.
VPWP = exp(24.4543 - 67.4509.*(100./TempK) - 4.8489.*log(TempK./100));
VPCorrWP = exp(-0.000544.*Sal);
VPSWWP = VPWP.*VPCorrWP;
VPFac = 1 - VPSWWP; % this assumes 1 atmosphere
end % end nested function (From CO2SYS, add dye constants)

function varargout=CalculatepHfCO2fromTATC(TAx, TCx)
global FugFac F;
% Outputs pH fCO2, in that order
% SUB FindpHfCO2fromTATC, version 01.02, 10-10-97, written by Ernie Lewis.
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, TA, TC, Sal, K(), T(), TempC, Pdbar
% Outputs: pH, fCO2
% This calculates pH and fCO2 from TA and TC at output conditions.
pHx   = CalculatepHfromTATC(TAx, TCx); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
fCO2x = CalculatefCO2fromTCpH(TCx, pHx);
varargout{1} = pHx;
varargout{2} = fCO2x;
end % end nested function (From CO2SYS)

function varargout=CalculatepHfromTATC(TAx, TCx)
global pHScale  WhichKs  WhoseKSO4 sqrSal Pbar RT;
global K0 fH FugFac VPFac ntps TempK logTempK;
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KI1 KI2;
global TB TF TS TP TSi F TI;
%Outputs pH
% SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis.
% Inputs: TA, TC, K(), T()
% Output: pH
% This calculates pH from TA and TC using K1 and K2 by Newton's method.
% It tries to solve for the pH at which Residual = 0.
% The starting guess is pH = 8.
% Though it is coded for H on the total pH scale, for the pH values occuring
% in seawater (pH > 6) it will be equally valid on any pH scale (H terms
% negligible) as long as the K Constants are on that scale.
%
% Made this to accept vectors. It will continue iterating until all
% values in the vector are "abs(deltapH) < pHTol". SVH2007
K1F=K1(F);   K2F=K2(F);   KWF =KW(F);
KP1F=KP1(F); KP2F=KP2(F); KP3F=KP3(F);  TPF=TP(F); KI2F = KI2(F);
TSiF=TSi(F); KSiF=KSi(F); TBF =TB(F);   KBF=KB(F); KI1F = KI1(F);
TSF =TS(F);  KSF =KS(F);  TFF =TF(F);   KFF=KF(F); TIF = TI(F);
vl          = sum(F);  % VectorLength
pHGuess     = 8;       % this is the first guess
pHTol       = 0.0001;  % tolerance for iterations end
ln10        = log(10); %
pHx(1:vl,1) = pHGuess; % creates a vector holding the first guess for all samples
deltapH     = pHTol+1;
while any(abs(deltapH) > pHTol)
    H         = 10.^(-pHx);
    Denom     = (H.*H + K1F.*H + K1F.*K2F);
    CAlk      = TCx.*K1F.*(H + 2.*K2F)./Denom;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    DyeTop    = KI1F.* KI2F - H.*H;
    DyeBot    = H.*H + KI1F.*H + KI1F.*KI2F;
    DAlk      = TIF.* DyeTop./DyeBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    FREEtoTOT = (1 + TSF./KSF); % pH scale conversion factor
    Hfree     = H./FREEtoTOT; % for H on the total scale
    HSO4      = TSF./(1 + KSF./Hfree); % since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); % since KF is on the free scale
    Residual  = TAx - CAlk - BAlk - OH - PAlk - SiAlk - DAlk + Hfree + HSO4 + HF;
    % find Slope dTA/dpH;
    % (this is not exact, but keeps all important terms);
    Slope     = ln10.*(TCx.*K1F.*H.*(H.*H + K1F.*K2F + 4.*H.*K2F)./Denom./Denom + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; % this is Newton's method
    % to keep the jump from being too big;
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pHx       = pHx + deltapH; % Is on the same scale as K1 and K2 were calculated...
end
varargout{1}=pHx;
end % end nested function

function varargout=CalculatefCO2fromTCpH(TCx, pHx)
global K0 K1 K2 F
% ' SUB CalculatefCO2fromTCpH, version 02.02, 12-13-96, written by Ernie Lewis.
% ' Inputs: TC, pH, K0, K1, K2
% ' Output: fCO2
% ' This calculates fCO2 from TC and pH, using K0, K1, and K2.
H            = 10.^(-pHx);
fCO2x        = TCx.*H.*H./(H.*H + K1(F).*H + K1(F).*K2(F))./K0(F);
varargout{1} = fCO2x;
end % end nested function (From CO2SYS)

function varargout=RevelleFactor(TAi, TCi)
% global WhichKs;
% ' SUB RevelleFactor, version 01.03, 01-07-97, written by Ernie Lewis.
% ' Inputs: WhichKs%, TA, TC, K0, K(), T()
% ' Outputs: Revelle
% ' This calculates the Revelle factor (dfCO2/dTC)|TA/(fCO2/TC).
% ' It only makes sense to talk about it at pTot = 1 atm, but it is computed
% '       here at the given K(), which may be at pressure <> 1 atm. Care must
% '       thus be used to see if there is any validity to the number computed.
TC0 = TCi;
dTC = 0.000001;% ' 1 umol/kg-SW
% ' Find fCO2 at TA, TC + dTC
TCi = TC0 + dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2plus = fCO2c;
% ' Find fCO2 at TA, TC - dTC
TCi = TC0 - dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2minus = fCO2c;
% CalculateRevelleFactor:
Revelle = (fCO2plus - fCO2minus)./dTC./((fCO2plus + fCO2minus)./TCi);
varargout{1}=Revelle;
end % end nested function  (From CO2SYS)

function varargout = DyeSalinity(T,S,pHd,TII,Equation)
%%pH and H
HIon = 10.^(-1.*pHd);
%%%Temperature convert unit from C to K
T = T + 273.15;

if (Equation ==1);
%%Equation=1 is Muller's fresh water equation, for wider salinity range
    a0 = 1.08071477e3;
    a1 = -1.35394946e-1;
    a2 = -1.98063716e2;
    a3 = 6.31924397e1;
    a4 = -5.18141866;
    b1 = -2.66457425e4;
    b2 = 5.08796578e3;
    b3 = -1.62454827e3;
    b4 = 1.33276788e2;
    c1 = -1.89671212e2;
    c2 = 3.49038762e1;
    c3 = -1.11336508e1;
    c4 = 9.12761930e-1;
    d1 = 3.27430677e-1;
    d2 = -7.51448528e-4;
    d3 = 3.94838229e-4;
    d4 = -6.00237846e-2;
    d5 = 1.90997693e-2;
    d6 = -1.56396488e-3;
PKI2e2(:,1) = a0 + (a1.*S.^0.5) + (a2.*S.^1.5) +(a3.*S.^2)+(a4.*S.^2.5)...
       + b1.*T.^-1+b2.*S.^1.5.*T.^-1+ b3.*S.^2.*T.^-1+b4.*S.^2.5.*T.^-1....
       + c1.*log(T)+c2.*S.^1.5.*log(T)+c3.*S.^2.*log(T)+c4.*S.^2.5.*log(T)...
       + d1.*T+d2.*S.^0.5.*T + d3.*S.*T + d4.*S.^1.5.*T + d5.*S.^2.*T + d6.*S.^2.5.*T;
%%%%snesitivity test for e2 value
e2 = ((2.22-2.306)/35).*S+2.306;
%%%get K2 term from PK2e2 term
lne2 = log10(e2);
PKI2 = PKI2e2 + lne2;
KI2 = 10.^(-1*PKI2);
%%%%PKI1 is the equation under 0.7 MNacl
PKI1 = - 782.62./T + 1.1131;
KI1 = 10.^PKI1;
Bottom = KI1.* KI2 + KI1.* HIon + HIon.* HIon;
H2I = TII.* HIon.* HIon./Bottom;
HI = TII.* KI1.* HIon./Bottom;
I2 = TII.* KI1.* KI2./Bottom;
end
% 
%%%Equation =2 is Liu's equation for purified dye
if (Equation ==2);
%%%constants
    a = -246.64209+0.315971.*S+2.8855.*10^-4.*S.^2;
    b = 7229.23864 - 7.098137.*S - 0.057034.*S.^2;
    c = 44.493382 - 0.052711.*S;
    d = 0.0781344;
    PKI2e2 = a+b./T+c.*log(T)-d.*T;
  
%%%%snesitivity test for e2 value
e2 = ((2.22-2.306)/35).*S+2.306;

%%%get K2 term from PK2e2 term
lne2 = log10(e2);
PKI2 = PKI2e2 + lne2  ;
KI2 = 10.^(-1*PKI2);

%%%%PKI1 is the equation under 0.7 MNacl
PKI1 = - 782.62./T + 1.1131;
KI1 = 10.^PKI1;
Bottom = KI1.* KI2 + KI1.* HIon + HIon.* HIon;
H2I = TII.* HIon.* HIon./Bottom;
HI = TII.* KI1.* HIon./Bottom;
I2 = TII.* KI1.* KI2./Bottom;
end

if (Equation ==3);  %Lai's equation for pure water only
    a = 2.0129398*1e3;
    b = 6.1409*1e-1;
    c = -5.024240*1e4;
    d = -3.543347*1e2;
    PKI2 = a + b.*T + c./T + d.*log(T);
  
%%%%snesitivity test for e2 value
e2 = 3.1400*1e-5.*(T.^2)-2.0527*1e-2.*T+5.6349;
KI2 = 10.^(-1*PKI2);
% 
%%%%PKI1 is the equation under 0.7 MNacl
PKI1 = - 782.62./T + 1.1131;                                   %Thsi is still from Liu et al
KI1 = 10.^PKI1;
Bottom = KI1.* KI2 + KI1.* HIon + HIon.* HIon;
H2I = TII.* HIon.* HIon./Bottom;
HI = TII.* KI1.* HIon./Bottom;
I2 = TII.* KI1.* KI2./Bottom;
end
OH_  =    10.^-(14-pHd);
TAd  =   (I2 - H2I - HIon + OH_);   %%(mol/kg)

%%output
varargout{1}= H2I;
varargout{2}= HI;
varargout{3}= I2;
varargout{4}= TAd;
end % end nested function

function varargout=CalculateAlkParts(pHx, TCx)
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KI1 KI2;
global TB TF TS TP TSi F TI;
% ' SUB CalculateAlkParts, version 01.03, 10-10-97, written by Ernie Lewis.
% ' Inputs: pH, TC, K(), T()
% ' Outputs: HCO3, CO3, BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF
% ' This calculates the various contributions to the alkalinity.
% ' Though it is coded for H on the total pH scale, for the pH values occuring
% ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
% ' negligible) as long as the K Constants are on that scale.
H         = 10.^(-pHx);
HCO3      = TCx.*K1.*H  ./(K1.*H + H.*H + K1.*K2);
CO3       = TCx.*K1.*K2./(K1.*H + H.*H + K1.*K2);
BAlk      = TB.*KB./(KB + H);
OH        = KW./H;
PhosTop   = KP1.*KP2.*H + 2.*KP1.*KP2.*KP3 - H.*H.*H;
PhosBot   = H.*H.*H + KP1.*H.*H + KP1.*KP2.*H + KP1.*KP2.*KP3;
PAlk      = TP.*PhosTop./PhosBot;
% this is good to better than .0006*TP:
% PAlk = TP*(-H/(KP1+H) + KP2/(KP2+H) + KP3/(KP3+H))
DyeTop    = KI1.* KI2 -H.*H;
DyeBot    = H.*H + KI1.*H + KI2.*KI1;
DAlk      = TI.* DyeTop./DyeBot;
SiAlk     = TSi.*KSi./(KSi + H);
FREEtoTOT = (1 + TS./KS);        %' pH scale conversion factor
Hfree     = H./FREEtoTOT;          %' for H on the total scale
HSO4      = TS./(1 + KS./Hfree); %' since KS is on the free scale
HF        = TF./(1 + KF./Hfree); %' since KF is on the free scale

varargout{1} = HCO3; varargout{2} = CO3;   varargout{3} = BAlk;  varargout{4} = OH;
varargout{5} = PAlk; varargout{6} = SiAlk; varargout{7} = Hfree; varargout{8} = HSO4;
varargout{9} = HF; varargout{10} = DAlk;
end % end nested function (From CO2SYS, add dye alkalinity term)

function varargout=FindpHOnAllScales(pH)
global pHScale K T TS KS TF KF fH ntps;
% ' SUB FindpHOnAllScales, version 01.02, 01-08-97, written by Ernie Lewis.
% ' Inputs: pHScale%, pH, K(), T(), fH
% ' Outputs: pHNBS, pHfree, pHTot, pHSWS
% ' This takes the pH on the given scale and finds the pH on all scales.
%  TS = T(3); TF = T(2);
%  KS = K(6); KF = K(5);% 'these are at the given T, S, P
FREEtoTOT = (1 + TS./KS);% ' pH scale conversion factor
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);% ' pH scale conversion factor
factor=nan(ntps,1);
F=pHScale==1;  %'"pHtot"
factor(F) = 0;
F=pHScale==2; % '"pHsws"
factor(F) = -log(SWStoTOT(F))./log(0.1);
F=pHScale==3; % '"pHfree"
factor(F) = -log(FREEtoTOT(F))./log(0.1);
F=pHScale==4;  %'"pHNBS"
factor(F) = -log(SWStoTOT(F))./log(0.1) + log(fH(F))./log(0.1);
pHtot  = pH    - factor;    % ' pH comes into this sub on the given scale
pHNBS  = pHtot - log(SWStoTOT) ./log(0.1) + log(fH)./log(0.1);
pHfree = pHtot - log(FREEtoTOT)./log(0.1);
pHsws  = pHtot - log(SWStoTOT) ./log(0.1);
varargout{1}=pHtot;
varargout{2}=pHsws;
varargout{3}=pHfree;
varargout{4}=pHNBS;
end % end nested function  (From CO2SYS)

function varagout =calculateRfrompH(spH,inputT,mS,Equation);

 mT = inputT + 273.15;

 if (Equation==1);
%Temperature 
    a0 = 1.08071477e3;
    a1 = -1.35394946e-1;
    a2 = -1.98063716e2;
    a3 = 6.31924397e1;
    a4 = -5.18141866;
    b1 = -2.66457425e4;
    b2 = 5.08796578e3;
    b3 = -1.62454827e3;
    b4 = 1.33276788e2;
    c1 = -1.89671212e2;
    c2 = 3.49038762e1;
    c3 = -1.11336508e1;
    c4 = 9.12761930e-1;
    d1 = 3.27430677e-1;
    d2 = -7.51448528e-4;
    d3 = 3.94838229e-4;
    d4 = -6.00237846e-2;
    d5 = 1.90997693e-2;
    d6 = -1.56396488e-3;
S_PKI2e2= a0 + (a1.*mS.^0.5) + (a2.*mS.^1.5) +(a3.*mS.^2)+(a4.*mS.^2.5)...
       + b1.*mT.^-1+b2.*mS.^1.5.*mT.^-1+ b3.*mS.^2.*mT.^-1+b4.*mS.^2.5.*mT.^-1....
       + c1.*log(mT)+c2.*mS.^1.5.*log(mT)+c3.*mS.^2.*log(mT)+c4.*mS.^2.5.*log(mT)...
       + d1.*mT+d2.*mS.^0.5.*mT + d3.*mS.*mT + d4.*mS.^1.5.*mT + d5.*mS.^2.*mT + d6.*mS.^2.5.*mT;
%%%Ratio with single volume dye
M = 10.^(spH - S_PKI2e2);
E1 = -0.007762 + 4.5174 * 1e-5.*mT;
E3_E2 = -0.020813 + 2.60262.*1e-4.*mT + 1.0436.*1e-4.*(mS -35);
R = (M + E1)./(1+M.*E3_E2);
varagout = R;
end

if (Equation==2);
    a = -246.64209+0.315971.*mS+2.8855.*10^-4.*mS.^2;
    b = 7229.23864 - 7.098137.*mS - 0.057034.*mS.^2;
    c = 44.493382 - 0.052711.*mS;
    d = 0.0781344;
    PKI2e2 = a+b./mT+c.*log(mT)-d.*mT;
    M = 10.^(spH - PKI2e2);
    E1 = -0.007762 + 4.5174 * 1e-5.*mT;
    E3_E2 = -0.020813 + 2.60262.*1e-4.*mT + 1.0436.*1e-4.*(mS -35);
    R = (M + E1)./(1+M.*E3_E2);
    varagout = R;
end

if (Equation ==3);
    a = 2.0129398*1e3;
    b = 6.1409*1e-1;
    c = -5.024240*1e4;
    d = -3.543347*1e2;
    PKI2 = a + b.*mT + c./mT + d.*log(mT);
    %%%%snesitivity test for e2 value
    E2 = 3.1400*1e-5.*(mT.^2)-2.0527*1e-2.*mT+5.6349;
    lne2 = log10(E2);
    E1 = 1.5036*1e-7.*(mT.^2)-5.8331*1e-5.*mT+1.0044*1e-2;
    E3 = -2.8764.*1e-6.*(mT.^2)+2.2697.*1e-3.*mT-2.9948.*1e-1;
    E3_E2 = E3./E2;
    PKI2e2 = PKI2-lne2;
    M = 10.^(spH - PKI2e2);
    R = (M + E1)./(1+M.*E3_E2);
    varagout = R;
end

end % end nested function

