%%%Example of how to use dyeperturbation function
%Details and descriptions please refer to the paper.
clc
clear

disp(' ')
disp('This is an example of the use of dyeperturbation.m')
disp(' ')
disp('We will generate a figure that shows the dye perturbation on pH')
disp(' ')

%%Approximate input conditions
Sal1 = [0:0.01:40];
Sal  = datasample(Sal1,10000);
TA   =  46.8537* Sal + 679.8349;
for n = [1:1:10000]
DIC_case1 = [TA(n)*0.8:1:TA(n)*1.5];
DIC(n)  = datasample(DIC_case1,1);
end

%sample properties
par1     =   TA;          % TA of the sample (in ummol/kg)
par2     =   DIC;          % DIC of the sample (in ummol/kg)
sal      =   Sal;              % Salinity of the sample
tempin   =    25;           % Temperature at input conditions, usually is the measurement temperature (in Celsius)
sil      =    0;          % Concentration of silicate  in the sample (in umol/kg)
po4      =    0;          % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1;          % pH scale at which the input pH is reported ("1" means "Total Scale") 
k1k2c    =    10;          % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c    =    1;          % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

%dye properties
dyec     = 2500;          % Concentration of the purified mCP dye stock (in umol/kg)
dyepH    =   8.0;          % pH of the purified mCP dye
dyeS     =   0;          % Salinity of the purifed mCP
SampM    =0.015;          % Mass of the sample (in unit kg)
DyeM     =2*1e-5;         % Mass of the dye(in unit kg)
mCPk     = 1;             %Choice of the mCP dissociation constants ("1" means Liu et al., 2011)

[DATA,HEADERS,NICEHEADERS]=dyeperturbation(par1,par2,sal,tempin,sil,po4,pHscale,k1k2c,kso4c,dyec,dyepH,dyeS,SampM,DyeM,mCPk);

IntpH = DATA(:,6);        % The intial calculated pH in total scale is in the 5 colum
dpH   = DATA(:,1);       % The pH differeneces between initial pH and pH with the dye addition
R     = DATA(:,36);       % The theoretical absorbance ratio from the pH with the dye addition


%%%%%Analyze
figure; clf
hold on
scatter(IntpH,dpH,15,TA,'filled') 
set(gca,'Fontsize',15)
set(gcf, 'Position',  [200, 200,800, 650])
xlim([7.0 8.5])
xlabel('seawater pH'); ylabel('  \Delta pH (seawater pH with dye addition - seawater pH)')
c = colorbar('northoutside')
c.Label.String = 'TA (\mumol/kg)';

MaxdpH = max(dpH);
MindpH = min(dpH);
ABS    = abs(dpH);
No_perturbation = IntpH(find(ABS == min(ABS)));

titles = [" Max dpH " " Min dpH " " 0 dpH "];
values = [MaxdpH MindpH No_perturbation];
result = [titles; values];
result = pad(result,'left');
result = join(result);
fprintf('%s\n', result);
text(8.0,MaxdpH-0.0005,['Dye S =' ,num2str(dyeS),' , ','Dye pH =', num2str(dyepH)],'Fontsize',15)
text(7.05,MindpH+0.0005,result,'Fontsize',15)
hold off


