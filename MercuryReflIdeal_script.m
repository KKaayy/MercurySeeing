% Supplementary material: Wohlfarth and Wöhler 2022
% Wavelength-dependent seeing systematically changes the normalized slope of telescopic reflectance spectra of Mercury.
% 13.01.22 K. Wohlfarth 
%% 1) Model Setup
fvx = linspace(-1.5,1.5,300);
fvy = linspace(-1.5,1.5,300);
fx = [-1.5 1.5];
fy = [-1.5 1.5];

selo = 0; % selo: sub-earth longitude
sela = 0; % sela: sub-earth latitude
ad = 7;   % ideal angular diameter
fvCenter = [sela 360 - selo];

fWLD = 0.570; % wavelength-difference 1.000 mum - 0.430  mum
reflmodel.model = 'amsa'; % Hapke 'amsa' or Kaasalainen-Shkuratov 'KS3'
switch reflmodel.model
    case 'amsa'
        load('fmRamsa.mat')
        fmR = fmRamsa;
    case 'KS3'
        load('fmRKS3.mat')
        fmR = fmRKS3;
end

%% Case 1: No seeing
% Reflectance
figure; 
himIdeal = imshow(fmR(:,:,1),[0 0.1],'XData',fx,'YData',fy);
colormap gray
hold on;
axesm ('ortho', 'Frame', 'off', 'Grid', 'on','GColor','white');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','white','LineWidth',2); 
axis equal
colorbar;
set(himIdeal, 'AlphaData', ~isnan(fmR(:,:,3)))
text(-1.4,-1.4,'$1 / sr$','Interpreter','Latex','Color','white')
text(-0.75,0,sprintf('X R0'));
[~,nxSample] = min(abs(fvx- (-0.7)));
[~,nySample] = min(abs(fvy- (0)));
R0 = fmR(nySample,nxSample,1);
disp(R0);

% Slope
fmIdealSlope = (fmR(:,:,1)-fmR(:,:,3))/fWLD;
figure; himIdealSlope = imshow(fmIdealSlope,[0 0.1],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'off', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
set(himIdealSlope, 'AlphaData', fmIdealSlope(:,:)>0)
text(-1.4,-1.4,'$1 / (sr\ \mu m)$','Interpreter','Latex')
text(-0.75,0,sprintf('X S0'));
[~,nxSample] = min(abs(fvx- (-0.70)));
[~,nySample] = min(abs(fvy- (0)));
S0 = fmIdealSlope(nySample,nxSample);
disp(S0);

% Normalized slope
fmIdealNormSlope = (fmR(:,:,1)-fmR(:,:,3))./fmR(:,:,3)/fWLD*100;
figure; himIdealNormSlope = imshow(fmIdealNormSlope,[200 300],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
set(himIdealNormSlope, 'AlphaData', fmIdealNormSlope>0)
text(-0.1,-0.90,sprintf('X\n A'));
[~,nxSample] = min(abs(fvx- (-0.04)));
[~,nySample] = min(abs(fvy- (-0.98)));
A = fmIdealNormSlope(nySample,nxSample);
disp(A);

text(-0.78,0,'X B');
[~,nxSample] = min(abs(fvx- (-0.72)));
[~,nySample] = min(abs(fvy- (0)));
B = fmIdealNormSlope(nySample,nxSample);
disp(B);
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

%% Case 2: Constant seeing
fwhm_to_sigma = 2*sqrt(2*log(2));
arcsec_ = 200/ad*1; % 1 arcsec in pixel 
seeing_factor = arcsec_/fwhm_to_sigma;

FWHM_seeing = 1.5;

sigma_1 = seeing_factor*FWHM_seeing; %seeing = 1.0 arcseconds
fmI_seeing = (imgaussfilt((fmR(:,1:300,:)),sigma_1));

%% 
% Reflectance
figure; 
himConstant = imshow(fmI_seeing(:,:,1),[0 0.1],'XData',fx,'YData',fy);
colormap gray
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','white');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','white','LineWidth',2); 
axis equal
colorbar; 
text(-1.4,-1.4,'$1/sr$','Interpreter','Latex','Color','white')
%set(himConstant, 'AlphaData', fmI_seeing1(:,:,2)>0)

% Slope
fmConstSlope = (fmI_seeing(:,:,1)-fmI_seeing(:,:,3))/fWLD;
figure; himConstantSlope = imshow(fmConstSlope,[0 0.1],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
%set(himConstantSlope, 'AlphaData', fmConstSlope>0)
set(himConstantSlope, 'AlphaData', fmI_seeing(:,:,3)>0.01*max(max(fmI_seeing(:,:,3))));
text(-1.4,-1.4,'$1 / (sr\ \mu m)$','Interpreter','Latex')
text(-0.8,0,sprintf('X S1'));
[~,nxSample] = min(abs(fvx- (-0.73)));
[~,nySample] = min(abs(fvy- (0)));
S1 = fmConstSlope(nySample,nxSample);
disp(S1);

% Normalized Slope
figure; 
fmConstNormSlope = (fmI_seeing(:,:,1)-fmI_seeing(:,:,3))./fmI_seeing(:,:,3)/fWLD*100;
himConstantNormSlope = imshow(fmConstNormSlope,[200 300],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
%set(himConstantNormSlope, 'AlphaData', ~isnan(fmConstNormSlope))
set(himConstantNormSlope, 'AlphaData', fmI_seeing(:,:,3)>0.01*max(max(fmI_seeing(:,:,3))));
text(0.03,-0.92,'X C');
[~,nxSample] = min(abs(fvx- (0.09)));
[~,nySample] = min(abs(fvy- (-0.9)));
C = fmConstNormSlope(nySample,nxSample);
disp(C);
text(-0.66,0,'X D');
[~,nxSample] = min(abs(fvx- (-0.60)));
[~,nySample] = min(abs(fvy- (0)));
D = fmConstNormSlope(nySample,nxSample);
disp(D);
text(-0.25,1.26,'X E');
[~,nxSample] = min(abs(fvx- (-0.19)));
[~,nySample] = min(abs(fvy- (1.28)));
E = fmConstNormSlope(nySample,nxSample);
disp(E);
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex');

%% Case 3: Wavelength-dependent seeing

fvMDISwl = [1000 750 430]*1e-6;
fwhm_to_sigma = 2*sqrt(2*log(2));
arcsec_ = 200/ad*1; % 1 arcsec in pixel 
seeing_factor = arcsec_/fwhm_to_sigma;
sigma_1000 = seeing_factor*0.5; %seeing = 05.arcseconds
sigma_750 = seeing_factor*0.5; %seeing = 05.arcseconds
sigma_430 = seeing_factor*0.5;

FWHM_seeing = 1.5;

% wavelength dependent FWHM (Boyd et al. 1978)
sigma_min = seeing_factor*FWHM_seeing; 
aa = 1/5; % power law
bb = sigma_min*fvMDISwl(3)^(aa); % coefficient to realize sigma_min
sigma = bb*fvMDISwl.^(-aa);
fmI_seeing(:,:,1) = (imgaussfilt((fmR(:,:,1)),sigma(1)));
fmI_seeing(:,:,2) = (imgaussfilt((fmR(:,:,2)),sigma(2)));
fmI_seeing(:,:,3) = (imgaussfilt((fmR(:,:,3)),sigma(3)));

%%
% Reflectance
figure; 
himWL = imshow(fmI_seeing(:,:,1),[0 0.1],'XData',fx,'YData',fy);
colormap gray
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','white');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','white','LineWidth',2); 
axis equal
colorbar;
text(-1.4,-1.4,'$1/sr$','Interpreter','Latex','Color','White')
%set(himConstant, 'AlphaData', fmI_seeing1(:,:,2)>0)

% Slope
fmWLSlope = (fmI_seeing(:,:,1)-fmI_seeing(:,:,3))/fWLD;
figure; himWLSlope = imshow(fmWLSlope,[0 0.1],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
%set(himWLSlope, 'AlphaData', fmWLSlope>0)
set(himWLSlope, 'AlphaData', fmI_seeing(:,:,3)>0.01*max(max(fmI_seeing(:,:,3))))
text(-1.4,-1.4,'$1 / (sr\ \mu m)$','Interpreter','Latex')
text(-0.8,0,sprintf('X S2'));
[~,nxSample] = min(abs(fvx- (-0.73)));
[~,nySample] = min(abs(fvy- (0)));
S2 = fmWLSlope(nySample,nxSample);
disp(S2);

% Normalized Slope
figure; 
fmWLNormSlope = (fmI_seeing(:,:,1)-fmI_seeing(:,:,3))./fmI_seeing(:,:,3)/fWLD*100;
himConstantNormSlope = imshow(fmWLNormSlope,[0 350],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
%colorbar;
cH = colorbar;
cH.TickLabelInterpreter = 'tex';
cH.TickLabels{1} = [sprintf('\\color[rgb]{%f,%f,%f} ', [0 0 0]), '<0'];
set(himConstantNormSlope, 'AlphaData', fmI_seeing(:,:,3)>0.01*max(max(fmI_seeing(:,:,3))))
text(-1.35,0,'X F','Color','w');
[~,nxSample] = min(abs(fvx- (-1.30)));
[~,nySample] = min(abs(fvy- (0)));
F = fmWLNormSlope(nySample,nxSample);
disp(F);
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

%% Ratio Normalized Slope wavelength-dependent seeing / Normalized Slope constant Seeing
figure; 
fmRatio = fmWLNormSlope./fmConstNormSlope;
himRatio = imshow(fmRatio,[0 1.20],'XData',fx,'YData',fy);
colormap jet
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on','GColor','black');
getm(gca);
setm(gca,'origin',[0 selo 180]);
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0,'EdgeColor','black','LineWidth',2); 
axis equal
colorbar;
set(himRatio, 'AlphaData', fmI_seeing(:,:,3)>0.01*max(max(fmI_seeing(:,:,3))))
cH = colorbar;
cH.TickLabelInterpreter = 'tex';
cH.TickLabels{1} = [sprintf('\\color[rgb]{%f,%f,%f} ', [0 0 0]), '<0'];
text(-0.4,-0.83,'X G');
[~,nxSample] = min(abs(fvx- (-0.35)));
[~,nySample] = min(abs(fvy- (-0.82)));
G = fmRatio(nySample,nxSample);
disp(G);
text(-0.13,0,'X H');
[~,nxSample] = min(abs(fvx- (-0.06)));
[~,nySample] = min(abs(fvy- (0)));
H = fmRatio(nySample,nxSample);
disp(H);

%% results at points
fvRes = [R0 S0 S1 S2 A B C D E F G H]

