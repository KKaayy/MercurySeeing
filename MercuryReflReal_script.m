% Supplementary material: Wohlfarth and Wöhler 2022
% Wavelength-dependent seeing systematically changes the normalized slope of telescopic reflectance spectra of Mercury.
% 13.01.22 K. Wohlfarth
%% 1) Model Setup
fvx = linspace(-1.5,1.5,300);
fvy = linspace(-1.5,1.5,300);
fx = [-1.5 1.5];
fy = [-1.5 1.5];

fWLD = 0.570; % wavelength-difference 1.000 mum - 0.430 mum
nObservation = 12; 
reflmodel = [];
reflmodel.model = 'amsa'; % Hapke 'amsa' or Kaasalainen-Shkuratov 'KS3'

% 05 Warell 1997
% 08 Warell 1999 <--
% 11 Warell 2002 North <--
% 12 Warell 2003 North <--
% 19 Vernazza 2008 <--
% 20 Erard 2006
% 21 Varatharajan 2018

switch nObservation
    case 1  % Warell and Limaye 2001, 20.10.1995
        selo = 273.71;  % selo: sub-earth longitude
        sela = 2.68;    % sela: sub-earth latitude
        sslo = 359.43;  % sslo: sub-solar longitude
        ssla = 0;       % ssla: sub-solar latitude
        al = 85.7;      % al: phase angle
        ad = 6.94;      % ad: angular diameter
        fvCenter = [sela 360 - selo];
    case 2  % Warell and Limaye 2001, 22.10.1995
        selo = 283.74;
        sela = 2.32;
        sslo = 359.64;
        ssla = 0;
        al = 75.9;
        ad = 6.57;
        fvCenter = [sela 360 - selo];
    case 3  % Warell and Limaye 2001, 19.04.1996
        selo = 91.43;
        sela = -2.35;
        sslo =2.01;
        ssla = 0;
        al = -89.4;
        ad = 7.09;
        fvCenter = [sela 360 - selo];
    case 4  % Warell and Limaye 2001, 22.11.1997
        selo = 213.22;
        sela = -1.27;
        sslo = 154.58;
        ssla = 0;
        al = -58.6;
        ad = 5.92;
        fvCenter = [sela 360 - selo];
    case 5  % Warell and Limaye 2001, 24.11.1997
        selo = 222.82;
        sela = -1.53;
        sslo = 159.42;
        ssla = 0;
        al = -63.4;
        ad = 6.12;
        fvCenter = [sela 360 - selo];
    case 6  % Warell and Limaye 2001, 09.07.1998
        selo = 307.68;
        sela = 6.34;
        sslo = 223.69;
        ssla = 0;
        al = -84;
        ad = 6.98;
        fvCenter = [sela 360 - selo];
    case 7  % Warell and Limaye 2001, 27.04.1999
        selo = 58.35;
        sela = -1.41;
        sslo = 143;
        ssla = 0;
        al = 75.6;
        ad = 6.59;
        fvCenter = [sela 360 - selo];
    case 8  % Warell 2003, 20.06.1999
        selo = 291.8;
        sela = 4.8;
        sslo = 207;
        ssla = 0;
        al = -84.8;
        ad = 7.0;
        fvCenter = [sela 360 - selo];
    case 9  % Warell 2003, 22.06.1999
        selo = 301.5;
        sela = 5.2;
        sslo = 212.5;
        ssla = 0;
        al = 89.1;
        ad = 7.2;
        fvCenter = [sela 360 - selo];
    case 10 % Warell and Blewett 2004, 01.07.2002
        selo = 277.1;
        sela = 5.3;
        sslo = 354.7;
        ssla = 0;
        al = 60.7;
        ad = 6.5;
        fvCenter = [sela 360 - selo];
    case 11 % Warell et al. 2006, 26.06.2002
        selo = 241;
        sela = 5.6;
        sslo = 341;
        ssla = 0;
        al = 99.4;
        ad = 7.7;
        fvCenter = [sela 360 - selo];
    case 12 % Warell et al. 2006, 17.08.2002 North
        selo = 197;
        sela = 8.2;
        sslo = 102;
        ssla = 0;
        al = -95.3;
        ad = 7.8;
        fvCenter = [sela 360 - selo];
    case 13 % Warell et al. 2006, 17.08.2002 South
        selo = 197;
        sela = 8.2;
        sslo = 102;
        ssla = 0;
        al = 95.3;
        ad = 7.8;
        fvCenter = [sela 360 - selo];
    case 14 % Vernazza et al. 2010, 28.02.2008
        selo = 141;
        sela = -7;
        sslo = 27;
        ssla = 0;
        al = 89;
        ad = 7.5;
        fvCenter = [sela 360 - selo];
    case 15 % Vernazza et al. 2010, 29.02.2008
        selo = 147;
        sela = -7;
        sslo = 27;
        ssla = 0;
        al = 87;
        ad = 7.4;
        fvCenter = [sela 360 - selo];
    case 16 % Vernazza et al. 2010, 22.03.2008
        selo = 254;
        sela = -4;
        sslo = 21;
        ssla = 0;
        al = 53;
        ad = 5.6;
        fvCenter = [sela 360 - selo];
    case 17 % Vernazza et al. 2010, 23.03.2008
        selo = 258;
        sela = -4;
        sslo = 21;
        ssla = 0;
        al = 51;
        ad = 5.6;
        fvCenter = [sela 360 - selo];
    case 18 % Vernazza et al. 2010, 13.05.2008
        selo = 120;
        sela = 1;
        sslo = 22;
        ssla = 0;
        al = 105;
        ad = 8.0;
        fvCenter = [sela 360 - selo];
    case 19 % Vernazza et al. 2010, 14.05.2008
        selo = 125;
        sela = 1;
        sslo = 22;
        ssla = 0;
        al = -107;
        ad = 8.2;
        fvCenter = [sela 360 - selo];
    case 20 % Erard et al. 2011, 16.06.2006
        selo = 233.24;
        sela = 4.82;
        sslo = 327.7;
        ssla = 0;
        al = -99.44;
        ad = 7.5;
        fvCenter = [sela 360 - selo];
    case 21 % Varatharajan et al. 2018, 16.12.2018
        selo = 301.49;
        sela = -3.92;
        sslo = 11.94;
        ssla = 0;
        al = 70.50;
        ad = 6.4;
        fvCenter = [sela 360 - selo];
end

switch reflmodel.model
    case 'amsa'
        load('cRamsa.mat');
        fmR = cRamsa{nObservation};
    case 'KS3'
        load('cRKS3.mat');
        fmR = cRKS3{nObservation};
end

%% Normalized Slope without seeing
fmNorm00 = (fmR(:,:,1)-fmR(:,:,3))./fmR(:,:,3)/fWLD*100;

figure; 
him00 = imshow(fmNorm00,[0 350],'XData',fx,'YData',fy);
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0); 
axis equal
colormap jet; 
colorbar;
set(him00, 'AlphaData', ~isnan(fmNorm00))

if(sela>0)
text(-0.01,0,strcat('X ',num2str(sela),'°N, ',num2str(360-selo),'°E'));
else
text(-0.01,0,strcat('X ',num2str(sela),'°S, ',num2str(360-selo),'°E'));
end
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')


%% Case 2: Constant seeing
fwhm_to_sigma = 2*sqrt(2*log(2));
arcsec_ = 200/ad*1; % 1 arcsec in pixel 
seeing_factor = arcsec_/fwhm_to_sigma;

FWHM_seeing = 1.5;

sigma_1 = seeing_factor*FWHM_seeing; %seeing = 1.0 arcseconds
fmI_seeingConst = (imgaussfilt((fmR(:,1:300,:)),sigma_1));

%% With constant seeing
fmNormSeeingConst = (fmI_seeingConst(:,:,1)-fmI_seeingConst(:,:,3))./fmI_seeingConst(:,:,3)/fWLD*100;
fmNormSeeingConst(fmNormSeeingConst<0) = NaN;
if(al<0)
    fmNormSeeingConst = fmNormSeeingConst(:,40:end);
    fvx = linspace(-1.5,1.5,300);
    fvy = linspace(-1.5,1.5,300);
end

figure;
him05 = imshow(fmNormSeeingConst,[0 350],'XData',fx,'YData',fy);
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0);
axis equal
colormap jet;
colorbar;
set(him05, 'AlphaData', ~isnan(fmNormSeeingConst))
text(-0.01,0,'X');
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

%% Case 3: Wavelength-dependent seeing

fvMDISwl = [1000 750 430]*1e-6;
fwhm_to_sigma = 2*sqrt(2*log(2));
arcsec_ = 200/ad*1; % 1 arcsec in pixel 
seeing_factor = arcsec_/fwhm_to_sigma;

FWHM_seeing = 1.5;

% wavelength dependent FWHM (Boyd et al. 1978)
sigma_min = seeing_factor*FWHM_seeing; 
aa = 1/5; % power law
bb = sigma_min*fvMDISwl(3)^(aa); % coefficient to realize sigma_min
sigma = bb*fvMDISwl.^(-aa);
fmI_seeingWL(:,:,1) = (imgaussfilt((fmR(:,:,1)),sigma(1)));
fmI_seeingWL(:,:,2) = (imgaussfilt((fmR(:,:,2)),sigma(2)));
fmI_seeingWL(:,:,3) = (imgaussfilt((fmR(:,:,3)),sigma(3)));

%% With seeing

fmNormSeeingWL = (fmI_seeingWL(:,:,1)-fmI_seeingWL(:,:,3))./fmI_seeingWL(:,:,3)/fWLD*100;
fmNormSeeingWL(fmNormSeeingWL<0) = NaN;
if(al<0)
   % fmNormSeeingWL = fmNormSeeingWL(:,40:end);
    fvx = linspace(-1.5,1.5,300);
    fvy = linspace(-1.5,1.5,300);
end

figure;
him05 = imshow(fmNormSeeingWL,[0 350],'XData',fx,'YData',fy);
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0);
axis equal
colormap jet;
colorbar;
set(him05, 'AlphaData', ~isnan(fmNormSeeingWL))

text(-0.01,0,'X');
switch nObservation
    case 11% warell 2002 north
        fvx = linspace(-1.5,1.5,300);
        nl = 122;
        nr = 127;
        nu = 55;
        nd = 85;
        rectangle('Position',[fvx(nl) fvx(nu) abs(fvx(nr)-fvx(nl)) abs(fvx(nu)-fvx(nd))]);
        fMean11 = mean(fmNormSeeingWL(nu:nd,nl:nr),'all')
        nl = 55;
        nr = 108;
        nu = 75;
        nd = 112;
        mean(fmNormSeeingWL(nu:nd,nl:nr),'all')
    case 12% warell 2003 north
        fvx = linspace(-1.5,1.5,300);
        nl = 30;
        nr = 280;
        nu = 80;
        nd = 86;
        rectangle('Position',[fvx(nl) fvx(nu) abs(fvx(nr)-fvx(nl)) abs(fvx(nu)-fvx(nd))]);
        fMean11 = mean(fmNormSeeingWL(nu:nd,nl:nr),'all','omitnan')
    case 19% Vernazza 2008
        line([fvx(200),fvx(280)],[0 0],'Color','black')
        mean(fmNormSeeingWL(149:151,200:280),'all','omitnan')
        line([fvx(210),fvx(210)],[fvx(150),fvx(280)],'Color','black','LineStyle','--')
        mean(fmNormSeeingWL(150:280,209:211),'all','omitnan')

end
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex');

%% Vernazza et al. 2010 translation
if(nObservation == 19)
sz = 300;  % size of image [pixels]
r0 = 100;   % disk radius [pixels]
u0 = round(sz/2);  % centre coordinates of disk [pixels]
v0 = round(sz/2);

% determine normal vector for each pixel of the disk
nv=nan(sz,sz,3);
for l=1:sz
    dv=l-v0;
    for k=1:sz
        du=k-u0;
        if(sqrt(du*du+dv*dv)<r0)
            un=du/r0;
            vn=dv/r0;
            wn=sqrt(1.0-un*un-vn*vn);
            nv(l,k,1)=un;
            nv(l,k,2)=vn;
            nv(l,k,3)=wn;
        end
    end
end

nTranslation = 15; % pixel of translation: disk-dianmeter is 2*r0 = 200 px

% compute illumination vector
vi=[-sind(al) 0.0 cosd(al)];

% compute cosine of emission angle for each pixel
vo = [0 0 1];
mu=nan(sz,sz);
for l=1:sz
    for k=1:sz
        if(~isnan(nv(l,k,1)))
            n=squeeze(nv(l,k,:));
            mu(l,k)=dot(n,vo);
        end
    end
end

fmEms = acosd(mu(:,:,1));

fmNormSeeingWL = (fmI_seeingWL(:,1:300,1)-fmI_seeingWL(:,1:300,3))./fmI_seeingWL(:,1:300,3)/fWLD*100;
fmNormSeeingWL_translated = [fmNormSeeingWL(:,nTranslation+1:end) nan(sz,nTranslation)];

fmNormSeeingWL_translated(fmNormSeeingWL_translated<0) = NaN;
figure; 
him15 = imshow(fmNormSeeingWL_translated,[0 350],'XData',fx,'YData',fy);
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180])
setm(gca,'scalefactor',r0/100)
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', r0/100);
plot(p,'FaceAlpha',0); 
axis equal
colormap jet; 
colorbar;
set(him15, 'AlphaData', ~isnan(fmNormSeeingWL_translated))
text(-0.01,0,'X');
line([fvx(200),fvx(280)],[0 0],'Color','black')
line([fvx(210),fvx(210)],[fvx(150),fvx(280)],'Color','black','LineStyle','--')
mean(fmNormSeeingWL(149:151,200:280),'all','omitnan')
mean(fmNormSeeingWL(150:280,209:211),'all','omitnan')
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

figure; 
hold on;
scatter(fmEms(150,200:280),fmNormSeeingWL(150,200:280),100,'o','k')
scatter(fmEms(150:280,210),fmNormSeeingWL(150:280,210),100,'x','k')
scatter(fmEms(150,200:280),fmNormSeeingWL_translated(150,200:280),100,'o','r')
scatter(fmEms(150:280,210),fmNormSeeingWL_translated(150:280,210),100,'x','r')
leg = legend('Normal disk solid line','Normal disk dashed line','Translated disk solid line','Translated disk dashed line',...
'Interpreter','Latex','Location','SouthWest');
leg.FontSize = 12;
xlabel('Emission Angle $(^\circ)$','Interpreter','Latex')
ylabel('Normalized Slope $[\%/\mu m]$','Interpreter','Latex')
grid on;
ylim([80 350])
end

%% Warell et al. 2002 rotate
if(nObservation == 8)
% Const
fx = [-1.5 1.5];
fy = [-1.5 1.5];
fmR_r = imrotate(fmR,-60,'crop');
fmR_rot = imrotate((fmR(:,:,1)-fmR(:,:,3))./fmR(:,:,3)/fWLD*100,-60,'crop');
fmR_rot(fmR_rot==0) = NaN;
figure; himRot = imshow(fmR_rot,[0 350],'Xdata',fx,'YData',fy); 
colormap jet; colorbar;
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180-60])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0); 
axis equal
colormap jet; 
set(himRot, 'AlphaData', ~isnan(fmR_rot))
set(himRot, 'AlphaData', fmR_r(:,:,3)>0.01*max(max(fmR_r(:,:,3))))
text(-0.01,0,'X');

if(sela>0)
text(-0.01,0,strcat('X ',num2str(sela),'°N, ',num2str(360-selo),'°E'));
else
text(-0.01,0,strcat('X ',num2str(sela),'°S, ',num2str(360-selo),'°E'));
end
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

%%
fmI_seeing_rot_const = imrotate(fmI_seeingConst,-60,'crop');
fmRotConst = imrotate((fmI_seeingConst(:,:,1)-fmI_seeingConst(:,:,3))./fmI_seeingConst(:,:,3)/fWLD*100,-60,'crop');
fmRotConst(fmRotConst==0) = NaN;
figure; himRot = imshow(fmRotConst,[0 350],'Xdata',fx,'YData',fy); 
colormap jet; 
colorbar;
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180-60])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0); 
axis equal
set(himRot, 'AlphaData', ~isnan(fmRotConst))
%set(himRot, 'AlphaData', fmR_r(:,:,3)>0.01*max(max(fmR_r(:,:,3))))

text(-0.01,0,'X');

nNLine = 8;
fvLine = linspace(-0.9,1.0,nNLine);
for k=1:nNLine
xline(fvLine(k))
end

fvLine = linspace(-1,1.1,nNLine);
for k=1:nNLine
xline(fvLine(k),'--')
end
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')


%%
fmI_seeing_rot_WL = imrotate(fmI_seeingWL,-60,'crop');
fmRotWL = imrotate((fmI_seeingWL(:,:,1)-fmI_seeingWL(:,:,3))./fmI_seeingWL(:,:,3)/fWLD*100,-60,'crop');
fmRotWL(fmRotWL==0) = NaN;
figure; 
himRot = imshow(fmRotWL,[0 350],'Xdata',fx,'YData',fy); 
hold on;
axesm ('ortho', 'Frame', 'on', 'Grid', 'on');
getm(gca)
setm(gca,'origin',[0 selo 180-60])
p = nsidedpoly(1000, 'Center', [0 0], 'Radius', 1);
plot(p,'FaceAlpha',0); 
axis equal
set(himRot, 'AlphaData', ~isnan(fmRotWL))

text(-0.01,0,'X');

nNLine = 8;
fvLine = linspace(-0.9,1.0,nNLine);
for k=1:nNLine
xline(fvLine(k))
end

fvLine = linspace(-1,1.1,nNLine);
for k=1:nNLine
xline(fvLine(k),'--')
end
cH2 = colorbar;
cH2.TickLabelInterpreter = 'tex';
cH2.TickLabels{1} = [sprintf('\\color[rgb]{%f,%f,%f} ', [0 0 0]), '<=0'];
colormap jet;
text(-1.4,-1.4,'$\% / \mu m$','Interpreter','Latex')

%%
fvx = linspace(-1.5,1.5,300);
fvy = linspace(-1.5,1.5,300);
fmI_seeing2 = fmI_seeingWL;
fmI_seeing2(fmI_seeing2<0)=0;
fvX1 = round(linspace(61,250,8))
fvX2 = round(linspace(51,260,8))

figure; 
hold on;
fmSlitAverage = zeros(4,8,3);
for k=1:8
    fvX1(k)
    %scatter(k,mean(fmI_seeing_rot(:,fvX1(k),2),'omitnan'),'.','k')
    %scatter(k,mean(fmI_seeing_rot(:,fvX2(k),2),'omitnan'),'.','r')
    fmSlitAverage(1,k,:) = mean(fmI_seeing_rot_const(:,fvX1(k),:),'omitnan');
    fmSlitAverage(2,k,:) = mean(fmI_seeing_rot_const(:,fvX2(k),:),'omitnan'); 
    fmSlitAverage(3,k,:) = mean(fmI_seeing_rot_WL(:,fvX1(k),:),'omitnan');
    fmSlitAverage(4,k,:) = mean(fmI_seeing_rot_WL(:,fvX2(k),:),'omitnan'); 
end

scatter(1:8,(fmSlitAverage(1,:,1)-fmSlitAverage(1,:,3))./fmSlitAverage(1,:,3)/fWLD*100,100,'k','o');
scatter(1:8,(fmSlitAverage(2,:,1)-fmSlitAverage(2,:,3))./fmSlitAverage(2,:,3)/fWLD*100,100,'k','x');

scatter(1:8,(fmSlitAverage(3,:,1)-fmSlitAverage(3,:,3))./fmSlitAverage(3,:,3)/fWLD*100,100,'r','o');
scatter(1:8,(fmSlitAverage(4,:,1)-fmSlitAverage(4,:,3))./fmSlitAverage(4,:,3)/fWLD*100,100,'r','x');

leg = legend('Constant Seeing solid slits','Constant seeing dashed slits','WL-dependent seeing solid slits','WL-dependent seeing dashed slits',...
'Interpreter','Latex','Location','SouthEast');
leg.FontSize = 12;
xlabel('Slit Position','Interpreter','Latex')
ylabel('Normalized Slope','Interpreter','Latex')
ylim([60 350])
grid on;
end

