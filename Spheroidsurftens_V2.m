%  <This program was written to measure the surface tension and the Young's modulus of magnetic spheroids>
%     Copyright (C) 2021  Irène Nagle, Myriam Reffay
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


Opt.Interpreter = 'tex';
Opt.WindowStyle = 'normal';    % to be able to change the font characteristics in message boxes with \TeX commands

uiwait(msgbox('\fontsize{12}TensioX Copyright (C) 2021 Irène Nagle and Myriam Reffay.                          This application comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions, for details see the license, GNU General Public License as published by the Free Software Foundation, version 3 or any later version.                    A tutorial can be found at http://www.msc.univ-paris-diderot.fr/spip.php?rubrique274&lang=en','', 'none', Opt));

[image_t0, path_t0] = uigetfile('*.jpg;*.png',' Select Spheroid Image at t=0');
image_t0 = strcat(path_t0,image_t0);
image_t0 = imread(image_t0);
% image_t0 = imread('ag9_t=0.jpg');

[image_tf,path_tf] = uigetfile('*.jpg;*.png','Select Spheroid Image at t=tf');
image_tf = strcat(path_tf,image_tf);
image_tf = imread(image_tf);

selpath = uigetdir('C:\','Select folder where results will be saved') ;

imshow(image_t0)


% prompt = 'Select left of initial spheroid (t=0)'
uiwait(msgbox('\fontsize{14} Select \bfleft\rm of initial spheroid (t=0) then press enter','', 'none', Opt));
[x_l0,z_l0] = ginput ;

% prompt = 'Select right of initial spheroid (t=0)'
uiwait(msgbox('\fontsize{14} Select \bfright\rm of initial spheroid (t=0) then press enter','', 'none', Opt));
[x_r0,z_r0] = ginput ;

% prompt = 'Select apex of initial spheroid (t=0)'
uiwait(msgbox('\fontsize{14} Select \bfapex\rm of initial spheroid (t=0) then press enter','', 'none', Opt));
[x_a0,z_a0] = ginput ;

% prompt = 'Select bottom of initial spheroid (t=0)'
uiwait(msgbox('\fontsize{14} Select \bfbottom\rm of initial spheroid (t=0) then press enter','', 'none', Opt));
[x_b0,z_b0] = ginput;


imshow(image_tf) ;

% prompt = 'Select left of flattened spheroid'
uiwait(msgbox('\fontsize{14} Select \bfleft\rm of flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_lf,z_lf] = ginput ;

% prompt = 'Select right of flattened spheroid'
uiwait(msgbox('\fontsize{14} Select \bfright\rm of flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_rf,z_rf] = ginput ;

% prompt = 'Select apex of flattened spheroid'
uiwait(msgbox('\fontsize{14} Select \bfapex\rm of flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_af,z_af] = ginput ;

% prompt = 'Select bottom of flattened spheroid'
uiwait(msgbox('\fontsize{14} Select \bfbottom\rm of flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_bf,z_bf] = ginput ;

% prompt = 'Select left of contact area'
uiwait(msgbox('\fontsize{14} Select \bfleft\rm of \bfcontact area\rm of the flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_cl,z_cl] = ginput ;

% prompt = 'Select right of contact area'
uiwait(msgbox('\fontsize{14} Select \bfright\rm of \bfcontact area\rm of the flattened spheroid (t=tf) then press enter','', 'none', Opt));
[x_cr,z_cr] = ginput ;

close

prompt = {'\fontsize{12} Enter the image scale in m/pixel','\fontsize{12} Enter the measured value of fv (grad(B)*Mv) in N/m^3','\fontsize{12} Enter an estimated value of gamma in mN/m';};
dlgtitle = 'Input';
dims = [1 75];
definput = {'2.21e-06','50000','10'};
answer= inputdlg(prompt,dlgtitle,dims,definput,Opt);
scale =str2double(answer(1));
fv=str2double(answer(2));
gamma=str2double(answer(3))*10^-3;

global hauteurexp
hauteurexp= (abs(z_af-z_bf))*scale; % in meters

global largeurexp
largeurexp= (abs(x_lf-x_rf))*scale; % in meters


global volumeexp
R0 = (abs(z_a0-z_b0)+abs(x_l0-x_r0))/4*scale;  % averaged radius of the initial spheroid at t=0 in meters
volumeexp = 4/3*pi*(R0)^3 ;      % in meters^3

contactexp=abs(x_cl-x_cr)*scale;    % in meters

%Starting values of b and c for the function fminimiser
b=1/R0;
c=fv/gamma;

X0=[b,c];
Xmin=fminsearch(@fminimiser,X0);

gamma_fitt = fv/Xmin(2)*10^3  ; % gamma in mN/m
E_fitt = (1-1/4)*pi/(contactexp/2)^3*fv*R0^4 ; %in Pa



if size(image_tf, 3) > 1
  image_tf = rgb2gray(image_tf);
end
[Rows,numCols] = size(image_tf);
x0=0;
y0=0;
width=numCols;
height=Rows;

% scale=2.212389381e-06; % echelle ici 2 mm=904 pixels

width_um = numCols * scale;

f = figure('visible','off'); 
set(gcf,'position')
hold on
X= [-flip(x(1:I));x(1:I)];
Z= [flip(z(1:I));z(1:I)];
largeur_fitt = max(X)-min(X);
hauteur_fitt = max(Z)-min(Z);


plot(X,Z,'k','LineWidth',10);
largeur_plot=1.1*largeur_fitt;
hauteur_plot=1.1*hauteur_fitt;
xlim([-(largeur_plot)/2 (largeur_plot)/2]);
ylim([-(hauteur_plot-hauteur_fitt)/2 hauteur_plot-(hauteur_plot-hauteur_fitt)/2]) ;
set(gca,'XTick',[], 'YTick', []);
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                 %[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);

set(gca,'visible','off');
set(gca,'xtick',[]);


Plot = fullfile(selpath,'Plot.jpg');
saveas(gca,Plot);

close 

P1=imread(Plot);
delete(Plot)
P2=P1(:,:,1);
[Rows1,numCols1] = size(P2);
P3 = imresize(P1,[hauteur_plot/scale,largeur_plot/scale],'bicubic');
B_rgb = imtranslate(P3,[(x_lf-(largeur_plot-largeur_fitt)/(2*scale)) (z_af-(hauteur_plot-hauteur_fitt)/(2*scale))],'FillValues',255,'OutputView','full');

if size(B_rgb, 3) > 1
  B = rgb2gray(B_rgb);
end

% transforms gray scale image_tf to RGB
[m n]=size(image_tf);
image_tf_rgb=zeros(m,n,3);
image_tf_rgb(:,:,1)=image_tf;
image_tf_rgb(:,:,2)=image_tf_rgb(:,:,1);
image_tf_rgb(:,:,3)=image_tf_rgb(:,:,1);
image_tf_rgb=image_tf_rgb/255;

for i=1:size(B_rgb,1)
    for j=1:size(B_rgb,2)
        if (B_rgb(i,j,1)>8 || B_rgb(i,j,2)>8 || B_rgb(i,j,3)>8)
            for h= 1:size(B_rgb,3)
                B_rgb(i,j,h)=image_tf_rgb(i,j,h)*255;
            end
        else
            B_rgb(i,j,:)=[255 0 0];
        end
    end
end


figure
imshowpair(image_tf_rgb,B_rgb,'blend','Scaling','joint')
profile_fitt = fullfile(selpath,'profile_fitt.jpg');
saveas(gca,profile_fitt);

results = fullfile(selpath,'results.txt');
fid = fopen(results, 'wt' );
fprintf( fid, 'gamma = %.2f mN/m \n E = %.0f Pa', gamma_fitt,E_fitt);
fclose(fid);

msgbox(sprintf('gamma = %.2f mN/m \n E = %.0f Pa', gamma_fitt,E_fitt))

