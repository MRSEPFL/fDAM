clear B1_fDAM
% load images data
Image1=load('ImageData1'); %size(slice,inplane1,inplane2)
Image2=load('ImageData2');

% only necessary if inverted values for GRE2 is expected CB1*alpha2>180°
% if angle between GRE1 and GRE2 is larger than 90° it's detected as inverted 
Inversion_detection= abs(angle(Image1)-angle(Image2));
Inversion_detection(Inversion_detection<=pi/2)=1;
Inversion_detection(Inversion_detection>pi/2)=-1;

% get magnitude images with inversion information
Image1=abs(Image1);
Image2=abs(Image2).*Inversion_detection;

nb_Slices=size(Image1,1);

% set T1 assumption
T1 =3.4;      % T1 in seconds

%Protocol parameter
alpha_1_nom = 59; %Nominal Flipangel set for GRE1 in degree
TR =5.4;    %repetition time in seconds; should be approx 1.6*T1 

% create look up table
res_lupt=1; %resolution of the look-up table in degree 
rfB1=1:res_lupt:alpha_1_nom*2; % Range of look-up table;
lib=  (  sind(2*rfB1).*(1-exp(-TR/T1))./(1-exp(-TR/T1)*cosd(2*rfB1)))...
    ./(2*sind(  rfB1).*(1-exp(-TR/T1))./(1-exp(-TR/T1)*cosd(  rfB1)));
%Loop over slices
for i=1:nb_Slices

% Signal ratio
S_ratio_DAM=abs(squeeze(Image2(i,:,:))./(2*squeeze(Image1(i,:,:))));
S_ratio_DAM(S_ratio_DAM>1.0)=nan; % exclude values which are larger then 1 due to noise

% Look up table fitting
[M,I]=min(abs(S_ratio_DAM(:)-lib),[],2);  
B1_fDAM(i,:,:)=reshape(rfB1(I),size(S_ratio_DAM)); 
end
B1_fDAM = B1_fDAM/alpha_1_nom; %convert actual FA map to nominal scaling factor


%% plot B1 maps
figure(11)
supblot_dim=ceil(sqrt(nb_Slices));

for k=1:nb_Slices

subplot(supblot_dim,supblot_dim,k)
imagesc(squeeze(B1_fDAM(k,:,:)));
ax=gca;
Colormap=colormap('hot');
ax.CLim=[0.5 1.5];

end

%%
save('B1_map.mat','B1_fDAM','TR','T1','alpha_1_nom')
