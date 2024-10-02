%   Creating mesh using iso2mesh for microstructure with three phase
%   Phase : Pore, Ionomer, Catalyst, Membrane, Current collector 
%   Updated : May 7, 2024
%   Creating mesh file from iso2mesh 
function xsave_mphtxt_file_wo_pore
% Main function
clc
clearvars
close all
ximage = xreadimagefile;
[node,elem,face] = xiso2mesh(ximage);
[face2,faceID] = xfaceID_tags(node,face);
xwriting_mphtxt_file(node,elem,face2,faceID)%,TPB_edge)
end
function [ximage] = xreadimagefile
% Reading Tiff image file
% Convert the tiff file from RGB color to 8 bit color in imageJ 
% >>imageJ -> Image -> Type -> 8bit color<<
[xfile,xpath] = uigetfile('*.tif');
ximage = tiffreadVolume(strcat(xpath,xfile));
% Converting to uint8 class variables
% adding additional layer of membrane and current collector
size_x = max(size(ximage(:,1,1)));
size_y = max(size(ximage(1,:,1)));
size_z = max(size(ximage(1,1,:)));

ximage_1 = uint8(zeros(size_x,size_y,size_z+220));
ximage_1(:,:,21:size_z+20) = uint8(ximage);
ximage_1(:,:,1:20) = uint8(4); % Current collector
ximage_1(:,:,size_z+21:end) = uint8(3); % Membrane

clear ximage
ximage = ximage_1;
figure
histogram(ximage)
% iso2mesh ignores any phase with value 0
% For simplicity, phases are arranged in order - Pore, Ionomer, Catalyst
for i = 1:size(ximage,1)            
    for j = 1:size(ximage,2)
        for k = 1:size(ximage,3)
            if(ximage(i,j,k) == 0)               % Pore -- Default (-)
                ximage(i,j,k) = 0;
            elseif(ximage(i,j,k) == 171)         % Ionomer -- Verify the pixel value (I)
                ximage(i,j,k) = 100;
            elseif(ximage(i,j,k) == 85)          % Catalyst -- Verify the pixel value (II)
                ximage(i,j,k) = 150;
            elseif(ximage(i,j,k) == 4)           % Current collector (IV)
                ximage(i,j,k) = 250;           
            elseif(ximage(i,j,k) == 3)           % Membrane (III)
                ximage(i,j,k) = 200;           
            end
        end
    end
end
end
function [node,elem,face] = xiso2mesh(ximage)
% by default, vol2mesh uses 'cgalsurf' method, which requires the following
opt(1).radbound=5;    % Pore        - element size bound
opt(2).radbound=1;    % Ionomer     - element size bound
opt(3).radbound=1;    % Catalyst    - element size bound
opt(4).radbound=10;   % Membrane
opt(5).radbound=10;   % Current collector
opt(1).side='upper';  %
opt(2).side='upper';  %
opt(3).side='upper';  %
opt(4).side='upper';  %
opt(5).side='upper';  %
%##################################################################
%#####################@@@@ VOL2MESH @@@@###########################
%##################################################################
%########### Add path to where iso2mesh file is placed ############
%##################################################################
tstart = tic;
%addpath('/Users/pjayapra/Documents/MATLAB/iso2mesh');
[node,elem,face]=v2m(ximage,[],opt,5,'cgalmesh');
% [node,elem,face] = vol2mesh(ximage,[],opt,10,'cgalmesh');
% Scaling each vertex as per voxel size on scanned image
%xinp = inputdlg({'Voxel value in micrometers'},'Voxel value', [1 50]);
% node = node.*(str2double(cell2mat(xinp)));
fprintf('Time : %f Minute\n',toc(tstart)/60);
xinp = 1e-8;    % Voxel size in um
node = node.*xinp;
end
function [face2,faceID] = xfaceID_tags(node,face)
% Face group
% Face 12 : Ionomer/Catalyst
% Face 13 : Ionomer/Membrane
% Face 14 : Ionomer/Current collector
% Face 23 : Catalyst/Membrane
% Face 24 : Catalyst/Current collector
global in12 in13 in14 in23 in24
global inter12 inter13 inter14 inter23 inter24 
global xf_mm xf_CC
n_faces = size(face,1);
% Identifying the interface mesh elements
% Update the face mesh information from Elements
face_temp = zeros(size(face(:,1:3)));                % Avoiding repeatation of faces
face_temp(1,:) = face(1,1:3);

% Check for repeation 
for nf1 = 2:n_faces
    % for same face ID
    if(mean(face(nf1,1:3) == face(nf1-1,1:3)) == 1)  % Check for immedate face Id
        face_temp(nf1,:) = 0;
        if(face(nf1,4) == face(nf1-1,4))
            face_temp(nf1,:) = 0;
        end
    else
        face_temp(nf1,:) = face(nf1,1:3);
    end
end

xindex = find(sum(face_temp,2) ~= 0);
face2 = face_temp(xindex,:);
faceID = zeros(size(face2(:,1)));
n_faces_2 = numel(face2(:,1));

% Tag for 2PB face elements
in12 = 1;       % Ionomer/Catalyst
in13 = 1;       % Ionomer/Membrane
in14 = 1;       % Ionomer/Current collector
in23 = 1;       % Catalyst/Membrane    
in24 = 1;       % Catalyst/Current collector
in11 = 1;
in22 = 1;
in33 = 1;
in44 = 1;

clear inter12 inter13 inter14 inter23 inter24 inter11 inter22 inter33 inter44
global inter12 inter13 inter14 inter23 inter24 inter11 inter22 inter33 inter44

% for nf1 = 1:n_faces_2
%     nf2  = nf1*2;
%     if (face(nf2,4) == face(nf2-1,4))
%         if(face(nf2,4) == 1)                % Pore
%             faceID(nf1) = n_faces_2*2 + nf1;
%         elseif(face(nf2,4)==2)              % Ionomer
%             faceID(nf1) = n_faces_2*3 + nf1;
%         elseif(face(nf2,4)==3)              % Catalyst
%             faceID(nf1) = n_faces_2*4 + nf1;
%         elseif(face(nf2,4)==4)              % Membrane
%             faceID(nf1) = n_faces_2*4 + nf1;
%         elseif(face(nf2,4)==5)              % Current collector
%             faceID(nf1) = n_faces_2*4 + nf1;
%         end
% %         faceID(nf1) = face(nf1,4); 
%     else
%         faceID(nf1) = nf2;
%         xt = 10*face(nf2,4) + face(nf2-1,4);
%         if(xt == 12 || xt == 21)           % Pore/Ionomer
%             inter12(in12,1) = nf2;
%             in12 = in12 + 1;    
%         elseif(xt == 23 || xt == 32)       % Ionomer/Catalyst
%             inter23(in23,1) = nf2;
%             in23 = in23 + 1;
%         elseif(xt == 13 || xt == 31)       % Pore/Catalyst
%             inter13(in13,1) = nf2;
%             in13 = in13 + 1;
%         elseif(xt == 14 || xt == 41)       % Pore/Membrane
%             inter14(in14,1) = nf2;
%             in14 = in14 + 1;
%         elseif(xt == 15 || xt == 51)       % Pore/CC
%             inter15(in15,1) = nf2;
%             in15 = in15 + 1;
%         elseif(xt == 24 || xt == 42)       % Catalyst/Membrane
%             inter24(in24,1) = nf2;
%             in24 = in24 + 1;
%         elseif(xt == 34 || xt == 43)       % Ionomer/Membrane
%             inter34(in34,1) = nf2;
%             in34 = in34 + 1;
%         elseif(xt == 25 || xt == 52)       % Catalyst/Current collector
%             inter25(in25,1) = nf2;
%             in25 = in25 + 1;
%         elseif(xt == 35 || xt == 53)       % Ionomer/Current collector
%             inter35(in35,1) = nf2;
%             in35 = in35 + 1;
%         end
%     end
% end

nf3 = 1;

for nf1 = 1:2:n_faces
    nf2  = nf1 + 1;
    if (face(nf1,4) == face(nf2,4))
        if(face(nf1,4) == 1)                % Ionomer
            faceID(nf3) = n_faces_2*2 + nf3;
            inter11(in11,1) = faceID(nf3);
            in11 = in11 + 1;
        elseif(face(nf1,4)==2)              % Catalyst
            faceID(nf3) = n_faces_2*3 + nf3;
            inter22(in22,1) = faceID(nf3);
            in22 = in22 + 1;
        elseif(face(nf1,4)==3)              % Membrane
            faceID(nf3) = n_faces_2*4 + nf3;
            inter33(in33,1) = faceID(nf3);
            in33 = in33 + 1;
        elseif(face(nf1,4)==4)              % Current collector
            faceID(nf3) = n_faces_2*5 + nf3;
            inter44(in44,1) = faceID(nf3);
            in44 = in44 + 1;
        end
    else
        faceID(nf3) = nf3;
        xt = 10*face(nf1,4) + face(nf2,4);
        if(xt == 12 || xt == 21)           % Ionomer/Catalyst
            inter12(in12,1) = faceID(nf3); 
            in12 = in12 + 1;    
        elseif(xt == 23 || xt == 32)       % Catalyst/Membrane
            inter23(in23,1) = faceID(nf3); 
            in23 = in23 + 1;
        elseif(xt == 13 || xt == 31)       % Ionomer/Membrane
            inter13(in13,1) = faceID(nf3); 
            in13 = in13 + 1;
        elseif(xt == 14 || xt == 41)       % Ionomer/Current collector
            inter14(in14,1) = faceID(nf3); 
            in14 = in14 + 1;
        elseif(xt == 24 || xt == 42)       % Catalyst/Current collector
            inter24(in24,1) = faceID(nf3); 
            in24 = in24 + 1;
        end
    end
    nf3 = nf3 + 1;
end

% Membrane interface end
z_max = max(node(:,3));
z1 = z_max - 0.005*z_max;
z2 = z_max;
xtemp1 = [];
for i = 1:numel(node(:,3))
    if(node(i,3) >= z1)
        xtemp1 = [xtemp1; i];
    end
end
% search for the faces and IDs
% xf_mm = max(faceID) + 1;
xftags = [];
for i = 1:numel(xtemp1)
    xin = find(face2(:,1)==xtemp1(i));
    % faceID(xin) = xf_mm;
    xftags = [xftags; faceID(xin)];
    xin = find(face2(:,2)==xtemp1(i));
    % faceID(xin) = xf_mm;
    xftags = [xftags; faceID(xin)];
    xin = find(face2(:,3)==xtemp1(i));
    % faceID(xin) = xf_mm;
    xftags = [xftags; faceID(xin)];
end
xftags = sort(xftags);
xf_mm = unique(xftags);

% Current collector end
z_min = min(node(:,3));
z1 = z_min + 1.2*z_min;
z2 = z_min;
xtemp2 = [];
for i = 1:numel(node(:,3))
    if(node(i,3) <= z1)
        xtemp2 = [xtemp2; i];
    end
end
% search for the faces and IDs
% xf_CC = max(faceID) + 1;
xftags = [];
for i = 1:numel(xtemp2)
    xin = find(face2(:,1)==xtemp2(i));
    % faceID(xin) = xf_CC;
    xftags = [xftags; faceID(xin)];
    xin = find(face2(:,2)==xtemp2(i));
    % faceID(xin) = xf_CC;
    xftags = [xftags; faceID(xin)];
    xin = find(face2(:,3)==xtemp2(i));
    % faceID(xin) = xf_CC;
    xftags = [xftags; faceID(xin)];
end
xftags = sort(xftags);
xf_CC = unique(xftags);

% Plotting the mesh file
figure
plotmesh(node(:,[2 1 3]),face2,'facealpha',0.7);
end
function xwriting_mphtxt_file(node,elem,face2,faceID)%,TPB_edge)
global inter12 inter13 inter14 
global inter23 inter24 xf_mm xf_CC
n_faces_2 = size(face2,1);
n_elements = size(elem,1);
n_nodes = size(node,1);
% Exporting the mesh infromation
% Writing the mesh file for COMSOL
% with selection of boundaries at interface
xfilename = strcat('PEM-without-Pore-',date,'.mphtxt');
fid = fopen(xfilename,'wt');
% Header
fprintf(fid,'# COMSOL mphtxt file creation - Synthetic_microstructure\n');
fprintf(fid,'# Major and minor version\n');
fprintf(fid,'0  1\n');
fprintf(fid,'12                      #number of tags\n');       
fprintf(fid,'5 mesh1                #tag 1\n'); % Mesh
fprintf(fid,'12 meshselcdom1        #tag 2\n'); % Ionomer - domain
fprintf(fid,'12 meshselcdom2        #tag 3\n'); % Catalyst - domain
fprintf(fid,'12 meshselcdom3        #tag 4\n'); % Membrane - domain
fprintf(fid,'12 meshselcdom4        #tag 5\n'); % Current collector - domain
fprintf(fid,'12 meshselcfac1        #tag 6\n'); 
fprintf(fid,'12 meshselcfac2        #tag 7\n'); 
fprintf(fid,'12 meshselcfac3        #tag 8\n'); 
fprintf(fid,'12 meshselcfac4        #tag 9\n');
fprintf(fid,'12 meshselcfac5        #tag 10\n'); 
fprintf(fid,'12 meshselcfac6        #tag 11\n');
fprintf(fid,'12 meshselcfac7        #tag 12\n');
%fprintf(fid,'9 edgeselc1            #tag 8\n');% 3PB          
fprintf(fid,'12                     #number of objects\n'); % Equal number of tags
fprintf(fid,'4 obj1                 #object 1\n');
fprintf(fid,'4 obj2                 #object 2\n');
fprintf(fid,'4 obj3                 #object 3\n');
fprintf(fid,'4 obj4                 #object 4\n');
fprintf(fid,'4 obj5                 #object 5\n');
fprintf(fid,'4 obj6                 #object 6\n');
fprintf(fid,'4 obj7                 #object 7\n');
fprintf(fid,'4 obj8                 #object 8\n');
fprintf(fid,'4 obj9                 #object 9\n');
fprintf(fid,'5 obj10                #object 10\n');
fprintf(fid,'5 obj11                #object 11\n');
fprintf(fid,'5 obj11                #object 12\n');
% Inclusion of information for object 0 ------------- Mesh
fprintf(fid,'0 0 1                  #Default condition\n');
fprintf(fid,'4 Mesh                 #Class type\n');
fprintf(fid,'4                      #version\n');
fprintf(fid,'3                      #dimension\n');
fprintf(fid,'%d                     #number of mesh vertices\n',n_nodes);
fprintf(fid,'1                      #Lowest index\n\n');
fprintf(fid,'# Mesh vertex coordinates\n');
% Mesh vertex coordinates
for ix = 1:n_nodes
    fprintf(fid,'%.16f %.16f %.16f               #%d\n',node(ix,1),node(ix,2),...
                                                  node(ix,3),ix);
end
fprintf(fid,'\n');
fprintf(fid,'2                      #number of mesh elemtent types\n');     % Change 3 from 2 if egde elements are included
%--------------------------------------------------------------------------
% Mesh tetrahedral  element
fprintf(fid,'# Type 0 --- Tetrahedral ---\n');
fprintf(fid,'3 tet                 #tetra\n\n\n');
fprintf(fid,'4                     #number of vertices per elements\n');
fprintf(fid,'%d                    #number of elements\n',n_elements);
fprintf(fid,'# Elements\n');
for it = 1 : n_elements
    fprintf(fid,'%d %d %d %d    #%d\n',elem(it,1),elem(it,2),...
                                                 elem(it,3),elem(it,4),it);
end
fprintf(fid,'\n\n');
% Tetrahedral element indices
fprintf(fid,'\n');
fprintf(fid,'%d                             #number enetity indices\n',n_elements);
% fprintf(fid,'%d\n',elem(:,5)');
for it = 1:n_elements
    fprintf(fid,'%d              #%d\n',elem(it,5),it);
end
fprintf(fid,'\n\n');
%--------------------------------------------------------------------------
% Mesh triangular element
fprintf(fid,'# Type 1 --- Triangular --- \n');
fprintf(fid,'3 tri                  #triangular\n');
fprintf(fid,'3                      #number of vertices per element\n');
fprintf(fid,'%d                     #number of elements\n',n_faces_2);
fprintf(fid,'# Elements\n');
for iq = 1:n_faces_2
    fprintf(fid,'%d %d %d       #%d\n',face2(iq,1:3),iq);
end
fprintf(fid,'\n\n');
% Triangular element indices
fprintf(fid,'\n');
fprintf(fid,'%d                             #number enetity indices\n',n_faces_2);
for iq = 1:n_faces_2
    fprintf(fid,'%d              #%d\n',faceID(iq),iq);
end
fprintf(fid,'\n\n');
%--------------------------------------------------------------------------
fprintf(fid,'# ---------------- Object 1 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'7 Ionomer          #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'3                  #Dimension\n');
fprintf(fid,'1                  #number of entities\n');
fprintf(fid,'#Entity tags\n');
fprintf(fid,'1\n');
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 2 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'8 Catalyst         #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'3                  #Dimension\n');
fprintf(fid,'1                  #number of entities\n');
fprintf(fid,'#Entity tags\n');
fprintf(fid,'2\n');
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 3 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'8 Membrane          #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'3                  #Dimension\n');
fprintf(fid,'1                  #number of entities\n');
fprintf(fid,'#Entity tags\n');
fprintf(fid,'3\n');
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 4 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'17 Current_collector   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'3                  #Dimension\n');
fprintf(fid,'1                  #number of entities\n');
fprintf(fid,'#Entity tags\n');
fprintf(fid,'4\n');
fprintf(fid,'\n\n');
%---------------------------- Edges - Interfaces --------------------------
fprintf(fid,'# ---------------- Object 5 -----------------\n');
fprintf(fid,'0 0 1               #Default version\n');
fprintf(fid,'9 Selection         #Class\n');
fprintf(fid,'0                   #Version\n');
fprintf(fid,'16 Ionomer/Catalyst #Label\n');
fprintf(fid,'5 mesh1             #Mesh tag\n');
fprintf(fid,'2                   #Dimension\n');
fprintf(fid,'%d                  #number of entities\n',numel(inter12));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(inter12)
    fprintf(fid,'%d\n',inter12(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 6 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'16 Ionomer/Membrane #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(inter13));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(inter13)
    fprintf(fid,'%d\n',inter13(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 7 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'25 Ionomer/Current_collector   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(inter14));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(inter14)
    fprintf(fid,'%d\n',inter14(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 8 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'17 Catalyst/Membrane   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(inter23));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(inter23)
    fprintf(fid,'%d\n',inter23(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 9 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'26 Catalyst/Current_collector   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(inter24));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(inter24)
    fprintf(fid,'%d\n',inter24(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 10 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'18 Membrane_interface   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(xf_mm));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(xf_mm)
    fprintf(fid,'%d\n',xf_mm(iq));
end
fprintf(fid,'\n\n');
fprintf(fid,'# ---------------- Object 11 -----------------\n');
fprintf(fid,'0 0 1              #Default version\n');
fprintf(fid,'9 Selection        #Class\n');
fprintf(fid,'0                  #Version\n');
fprintf(fid,'17 Current_collector   #Label\n');
fprintf(fid,'5 mesh1            #Mesh tag\n');
fprintf(fid,'2                  #Dimension\n');
fprintf(fid,'%d                 #number of entities\n',numel(xf_CC));
fprintf(fid,'#Entity tags\n');
for iq = 1:numel(xf_CC)
    fprintf(fid,'%d\n',xf_CC(iq));
end
fprintf(fid,'\n\n');
fclose(fid);
fprintf('Finished writing mphtxt file\n');
disp(datetime('now','Format','HH:mm:ss'));
end