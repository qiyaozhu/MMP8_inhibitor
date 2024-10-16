function Combine_extension

addpath("/mnt/home/qzhu1/my_functions");
addpath("/mnt/home/qzhu1/RamaMap");

rama_phe = load('ramamap_phe.mat').Map;
rama_dphe = load('ramamap_dphe.mat').Map;

n = 18;

% Grid for preferred backbone regions
MMP8 = load("MMP8_grid.mat");
progrid = MMP8.progrid;
rewardgrid = MMP8.rewardgrid;
minX = MMP8.minX;
minY = MMP8.minY;
minZ = MMP8.minZ;
grid_dx = MMP8.grid_dx;
nX = MMP8.nX;
nY = MMP8.nY;
nZ = MMP8.nZ;
min_s = MMP8.min_s;
target_grid = MMP8.target_grid;

MMP8 = pdb2mat("MMP8.pdb");
proX = MMP8.X;
proY = MMP8.Y;
proZ = MMP8.Z;
proatom_names = MMP8.atomName;
proatom_number = MMP8.atomNum;
prores_names = MMP8.resName;
prores_number = MMP8.resNum;
proelements = MMP8.element;

solname = "DASP_DNLU_backbone_sol2";
stub = pdb2mat(solname+"_0001.pdb");
stubX = stub.X;
stubY = stub.Y;
stubZ = stub.Z;
stubatom_names = stub.atomName;
stubatom_number = stub.atomNum;
stubres_names = stub.resName;
stubres_number = stub.resNum;
stubelements = stub.element;

fake_H_ind = [];
for i = 1 : length(stubX)
    if stubatom_names(i) == "N"
        if stubres_number(i) == 1
            N1 = [stubX(i); stubY(i); stubZ(i)];
        else
            N2 = [stubX(i); stubY(i); stubZ(i)];
        end
    elseif stubatom_names(i) == "CA"
        if stubres_number(i) == 1
            CA1 = [stubX(i); stubY(i); stubZ(i)];
        else
            CA2 = [stubX(i); stubY(i); stubZ(i)];
        end
    elseif stubatom_names(i) == "C"
        if stubres_number(i) == 1
            C1 = [stubX(i); stubY(i); stubZ(i)];
        else
            C2 = [stubX(i); stubY(i); stubZ(i)];
        end
    elseif stubatom_names(i) == "O" && stubres_number(i) == 1
        O1 = [stubX(i); stubY(i); stubZ(i)];
    elseif stubatom_names(i) == "H" && stubres_number(i) == 2
        H2 = [stubX(i); stubY(i); stubZ(i)];
    elseif stubatom_names(i) == "O" && stubres_number(i) == 2
        fake_O_ind = i;
    elseif ismember(stubatom_names(i), ["1H","2H","3H"]) && stubres_number(i) == 1
        fake_H_ind = [fake_H_ind, i];
    end
end
fake_H_ind = sort(fake_H_ind);

pepX = [];
pepY = [];
pepZ = [];
pepatom_names = [];
pepatom_number = [];
pepres_names = [];
pepres_number = [];
pepelements = [];
for fi = 1 : 7
    pep = pdb2mat(fi+"Extension_DASP_DNLU_backbone_sol2_0001.pdb");
    pepX = [pepX, pep.X];
    pepY = [pepY, pep.Y];
    pepZ = [pepZ, pep.Z];
    pepatom_names = [pepatom_names, pep.atomName];
    pepatom_number = [pepatom_number, pep.atomNum];
    pepres_names = [pepres_names, pep.resName];
    pepres_number = [pepres_number, pep.resNum];
    pepelements = [pepelements, pep.element];
end

candEnergyRamp = [];
candTargetRamp = [];
for fi = 1 : 7
    matfile = load(fi+"Extension_DASP_DNLU_backbone_sol2_0001.mat");
    candEnergyRamp = [candEnergyRamp, matfile.candEnergyRamp];
    candTargetRamp = [candTargetRamp, matfile.candTargetRamp];
end

[M,I] = sort(candEnergyRamp);
Coordinates = zeros(7*n,3*length(I));
for i = 1 : length(I)
    Coordinates(:,3*i-2:3*i) = [pepX((I(i)-1)*7*n+1:I(i)*7*n).', pepY((I(i)-1)*7*n+1:I(i)*7*n).', pepZ((I(i)-1)*7*n+1:I(i)*7*n).'];
end
candEnergyRamp = candEnergyRamp(I);
candTargetRamp = candTargetRamp(I);

% indices = clustering_RMSD(Coordinates(repelem([0:n-1]*7,4)+repmat([1,3,6,7],1,n),:));
% Coordinates = Coordinates(:,repelem((indices-1)*3,3)+repmat([1,2,3],1,length(indices)));
% candEnergyRamp = candEnergyRamp(indices);
% candTargetRamp = candTargetRamp(indices);

target_threshold = min_s*size(target_grid,1)*0.8;
npep = size(Coordinates,2) / 3;

NF = 8;
NS = ceil(npep/NF);
peptide = cell(1,npep);
for nf = 1 : NF
    file = fopen("Backbone_list_"+solname+"_"+nf+".txt", "w");
    fclose(file);

    startS = (nf-1)*NS+1;
    endS = nf*NS;
    if nf == NF
        endS = npep;
    end

    parfor np = startS : endS
        coordinates = Coordinates(:,3*npep-2:3*npep);
        E_target = target_score(coordinates, progrid, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid);

        if E_target <= target_threshold
            % Stub1 Rama map check after loop closure
            stubnames = solname.split("_");
            s1 = stubnames(1);
            s1resname = char(s1);
            conf = "";
            if s1resname(1) == 'D'
                s1resname = string(s1resname(1:4));
                conf = "D";
            else
                s1resname = string(s1resname(1:3));
            end

            if ismember(s1resname, ["BP5","DBP5"])
                rama_s1 = load("ramamap_"+lower(conf+"HIS")+".mat").Map;
            elseif ismember(s1resname, ["SCC","DSCC"])
                rama_s1 = load("ramamap_"+lower(conf+"CYS")+".mat").Map;
            else
                rama_s1 = load("ramamap_"+lower(s1resname)+".mat").Map;
            end

            lastC = coordinates(end-1,:);
            phi_s1 = torsion(lastC.',N1,CA1,C1);
            psi_s1 = torsion(N1,CA1,C1,N2);
            phi_s2_ind = floor((phi_s1+pi)/deg2rad(10))+1;
            psi_ind = floor((psi_s1+pi)/deg2rad(10))+1;

            if rama_s1(phi_s2_ind,psi_ind) == 1
                peptide{np}.X = [proX, pepX((np-1)*n*7+10*7+1:np*n*7), pepX((np-1)*n*7+2), stubX(1:fake_H_ind(1)-1), stubX(fake_H_ind(end)+1:fake_O_ind-1), stubX(fake_O_ind+1:end), pepX((np-1)*n*7+2*7:(np-1)*n*7+10*7)];
                peptide{np}.Y = [proY, pepY((np-1)*n*7+10*7+1:np*n*7), pepY((np-1)*n*7+2), stubY(1:fake_H_ind(1)-1), stubY(fake_H_ind(end)+1:fake_O_ind-1), stubY(fake_O_ind+1:end), pepY((np-1)*n*7+2*7:(np-1)*n*7+10*7)];
                peptide{np}.Z = [proZ, pepZ((np-1)*n*7+10*7+1:np*n*7), pepZ((np-1)*n*7+2), stubZ(1:fake_H_ind(1)-1), stubZ(fake_H_ind(end)+1:fake_O_ind-1), stubZ(fake_O_ind+1:end), pepZ((np-1)*n*7+2*7:(np-1)*n*7+10*7)];

                filename = "pdb_input/Backbone_"+solname+"_"+np+".pdb";
                file = fopen(filename, "w");
                fclose(file);

                peptide{np}.outfile = filename;
                peptide{np}.atomName = [proatom_names, pepatom_names((np-1)*n*7+10*7+1:np*n*7), "H", stubatom_names(1:fake_H_ind(1)-1), stubatom_names(fake_H_ind(end)+1:fake_O_ind-1), stubatom_names(fake_O_ind+1:end), "O", pepatom_names((np-1)*n*7+2*7+1:(np-1)*n*7+10*7)];
                peptide{np}.atomNum = 1 : length(peptide{np}.X);
                peptide{np}.resName = [prores_names, pepres_names((np-1)*n*7+10*7+1:np*n*7), stubres_names(3:end), pepres_names((np-1)*n*7+2*7+1:(np-1)*n*7+10*7)];
                peptide{np}.resNum = [prores_number, prores_number(end)+pepres_number((np-1)*n*7+10*7+1:np*n*7)-10, ...
                    prores_number(end)+stubres_number(3:end)+8, prores_number(end)+pepres_number((np-1)*n*7+2*7+1:(np-1)*n*7+10*7)+8];
                peptide{np}.element = [proelements, pepelements((np-1)*n*7+10*7+1:np*n*7), "H", stubelements(1:fake_H_ind(1)-1), stubelements(fake_H_ind(end)+1:fake_O_ind-1), stubelements(fake_O_ind+1:end), "O", pepelements((np-1)*n*7+2*7+1:(np-1)*n*7+10*7)];
                mat2pdb(peptide{np});
                fprintf("\nExtension "+np+" done!\n");

                file = fopen("Backbone_list_"+solname+"_"+nf+".txt", "a");
                fprintf(file, filename+"\n");
                fclose(file);
            else
                fprintf("\nExtension "+np+" phi of "+s1+" not satisfied!\n");
            end
        else
            fprintf("\nExtension "+np+" target threshold not satisfied!\n");
        end
    end
end
end


function chi = torsion(p1, p2, p3, p4)
b1 = p2-p1;
b2 = p3-p2;
b3 = p4-p3;
n1 = cross(b1, b2)/norm(cross(b1, b2));
n2 = cross(b2, b3)/norm(cross(b2, b3));
x = dot(n1, n2);
y = dot(cross(n1, n2), b2/norm(b2));
chi = atan2(y, x);
end

function t = T(angle)
t = [cos(pi-angle), -sin(pi-angle), 0; sin(pi-angle), cos(pi-angle), 0; 0, 0, 1];
end

function r = R(angle)
r = [1, 0, 0; 0, cos(angle), -sin(angle); 0, sin(angle), cos(angle)];
end

function plot_peptide(coordinates, filename)
n = size(coordinates,1)/7;
m = size(coordinates,2)/3;

atom_names = repmat(["N", "H", "CA", "1HA", "2HA", "C", "O"], 1, n);
atom_number = 1:7*n;
res_names = repmat(["GLY"], 1, 7*n);
res_number = repelem([1:n], 7);
elements = repmat(["N", "H", "C", "H", "H", "C", "O"], 1, n);

file = fopen(filename, 'w');
fclose(file);

for pep = 1 : m
    peptide.X = coordinates(:,pep*3-2);
    peptide.Y = coordinates(:,pep*3-1);
    peptide.Z = coordinates(:,pep*3);

    peptide.outfile = filename;
    peptide.atomName = atom_names;
    peptide.atomNum = atom_number;
    peptide.resName = res_names;
    peptide.resNum = res_number;
    peptide.element = elements;
    mat2pdb(peptide);
end
end

% Function to make the angles within [-pi, pi]
function angles_periodic = periodic(angles)

angles_periodic = angles;
for i = 1 : length(angles)
    if angles(i) < -pi
        angles_periodic(i) = angles(i) + 2*pi;
    elseif angles(i) > pi
        angles_periodic(i) = angles(i) - 2*pi;
    end
end
end

% Function to calculate scores for the protein grid
function score = target_score(coordinates, progrid, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid)
score = 0;
target_score = zeros(size(target_grid,1),1);

for i = 1 : size(coordinates,1)
    xpos = floor((coordinates(i,1)-minX)/grid_dx)+1;
    ypos = floor((coordinates(i,2)-minY)/grid_dx)+1;
    zpos = floor((coordinates(i,3)-minZ)/grid_dx)+1;

    s = 0;
    if xpos>=1 && xpos<=nX && ypos>=1 && ypos<=nY && zpos>=1 && zpos<=nZ
        s = progrid(xpos,ypos,zpos);
    end

    if s < 0
        % Only CA atoms contribute to rewards
        if mod(i,7) == 3
            rewards = rewardgrid(xpos,ypos,zpos,:);
            for t = 1 : length(target_score)
                if rewards(t) < target_score(t)
                    target_score(t) = rewards(t);
                end
            end
        end
        s = 0;
    end

    score = score + s;
end

score = score + sum(target_score);
% target_score
end

% Clustering the good candidates
function indices = clustering(Angles)
N = size(Angles,2);
lib = 1 : N;
indices = [];

while ~isempty(lib)
    members = [];
    center = lib(1);
    centerang = Angles(:,center);
    for i = 1 : length(lib)
        ang = Angles(:,lib(i));
        dev = sum(abs(ang-centerang));
        if dev < pi
            members = [members, lib(i)];
        end
    end
    lib = setdiff(lib, members);
    indices = [indices, center];
end
end


function indices = clustering_RMSD(Coordinates)
n = size(Coordinates,1)/4;
N = size(Coordinates,2)/3;
lib = 1 : N;
indices = [];

while ~isempty(lib)
    members = [];
    center = lib(1);
    centercoor = Coordinates(:,center*3-2:center*3);
    for i = 1 : length(lib)
        coor = Coordinates(:,lib(i)*3-2:lib(i)*3);
        [~, ~, rmsd] = kabsch_algorithm(coor, centercoor);
        if rmsd < 0.5
            members = [members, lib(i)];
        end
    end
    lib = setdiff(lib, members);
    indices = [indices, center];
end
end

% KABSCH_ALGORITHM calculates the optimal rigid body transformation
% that aligns two sets of 3D points (P and Q) using the Kabsch algorithm.
%
% Inputs:
%   P: Nx3 array of points to be aligned
%   Q: Nx3 array of target points
%
% Outputs:
%   R: 3x3 rotation matrix
%   t: 3x1 translation vector
%   rmsd: root-mean-square deviation between the aligned points
function [R, t, rmsd] = kabsch_algorithm(P, Q)

% Calculate the centroids of the two sets of points
centroid_P = mean(P, 1);
centroid_Q = mean(Q, 1);

% Center the points by subtracting their centroids
P_centered = P - centroid_P;
Q_centered = Q - centroid_Q;

% Calculate the covariance matrix of the centered points
covariance_matrix = P_centered' * Q_centered;

% Calculate the optimal rotation matrix using singular value decomposition (SVD)
[U, ~, V] = svd(covariance_matrix);
rotation_matrix = V * U';

% If the determinant of the rotation matrix is negative, we need to flip one axis
if det(rotation_matrix) < 0
    V(:, 3) = -V(:, 3);
    rotation_matrix = V * U';
end

% Calculate the translation vector
translation_vector = centroid_Q' - rotation_matrix * centroid_P';

% Apply the rotation and translation to the original set of points
P_aligned = (rotation_matrix * P')' + translation_vector';

% Calculate the root-mean-square deviation (RMSD) between the aligned points
rmsd = sqrt(sum(sum((Q - P_aligned).^2)) / size(P, 1));

% Output the rotation matrix, translation vector, and RMSD
R = rotation_matrix;
t = translation_vector;
end