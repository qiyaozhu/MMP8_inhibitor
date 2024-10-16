function DASP_DNLU_Backbone_extension_1

addpath("/mnt/home/qzhu1/my_functions");
addpath("/mnt/home/qzhu1/RamaMap");

rama_phe = load('ramamap_phe.mat').Map;
rama_dphe = load('ramamap_dphe.mat').Map;

filename = "DASP_DNLU_backbone_sol2_0001.pdb";
peptide = pdb2mat(filename);
X = peptide.X;
Y = peptide.Y;
Z = peptide.Z;
atom_names = peptide.atomName;
res_num = peptide.resNum;

Stubs = [];
for i = 1 : length(X)
    if ~((res_num(i)==1 && ismember(atom_names(i), ["1H", "2H", "3H"])) || (res_num(i)==2 && atom_names(i)=="O"))
        Stubs = [Stubs; X(i), Y(i), Z(i)];
    end
end

for i = 1 : length(X)
    if atom_names(i) == "N"
        if res_num(i) == 1
            N1 = [X(i); Y(i); Z(i)];
        else
            N2 = [X(i); Y(i); Z(i)];
        end
    elseif atom_names(i) == "CA"
        if res_num(i) == 1
            CA1 = [X(i); Y(i); Z(i)];
        else
            CA2 = [X(i); Y(i); Z(i)];
        end
    elseif atom_names(i) == "C"
        if res_num(i) == 1
            C1 = [X(i); Y(i); Z(i)];
        else
            C2 = [X(i); Y(i); Z(i)];
        end
    elseif atom_names(i) == "O" && res_num(i) == 1
        O1 = [X(i); Y(i); Z(i)];
    elseif atom_names(i) == "H" && res_num(i) == 2
        H2 = [X(i); Y(i); Z(i)];
    end
end

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

% rotation matrices and translation vector
T_NR = [0.5255, -0.8508, 0; 0.8508, 0.5255, 0; 0, 0, 1];
T_CaR = [0.3616, -0.9323, 0; 0.9323, 0.3616, 0; 0, 0, 1];
T_CR = [0.4415, -0.8973, 0; 0.8973, 0.4415, 0; 0, 0, 1];
R_omega = [1, 0, 0; 0, -1, 0; 0, 0, -1];

x = CA1 - N1;
x = x/norm(x);
u = C1 - N1;
y = u - dot(u,x)*x;
y = y/norm(y);
z = cross(x, y);

phi_s2 = torsion(C1,N2,CA2,C2);
CAS1R = acos(dot((C1 - CA1), (N1 - CA1)) / norm(C1 - CA1) / norm(N1 - CA1));
CAS2R = acos(dot((C2 - CA2), (N2 - CA2)) / norm(C2 - CA2) / norm(N2 - CA2));
rotation = T(CAS1R) * R(torsion(N1,CA1,C1,N2)) * T_CR * R_omega * T_NR * R(phi_s2) * T(CAS2R);

% Check Rama map
stubnames = filename.split("_");
s1 = stubnames(1);
s1resname = char(s1);
conf = "";
if s1resname(1) == 'D'
    s1resname = string(s1resname(1:4));
    conf = "D";
else
    s1resname = string(s1resname(1:3));
end
psi_s1 = torsion(N1,CA1,C1,N2);
psi_s1_ind = floor((psi_s1+pi)/deg2rad(10))+1;
if ismember(s1resname, ["BP5","DBP5"])
    rama_s1 = load("ramamap_"+lower(conf+"HIS")+".mat").Map;
elseif ismember(s1resname, ["SCC","DSCC"])
    rama_s1 = load("ramamap_"+lower(conf+"CYS")+".mat").Map;
else
    rama_s1 = load("ramamap_"+lower(s1resname)+".mat").Map;
end

% Stub 2 Rama
phi_s2_ind = floor((phi_s2+pi)/deg2rad(10))+1;
s2 = stubnames(2);
s2resname = char(s2);
conf = "";
if s2resname(1) == 'D'
    s2resname = string(s2resname(1:4));
    conf = "D";
else
    s2resname = string(s2resname(1:3));
end

if ismember(s2resname, ["NLU","DNLU"])
    rama_s2 = load("ramamap_"+lower(conf+"LEU")+".mat").Map;
    display("Loading RamaMap of "+lower(conf+"LEU"));
elseif conf == ""
    rama_s2 = rama_phe;
else
    rama_s2 = rama_dphe;
end
psimap_s2 = find(rama_s2(phi_s2_ind,:)== 1);

% Total size
n = 18;

% Ramachandran plot for symmetric glycine
Boundaries = load('ramabin_glycine.mat').Boundaries.';

% atom properties
[LJ_radius, LJ_well, Coulumb, LK_DG, LK_Lambda, LK_V] = FA_parameter(n); % parameters
D = n_bonds(n); % number of bonds between any atom pair

% Simulated annealing parameters
w_rep = 0.5;
w_hbd_ramping = 10;
w_cyc = 10;
w_target = 2;

t0_ramping = 4*n;
k0 = 0.6;
b = 18;
c = 20;

M = 20000;
ramp0 = 1;

% criteria for good candidates
rep_cri = polyval([0.5932    0.9627], n);
cyc_cri = 0.7;
count_cri = ceil(n/3);
target_cri = min_s*size(target_grid,1)*0.8;
touch_threshold = 1.5;

% Different initial angles
InitAngles = [-2, -1.7, -1.9, 1.9, 1.7, 2; ...
    -2.4, 0.17, 2.6, -2.6, -0.17, 2.4];

% Simulated annealing to find cyclic peptides with at least n/3 hydrogen
% bonds and low energies, espeically the repulsive VDW energy
GoodAnglesRamp = [];
GoodEnergyRamp = [];
GoodTargetRamp = [];

IP = load("initialpoints_combdesign.mat").InitialPoints;
Repeat = size(IP,2)/8;

Angles_Initial = zeros(2*n-3, Repeat);
for repeat = 1 : Repeat
    psi_s2 = (IP(1,repeat)-1)*deg2rad(10)-pi;
    bin = IP(2:end,repeat);
    angles_initial = [psi_s2; reshape(InitAngles(:,bin),[],1)];
    Angles_Initial(:,repeat) = angles_initial;
end

parfor repeat = 1 : Repeat
    % initial coordinates and energies
    angles_initial = Angles_Initial(:,repeat);
    [coordinates, cyc] = peptide_cyc_MMP8(angles_initial, N1, CA1, C1, O1, N2, H2, CA2, C2, rotation, x, y, z);
    rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
    [hbond_ramping, count, oversat] = E_hbond_ramping(coordinates, ramp0);
    E_target = target_score(coordinates, progrid, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid);
    E_total_ramping = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping + w_target*E_target;

    accept_ramping = 0;
    angles = angles_initial;
    for i = 1 : M
        % generate new random move
        k = k0/(1+b*i/M);
        psi_p = k*(2*rand-1);
        if rama_s2(phi_s2_ind, floor((periodic(psi_s2+psi_p) + pi) / deg2rad(10)) + 1) == 1
            psi_s2_new = psi_s2 + psi_p;
        else
            psi_s2_new = psi_s2;
        end        

        p = random_move(k,n-2);
        [~, extension_new] = feasible(p, angles(2:end), Boundaries);
        angles_new = [psi_s2_new; extension_new];

        ramp = ramp0*(M-i)/M;
        [coordinates, cyc] = peptide_cyc_MMP8(angles_new, N1, CA1, C1, O1, N2, H2, CA2, C2, rotation, x, y, z);
        rep = fa_rep_manhattan(coordinates, D, LJ_radius, LJ_well);
        [hbond_ramping, count, oversat] = E_hbond_ramping(coordinates, ramp);
        E_target = target_score(coordinates, progrid, rewardgrid, minX, minY, minZ, grid_dx, nX, nY, nZ, target_grid);
        E_total_ramping_new = w_cyc*cyc + w_rep*rep + w_hbd_ramping*hbond_ramping + w_target*E_target;

        pass = false;
        temp = t0_ramping/(1+c*i/M);
        if E_total_ramping_new <= E_total_ramping
            pass = true;
        else
            prob = exp(1)^((E_total_ramping-E_total_ramping_new)/temp);
            if rand <= prob
                pass = true;
            end
        end

        % pass all layers, this random move is accepted
        if pass
            angles = angles_new;
            E_total_ramping = E_total_ramping_new;
            accept_ramping = accept_ramping + 1;

            if rep <= rep_cri && cyc <= cyc_cri && count >= count_cri && oversat <= 0 && E_target <= target_cri
                lastC = coordinates(end-1,:);
                phi_s1 = torsion(lastC.',N1,CA1,C1);
                phi_s1_ind = floor((phi_s1+pi)/deg2rad(10))+1;
                if rama_s1(phi_s1_ind,psi_s1_ind) == 1
                    % Check if steric clashes with the stubs
                    stubs = [Stubs; coordinates(2,:); coordinates(14,:)];
                    backbone = [coordinates(16:21,:); coordinates(end-6:end-2,:); coordinates(end,:)];
                    Dis = pdist2(stubs, backbone);
                    Dis_score = sum(Dis <= touch_threshold, "all");
                    display("Steric clash with the stubs "+Dis_score);

                    if Dis_score == 0
                        GoodAnglesRamp = [GoodAnglesRamp, angles];
                        GoodEnergyRamp = [GoodEnergyRamp, E_total_ramping];
                        GoodTargetRamp = [GoodTargetRamp, E_target];
                    end
                end
            end
        end
    end

    display("Ramp accept = "+accept_ramping);
end

CandAnglesRamp = [];
candEnergyRamp = [];
candTargetRamp = [];

if ~isempty(GoodEnergyRamp)
    [M,I] = sort(GoodEnergyRamp);
    GoodAnglesRamp = GoodAnglesRamp(:,I);
    indices = clustering(GoodAnglesRamp);
    CandAnglesRamp = GoodAnglesRamp(:,indices);
    candEnergyRamp = GoodEnergyRamp(I(indices));
    candTargetRamp = GoodTargetRamp(I(indices));

    Coordinates = zeros(7*n,3*length(indices));
    for i = 1 : length(indices)
        [coordinates, cyc] = peptide_cyc_MMP8(CandAnglesRamp(:,i), N1, CA1, C1, O1, N2, H2, CA2, C2, rotation, x, y, z);
        Coordinates(:,3*i-2:3*i) = coordinates;
    end
    
    indices = clustering_RMSD(Coordinates(repelem([0:n-1]*7,4)+repmat([1,3,6,7],1,n),:));
    Coordinates = Coordinates(:,repelem((indices-1)*3,3)+repmat([1,2,3],1,length(indices)));
    candEnergyRamp = candEnergyRamp(indices)
    candTargetRamp = candTargetRamp(indices)

    plot_peptide(Coordinates, "1Extension_"+filename);
    prefix = filename.split(".");
    save("1Extension_"+prefix(1)+".mat", "candEnergyRamp", "candTargetRamp");
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