function [Cl_a, Cl_a_sec_fn, Cd_a, Cd_a_sec_fn, ...
          chord, pitch] = Import_Data(D)
%IMPORT_DATA

%% Import lift data
fid = csvread('CL_alpha.csv', 2);

Cl_a.sec_1 = [-4 6; -0.0001 0.0001];
Cl_a.sec_2 = fid(1:2,1:2)';
Cl_a.sec_3 = fid(1:2,3:4)';
Cl_a.sec_4 = fid(1:2,5:6)';
Cl_a.sec_5 = fid(1:2,7:8)';
Cl_a.sec_6 = fid(1:2,9:10)';
Cl_a.sec_7 = fid(1:2,11:12)';
Cl_a.sec_8 = fid(1:2,13:14)';
Cl_a_sec_fn = fieldnames(Cl_a);

%% Import L/D data
fid = csvread('L_D_alpha.csv', 2);

Cd_a.sec_1 = [-6 6; 0.00025 0.00025];
Cd_a.sec_2 = [-6 6; 0.001 0.001];

a = nonzeros(fid(:,1)');
b = nonzeros(fid(:,2)');
Cd_a.sec_3 = [a';b'];

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
Cd_a.sec_4 = [a';b'];

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
Cd_a.sec_5 = [a';b'];

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
Cd_a.sec_6 = [a';b'];

a = nonzeros(fid(:,9)');
b = nonzeros(fid(:,10)');
Cd_a.sec_7 = [a';b'];

a = nonzeros(fid(:,11)');
b = nonzeros(fid(:,12)');
Cd_a.sec_8 = [a';b'];

Cd_a_sec_fn = fieldnames(Cd_a);

%% Import chord length distribution
fid = csvread('Chord_Pitch_radial.csv',2);

chord = fid(:,1:2)';
chord(2,:) = chord(2,:)*D;

%% Import pitch distribution
fid = csvread('Chord_Pitch_radial.csv',2);

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
pitch.deg_25 = [a';b'];
pitch.deg_25(2,:) = pitch.deg_25(2,:)*D;

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
pitch.deg_35 = [a';b'];
pitch.deg_35(2,:) = pitch.deg_35(2,:)*D;

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
pitch.deg_45 = [a';b'];
pitch.deg_45(2,:) = pitch.deg_45(2,:)*D;


end

