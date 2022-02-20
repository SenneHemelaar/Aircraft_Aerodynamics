function [ forces ] = forces(filename)

avlfile = strcat(filename,'.forces');
file = textread(avlfile, '%s', 'delimiter', '\n','whitespace', '');

line = split(file{26},' ');
forces.CDind = str2double(line{end});

line = split(file{27},' ');
forces.CDff = str2double(line{16});

forces.strip1 = importdata(avlfile,' ',72);
forces.strip1 = forces.strip1.data;
% try
forces.strip2 = importdata(avlfile, ' ',102+2*size(forces.strip1,1));
forces.strip2 = forces.strip2.data;
% end
end