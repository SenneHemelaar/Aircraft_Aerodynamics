function writeAVL(avlfile, lw, cwr, lamw, epsr, epst, Lamw,...
                  phiw, Sref, cref, bref, Ncm, Nsm, Ncw, Nsw)
% Create AVL file

% Open the file with write permission
fid = fopen(strcat(avlfile,'.avl'), 'w');

fprintf(fid, 'Plane\n');
Mach = 0.7;
fprintf(fid, '%f\n', Mach);     
fprintf(fid, '0 0 0.0\n');
fprintf(fid, '%f %f %f\n', Sref, cref, bref);
fprintf(fid, '0.50 0.0 0.0\n\n');

%%%============================== WING =================================%%%
fprintf(fid, 'SURFACE\nMainWing\n');
fprintf(fid, '%i 1.0 %i -2.0\n', Ncm, Nsm);
fprintf(fid, 'COMPONENT\n1\nYDUPLICATE\n0.0\nANGLE\n0.0\n');
% Sections
fprintf(fid, 'SECTION\n0. 0. 0. 5.5 0.0 0 0\nCLAF\n1.0\n');
fprintf(fid, 'SECTION\n3.5 14.0 0.0 2.0 0.0 0 0\nCLAF\n1.0\n\n');

%%%============================ WINGLET ================================%%%
fprintf(fid, 'SURFACE\nWinglet\n');
fprintf(fid, '%i 1.0 %i 1.0\n', Ncw, Nsw);
fprintf(fid, 'COMPONENT\n1\nYDUPLICATE\n0.0\nANGLE\n0.0\n');
% Sections
fprintf(fid, 'SECTION\n%f 14.0 0.0 %f %f 0 0\nCLAF\n1.0\n',...
        5.5-cwr, cwr, epsr);
fprintf(fid, 'SECTION\n%f %f %f %f %f 0.0 0\nCLAF\n1.0\n',...
        5.5-cwr+lw/cos(Lamw), 14.0+lw*sin(phiw), lw*cos(phiw),...
        lamw*cwr, epst);

fclose(fid);
end
