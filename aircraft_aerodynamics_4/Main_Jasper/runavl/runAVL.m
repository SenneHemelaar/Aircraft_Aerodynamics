function a = runAVL(avlfile, rho, Mach)

% Create run file
% Open the file with write permission
fid = fopen('runavl/bin/runfile.x', 'w');

%Load the AVL definition of the aircraft
fprintf(fid, 'LOAD %s\n', strcat(avlfile,'.avl'));

%Load mass parameters
% fprintf(fid, 'MASS %s\n', strcat(filename,'.mass'));
% fprintf(fid, 'MSET\n');
%Change this parameter to set which run cases to apply 
% fprintf(fid, '%i\n',   0); 

%Disable Graphics
fprintf(fid, 'PLOP\ng\n\n'); 

%Open the OPER menu
fprintf(fid, '%s\n',   'OPER');   

% Change forces output
fprintf(fid, 'o\n');  
fprintf(fid, 'p\n');  
fprintf(fid, 't,t,t,f\n');
fprintf(fid, '\n');  

%Define the run case
fprintf(fid, 'a\n'); 
fprintf(fid, 'c 0.6\n');       
% fprintf(fid, 'v %6.4f\n',velocity);
fprintf(fid, 'M\n');   
fprintf(fid, 'MN %g\n', Mach);   
fprintf(fid, 'D %g\n', rho);   
fprintf(fid, '\n');

%Options for trimming
%fprintf(fid, '%s\n',   'd1 rm 0'); %Set surface 1 so rolling moment is 0
%fprintf(fid, '%s\n',   'd2 pm 0'); %Set surface 2 so pitching moment is 0

%Run the Case
fprintf(fid, 'x\n'); 

%Save the force data
fprintf(fid, '%s\n',   'w'); 
fprintf(fid, 'runavl/res/%s%s\n',avlfile,'.forces');   

%Drop out of OPER menu
fprintf(fid, '\n');
 
%Quit Program
fprintf(fid, 'Quit\n'); 

%Close File
fclose(fid);

%Run AVL using 
%dos(strcat('runavl\avl.exe',' bin\',avlfile,'.avl'));
dos(strcat('runavl\avl.exe',' <  ','runavl\bin\runfile.x'));
% [status,result] = dos(strcat('bin/runfile.x'));
% [status,result] = dos(strcat('.\avl',' bin/',avlfile,' < ',' bin/runfile.x'));
end
