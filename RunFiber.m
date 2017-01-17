% This script run Fiber class and output modes refractive index and U
% values
clear
clc
format long
% wavelength of source in meter
lambda = 650e-9; 
% core radius in meter
%r = 5.5e-6; 
r = 4.5e-6; 
% core refrative index
%Nc = 1.456  1.5699;
Nc = 1.448918;
% cladding refractive index
Ng = 1.444418;

fi = Fiber(lambda,r,Nc,Ng);
fi.eigenValuesTE;
disp('Fiber TE U values');
disp(fi.uTE);
fi.eigenValuesTM;
disp('Fiber TM U values');
disp(fi.uTM);
fi.eigenValuesHyb;
disp('Fiber HE U values');
disp(fi.uHE);
disp('Fiber EH U values');
disp(fi.uEH);
% ==================================
disp('Fiber TE Neff values');
disp(fi.NeffTE);
disp('Fiber TM Neff values');
disp(fi.NeffTM);
disp('Fiber HE Neff values');
disp(fi.NeffHE);
disp('Fiber EH Neff values');
disp(fi.NeffEH);

