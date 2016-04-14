fprintf('Compiling MEX files\n');
cd MEX
make
fprintf('Compiling METIS MEX files\n');
cd ../clustering/metis-5.0.2/metismex/
make
cd ../../../