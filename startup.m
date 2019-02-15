function startup
% Local startup file. Modifies the path.
disp('--------------------------------------------------------------------------------------')
disp('Loading local startup.m ...');
disp('--------------------------------------------------------------------------------------')

fprintf('\tAdding current directory (including subfolders) to the path ...');

disp(' ')
disp(pwd)
addpath(genpath(pwd));

disp('--------------------------------------------------------------------------------------')
disp('Local startup finished.');
disp('--------------------------------------------------------------------------------------')