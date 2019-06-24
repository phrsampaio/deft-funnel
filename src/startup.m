function startup
% Local startup file. Modifies the path.
disp('--------------------------------------------------------------------------------------')
disp('Loading local startup.m ...');
disp('--------------------------------------------------------------------------------------')

cd ..
fprintf('\tAdding deft-funnel directory (including subfolders) to the path ...');

disp(' ')
disp(pwd)
addpath(genpath(pwd));

disp('--------------------------------------------------------------------------------------')
disp('Local startup finished.');
disp('--------------------------------------------------------------------------------------')
