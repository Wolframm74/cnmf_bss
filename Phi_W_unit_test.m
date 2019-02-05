%, outtensor_2nd_real_FOM_mex, outtensor_2nd_cx_FOM_mex)
function Phi_W_unit_test(outtensor_1st_real_NOF_mex, outtensor_1st_cx_NOF_mex)

display('Hello World');

save('/home/tung/Documents/School/thesis2/mycode/nguyen_2015_2/code_vectorised_16/bottom_up/mex/W_outtensors_mex_file.mat', 'outtensor_1st_real_NOF_mex', 'outtensor_1st_cx_NOF_mex');
% save('/home/tung/Documents/School/thesis2/mycode/matlab/nguyen_2015_2/code_vectorised_16/bottom_up/mex/W_outtensors_mex_file.mat', 'outtensor_1st_real_NOF_mex', 'outtensor_1st_cx_NOF_mex', 'outtensor_2nd_real_FOM_mex', 'outtensor_2nd_cx_FOM_mex');

% load('/home/tung/Documents/School/thesis2/mycode/matlab/nguyen_2015_2/code_vectorised_16/bottom_up/mex/W_outtensors_file.mat');
% 
% outtensor_1st_real_NOF_diff=sum(sum(sum(outtensor_1st_real_NOF-outtensor_1st_real_NOF_mex)))
% outtensor_1st_cx_NOF_diff=sum(sum(sum(outtensor_1st_cx_NOF-outtensor_1st_cx_NOF_mex)))
% outtensor_2nd_real_FOM_diff=sum(sum(sum(outtensor_2nd_real_FOM-outtensor_2nd_real_FOM_mex)))
% outtensor_2nd_cx_FOM_diff=sum(sum(sum(outtensor_2nd_cx_FOM-outtensor_2nd_cx_FOM_mex)))

end