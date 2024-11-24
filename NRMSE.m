function error = NRMSE(obj,exa_obj)
%NRMSE Calculates the normalized root-mean-square error.
%==========================================================================

norm_obj = abs(abs(obj) - abs(exa_obj));
error = sqrt(mean(norm_obj.^2,[1 2])) / sqrt(mean(abs(exa_obj).^2,[1 2]));

end

