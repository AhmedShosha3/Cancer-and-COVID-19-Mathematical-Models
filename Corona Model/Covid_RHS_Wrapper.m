function f = Covid_RHS_Wrapper(t, y, params)
% returns column vector f for cancer RHS
f = Covid_Model(t, y, params); % reuse cancerSystem (it returns column)
end
