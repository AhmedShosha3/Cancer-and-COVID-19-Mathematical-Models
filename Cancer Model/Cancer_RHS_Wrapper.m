function f = Cancer_RHS_Wrapper(t, y, C)
% returns column vector f for cancer RHS
f = Cancer_Model(t, y, C); % reuse cancerSystem (it returns column)
end
