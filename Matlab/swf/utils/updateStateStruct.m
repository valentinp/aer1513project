function [updatedStates] = updateStateStruct( currentStates, dx)
%UPDATESTATESTRUCT Updates states by applying a step dx

updatedStates = currentStates;

for stIdx = 1:length(currentStates)
    dr = dx(1+(stIdx-1)*6:3+(stIdx-1)*6);
    phi = dx(4+(stIdx-1)*6:6+(stIdx-1)*6);
    
    updatedStates{stIdx}.r_vi_i = currentStates{stIdx}.r_vi_i + dr;
    updatedStates{stIdx}.C_vi = Cfrompsi(phi)*currentStates{stIdx}.C_vi;

end


end

