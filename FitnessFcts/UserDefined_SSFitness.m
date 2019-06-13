function [Fitness,output] = UserDefined_SSFitness(SS) % function Fitness = UserDefined_SSFitness(SS,INPUT1,...,INPUTN)

keyboard

% Add your fitness calculation as a function of the input Stacking Sequence
% returned by the GA
% Fitness = TO COMPLETE HERE                        %(this field is compulsory)


% The output structure can contain any data you whish the code to output.
% Addtionally, if you have constraints to enforce indirectly you must
% somehow relate the number of violated constraints with fitness
% calculations, for instance:
% Fitness =  Fitness * NviolatedConstraints

% output.NViolatedConst =  TO COMPLETE HERE        % set to 0 if no constraint is calculated (this field is compulsory) - used during inipop generation 
% output.YourFIELD  =                              % optional output


% If more input are required (e.g. Input1), you will have to use a wrapper
% function because the default Matlab GA only accepted fitness functions 
% with 1 vector input of design variables. 
% Wrapper function example:
% WrappedFitnessFct = @(SS) UserDefined_SSFitness(SS,Input2);


end