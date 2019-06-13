function [Fitness,output] = UserDefined_LPFitness (LP)

keyboard

%  Add your fitness calculation as a function of the input Lamination
% Parameters returned by the GA.


% ---
% Fitness =  To COMPLETE HERE                     %(this field is compulsory)
% ---


% The output structure can contain any data you whish the code to output.
% Addtionally, if you have constraints to enforce indirectly you must
% somehow relate the number of violated constraints with fitness
% calculations, for instance:
% Fitness =  Fitness * NviolatedConstraints


% ---               
% output.NViolatedConst =  To COMPLETE HERE      % set to 0 if no constraint is calculated (this field is compulsory) - used during inipop generation 
% output.YourFIELD  =                            % optional output you which to display at the end of the optimisation
% ---


% If more input to the fitness function are required (e.g. Input2), you will have to use a wrapper
% function because the default Matlab GA only accepted fitness functions 
% with 1 vector input of design variables. 
% Wrapper function example:
% WrappedFitnessFct = @(LP) UserDefine_LPFitness(LP,Input2);


end