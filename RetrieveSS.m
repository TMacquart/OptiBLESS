% ----------------------------------------------------------------------- %
% Copyright (c) <2015>, <Terence Macquart>
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those
% of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of the FreeBSD Project.
% ----------------------------------------------------------------------- %


% =====                                                              ==== %
%    Use GA to find stacking sequence ply angles matching Lam. Param.     %
%          and varying thickness employing Stacking Sequence Table        %
%                       Terence Macquart (09/11/2015)                     %
%
%  Blended laminate
%  ----------------
%  objective:   match other LP by removing plies while minimising thickness 
%  constraints implemented: covering, max stoping, continuity, alternation
%
% 
%  The individual of GA is composed of 3 parts:
%  --------------------------------------------
%  [ [Nply(1) ... Nply(Nsec)]                   -- the Number of plies
%  [ Theta(1) ... Theta(N)   nan ... nan]       -- N is the guide Nply
%  [ Drop(1)  ... Drop(M)    nan ... nan]       -- M is the Delta Nply
% =====                                                              ==== %

function [output] = RetrieveSS(Objectives,Constraints,GAoptions)


% ---
if 1    % Define objectives (i.e. LP to match)

   NplyIni = cell2mat(Objectives.Table(2:end,2));

    if Constraints.Sym    && Constraints.Balanced,
        LamType = 'Balanced_Sym';
        NpliesUpBound   = round(NplyIni(:,2)/4)*4;
        NpliesLowBound  = round(NplyIni(:,1)/4)*4;
        Nmax = max(NpliesUpBound)/4 ;
        Nmin = min(NpliesLowBound)/4;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):4:NpliesUpBound(i))/4;
        end
    end
    if Constraints.Sym    && ~Constraints.Balanced,
        LamType = 'Sym';
        NpliesUpBound  = round(NplyIni(:,2)/2)*2;
        NpliesLowBound  = round(NplyIni(:,1)/2)*2;
        Nmax = max(NpliesUpBound)/2 ;
        Nmin = min(NpliesLowBound)/2;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):2:NpliesUpBound(i))/2;
        end
    end
    if ~Constraints.Sym   && Constraints.Balanced,
        LamType = 'Balanced';
        NpliesUpBound   = round(NplyIni(:,2)/2)*2;
        NpliesLowBound  = round(NplyIni(:,1)/2)*2;
        Nmax = max(NpliesUpBound)/2 ;
        Nmin = min(NpliesLowBound)/2;
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = (NpliesLowBound(i):2:NpliesUpBound(i))/2;
        end
    end
    if  ~Constraints.Sym  && ~Constraints.Balanced,
        LamType = 'Generic';
        NpliesUpBound   = round(NplyIni(:,2));
        NpliesLowBound  = round(NplyIni(:,1));
        Nmax = max(NpliesUpBound) ;
        Nmin = min(NpliesLowBound);
        for i=1:length(NpliesUpBound)
            AllowedNplies{i} = NpliesLowBound(i):NpliesUpBound(i);
        end
    end

    [~,sortIndex] = sort(NplyIni(:,1),'descend');
    AllowedNplies = AllowedNplies(sortIndex);
    SortedTable   =  [Objectives.Table(1,:); Objectives.Table(1+sortIndex,:)];
    
    if size(NplyIni,1) == 1, Nmin = Nmax; end % particular case for Unique SS retrieval
    

   Nrange = cellfun(@max,AllowedNplies,'UniformOutput', true) - cellfun(@min,AllowedNplies,'UniformOutput', true); % max ply - min ply per lam.
   Nsec = zeros(length(Nrange),1);
   Nsec(Nrange>0)=1; % number of variable thickness lam.

end
% keyboard


% ---
if 1    % Set GA options
    options  = gaoptimset('PopulationSize' ,GAoptions.Npop,'Generation',GAoptions.Ngen, ...
        'StallGenLimit',GAoptions.NgenMin,'EliteCount',ceil(GAoptions.Elitism*GAoptions.Npop),...
        'FitnessLimit' ,1e-5,'TolFun' ,1e-10,'PlotInterval',1,'CrossoverFraction',GAoptions.PC);
    if GAoptions.Plot
        options  = gaoptimset(options,'PlotFcns' ,{@GACustomPlot}); %gaplotbestf
    end
    
    ConstraintVector = Constraints.Vector;                                  % [Damtol  Rule10percent  Disorientation  Contiguity  DiscreteAngle  InernalContinuity  Covering];
    DeltaAngle       = Constraints.DeltaAngle;
    
    if ~Constraints.Balanced
        Nvar = sum(Nsec) + Nmax +(Nmax-Nmin) ;   % Nsec (thickness) + Nmax Guide plies + (Nmax-Nmin) potential drops 
    else
        Nvar = sum(Nsec) + Nmax*2 +(Nmax-Nmin) ; % if balanced, need ply shuffling = + Nmax var
    end
    
    fct_handle = @(x)Eval_Fitness(x,Objectives,Constraints,Nsec,AllowedNplies,SortedTable,LamType);  % handle of the fitness function becomes constraint
    
    % ---
    if 1  % design variable boundaries
        if ConstraintVector(5) % if DiscreteAngle
            IntegerDV = 1:Nvar;                                             % discrete
            Nd_state  = length(-90:DeltaAngle:90);                           % number of discrete state
            LB = 0*ones(Nmax,1)            ;
            UB = (Nd_state-1)*ones(Nmax,1);
        else
            IntegerDV = [[1:sum(Nsec)] [(sum(Nsec)+Nmax):(sum(Nsec)+2*Nmax-Nmin)]];             % continuous
            LB = -90*ones(Nmax,1);
            UB = 90*ones(Nmax,1);
        end
        
        if Constraints.Balanced % add suffling Bounds
            LB = [LB; ones(Nmax,1)];
            UB = [UB; 2*Nmax*ones(Nmax,1)];
        end
        
         if ConstraintVector(5) % if DiscreteAngle
            LB = [LB            ; ones(Nmax-Nmin,1)];
            UB = [UB ; Nmax*ones(Nmax-Nmin,1)];
         else
            LB = [LB; ones(Nmax-Nmin,1)];
            UB = [UB;  Nmax*ones(Nmax-Nmin,1)];
         end
        
        
        for iply=sum(Nsec):-1:1
            LB = [AllowedNplies{iply}(1)   ;LB];
            UB = [AllowedNplies{iply}(end) ;UB];
        end
        
        if ConstraintVector(1) % if Damtol, make the first ply +- 45
            if ~ConstraintVector(5),    
                IntegerDV = [sum(Nsec)+1 IntegerDV];
            end
            LB(sum(Nsec)+1) = 1;
            UB(sum(Nsec)+1) = 2;
        end
        

    end

end

% keyboard
% --- Generate Ini. Pop.
for i = 1:5
    try
        [IniPop] = Generate_IniPop (Nvar,GAoptions.Npop,Nmax,Nmin,Constraints,Nsec,AllowedNplies,LamType);
        break; 
    catch
        fprintf('Retrying IniPop Gen. ...\n');
        if i == 5
            error(' Did not manage to generate a feasible initial population.')
%             
%              UserResponse = input('Non-feasible Elmts in the IniPop, Do you want to continue (y|n)','s');
%              if strcmp(UserResponse,'y')
%                  temp = Constraints;
%                  temp.Vector = 0* temp.Vector;
%                  [IniPop] = Generate_IniPop (Nvar,GAoptions.Npop,Nsec,Nmax,Nmin,temp,LamType);
%              else 
%                  error('Code stopped by the user')
%              end
        end
    end
end
% IniPop(1,:) = [3 1 4 2 5 6 7 8 2]
% IniPop(1,:) = [[3 1 4 2 1 3 4 2] fliplr([3 1 4 2 1 3 4 2])]
% IniPop(1,:) = [3 1 4 2 , 2 4 1 3 , 5 6 7 8 , 9 10 11 12] 
% IniPop(1,:) = [3 4 1 2 1 0 3 2]

%  IniPop(1,:) = [(90+[45   -15    30   -55    40   -80   -40   -80   -70    60])/Constraints.DeltaAngle [2 3 8]]
options = gaoptimset(options,'InitialPopulation' ,IniPop);



%  keyboard

% --- run GA
fprintf(strcat('Running GA \n'))
[xOpt,fval,~,StandardOutputGA] = ga(fct_handle,Nvar,[],[],[],[],LB,UB,[],IntegerDV,options);
display('GA(s) Terminated Successfully')



% --- Results
[~,output] = fct_handle(xOpt);

SS = output.SS;

if strcmp(Objectives.Type,'LP')
    LPMatched = output.LP;
    Table = [{'Lam #'} {'Nplies SST'} {'Ply Angles'} {'LP2Match'} {'LPOpt'} {'NormE'} {'RMSE'} {'MAE'} {'MaxAE'}];
    for j = 1:length(AllowedNplies)
        LP2Match    = Objectives.Table{j+1,3};
        ScalingCoef = Objectives.Table{j+1,4};
        
        
        QualIndex1 = norm( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex2 = rms( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex3 = mae( (LPMatched(:,j) - LP2Match(:)).*ScalingCoef );
        QualIndex4 = max( abs((LPMatched(:,j) - LP2Match(:)).*ScalingCoef) );
        Table      = [Table ;  {j} {length(SS{j})} SS(j) {LP2Match} {LPMatched(:,j)} {QualIndex1} {QualIndex2} {QualIndex3} {QualIndex4}];
    end
end


if strcmp(Objectives.Type,'ABD')

    Table = [{'Lam #'} {'Nplies SST'} {'Ply Angles'} {'A2Match'} {'AOpt'} {'Error % A'} {'Error Norm A'} {'Error RMS A'} ...
            {'B2Match'} {'BOpt'} {'Error % B'} {'Error Norm B'} {'Error RMS B'} ...
            {'D2Match'} {'DOpt'} {'Error % D'} {'Error Norm D'} {'Error RMS D'}];
    for j = 1:length(AllowedNplies)
        A_Matched = output.A{j};
        B_Matched = output.B{j};
        D_Matched = output.D{j};
        A2Match   = Objectives.Table{j+1,3};
        B2Match   = Objectives.Table{j+1,4};
        D2Match   = Objectives.Table{j+1,5};
        
        AScaling = Objectives.Table{j+1,6};
        BScaling = Objectives.Table{j+1,7};
        DScaling = Objectives.Table{j+1,8};
    
        QualIndex1A = 100*sum(abs(  AScaling(:).*((A_Matched(:) - A2Match(:))./A2Match(:)) ));
        QualIndex2A = norm( AScaling(:).*(A_Matched(:) - A2Match(:)) );
        QualIndex3A = rms(  AScaling(:).*(A_Matched(:) - A2Match(:)) );
        
        QualIndex1B = 100*sum(abs(  BScaling(:).*((B_Matched(:) - B2Match(:))./B2Match(:)) ));
        QualIndex2B = norm( BScaling(:).*(B_Matched(:) - B2Match(:)) );
        QualIndex3B = rms(  BScaling(:).*(B_Matched(:) - B2Match(:)) );
        
        QualIndex1D = 100*sum(abs(  DScaling(:).*((D_Matched(:) - D2Match(:))./D2Match(:)) ));
        QualIndex2D = norm( DScaling(:).*(D_Matched(:) - D2Match(:)) );
        QualIndex3D = rms(  DScaling(:).*(D_Matched(:) - D2Match(:)) );
        
        Table = [Table ;  {j} {length(SS{j})} SS(j) ...
                    {A2Match} {A_Matched} {QualIndex1A} {QualIndex2A} {QualIndex3A}...
                    {B2Match} {B_Matched} {QualIndex1B} {QualIndex2B} {QualIndex3B} ...
                    {D2Match} {D_Matched} {QualIndex1D} {QualIndex2D} {QualIndex3D}];
    end

end

output.NfctEval  = StandardOutputGA.funccount;
output.NGen      = StandardOutputGA.generations;
output.Table     = Table;
output.xOpt      = xOpt;
output.fval      = fval;

if ~output.FEASIBLE
    warning('The SST GA could not find a single feasible solution!')
end
end
