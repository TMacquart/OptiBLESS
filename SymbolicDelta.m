clear all;close all;clc; format short g; format compact;


syms MT1 MT2 MT3 MT4 MT5 MT6 MT7 MT8 MT9 MT10 DT1 DT2 DT3 DT4 DT5 DT6 DT7 DT8 DT9 DT10
BalancedVector_S = [MT2 MT3 MT4 MT5 MT6 MT7 MT8 MT9 MT10];
DeltaVector_S    = [DT2 DT3 DT4 DT5 DT6 DT7 DT8 DT9 DT10];


Nrun = 100;
NInfeasible = 0;
Ndisorientation = zeros(Nrun,1);
for iRun = 1:Nrun
    display(iRun)
    Solvable = true;
    %% Implicit Approach
    if 0
        Thetas      = randi(37,1,10)*5 -5 -90;
        BalancedLoc = randi(20,1,10);
        
        ThetasBalanced = [Thetas; [1:length(Thetas)]; [1:length(Thetas)]];
        
        BalancedAngles = [-Thetas' BalancedLoc' [1:length(BalancedLoc)]'];
        BalancedAngles = sortrows(BalancedAngles,2)';
        
        for j = 1:length(BalancedAngles(1,:))
            if BalancedAngles(2,j)>size(ThetasBalanced,2)
                ThetasBalanced = [ThetasBalanced BalancedAngles(:,j)];
            else
                ThetasBalanced = [ThetasBalanced(:,1:BalancedAngles(2,j)-1) BalancedAngles(:,j)  ThetasBalanced(:,BalancedAngles(2,j):end)];
            end
        end
        Thetas = ThetasBalanced(1,:);
    end
    
    
    
    %% New Approach
    if 1
        % syms MT1 MT2 MT3 MT4 MT5 DT1 DT2 DT3 DT4 DT5
        % BalancedVector_S = [MT2 MT3 MT4 MT5];
        % DeltaVector_S    = [DT2 DT3 DT4 DT5];
        % DeltaVector = [10 -50 30 40];
        % ThetaValues = cell(5,1);
        
        Nunknowns = length(DeltaVector_S);
        NonSolved = ones(Nunknowns,1);
        ThetaValues = cell(Nunknowns,1);
        ThetaValues{1} = 45;
        DeltaVector = randi(19,1,Nunknowns)*5 -5 -45;
        
        % V = [ThetaValues{1} MT3 DT2 MT4 DT3 DT4 MT2 MT5 -ThetaValues{1} DT5]
        % V = [ThetaValues{1} DT2 MT3 MT4 DT3 DT4 MT2 MT5 -ThetaValues{1} DT5]
        % V = [ThetaValues{1}    MT3      DT2     MT4      DT4     MT2      MT5      -ThetaValues{1}     DT5     DT3 ];
        % V = [ThetaValues{1}   MT3 MT2 DT3 DT4 DT5 MT5 DT2  MT4      -ThetaValues{1}  ]; % non solvable
        % V = [ThetaValues{1}   MT2 MT3  DT3 DT4 DT5 MT5 DT2  MT4      -ThetaValues{1}  ];
       
        Combi = [BalancedVector_S DeltaVector_S];
        V = [ThetaValues{1} Combi(randperm(Nunknowns*2,Nunknowns*2)) -ThetaValues{1}];
        
        % V = [V(1:10) 90 V(11:end)]
        
        % --- Look for directly computable Theta Values (i.e. if numerical value before DT)
        Computable = true;
        while Computable
            % Update DeltaTheta and Numeric indexes
            MThetaMask      = ismember(V,BalancedVector_S);
            MThetaIndexes   = find(ismember(V,BalancedVector_S));
            DeltaMask       = ismember(V,DeltaVector_S);
            DeltaIndexes    = find(ismember(V,DeltaVector_S));
            NumericMask     = ~(DeltaMask+MThetaMask);
            NumericIndexes  = find(NumericMask);
            
            
            ComputIndex = DeltaIndexes(find(ismember((DeltaIndexes-1),NumericIndexes)));
            
            if isempty(ComputIndex)
                Computable = false;
            else
                for i=1:length(ComputIndex)
                    % compute Thetas and substitute into V
                    [~,j] = ismember(V(ComputIndex(i)),DeltaVector_S);
                    V(ComputIndex(i)) = V(ComputIndex(i)-1) + DeltaVector(j);
                    ThetaValues{j+1}  = eval(V(ComputIndex(i)));
                    
                    % substitute -Thetas and substitute into V as well
                    MTIndex = find(ismember(V,BalancedVector_S(j))); % Minus theta index in V
                    V(MTIndex) = -ThetaValues{j+1};
                    Nunknowns  = Nunknowns -1;
                end
            end
        end
        % ---
        
        
        % ---  Set Equation to solve the other unknowns
        EqCell     = cell(5,1);
        EqIndex    = 1;
        SolvedEq   = [];
        Solved     = false; % becomes true if at least 1 Eq. has been solved
        UnknownVar = [];
        UnknownInd = [];
        for i = 1:length(DeltaIndexes)
            EqNunknown      = 0;
            EqCell{EqIndex} = 0;
            
            
            IndexDT  = find(ismember(DeltaVector_S,V(DeltaIndexes(i))));
            
            IndexMT = find(ismember(V,BalancedVector_S(IndexDT)));
            
            
            UnknownInd = [UnknownInd; IndexDT];
            UnknownVar = [UnknownVar;V(IndexMT)];
            EqCell{EqIndex} =  EqCell{EqIndex} + V(DeltaIndexes(i)-1) + V(DeltaIndexes(i)) == -V(IndexMT);
            
            [Bool,Index] = ismember(V(DeltaIndexes(i)-1),DeltaVector_S);
            if Bool
                EqCell{EqIndex} = subs(EqCell{EqIndex},V(DeltaIndexes(i)-1),-BalancedVector_S(Index));
            end
            
            EqCell{EqIndex} = subs(EqCell{EqIndex},DeltaVector_S,DeltaVector);
            
            
            EqIndex = EqIndex+1;
        end
        
        % --- if no Equations have been solved, need to solve coupled Eq.
        if ~isempty(UnknownInd)
            [UnknownInd,SortIndex] = sort(UnknownInd);
            UnknownVar = UnknownVar(SortIndex);
            [A,b] = equationsToMatrix(EqCell{SortIndex},UnknownVar);
            
            if rank(A)==length(A)
                SolvedVar = A^-1*b;
                
                
                for i = 1:length(UnknownInd)
                    ThetaValues{1+UnknownInd(i)} = -double(SolvedVar(i));
                    V(find(ismember(V,DeltaVector_S(UnknownInd(i)))))    = ThetaValues{UnknownInd(i)+1};
                    V(find(ismember(V,BalancedVector_S(UnknownInd(i))))) = -ThetaValues{UnknownInd(i)+1};
                end
                Nunknowns = 0;
            else
                Solvable = false;
                NInfeasible = NInfeasible +1;
                %             error('Non-Solvable')
            end
        end
        
        
        if Solvable
            % update results above 90 and below -90 accordinlgy
            Thetas = double(V);
            
            
            Index90plus  = find(Thetas>90);
            Index90Minus = find(Thetas<-90);
            
            Index180plus  = find(Thetas>180);
            Index180minus = find(Thetas<-180);
            
            if ~isempty(Index180plus) || ~isempty(Index180minus)
                keyboard
            end
            
            Thetas(Index90plus)  = -180 + Thetas(Index90plus);
            Thetas(Index90Minus) = 180  + Thetas(Index90Minus);
        end
        %         sum(Thetas)
        %         display(Thetas)
    end
    
    
    
    %% check disorientation
    if Solvable
        DetlaAngle      = ComputeDeltaAngle(Thetas); % delta fibre angles between plies
        Ndisorientation(iRun) = length(find(DetlaAngle>45));
    end
    
end

hist(Ndisorientation)
NInfeasible