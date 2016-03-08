
clear all
close all
clc

syms T1 T2 T3 T4 T5 DT1 DT2 DT3 DT4 DT5

Vsym = [T1 T2 T3 T4 T5]
Dsym = [DT1 DT2 DT3 DT4 DT5]
Theta1 = 45;
DTheta = [0 10 -50 30 40];
ThetaCell = cell(5,1);


Nunknowns = length(

% V      = [Theta1 Dsym(2) T3 T4 Dsym(3) Dsym(4) T2 T5 -Theta1 Dsym(5)] %  work for this one
V      = [Theta1 T3 Dsym(2) T4 Dsym(4) T2 T5 -Theta1 Dsym(5) Dsym(3) ]

SolvedEq = [];
EqIndex  = 1;
EqCell   = cell(5,1);
for i=length(V):-1:2
    
    [Boolean,Index] = ismember(V(i),Dsym);
    
    if Boolean
        % the i index of V is unknown
        % keyboard
        
        Boolean1 = ismember(V(i-1),Dsym);
        Boolean2 = ismember(V(i-1),Vsym);
        if ~Boolean1 && ~Boolean2
            V(i)=V(i-1)+DTheta(Index);
            V(ismember(V,Vsym(Index))) = -V(i);
            
        else
            % rewind back from DT until you reach T
            IndexT = find(ismember(V,Vsym(Index)),1);
            EqCell{EqIndex} = 0;
            Nunknown = 1;
            for j = i:-1:IndexT
                if j == IndexT
                    EqCell{EqIndex} =  EqCell{EqIndex} - V(j);
                else
                    
                    [Boolean1,Index1] = ismember(V(j),Dsym);
                    [Boolean2,Index2] = ismember(V(j),Vsym);
                    if Boolean1
                        %                         if known value, replace
                        EqCell{EqIndex} =  EqCell{EqIndex} + DTheta(Index1);
                        %                         EqCell{EqIndex} =  EqCell{EqIndex} +  V(j);
                    elseif Boolean2
                        % leave unknown
                        EqCell{EqIndex} =  EqCell{EqIndex} - V(j);
                        Nunknown = Nunknown +1;
                    else
                        % replace by numeric
                        EqCell{EqIndex} =  EqCell{EqIndex} + V(j);
                    end
                end
            end
            EqCell{EqIndex} = ( EqCell{EqIndex} == Vsym(Index));
            
           
            if Nunknown == 1
                % solve directly
                Sol = solve(EqCell{EqIndex});
                V   = subs(V,V(j),-Sol);
                SolvedEq = [SolvedEq ; EqIndex];
                %
                ThetaCell{IndexT}= Sol;
            end
            
            EqIndex = EqIndex+1;
        end
    end
    
end

%  keyboard

 UnsolvedIndex = find(~cellfun(@isempty,EqCell));
 UnsolvedIndex(SolvedEq) = [];
 Nunsolved     = length(UnsolvedIndex);
 
 
 
if Nunsolved>0
    keyboard
    S = solve(EqCell{UnsolvedIndex(1)});
end
% V = subs(V,T3,-S);
% V = subs(V,T4,-S.T4);


RemIndex = find(ismember(V,Dsym));
for i = 1:length(RemIndex)
    DeltaIndex = find(Dsym==V(RemIndex(i)));
    V(RemIndex(i)) = V(RemIndex(i)-1) + DTheta(DeltaIndex);
end

display(V)


% individual checks are (might need to repeart until no more unknown)
% 1 - if the value before DeltaTheta is known you can calculate theta
% 2 - if theta values is known, replace the -thetas


