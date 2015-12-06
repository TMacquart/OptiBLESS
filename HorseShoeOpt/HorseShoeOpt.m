% ======================================================================= %
%                        HorseShoe Blending Example                       %
%
% fmincon(@EvaluationFct,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% ======================================================================= %
clear all; close all; clc; format short g; format compact;

global Parameters OptConstraint

% ---
if 1    % Problem definition
    Parameters.Dim{1}      = [18 24];                                % Panel dimension (inch)
    Parameters.Dim{2}      = [20 12];                                % Panel dimension (inch)
    Parameters.PanelDim    = [1 1 2 2 2 2 2 2 1 1 1 1 2 2 2 2 2 2]';
    Parameters.Nx          = [700 375 270 250 210 305 290 600 1100 900 375 400 330 190 300 815  320 300]';      % Paper Problem
    Parameters.Ny          = [400 360 325 200 100 360 195 480 600  400 525 320 330 205 610 1000 180 410]';      % Paper Problem
    Parameters.E1          = 20.5e6;
    Parameters.E2          = 1.31e6;
    Parameters.G12         = 0.62e6;
    Parameters.v12         = 0.32;
    Parameters.ply_t       = 0.0075;
    Parameters.ply_tMax    = 0.0075 *100;
    
    Parameters.connectivity = [
        0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
        1 0 1 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0
        0 1 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0
        0 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0
        0 0 1 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0
        1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0
        1 1 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 0 1 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 0 1 0 0 1 0 0
        0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1
        0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1 0
        0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 1
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0];
    
    
    % --- for gradients
    v21 = Parameters.v12*Parameters.E2/Parameters.E1;
    Q11 = Parameters.E1/(1-Parameters.v12*v21);
    Q22 = Parameters.E2/(1-Parameters.v12*v21);
    Q12 = Parameters.v12*Parameters.E2/(1-Parameters.v12*v21);
    Q66 = Parameters.G12;
    
    % Material invariants
    U1 = 1/8*(3*Q11+3*Q22+2*Q12+4*Q66);
    U2 = 1/2*(Q11-Q22);
    U3 = 1/8*(Q11+Q22-2*Q12-4*Q66);
    U4 = 1/8*(Q11+Q22+6*Q12-4*Q66);
    U5 = 1/8*(Q11+Q22-2*Q12+4*Q66);
    
    % Stiffness matrix
    Parameters.T0 = [U1, U4, 0 ;
        U4, U1, 0 ;
        0,  0,  U5];
    
    Parameters.T1 = [U2, 0,  0;
        0, -U2, 0;
        0,  0,  0];
    
    Parameters.T2 = [0,    0,    U2/2;
        0,    0,    U2/2;
        U2/2, U2/2,  0 ];
    
    Parameters.T3 = [U3, -U3, 0;
        -U3,  U3, 0;
        0,   0,  -U3];
    
    Parameters.T4 = [0,   0,  U3;
        0,   0, -U3;
        U3, -U3, 0 ];
    
    %     D   = t^3/12 * (T0 +W1*T1 +W2*T2 +W3*T3 +W4*T4)
    %     D11 = t^3/12 * (T0(1,1) + W1*T1(1,1) + W3*T3(1,1))
    %     D12 = t^3/12 * (T0(1,2)              + W3*T3(1,2))
    %     D22 = t^3/12 * (T0(2,2) + W1*T1(2,2) + W3*T3(2,2))
    %     D66 = t^3/12 * (T0(3,3)              + W3*T3(3,3))
    
    if 0
        ind = 1;
        m   = 1;
        n   = 1;
        a   = 24;
        b   = 18;
        t_array = 0.1:0.01:0.15; %0.0075:0.0075:0.0075*25;
        w1_array = -1:0.02:1;
        w3_array = -1:0.02:1;
        for tindex = 1:length(t_array)
            for w1index = 1:length(w1_array)
                for w3index = 1:length(w3_array)
                    
                    t(ind)  = t_array(tindex);
                    w1(ind) = w1_array(w1index);
                    w3(ind) = w3_array(w3index);
                    Nx      = -250;
                    Ny      = -100;
                    
                    %                     den   = (m/a)^2 +(n/b)^2*(Ny/Nx);
                    den   = (Nx)*(m/a)^2 +(Ny)*(n/b)^2;
                    acoef = pi^2*(m/a)^4/den;
                    bcoef = pi^2*2*(m/a)^2*(n/b)^2/den;
                    ccoef = pi^2*(n/b)^4/den;
                    dcoef = pi^2*4*(m/a)^2*(n/b)^2/den;
                    
                    alpha = acoef * Parameters.T0(1,1) + bcoef*Parameters.T0(1,2) + ccoef*Parameters.T0(2,2) + dcoef*Parameters.T0(3,3);
                    beta  = acoef * Parameters.T1(1,1)                            + ccoef*Parameters.T1(2,2) ;
                    gamma = acoef * Parameters.T3(1,1) + bcoef*Parameters.T3(1,2) + ccoef*Parameters.T3(2,2) + dcoef*Parameters.T3(3,3);
                    
                    HessianMatrix = -[t(ind)/2*(alpha + w1(ind)*beta + w3(ind)*gamma)    t(ind)^2/4*beta     t(ind)^2/4*gamma
                        t(ind)^2/4*beta                             0               0
                        t(ind)^2/4*gamma                            0               0];
                    
                    if  2*w1(ind)^2*(1-w3(ind)) + w3(ind)^2 - 1 < 0
                        tsave(ind)  = t(ind);
                        w1save(ind) = w1(ind);
                        w3save(ind) = w3(ind);
                        lambda(ind) = t(ind)^3/12*(alpha + w1(ind)*beta + w3(ind)*gamma);
                        
                        %                                                 if find(eig(HessianMatrix)<-0.001)
                        %                                                     fr
                        %                                                 end
                        
                        %                         LP      = [zeros(1,8) w1(ind) 0 w3(ind) 0 ]';
                        %                         [~,~,D] = LP2ABD (Parameters.E1,Parameters.E2,Parameters.v12,Parameters.G12,t(ind),LP,false);
                        %                         den2            = ((m/a)^2*Nx + Ny*(n/b)^2);
                        %                         Constraint(ind) = pi^2* ((D(1,1)*(m/a)^4 + 2*(D(1,2)+2*D(3,3))*(m/a)^2*(n/b)^2 + D(2,2)*(n/b)^4)) / den2;  % should be 1 by default instead of 1.1
                        
                        ind = ind + 1;
                    end
                end
            end
        end
        Plot3k({tsave(:) w1save(:) w3save(:)},                                      ...
            'Plottype','scatter','FontSize',12,                           ...
            'ColorData',lambda(:),'ColorRange',[min(lambda) max(lambda)],'Marker',{'o',2}, ...
            'Labels',{'','thickness','w1','w3','lambda'});
        fr
    end
    
end
% ---

% ---
if 1    % setup
    Tini     = rand(18,1);                                                   % initial thickness
    LPIni    = 2*(rand(18*2,1)-0.5);                                         % initial Lam. Param.
    dvIni    = [Tini;LPIni];                                                 % initial design variable vector
    lb       = [0.01*ones(18,1);   -1*ones(18*2,1)];                         % variables lower bound vector
    ub       = [ ones(18,1);       +1*ones(18*2,1)];                         % variables upper bound vector
    
    
    % --- User Setting
    OptConstraint.Analytical = true;                % use analytical gradient
    
    OptConstraint.Range      = 'global';            % ('global' or 'local') blending constraint must be satisfied between all panels
    OptConstraint.V1Donly    = true;
    OptConstraint.V3Donly    = true;
    OptConstraint.V1DV3D     = true;
    
    OptConstraint.ScoefVkd   = 0.1;                 % (alpha) Spherical Coeff for V1D, V3D individually
    OptConstraint.ScoefV13d  = 0.1;                 % (Beta) Spherical Coeff for V1D coupled V3D
    OptConstraint.FoS        = 1.00;                   % Factor of Safety for the Buckling Factors
    
    OptConstraint.OverRepair = true;
    OptConstraint.MinChange  = false;               % enforce a minimum change between different laminates
    OptConstraint.MinValue   = 0.01;
    Parameters.mMax          = 1;
    Parameters.nMax          = 1;
    % ---
    
    
    OptConstraint.FixedThick = [];
    if OptConstraint.Analytical
        options  = optimset('Algorithm','interior-point','GradConstr','on','GradObj','on','DerivativeCheck','on','DiffMaxChange',1e-5,...
            'DiffMinChange',1e-6,'TolFun',1e-6,'TolX',1e-6,'TolCon',1e-4,'PlotFcns',{@optimplotfval @optimplotconstrviolation},...
            'FinDiffType','central','MaxFunEvals',25000,'FinDiffRelStep',1e-7);
    else
        options  = optimset('Algorithm','interior-point','FinDiffRelStep',1e-4,'MaxFunEvals',50000,'DiffMaxChange',1e-8,'DiffMinChange',1e-9,'TolFun',1e-6,'TolX',1e-6,'TolCon',1e-7,'PlotFcns',{@optimplotfval @optimplotconstrviolation});
    end
end
% ---

% ----
if 1    % no Blending Constraint
    OptConstraint.Type   = '';
    [x0,fval0,~,output0] = fmincon(@EvaluationFct,dvIni,[],[],[],[],lb,ub,@MyConstraints,options);
    temp                 = MyConstraints(x0);
    BucklingFactors      = temp(19:36);
    XContinousOpt        = num2cell([[1:18]' [x0(1:18) x0(19:2:end) x0(20:2:end) x0(1:18,1)*Parameters.ply_tMax/Parameters.ply_t] BucklingFactors]);
    XContinousOpt        = [[{'Panel #'} {'t(inch)'} {'V1D'} {'V3D'} {'Nply'} {'Buckling Constraint'}] ; XContinousOpt];
    MaxViolation.NoConstraints = max(BucklingFactors);
    OptNply.NoConstraints      = sum(cell2mat(XContinousOpt(2:end,5)));
end
% ----

% ----
if 1    % non convex Blending Constraints
    OptConstraint.Type   = 'NonConvexConstraints';
    
    try
        [x2,fval2,~,output2] = fmincon(@EvaluationFct,rand(size(x0)),[],[],[],[],lb,ub,@MyConstraints,options);
        
        [x2,fval2,~,output2] = fmincon(@EvaluationFct,x0,[],[],[],[],lb,ub,@MyConstraints,options);
    catch exception
        warning('Derivative Check Failed.')
        UserInput = input(' Would you like to continue? Y/N:','s');
        if UserInput=='N'
            error('The Code has been stopped!')
        else
            options.DerivativeCheck='off';
            [x2,fval2,~,output2] = fmincon(@EvaluationFct,x0,[],[],[],[],lb,ub,@MyConstraints,options);
            options.DerivativeCheck='on';
        end
    end
    display(output2)
    
    % ---
    if 0    % MultiStart
        problem = createOptimProblem('fmincon','x0',x0,'objective',@EvaluationFct,'lb',lb,'ub',ub,'nonlcon',@MyConstraints,'options',options);
        [xmulti,errormulti] = run(MultiStart(),problem,20);
        
        EvaluationFct(xmulti)
        temp             = MyConstraints(xmulti);
        max(temp(19:36))                            % max buckling factor constraints (ideally negative)
        x2 = xmulti;
    end
    % ---
    
    temp             = MyConstraints(x2);
    BucklingFactors  = temp(19:36);
    XContinousOpt    = num2cell([[1:18]' [x2(1:18) x2(19:2:end) x2(20:2:end) x2(1:18,1)*Parameters.ply_tMax/Parameters.ply_t] BucklingFactors]);
    XContinousOpt    = [[{'Panel #'} {'t(inch)'} {'V1D'} {'V3D'} {'Nply'}  {'Buckling Constraint'}] ; XContinousOpt];
    
    MaxViolation.Blended = max(BucklingFactors);
    OptNply.Blended      = sum(cell2mat(XContinousOpt(2:end,5)));
    
    % ---
    if 0    % Circle plots
        N     = x2(1:18);
        V1D   = x2(19:2:end);
        V3D   = x2(20:2:end);
        array = [N V1D V3D [1:18]'] ;
        array = sortrows(array,-1);
        for i = 1:17
            X_N(i) = (array(i,1)-array(i+1,1))/array(i,1);
        end
        
        figure(3); hold on;
        FfunctionVkD    = 4;
        FfunctionV1DV3D = 5.1443;
        ColorVectors    = [{'s blue'} {'s magenta'} {'s cyan'} {'s red'} {'s green'}  {'s black'} {'s yellow'}];
        ColorVector     = [{'blue'} {'magenta'} {'cyan'} {'red'} {'green'}  {'black'} {'yellow'}];
        
        % ---
        if 0    % plot V1D line + circles ( Delta V1D constraint )
            xlabel('V1D')
            ylabel(' Empty ')
            for i = 1:7
                MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
                
                V1D = array(i,2);
                plot(V1D,0,ColorVectors{i})
                text(V1D,0,strcat(num2str(i)))
                PlotCircle(V1D,0,0,MaxSphereRadius,ColorVector{i},3)
            end
            
            % BETWEEN 1 and 7
            X_NExtra = (array(1,1)-array(3,1))/array(1,1);
            MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_NExtra^2 -36*X_NExtra^3 +60*X_NExtra^4 -48*X_NExtra^5 + 16*X_NExtra^6);
            MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
            
            PlotCircle(array(1,2),0,0,MaxSphereRadius,'-- red',3)
            
            
        end
        % ---
        
        % ---
        if 0    % plot V3D line + circles ( Delta V3D constraint )
            xlabel('Empty')
            ylabel(' V3D ')
            for i = 1:7
                MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
                
                V3D = array(i,3);
                plot(0,V3D,ColorVectors{i})
                text(0,V3D,strcat(num2str(i)))
                PlotCircle(0,V3D,0,MaxSphereRadius,ColorVector{i},3)
            end
        end
        % ---
        
        % ---
        if 1    % plot V1D V3D plot + circles ( Delta V1D,V3D constraint )
            PlotLPSpace('V1A','V3A',3)
            xlabel('V1D')
            ylabel('V3D')
            for i = 1:7
                
                MaxDeltaV13DSquare = FfunctionV1DV3D * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius    = sqrt(OptConstraint.ScoefV13d  * MaxDeltaV13DSquare);
                
                V1D = array(i,2);
                V3D = array(i,3);
                plot(V1D,V3D,ColorVectors{i})
                text(V1D,V3D,strcat(num2str(i)))
                PlotCircle(V1D,V3D,0,MaxSphereRadius,ColorVector{i},3)
            end
        end
        % ---
        
    end
    % ---
end
% ----


% ----
if 1    % Repair Step - Add plies if constraint not respected
    Xrepaired        = cell2mat(XContinousOpt(2:end,2:4));
    Xrepaired(:,1)   = round(Xrepaired(:,1)* Parameters.ply_tMax/Parameters.ply_t) /Parameters.ply_tMax*Parameters.ply_t;
    Xrepaired(:,1)   = 2*round(Xrepaired(:,1)*100/2)/100;                 % made symmetric
    temp             = MyConstraints(Array2dv(Xrepaired));
    BucklingFactors  = temp(19:36);
    
    if OptConstraint.OverRepair
        Xrepaired((BucklingFactors>0),1) = Xrepaired((BucklingFactors>0),1) + 2*Parameters.ply_t/Parameters.ply_tMax;
        temp             = MyConstraints(Array2dv(Xrepaired));
        BucklingFactors  = temp(19:36);
    end
    
    Xrepaired = num2cell([[1:18]' [Xrepaired(:,1) Xrepaired(:,2) Xrepaired(:,3) Xrepaired(:,1)*Parameters.ply_tMax/Parameters.ply_t] BucklingFactors]);
    Xrepaired = [[{'Panel #'} {'t(inch)'} {'V1D'} {'V3D'} {'Nply'} {'Buckling Constraint'}] ; Xrepaired];
    
    MaxViolation.Repaired = max(BucklingFactors);
    OptNply.Repaired      = sum(cell2mat(Xrepaired(2:end,5)));
end
% ----

% ----
if 1    % Run Opti after repaired fixed thickness
    OptConstraint.FixedThick = cell2mat(Xrepaired(2:end,2));
    dvRepaired               = Array2dv(cell2mat(Xrepaired(2:end,2:4)));
    LP_Only                  = dvRepaired(19:end);
    
    try
        [xt,fvalt,~,outputt] = fmincon(@EvaluationBucklingFactor,LP_Only,[],[],[],[],lb(19:end),ub(19:end),@AllConstraint_LPOnly,options);
    catch exception
        warning('Derivative Check Failed.')
        UserInput = input(' Would you like to continue? Y/N:','s');
        if UserInput=='N'
            error('The Code has been stopped!')
        else
            options.DerivativeCheck='off';
            [xt,fvalt,~,outputt] = fmincon(@EvaluationBucklingFactor,LP_Only,[],[],[],[],lb(19:end),ub(19:end),@AllConstraint_LPOnly,options);
            options.DerivativeCheck='on';
        end
    end
    display(outputt)
    
    % ---
    if 0    % MultiStart
        problem = createOptimProblem('fmincon','x0',LP_Only,'objective',@EvaluationFct_LPOnly,'lb',lb(19:end),'ub',ub(19:end),'nonlcon',@MyConstraints_LPOnly,'options',options);
        [xmulti,errormulti] = run(MultiStart(),problem,4);
        
        EvaluationFct_LPOnly(xmulti)
        temp             = MyConstraints_LPOnly(xmulti);
        max(temp(19:36))                            % max buckling factor constraints (ideally negative)
        xt = xmulti;
    end
    % ---
    
    xt               = [OptConstraint.FixedThick; xt];
    temp             = MyConstraints(xt);
    BucklingFactors  = temp(19:36);
    XContinousFixedT = num2cell([[1:18]' [xt(1:18) xt(19:2:end) xt(20:2:end) xt(1:18,1)*Parameters.ply_tMax/Parameters.ply_t] BucklingFactors]);
    XContinousFixedT = [[{'Panel #'} {'t(inch)'} {'V1D'} {'V3D'} {'Nply'}  {'Buckling Factor'}] ; XContinousFixedT];
    
    MaxViolation.RepairedBlended = max(BucklingFactors);
    OptNply.RepairedBlended      = sum(cell2mat(XContinousFixedT(2:end,5)));
    
    % ---
    if 0    % Circle plots
        N     = xt(1:18);
        V1D   = xt(19:2:end);
        V3D   = xt(20:2:end);
        array = [N V1D V3D [1:18]'] ;
        array = sortrows(array,-1);
        for i = 1:17
            X_N(i) = (array(i,1)-array(i+1,1))/array(i,1);
        end
        
        figure(3); hold on;
        FfunctionVkD    = 4;
        FfunctionV1DV3D = 5.1443;
        ColorVectors    = [{'s blue'} {'s magenta'} {'s cyan'} {'s red'} {'s green'}  {'s black'} {'s yellow'}];
        ColorVector     = [{'blue'} {'magenta'} {'cyan'} {'red'} {'green'}  {'black'} {'yellow'}];
        
        % ---
        if 0    % plot V1D line + circles ( Delta V1D constraint )
            xlabel('V1D')
            ylabel(' Empty ')
            for i = 1:7
                MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
                
                V1D = array(i,2);
                plot(V1D,0,ColorVectors{i})
                text(V1D,0,strcat(num2str(i)))
                PlotCircle(V1D,0,0,MaxSphereRadius,ColorVector{i},3)
            end
            
            % BETWEEN 1 and 7
            X_NExtra = (array(1,1)-array(3,1))/array(1,1);
            MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_NExtra^2 -36*X_NExtra^3 +60*X_NExtra^4 -48*X_NExtra^5 + 16*X_NExtra^6);
            MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
            
            PlotCircle(array(1,2),0,0,MaxSphereRadius,'-- red',3)
            
            
        end
        % ---
        
        % ---
        if 0    % plot V3D line + circles ( Delta V3D constraint )
            xlabel('Empty')
            ylabel(' V3D ')
            for i = 1:7
                MaxDeltaVkDSquare = FfunctionVkD * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius   = sqrt(OptConstraint.ScoefVkd * MaxDeltaVkDSquare);
                
                V3D = array(i,3);
                plot(0,V3D,ColorVectors{i})
                text(0,V3D,strcat(num2str(i)))
                PlotCircle(0,V3D,0,MaxSphereRadius,ColorVector{i},3)
            end
        end
        % ---
        
        % ---
        if 0    % plot V1D V3D plot + circles ( Delta V1D,V3D constraint )
            PlotLPSpace('V1A','V3A',3)
            xlabel('V1D')
            ylabel('V3D')
            for i = 1:7
                
                MaxDeltaV13DSquare = FfunctionV1DV3D * 4 * (9*X_N(i)^2 -36*X_N(i)^3 +60*X_N(i)^4 -48*X_N(i)^5 + 16*X_N(i)^6);
                MaxSphereRadius    = sqrt(OptConstraint.ScoefV13d  * MaxDeltaV13DSquare);
                
                V1D = array(i,2);
                V3D = array(i,3);
                plot(V1D,V3D,ColorVectors{i})
                text(V1D,V3D,strcat(num2str(i)))
                PlotCircle(V1D,V3D,0,MaxSphereRadius,ColorVector{i},3)
            end
        end
        % ---
        
    end
    % ---
    
end
% ----

% --- % Continous Results
if ~isempty(find(BucklingFactors>0,1)),
    warning('Buckling constraints are not satisfied by the Optimal Continuous Design');
    UserInput = input(' Would you like to continue? Y/N:','s');
    if UserInput=='N'
        error('The Code has been stopped!')
    end
end
display('Buckling constraints are satisfied')

% DesignVarToUseforGA = XContinousOpt;
DesignVarToUseforGA   = XContinousFixedT;

dvContinuous.Nply            = round(cell2mat(DesignVarToUseforGA(2:end,2))*Parameters.ply_tMax/Parameters.ply_t);
dvContinuous.V1D             = cell2mat(DesignVarToUseforGA(2:end,3));
dvContinuous.V3D             = cell2mat(DesignVarToUseforGA(2:end,4));
dvContinuous.BucklingFactors = cell2mat(DesignVarToUseforGA(2:end,6));

display(DesignVarToUseforGA)
display(MaxViolation)
display(OptNply)
% ---

OptConstraint.FoS    = 1.0;
cd GA_MatchingLP
run GA_LP2SS_MatchingLP


% =============================================================================
% checking final results
cd ..

NGuidePlies = max(dvContinuous.Nply)/2;
[laminates] = FromIndividual2laminate (xOpt{minIndex},NGuidePlies);
LPArray     = zeros(12,18);
for i = 1 : 18
    ply_candidate = [laminates{i} fliplr(laminates{i})]';
    t(i)          = length(laminates{i})*2*Parameters.ply_t;
    LPArray(:,i)  = SS2LP(Parameters.ply_t,ply_candidate);
end

V1D   = LPArray(9,:);
V3D   = LPArray(11,:);
tNorm = t/Parameters.ply_tMax;
dv = tNorm';
for i = 1 : 18
    dv = [dv; [V1D(i)  V3D(i)]'];
end

warning('remove Safety factor for final constraint Eval')

FinalConstraints     = MyConstraints(dv);

LowestLPfeasibility  = max(FinalConstraints(1:18));
LowestBucklingFactor = max(FinalConstraints(19:36));
LowestBlending       = max(FinalConstraints(37:end));

display(LowestLPfeasibility)
display(LowestBucklingFactor)
display(LowestBlending)

sum(FinalConstraints(FinalConstraints>0))
fval % sum of positive constraint violation


Discrete.OptimalAngles = nan*ones(18,max(cellfun(@length,laminates)));
for i = 1 : 18
    Discrete.OptimalAngles(i,1:length(laminates{i})) = laminates{i};
    
    Ply_candidate = [laminates{i} fliplr(laminates{i})]';
    LP(:,i) = SS2LP( Parameters.ply_t,Ply_candidate);
    
    Discrete.V1D(i) = LP(9,i);
    Discrete.V3D(i) = LP(11,i);
end

Discrete.BucklingFactors = FinalConstraints(19:36);

% --- Discrete Results
DiscreteResults = num2cell([[1:18]' dvContinuous.Nply Discrete.V1D' Discrete.V3D' Discrete.BucklingFactors]);
DiscreteResults = [[{'Panel #'} {'Nplies'} {'V1D'} {'V3D'} {'Buckling Constraint'}]; DiscreteResults];
display(DiscreteResults)

OptNply.MaxViolation = max(cell2mat(DiscreteResults(2:end,4)));
OptNply.discrete     = sum(cell2mat(DiscreteResults(2:end,1)));

if NGA>1
    display(GAPerfIndex)
    display(AvgMaxBucklingConstraints)
    display(BestDesignBucklingConstraints)
end


% ========================================================================
% ---
if 0    % Benchmark with lit. results
    
    % Continuous Optimal
    %     OptiT  = [0.238 0.203 0.150 0.132 0.112 0.155 0.134 0.177 0.276 0.254 0.214 0.204 0.153 0.130 0.179 0.217 0.133 0.160]';
    %     OptiV1 = [0.383 0.383 -0.832 -0.741 -0.581 -0.832 -0.692 -0.741 0.383 0.383 0.383 0.383 -0.799 -0.818 -0.832 -0.832 -0.635 -0.832]';
    %     OptiV3 = [-0.707 -0.707 0.383 0.105 -0.319 0.383 -0.039 0.105 -0.707 -0.707 -0.707 -0.707 0.283 0.340 0.383 0.383 -0.182 0.383]';
    
    % Discrete Optimal
    OptiT  = [0.240 0.210 0.150 0.135 0.120 0.150 0.135 0.180 0.285 0.255 0.210 0.210 0.150 0.135 0.180 0.225 0.135 0.165 ]';
    OptiV1 = [0.387 0.383 -0.825 -0.621 -0.5 -0.846 -0.568 -0.609 0.384 0.383 0.366 0.383 -0.668 -0.676 -0.834 -0.833 -0.517 -0.837]';
    OptiV3 = [-0.604 -0.617 0.425 -0.145 -0.5 0.444 -0.310 -0.201 -0.616 -0.612 -0.634 -0.617 -0.040 -0.017 0.412 0.409 -0.450 0.446]';
    
    OptiT_Normalised = OptiT / Parameters.ply_tMax;
    
    OptiSamuel = [OptiT_Normalised];
    for i = 1 : 18
        OptiSamuel = [OptiSamuel; OptiV1(i); OptiV3(i)];
    end
    
    EvaluationFct(OptiSamuel)
    MyConstraints(OptiSamuel)
end
% ---


% ---
if 0    % Adams results validation (460 ply, fully blended)
    OptimalAngles =[
        nan nan +45 -45  30 +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  nan nan  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  nan nan nan -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  nan nan  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan +30 nan  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        +45 -45 +45 -45  30 +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan -45 +45 -45  30 +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan  30 +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  nan nan  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan +30 nan  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        +45 -45 +45 -45  30 +30 -30  45  +45 -45  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  nan nan  60 -60 60 -60 60 -60 -75 60 75
        nan nan nan nan nan nan nan nan  +45 -45  60 -60 60 -60 60 -60 -75 60 75];
    
    for i = 1 : 18
        laminates{i}  = OptimalAngles(i,~isnan(OptimalAngles(i,:)));
        Ply_candidate = [laminates{i} fliplr(laminates{i})]';
        
        LP(:,i) = SS2LP(Parameters.ply_t,Ply_candidate);
        Nply(i) = length(laminates{i})*2;
        V1D(i)  = LP(9,i);
        V3D(i)  = LP(11,i);
    end
    
    tnorm = [Nply*Parameters.ply_t/Parameters.ply_tMax];
    
    dvArray= [tnorm' V1D' V3D']
    OptiAdams = Array2dv(dvArray)
    
    EvaluationFct(OptiAdams)
    MyConstraints(OptiAdams)
    
end
% ---

