classdef SDMAfit < handle
    
    properties
        
        fitobject
        gofstats
        output
        residualstats
        name
        formula
        coefficients
        r2adj
        
        
    end
    
    
    
    properties (Transient)
        parentSDMAanalysis
    end
    
    
    methods
        
        
        %---------------------------------------------
        function obj = SDMAfit(parent,x,y,equationstruct)    
            
            obj.parentSDMAanalysis = parent;
            
               fo = fitoptions('Method','NonlinearLeastSquares', ...
                'MaxFunEvals', 2000, 'MaxIter', 2000);
            
            numcoeffs = equationstruct.numcoeffs;
            startpoint = ones(numcoeffs,1);
            fo.StartPoint = startpoint;
            ft  = fittype(equationstruct.formula,'options',fo);
            
            [eqfit,gof,out] = fit(x,y,ft);
            
            fitstat = stats(out);
            
            obj.fitobject = eqfit;
            obj.gofstats = gof;
            obj.output = out;
            obj.residualstats = fitstat;
            obj.name = equationstruct.name;
            obj.formula = equationstruct.formula;
            obj.coefficients = coeffvalues(eqfit);
            obj.r2adj = obj.gofstats.adjrsquare;
            
        end
        
        
        
        
        

        
    end
    
    
    
end






%---------------------
% Stats function
function statresults = stats(output)

% Statistical tests for model fit
if ~isempty(output)
    ksres = (output.residuals - 0) / var(output.residuals); % normalize residuals...?
    [~,ksP] = kstest(ksres);            % Kolgorov-smirnov test
    [~,adP] = adtest(output.residuals); % Anderson-Darling test
    [~,liP] = lillietest(output.residuals); % Lillefors test
    [~,jbP] = jbtest(output.residuals);     % Jarque-Bare test
end

statresults.ks = ksP;
statresults.ad = adP;
statresults.li = liP;
statresults.jb = jbP;
end
