classdef SDMAanalysis < handle
    
    properties
        
        analyzedby
        datetimeanalyzed
        channelused
        fitdatax
        fitdatay
        bestfit
        chosenfit
        threshold
        fitstart
        fitstop
        type
        amprange
        sigmathreshold
        
        identifier
        fits
        globalthresh
        
        
        fitequationsfilename
        
        isweighted
        weights
        weightfunction
        
        
        
    end
    
    
    properties (Transient)
        
        parentSDMAseizure
        parentSDMApatient
        
    end
    
    
    methods
        
        
        
        %-----------------------------------------
        function obj = SDMAanalysis(parentseizure, parentpatient, type, identifier)
            
            obj.parentSDMAseizure = parentseizure;
            obj.parentSDMApatient = parentpatient;
            obj.type = type;
            obj.identifier = identifier;
            obj.sigmathreshold = 0.05;
            obj.fitequationsfilename = './Settings/fitequationslog.mat';
            
        end
        
        
        
        
        
        %---------------------------------------------
        function obj = beginisi(obj)
            
            load(obj.fitequationsfilename);
            
            [obj.fitdatax, obj.fitdatay, obj.threshold, obj.fitstart, obj.fitstop] = ...
                obj.getisioramp(obj.parentSDMAseizure.data, obj.threshold, obj.fitstart, obj.fitstop, ...
                obj.parentSDMAseizure.samplerate, 'isi', obj.amprange);
            
            
            if ~isempty(obj.fitdatax) && ~isempty(obj.fitdatay)
                
                warning('off','all');
                for i=1:length(fitequations)  %#ok<USENS>
                    nf = SDMAfit(obj,obj.fitdatax,obj.fitdatay,fitequations{i});
                    obj.fits = [obj.fits, nf];
                end
                warning('on','all');
                
                rsqadj = [];
                
                for i = 1:length(obj.fits)
                    rsqadj(i) = obj.fits(i).gofstats.adjrsquare; %#ok<AGROW>
                end
                
                [~,bestfitidx] = max(rsqadj);%#ok
                
                %%%%%%%%% For right now only - just to force best fit to log
                bestfitidx = 1;
                
                obj.bestfit = obj.fits(bestfitidx);
                obj.chosenfit = obj.fits(bestfitidx);
                
                statthresh = obj.sigmathreshold;
                
                fprintf('\nBest fitting model is: ''%s''.\n',obj.bestfit.name);
                
                if obj.fits(bestfitidx).residualstats.ks > statthresh
                    fprintf('\nKS Test: Non-normal residuals for best-fitting model!\n');
                elseif obj.fits(bestfitidx).residualstats.ad > statthresh;
                    fprintf('\nAD Test: Non-normal residuals for best-fitting model!\n');
                elseif obj.fits(bestfitidx).residualstats.li > statthresh;
                    fprintf('\nLI Test: Non-normal residuals for best-fitting model!\n');
                elseif obj.fits(bestfitidx).residualstats.jb > statthresh;
                    fprintf('\nJB Test: Non-normal residuals for best-fitting model!\n');
                end
                
                screensize = get(groot,'ScreenSize');
                f = figure('Position',[1 screensize(4) screensize(3)/2 screensize(4)/2]);
                
                hold on
                plot(obj.fits(bestfitidx).fitobject,obj.fitdatax,obj.fitdatay);
                
                
                titlestr = sprintf('%s isi fit- best fitting model: %s- %s, coeffs: %s', ...
                    obj.identifier,obj.bestfit.name,obj.bestfit.formula, obj.bestfit.coefficients);
                
                title(titlestr);
                xlabel('Time until end of seizure (s)');
                ylabel('ISI (s)');
                hold off
                
                figfile = sprintf('./Figures/%s', obj.identifier);
                
                savefigure(f, figfile);
                
                %                 filename = sprintf('%s.txt',obj.parentSDMApatient.patientkey);
                %                 writestring = sprintf('%s; %s; %s; %s; %d; %d; %d\n', ...
                %                     obj.identifier,obj.bestfit.name,obj.bestfit.formula,mat2str(obj.bestfit.coefficients),...
                %                     obj.threshold, obj.fitstart, obj.fitstop);
                %
                %                 fileID = fopen(filename,'at');
                %
                %                 fprintf(fileID, '%s', writestring);
                %                 fclose(fileID);
                
                close all
                
            end
            
            
            
        end
        
        
        
        
        
        
        
        
        
        % Get ISI or AMP function
        function [x,y,thresh_ret,isiStart_ret,isiStop_ret] = getisioramp(obj,seizureData, thresh, isiStart, isiStop, Fs, isioramp, amprange)
            
            addpath('./ISI Fitting');
            
            thresh_ret = thresh;
            isiStart_ret = isiStart;
            isiStop_ret = isiStop;
            
            try
                if strcmp(isioramp,'isi')
                    [x, y] = isi(seizureData, thresh, isiStart, isiStop);
                elseif strcmp(isioramp,'amp')
                    [x, y] = amp2range(seizureData, thresh, isiStart, isiStop, Fs, amprange);
                end
            catch error
                
                fprintf('Error was generated: %s\n',error.identifier);
                warning('Maybe ISI or AMP Parameters are incorrect.');
                
                fprintf('\nSz End with buffer: %.2f; Sz End without buffer: %.2f\n', isiStop_ret, isiStop_ret - 30);
                fprintf('Sz Start with buffer: %.2f; Sz Start without buffer: %.2f\n\n', isiStart_ret, isiStart_ret + 30);
                
                if ~isempty(obj.globalthresh)
                    threshstring = sprintf('Global threshold found: %d. Apply it and defaults above? (y/n) ',obj.globalthresh);
                    answer = input(threshstring,'s');
                    if strcmp(answer,'y')
                        thresh_ret = obj.globalthresh;
                        isiStart_ret = isiStart_ret + 30;
                            isiStop_ret  = isiStop_ret - 30 + 2;
                            if isiStop_ret - 30 > 120
                                isiStart_ret = isiStop_ret - 120;
                            end
                        
                    else
                        thresh_ret = input('Input new threshold: ');
                        
                        ask = input('Fit to default start and stop? (y/n)','s');
                        if strcmp(ask, 'y')
                            isiStart_ret = isiStart_ret + 30;
                            isiStop_ret  = isiStop_ret - 30 + 2;
                            if isiStop_ret - 30 > 120
                                isiStart_ret = isiStop_ret - 120;
                            end
                        else
                            isiStop_ret = input('Input new isiStop: ');
                            
                            if isiStop_ret - 30 > 120
                                isiStart_ret = isiStop_ret - 120;
                            else
                                isiStart_ret = input('Input new isiStart: ');
                            end
                        end
                    end
                else
                    thresh_ret = input('Input new threshold: ');
                    
                    ask = input('Fit to default start and stop? (y/n)','s');
                        if strcmp(ask, 'y')
                            isiStart_ret = isiStart_ret + 30;
                            isiStop_ret  = isiStop_ret - 30 + 2;
                            if isiStop_ret - 30 > 120
                                isiStart_ret = isiStop_ret - 120;
                            end
                        else
                            isiStop_ret = input('Input new isiStop: ');
                            
                            if isiStop_ret - 30 > 120
                                isiStart_ret = isiStop_ret - 120;
                            else
                                isiStart_ret = input('Input new isiStart: ');
                            end
                        end
                end
                
                
                [x,y,thresh_ret,isiStart_ret,isiStop_ret] = ...
                    obj.getisioramp(seizureData, thresh_ret, isiStart_ret, isiStop_ret, Fs, isioramp,amprange);
            end
            
        end
        
        
        
        
        
        
        
        
        
    end
    
end







%---------------------
% Save Figure function
function savefigure(handle,filename)
savefig(handle,filename);
end













