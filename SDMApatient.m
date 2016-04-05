
classdef SDMApatient < handle
    
    properties
        
        notes
        gender
        handedness
        ethnicity
        ageatonset
        seizurehistory
        datasource
        patientkey
        seizures
        analyses
        pathology
        
        globalchannel
        
%         parentSDMAprocess
        
        datadirectory
        datafile
        
        analysistype
        analysiscoeffs
        analysisr2adj
        
        startingsz
        maxsztoanalyze
        maxsztypetodownload
        willlimitnumszdownloads
        willlimitnumanalyses
        
        savebuffer
        globalthreshold
        
    end
    
    
    properties (Dependent)
        
        hasseizures
        numseizures
        hasanalyses
        numanalyses
        
        
        
    end
    
    
    properties (Transient)    
        parentSDMAprocess
    end
    
    
    methods
        
        
        
        
        %-----------------------------------------
        function obj = SDMApatient(parent, patientkey)
            
            obj.patientkey = patientkey;
            obj.parentSDMAprocess = parent;
            obj.seizures = {};
            obj.analyses = {};
            obj.analysistype = {};
            obj.analysiscoeffs = {};
            
            obj.datadirectory = sprintf('%s%s%s', obj.parentSDMAprocess.datadirectory, filesep, obj.patientkey);
            obj.datafile = sprintf('%s%s%s.mat', obj.datadirectory, filesep, obj.patientkey);
            
            obj.savebuffer = 30; % 30 second buffer on beginning and end of seizures
            
            obj.willlimitnumszdownloads = 0;
            obj.willlimitnumanalyses = 0;
            
            obj.maxsztoanalyze = 50;
            obj.maxsztypetodownload = 50;
            
        end
        
        
        
        
        
        
        
        
        %-----------------------------------------
        function hasseizures = get.hasseizures(obj)
            
            if size(obj.seizures) > 0
                hasseizures = 1;
            else
                hasseizures = 0;
            end
        end
        
        
        %-----------------------------------------
        function numseizures = get.numseizures(obj)
            
            numseizures = length(obj.seizures);
            
        end
        
        
        
        %-----------------------------------------
        function hasanalyses = get.hasanalyses(obj)
            
            if size(obj.analyses) > 0
                hasanalyses = 1;
            else
                hasanalyses = 0;
            end
        end
        
        
        %-----------------------------------------
        function numanalyses = get.numanalyses(obj)
            
            numanalyses = length(obj.analyses);
            
        end
        
        
        
        
        %---------------------------------------------
        function sz = addnewseizure(obj,varargin)
            
            
            id = sprintf('%ssz%d', obj.patientkey, obj.numseizures+1);
            
            sz = SDMAseizure(obj,id);
            sz.parentSDMApatient = obj;
            sz.parentSDMAprocess = obj.parentSDMAprocess;
            
            obj.seizures = [obj.seizures, sz];
            obj.parentSDMAprocess.seizures = ...
                [obj.parentSDMAprocess.seizures, sz];
            
            
            
        end
        
%         %---------------------------------------------
%         function delete(obj)
%             
%             prompt = sprintf('Do you wish to save patient %s before deletion? (y/n): ', obj.patientkey);
%             str = input(prompt,'s');
%             
%             if strcmp(str,'y')
%                 obj.save;
%             end
%             
%         end
%         
        
        %---------------------------------------------
        function save(obj)
            
            if ~exist(obj.datadirectory,'dir')
                mkdir(obj.datadirectory);
            end
            
            eval(sprintf('%s = obj',obj.patientkey));
            save(obj.datafile,obj.patientkey,'-append');
            
        end
        
        
        %---------------------------------------------
        function downloadszfromieeg(obj,varargin)
            
            isi_path = './ieeg-matlab-1.9.5/IEEGToolbox';
            
            chantouse = obj.globalchannel;
            
            addpath(isi_path);
            session = IEEGSession(obj.patientkey, obj.parentSDMAprocess.userName, obj.parentSDMAprocess.ieegPwdFile);
            
            % Store dataset info
            dataset         = session.data;
            numchannels     = length(dataset.channels);
            samplerate      = dataset.sampleRate;
            filter          = dataset.filter;
            resample        = dataset.resample;
            
            if ~isempty(varargin)
                
                % TODO: add in error checking for varargin
                % also ability to download multiple seizures
                fprintf('\nGetting data for seizure.\n');
                
                chantouse = varargin{1};
                start = varargin{2};
                stop = varargin{3};
        
                act_start = (start - obj.savebuffer) * samplerate;
                act_stop  = (stop + obj.savebuffer) * samplerate;
                data = dataset.getvalues(act_start:act_stop, chantouse);
                
                data = processdata(data);
                newsz = createszobject(data, chantouse, act_start, act_stop);
                
            else
            
            % Initialize
            antns = {};
            antnnames = {};
            annLayeridx_touse = [];
            numAnnLayers = length(dataset.annLayer);
            
            % If more than 1 annotation layer,
            %       look for layers with the word 'Seizure' in their name
            if numAnnLayers > 1
                for i = 1:numAnnLayers
                    annLayerName = dataset.annLayer(i).name;
                    if ~isempty(strfind(annLayerName,'Seizure'))
                        % When found add layer to list of layers to use in
                        %   getting annotations
                        annLayeridx_touse = [annLayeridx_touse i]; %#ok
                    end
                end
                % If only one layer found, seizure annotation events may still be in it
            elseif numAnnLayers == 1
                annLayeridx_touse = 1;
            end
            
            % Loop through all 'Seizure' containing layers and get seizure events
            for i=1:length(annLayeridx_touse)
                
                antns{i} = getEvents(dataset.annLayer(annLayeridx_touse(i)), 0);  %#ok<AGROW>
                antnnames{i} = dataset.annLayer(annLayeridx_touse(i)).name;  %#ok<AGROW>
                
                samplebuffer = obj.savebuffer * 1e6;
                startingtimes = [antns{i}.start] - samplebuffer;
                
                widths = ([antns{i}.stop] - [antns{i}.start]) + (2 * samplebuffer);
                
                
                % Loop through each seizure in each annotation layer and get data, then create new seizure object with that data
                for j=1:length(antns{i})
                    
                    fprintf('\nGetting data for seizure %d/%d for annlayer %d.\n',j,length(antns{i}),i);
                    
                    if (widths(j) - (2 * samplebuffer)) > 2 * 1e6
                        data = dataset.getvalues( startingtimes(j), widths(j), chantouse);
                        
                        data = processdata(data);
                        
                        start = startingtimes(j);
                        stop = startingtimes(j) + widths(j);
                        
                        newsz = createszobject(data, chantouse, start, stop);

%                         newsz.save;
                        
                    end
                    
                    if j > obj.maxsztypetodownload && obj.willlimitnumszdownloads == 1
                        break
                    end
                    
                end
                
            end
            end
           
            function processed_data = processdata(data)
                
                processed_data = data;
                
               nans = isnan(processed_data) == 1;
                        processed_data(nans) = 0;
                        
                        processed_data = filterdata(processed_data, samplerate);
                        processed_data = removeDC(processed_data);
                        
                        processed_data = addtseries(processed_data, samplerate); 
            end
            
            function szobject = createszobject(data, channel, start, stop)
                 szobject = addnewseizure(obj);
                        
                        szobject.rawdata = data;
                        szobject.data = data;
                        
                        szobject.samplerate = samplerate;
                        szobject.rawfiltered = filter;
                        szobject.rawresampled = resample;
                        
                        szobject.totalchannels = numchannels;
                        szobject.channelsaved = channel;
                        
                        szobject.start = start;
                        szobject.stop = stop;
                        
                        szobject.buffer = obj.savebuffer;
                        
                        if exist('antnnames','var')
                            szobject.annotationtype = antnnames{i};
                        else
                            szobject.annotationtype = '';
                        end
            end
        end
        
        
        
        
        %---------------------------------------------
        function analyzeszwithisi(obj,varargin)
            
            globalthresh = [];%#ok
            
            if ~isempty(varargin)
                if length(varargin) == 1
                    obj.startingsz = varargin{1};
                elseif length(varargin) == 2
                    obj.startingsz = varargin{1};
                    globalthresh = varargin{2}; %#ok
                end
            else
                obj.startingsz = 1;
            end
            
            for i=obj.startingsz:length(obj.seizures)
            
                fprintf('Analyzing seizure %d/%d\n', i, obj.numseizures);
                
                if ~isempty(globalthresh)
                    analysis = obj.seizures(i).analyzewithtype('ISI',globalthresh);
                else
                    analysis = obj.seizures(i).analyzewithtype('ISI');
                end
                obj.analyses = [obj.analyses, analysis];
                at = sprintf('%s: %s', analysis.bestfit.name, analysis.bestfit.formula);
    
                obj.analysistype = [obj.analysistype, at];
                obj.analysiscoeffs = [obj.analysiscoeffs, analysis.bestfit.coefficients];
                obj.analysisr2adj = [obj.analysisr2adj, analysis.bestfit.r2adj];
                
%                 obj.save;
                
                if i > obj.maxsztoanalyze && obj.willlimitnumanalyses;
                    break
                end
            
            end
            
        end
        
        
        
        
        
        
        
        
        
    end
    
    
    
end



%---------------------------------
% Filter (highpass) data function
function data = filterdata(data, sampleRate)
[b, a] = butter(1,1/(sampleRate/2),'high');
data = filtfilt(b,a,data);
end





%---------------------
% Remove DC function
function data = removeDC(data)
data = data - mean(data);
end






%---------------------
% Add time series function
function data = addtseries(data, sampleRate)
Fs = sampleRate;
t = (0:1/Fs:(length(data)/Fs - 1/Fs))';
data = [t data];
end
