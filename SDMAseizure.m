classdef SDMAseizure < handle
    
    properties
        
        start
        stop
        
        
        rawdata
        data
        rawfiltered
        rawresampled
        samplerate
        buffer
        
        channelsaved
        totalchannels
        annotationtype
        
        willdetrend
        willdecimate
        willhighpass
        
        notes
        
        decimatefactor
        hpcutoff
        
        analyses
        
        
        
        datadirectory
        datafile
        identifier
        
    end
    
    
    properties (Dependent)
        
        hasanalysis
        szlength
        
    end
    
    properties (Transient)
        
        parentSDMApatient
        parentSDMAprocess
    
    end
        
    
    
    methods
        
        
        
        %-----------------------------------------
        function obj = SDMAseizure(patient, identifier)
            
            
            obj.parentSDMApatient = patient;
            obj.identifier = identifier;
            
            obj.datadirectory = sprintf('%s%s%s', obj.parentSDMApatient.datadirectory, filesep, identifier);
            obj.datafile = sprintf('%s%s%s.mat', obj.datadirectory, filesep, identifier);
            
        end
        
        
        %-----------------------------------------
        function szlength = get.szlength(obj)
            szlength = length(obj.data) / obj.samplerate;
        end
        
        
        %-----------------------------------------
        function hasanalysis = get.hasanalysis(obj)
            
            if size(obj.analyses) > 0
                hasanalysis = 1;
            else
                hasanalysis = 0;
            end
        end
        
        
%         %-----------------------------------------
%         function delete(obj)
%             
%             %             prompt = sprintf('Do you wish to save seizure before deletion? (y/n): ');
%             %             str = input(prompt,'s');
%             %
%             %             if strcmp(str,'y')
%             obj.save;
%             %             end
%             
%         end
%         
        
        %---------------------------------------------
        function save(obj)
            
            if ~exist(obj.datadirectory,'dir')
                mkdir(obj.datadirectory);
            end
            
            eval(sprintf('%s = obj;',obj.identifier));
            save(obj.datafile,obj.identifier);
            
        end
        
        
        
        %---------------------------------------------
        function an = analyzewithtype(obj,type,varargin)
            
            id = sprintf('%s%s', obj.identifier, type);
            an = SDMAanalysis(obj, obj.parentSDMApatient, type, id);
            an.channelused = obj.channelsaved;
            an.fitstart = 1;
            an.fitstop = obj.szlength - 0.1;
            
            if ~isempty(varargin)
                an.globalthresh = varargin{1};
            end
            
            if strcmp(type,'ISI')
                an = an.beginisi;
            elseif strcmp(type,'AMP')
                an = an.beginamp;
            end
            
            obj.analyses = [obj.analyses, an];
            
        end
    end
    
    
    
end