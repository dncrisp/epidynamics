%{

SDMA - Seizure Data Management and Analysis
-------------------------------------------
%}

classdef SDMAprocess < handle
    
    
    %---------------------
    % Property definitions
    properties
        
        % Paths and filenames
%         patients
%         seizures
%         analyses
        name
        possibledatasources
        datadirectory
        
           patients
        seizures
        analyses
        
        % Paths and filenames
        seizureTrackerFilename
        defaultspath
        logspath
        ieegPwdFile
        savepath
        openpath
        savefilename
        toolbox_path
        isi_path
        figurepath
        um_datapath
        frbg_datapath
        
        
        % IEEG related
        userName
        portalKeys
        
        % Processing settings
        pathdelim
        writetotracker
        savefits
        highpass
        savebuffer
        robustfit
        lfc
        thresh
        todetrend
        todecimate
        decimatefactor
        samplingrate
        amprange
        multiseizure
        
        % Data settings
        normalizefits
        patientKeyNums
        defaults
        hassz
        start
        stop
        isiStart
        saveonly
        isiStop
        source
        annotations
        chanNums_vec
        endISI
        seizureData
        fitequations
        ansznum
        
        
        
    end
    
    
    properties (Dependent)
        numpatients
        numseizures
        numanalyses
        
     
    end
    
    
%     methods (Static)
%         
%             
%         function obj = loadobj(s)
%             
%            if strcmp(class(s),'SDMAprocess')
%                
%                obj = s;
%                
%                loadprocesspatients(obj);
%                
%                
%         end
%          
%         
%     end
    
    
    methods
        
        
        
        %%
        
        %---------------------------------------------
        % Main constructor
        function obj = SDMAprocess(name)
            
            obj.name = name;
            obj.datadirectory = sprintf('.%sData',filesep);
       
            obj.possibledatasources = {'umich', 'frbg', 'ieeg'};
            
            obj.ieegPwdFile = './ieeg-matlab-1.9.5/jar_ieeglogin';
            obj.userName = 'jaredmsc';
            
        end
   
        % Dependent constructors
        
        %---------------------------------------------
        function numpatients = get.numpatients(obj)
            numpatients = length(obj.patients);
        end
        
        
        %---------------------------------------------
        function numseizures = get.numseizures(obj)
            numseizures = length(obj.seizures);
        end
        
        %---------------------------------------------
        function numanalyses = get.numanalyses(obj)
            numanalyses = length(obj.analyses);
        end
        
        
%         %---------------------------------------------
%         function patients = get.patients(obj)
%            
%             patients = {};
%             
%             cd(obj.datadirectory);
%             ddir = dir;
%             
%             for i = 4:length(ddir)
%                 
%                 pname = ddir(i).name;
%                 cd(pname);
%                 
%                 pfile = sprintf('%s.mat', pname);
%                 load(pfile);
%                 %np = {};
%                 eval(sprintf('np = %s', pname));   
%                 patients = [patients, np]; %#ok
%                 
%                 cd ..
%             end
%             
%             cd ..
%             
%             
%         end
%         
%         function patients = set.patients(obj, value)
%             
%             patients = [obj.patients, value];
%             
%         end
        
        
        % Object management
        
%         %---------------------------------------------
%         function delete(obj)
%             
%             prompt = sprintf('Do you wish to save SDMAprocess %s before deletion? (y/n): ', obj.name);
%             str = input(prompt,'s');
%             
%             if strcmp(str,'y')
%                 obj.save;
%             end
%         end
        
        %---------------------------------------------
        function save(obj)
            
            filename = sprintf('%s.mat',obj.name);
            eval(sprintf('%s = obj',obj.name));
            save(filename,obj.name);
        end
        
        
        
        
        
        
        
        
        
        
        
        %---------------------------------------------
        function pt = addnewpatient(obj,patientkey)
            
            pt = SDMApatient(obj, patientkey);
            obj.patients = [obj.patients, pt];
            
            
            
        end
        
        
        
        
        
    end

end


