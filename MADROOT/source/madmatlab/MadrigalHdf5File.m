classdef MadrigalHdf5File
    % class MadrigalHdf5File allows the creation of Madrigal Hdf5 files via
    % Matlab.  The general idea of this class is to simply write out all
    % data for the file in a Matlab struct array.  Then the python script
    % createMadrigalHdf5FromMatlab.py is called from this scriptto create all 
    % metadata and alternate array layouts.  This keeps the amount of Matlab 
    % code here at a minimum.  If the output file extension is *.mat, then 
    % only the *.mat is created, and the user must call createMadrigalHdf5FromMatlab.py
    % at a later time themselves.  See file testMadrigalHdf5File.m for example
    % usage.
    %
    % $Id: MadrigalHdf5File.m 6538 2018-07-05 19:26:34Z brideout $
    properties
        filename
        extension % either .hdf5, .h5. .hdf, or .mat
        oneDParms
        independent2DParms
        twoDParms
        arraySplittingParms
        allParms % a combination of oneDParms, independent2DParms, and twoDParms
        recordCount % number of records in file so far
        lastRecord % index of last records added.  Starts with 0. -1 if no records
        data % structure array with fields = stdParms + oneDParms + independent2DParms + twoDParms
        principleInvestigator % the following are strings that used fill up catalog record
        expPurpose            % default for all is empty string
        expMode
        cycleTime
        correlativeExp
        sciRemarks
        instRemarks
        kindatDesc % the following are strings that used fill up header record
        analyst    % default for all is empty string
        comments
        history
        skipArray % bool that determines whether to skip array layout.
    end
    
    properties (Constant)
        stdParms = cellstr(char('recno', 'kinst', 'kindat', 'ut1_unix', 'ut2_unix'));
    end
    
    methods
        function madFile = MadrigalHdf5File(filename, oneDParms, ...
                independent2DParms, twoDParms, arraySplittingParms, ...
                skipArray)
            % Object constructor for MadrigalHdf5File
            % Inputs:
            %   filename - the filename to write to.  Must end *.hdf5, .h5,
            %   .hdf, or .mat.  If .mat, writes a Matlab file that must
            %   later be converted to Madrigal using createMadrigalHdf5FromMatlab.py
            %   oneDParms - a cell array of strings representing 1D parms.  May be
            %     empty.  Example:
            %     cellstr(char('azm', 'elm', 'sn', 'beamid'))
            %   independent2DParms - a cell array of strings representing independent 
            %     2D parms.  May be empty (ie, {}). Examples:
            %       cellstr(char('range'))
            %       cellstr(char('gdlat', 'glon'))
            %       cellstr(char())
            %   twoDParms - a cell array of strings representing dependent 
            %     2D parms.  May be empty (ie, {}). Examples:
            %       cellstr(char('ti', 'dti', 'ne', 'dne'))
            %       cellstr(char())
            %   arraySplittingParms - a cell array of strings representing  
            %     parameters whose values are used to split arrays.  May 
            %     be empty, in which case set to {}. Example:
            %       cellstr(char('beamid'))
            %  skipArray - optional argument.  If set to true, no array
            %     layout created.  If false or not passed in, array layout
            %     created if any 2D variables.
            madFile.filename = filename;
            % verify a valid file extension
            [pathstr,name,ext] = fileparts(filename);
            if (~(strcmp(ext, '.hdf5') | ...
                  strcmp(ext, '.h5') | ...
                  strcmp(ext, '.hdf') | ...
                  strcmp(ext, '.mat')))
              ME = MException('MadrigalHdf5File:invalidExtenstion', ...
                  'Illegal extension %s found', ext);
              throw(ME)
            end
            madFile.extension = ext;
            madFile.oneDParms = cellstr(oneDParms);
            % change all parameters with + to be __plus__
            madFile.oneDParms = strrep(madFile.oneDParms, '+', '__plus__');
            madFile.independent2DParms = cellstr(independent2DParms);
            % change all parameters with + to be __plus__
            madFile.independent2DParms = strrep(madFile.independent2DParms, '+', '__plus__');
            madFile.twoDParms = cellstr(twoDParms);
            % change all parameters with + to be __plus__
            madFile.twoDParms = strrep(madFile.twoDParms, '+', '__plus__');
            madFile.arraySplittingParms = cellstr(arraySplittingParms);
            % change all parameters with + to be __plus__
            madFile.arraySplittingParms = strrep(madFile.arraySplittingParms, '+', '__plus__');
            madFile.recordCount = 0; % no records added yet
            madFile.lastRecord = -1; % index of last record added.
            madFile.data = struct([]);
            
            madFile.allParms = cellstr(char(char(madFile.stdParms), ...
                char(madFile.oneDParms), char(madFile.independent2DParms), ...
                char(madFile.twoDParms)));
            % verify no overlapping
            if (length(madFile.allParms) ~= length(unique(madFile.allParms)))
                ME = MException('MadrigalHdf5File:invalidParameters', ...
                  'Illegal duplicate parameters found in inputs');
                throw(ME)
            end
            
            % set all catalog and header strings to default empty strings
            madFile.principleInvestigator = '';
            madFile.expPurpose = '';
            madFile.expMode = '';
            madFile.cycleTime = '';
            madFile.correlativeExp = '';
            madFile.sciRemarks = '';
            madFile.instRemarks = '';
            madFile.kindatDesc = '';
            madFile.analyst = '';
            madFile.comments = '';
            madFile.history = '';
            
            if (nargin > 5)
                madFile.skipArray = skipArray;
            else
                madFile.skipArray = false;
            end
            
        end % end MadrigalHdf5File constructor
        
        
        function madFile = appendRecord(madFile, ut1_unix, ut2_unix, kindat, ...
                kinst, numRows)
            % appendRecord adds a new record to MadrigalHdf5File.  It
            % returns the record number of the present row (first will be
            % 0)
            %  Inputs:
            %    madFile - the created MadrigalHdf5File object
            %    ut1_unix, ut2_unix - unix start and end time of record in
            %      float seconds since 1970-01-01
            %    kindat - integer kind of data code. See metadata.
            %    kinst - integer instrument code.  See metadata.
            %    numRows - number of rows of 2D data.  If all 1D data, set
            %    to 1
            %  Returns:
            %    the record number of the present row (first will be 0)
            %  Affects:
            %    Updates madFile.recordCount, appends to madFile.data the
            %    number of rows numRows with all data except stdParms set
            %    to NaN.  Use set1D and set2D to populate that record using
            %    recNum as index.
            if (ut1_unix > ut2_unix)
                ME = MException('MadrigalHdf5File:invalidTimes', ...
                  'ut1_unix > ut2_unix - illegal');
              throw(ME)
            end
            thisArr = NaN(numRows, length(madFile.allParms));
            thisTable = array2table(thisArr, 'VariableNames',madFile.allParms);
            
            % set all stdParms
            tmp_arr = ones(numRows,1);
            tmp_arr = ut1_unix;
            thisTable(:,'ut1_unix') = num2cell(tmp_arr);
            tmp_arr = ut2_unix;
            thisTable(:,'ut2_unix') = num2cell(tmp_arr);
            tmp_arr = kindat;
            thisTable(:,'kindat') = num2cell(tmp_arr);
            tmp_arr = kinst;
            thisTable(:,'kinst') = num2cell(tmp_arr);
            tmp_arr = madFile.recordCount;
            thisTable(:,'recno') = num2cell(tmp_arr);
            
            d = size(madFile.data);
            if (d(1) == 0)
                madFile.data = thisTable;
            else
                madFile.data = [madFile.data; thisTable];
            end
            
            madFile.recordCount = 1 + madFile.recordCount; 
            madFile.lastRecord = 1 + madFile.lastRecord;
            
        end % end appendRecord
         
        
        function lastRecord = get.lastRecord(madFile)
            lastRecord = madFile.lastRecord;
        end % lastRecord get function
        
        
        function data = get.data(madFile)
            data = madFile.data;
        end % lastRecord get function
        
        
        function madFile = set1DParm(madFile, parm, value, lastRec)
            % set1DParm sets the values of 1D parm parm to value value for
            % record with lastRecord value lastRec
            % change all parameters with + to be __plus__
            parm = strrep(parm, '+', '__plus__');
            if (~ismember(parm, madFile.oneDParms))
                ME = MException('MadrigalHdf5File:invalidparm', ...
                  'parm %s not in oneDParms', parm );
              throw(ME)
            end
            rows = madFile.data.recno == lastRec;
            tmpArr = ones(length(find(rows)),1);
            tmpArr = value;
            madFile.data(rows,parm) = num2cell(tmpArr);
        end % end set1DParm
        
        
        function madFile = set2DParm(madFile, parm, values, lastRec)
            % set2DParm sets the values of 2D parm parm to value values for
            % record with lastRecord value lastRec
            % change all parameters with + to be __plus__
            parm = strrep(parm, '+', '__plus__');
            if (~ismember(parm, madFile.twoDParms) & ...
                ~ismember(parm, madFile.independent2DParms))
                ME = MException('MadrigalHdf5File:invalidparm', ...
                  'parm %s not in twoDParms or independent2DParms', parm );
              throw(ME)
            end
            rows = madFile.data.recno == lastRec;
            newValues = reshape(values, [length(values) ,1]);
            madFile.data(rows,parm) = num2cell(newValues);
        end % end set2DParm
        
        
        
        function madFile = setCatalog(madFile, principleInvestigator, expPurpose, expMode, ...
                                      cycleTime, correlativeExp, sciRemarks, instRemarks)
            % setCatalog allows setting extra information in the catalog
            % record.  This method is optional.  Even if this method is not
            % called, the catalog record will contain a description of the
            % instrument (kinst code and name) and kind of data brief
            % description, along with a list of description of the
            % parameters in the file, and the first and last times of the
            % measurements.
            %
            % Inputs:
            %
            %    principleInvestigator - Names of responsible Principal Investigator(s) or 
            %       others knowledgeable about the experiment.
            %    expPurpose - Brief description of the experiment purpose
            %    expMode - Further elaboration of meaning of MODEXP; e.g. antenna patterns 
            %       and pulse sequences.
            %    cycleTime - Minutes for one full measurement cycle - must
            %       be numeric
            %    correlativeExp - Correlative experiments (experiments with related data)
            %    sciRemarks - scientific remarks
            %    instRemarks - instrument remarks
            %
            if nargin > 1
                madFile.principleInvestigator = principleInvestigator;
            end
            if nargin > 2
                madFile.expPurpose = expPurpose;
            end
            if nargin > 3
                madFile.expMode = expMode;
            end
            if nargin > 4
                if ~isnumeric(cycleTime)
                    ME = MException('MadrigalHdf5File:invalidArgument', ...
                        'cycleTime not numeric');
                    throw(ME)
                end
                madFile.cycleTime = cycleTime;
            end
            if nargin > 5
                madFile.correlativeExp = correlativeExp;
            end
            if nargin > 6
                madFile.sciRemarks = sciRemarks;
            end
            if nargin > 7
                madFile.instRemarks = instRemarks;
            end
        end % end setCatalog
        
        
        function madFile = setHeader(madFile, kindatDesc, analyst, comments, history)
            % setHeader allows setting extra information in the header
            % record.  This method is optional.
            %
            % Inputs:
            %
            %    kindatDesc - description of how this data was analyzed (the kind of data)
            %    analyst - name of person who analyzed this data
            %    comments - additional comments about data (describe any instrument-specific parameters)
            %    history - a description of the history of the processing of this file
            %
            if nargin > 1
                madFile.kindatDesc = kindatDesc;
            end
            if nargin > 2
                madFile.analyst = analyst;
            end
            if nargin > 3
                madFile.comments = comments;
            end
            if nargin > 4
                madFile.history = history;
            end
        end % end setHeader
        
        
        function write(madFile)
            % write writes out the complete Hdf5 file to madFile.filename,
            % or if the extension is *.mat, only writes out Matlab file
            % without conversion to Hdf5 (which user must do later with 
            % createMadrigalHdf5FromMatlab.py
            filename = madFile.filename;
            oneDParms = madFile.oneDParms;
            independent2DParms = madFile.independent2DParms;
            twoDParms = madFile.twoDParms;
            arraySplittingParms = madFile.arraySplittingParms;
            data = table2array(madFile.data);
            principleInvestigator = madFile.principleInvestigator;
            expPurpose = madFile.expPurpose;
            expMode = madFile.expMode;
            cycleTime = madFile.cycleTime;
            correlativeExp = madFile.correlativeExp;
            sciRemarks = madFile.sciRemarks;
            instRemarks = madFile.instRemarks;
            kindatDesc = madFile.kindatDesc;
            analyst = madFile.analyst;
            comments = madFile.comments;
            history = madFile.history;
            skipArray = madFile.skipArray;
            
            
            if ~strcmp(madFile.extension, '.mat')
                outputMatlabFile = strcat(madFile.filename, '.mat');
                save(outputMatlabFile, 'filename', 'oneDParms', ...
                    'independent2DParms', 'twoDParms', 'arraySplittingParms', ...
                    'data', 'principleInvestigator', 'expPurpose', 'expMode', ...
                    'cycleTime', 'correlativeExp', 'sciRemarks', 'instRemarks', ...
                    'madFile', 'kindatDesc', 'analyst', 'comments', 'history', ...
                    'skipArray');
            else
                outputMatlabFile = madFile.filename;
                save(madFile.filename, 'filename', 'oneDParms', ...
                    'independent2DParms', 'twoDParms', 'arraySplittingParms', ...
                    'data', 'principleInvestigator', 'expPurpose', 'expMode', ...
                    'cycleTime', 'correlativeExp', 'sciRemarks', 'instRemarks', ...
                    'madFile', 'kindatDesc', 'analyst', 'comments', 'history', ...
                    'skipArray');
                return % because not converting to Hdf5 yet
            end 
            
            % create python command to create Hdf5 file from *.mat file if
            % creating a Madrigal Hdf5 file
            madroot = getenv('MADROOT');
            if length(madroot) == 0
                ME = MException('MadrigalHdf5File:missingEnvVariable', ...
                        'MADROOT env variable not set - required');
                    throw(ME)
            end
            execScript = fullfile(madroot, 'bin', 'createMadrigalHdf5FromMatlab.py');
            cmd = sprintf('%s %s', execScript, outputMatlabFile);
            disp(cmd);
            [status,cmdout] = system(cmd);
            disp(cmdout);
        end
            
        
        
    end % end methods

end % classdef


function convertToMadrigal(matFile, madrigalFile)
    % convertToMadrigal converts a matlab mat file to Madrigal Hdf5 file
    %   Inputs:
    %      matFile - existing Matlab .mat file created earlier
    %      madrigalFile - madrigal file to create.  Must end *.hdf5, .h5,
    %         .hdf, or .mat.
    %      
    % create python command to create Hdf5 file from *.mat file
    madroot = getenv('MADROOT');
    if length(madroot) == 0
        ME = MException('MadrigalHdf5File:missingEnvVariable', ...
                'MADROOT env variable not set - required');
            throw(ME)
    end
    execScript = fullfile(madroot, 'bin', 'createMadrigalHdf5FromMatlab.py');
    cmd = sprintf('%s %s', execScript, matFile);
    disp(cmd);
    [status,cmdout] = system(cmd);
    disp(cmdout);
end % end convertToMadrigal