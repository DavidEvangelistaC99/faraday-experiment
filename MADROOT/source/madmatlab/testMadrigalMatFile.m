% test/example script to exercise MadrigalHdf5File
%
% $Id: testMadrigalMatFile.m 6477 2018-04-04 17:43:49Z brideout $
filename = '/Users/brideout/Documents/workspace/mad3_trunk/source/madmatlab/example.mat';
oneDParms = cellstr(char('azm', 'elm', 'sn', 'beamid'));
independent2DParms = cellstr(char('range'));
twoDParms = cellstr(char('ti', 'dti'));
arraySplittingParms = cellstr(char('beamid'));

% some hard-coded fake data
ut1_unix = 1.0E9;
ut2_unix = 1.0E9 + 100;
kindat = 3410;
kinst= 30;
numRows = 5;
azm = 100.0;
elm = 45.0;
sn = 0.5;
range = [100.0, 150.0, 200.0, 250.0, 300.0];
ti = [1000.0, 1100.0, 1200.0, 1300.0, 1400.0];
dti = [100.0, 150.0, 200.0, 250.0, 300.0];
beamids = [1,2,1,2,1,2,1,2,1,2];

% new file must not already exist
delete(filename);


madFile = MadrigalHdf5File(filename, oneDParms, ...
          independent2DParms, twoDParms, arraySplittingParms);
      
for rec = 1:10
    madFile = madFile.appendRecord(ut1_unix, ut2_unix, kindat, ...
        kinst, numRows);

    % set 1D values
    madFile = madFile.set1DParm('azm', azm, madFile.lastRecord);
    madFile = madFile.set1DParm('elm', elm, madFile.lastRecord);
    madFile = madFile.set1DParm('sn', sn, madFile.lastRecord);
    madFile = madFile.set1DParm('beamid', beamids(rec), madFile.lastRecord);
    
    % set 2D and independent variables
    madFile = madFile.set2DParm('range', range, madFile.lastRecord);
    madFile = madFile.set2DParm('ti', ti, madFile.lastRecord);
    madFile = madFile.set2DParm('dti', dti, madFile.lastRecord);
    
    % time must always increase in Madrigal 3
    ut1_unix = ut1_unix + 100;
    ut2_unix = ut1_unix + 200;
    
end

% add catalog and header info
principleInvestigator = 'Bill Rideout';
expPurpose = 'Measure the ionosphere';
expMode = 'Vector measurements';
cycleTime = 20.0;
correlativeExp = ''; 
sciRemarks = 'Big solar storm'; 
instRemarks = 'Working well!!!';
madFile = madFile.setCatalog(principleInvestigator, expPurpose, expMode, ...
                            cycleTime, correlativeExp, sciRemarks, instRemarks);

kindatDesc = 'Regular processing';
analyst = 'Phil Erickson'; 
comments = 'Include unusual parameter description here';
history = 'Reprocessed four times';
madFile = madFile.setHeader(kindatDesc, analyst, comments, history);

write(madFile);
    
