% analysis base
% spikes
global DIRS
DIRS.analysis = '/Users/Aniruddh/Documents/Projects/ResidualDynamics/toolboxes/Dataloader/';
DIRS.user = 'ag';

project = ProjectCreate('trial')
bool = ProjectDelete('trial');

project = ProjectCreate('trial')
project = ProjectCreate('trial2','description')
project = ProjectCreate('trial3','','testmodel')
project = ProjectCreate('trial3','','testmodel',{'par1','par2'})
bool = ProjectDelete('trial3');

% project = ProjectCreate('trial','','testmodel',{'par1','par2'})
% unit = UnitLoad(DIRS.spikes,'catz006',1,1);
% analysis = AnalysisAdd('trial', unit);
% analysis = AnalysisAdd('trial', unit,rand(1,2),rand(3,2));
% bool = AnalysisDelete('trial',unit);
% 
% analysis = AnalysisAdd('trial', unit,rand(1,2),rand(3,2));
% bool = AnalysisUserDataSet('trial',unit,'oripref',45);
% AnalysisUserDataGet('trial',unit,'oripref')
% bool = AnalysisUserDataDelete('trial',unit,'oripref');
% analysis = AnalysisLoad('trial',unit)
