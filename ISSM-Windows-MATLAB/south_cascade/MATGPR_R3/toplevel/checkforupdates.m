function checkforupdates(response)
%
% CHECKFORUPDATES: Checks the MATGPR home page for updates/upgrades and
%                  bug fixes - informs the user.
%
%          Usage : checkforupdates()
%
%         Author : Andreas Tzanis,
%                  Department of Geophysics, 
%                  University of Athens
%                  atzanis@geol.uoa.gr
%
% Copyright (C) 2007, 2013, Andreas Tzanis. All rights reserved.
%
%  This program is copyrighted part of the matGPR software. Please do not
%  distribute or reverse-engineer and/or modify this program under any
%  circumstances and without the written permission of the author and
%  copyright owner.  
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%

global ENVAR
try
     new_version = urlread('http://users.uoa.gr/~atzanis/matgpr/version.txt');
catch                                                              %#ok<CTCH>
% in case user cannot access our site
    new_version = ENVAR.CURRENT_VERSION;
end
if ~strcmp(ENVAR.CURRENT_VERSION,new_version)
    wwandw = which('web'); 
    if ~isempty(wwandw),                  % MATLAB Web browser is available
        msgtext = ['The new matGPR version ' new_version ' has been released!'];
        msg = questdlg(msgtext,'matGPR Updates', 'Learn More', 'Download', 'Cancel', 'Learn More');
        waitfor(msg);
        if strcmp(msg,'Learn More'),
            %web('file:///C:\Users\Erric/Desktop/WORK/GPR/MATGPR/Updates_and_bug_fixes.txt')
            web('http://users.uoa.gr/~atzanis/matgpr/Updates_and_bug_fixes.txt')
            newtext = 'Go to "Check for Updates" to Download!';
            newmsg = msgbox(newtext,'matGPR Updates','help');
            waitfor(newmsg);
        end
        if strcmp(msg, 'Download'),
            web('http://users.uoa.gr/~atzanis/matgpr/matgpr.html')
        end
    else                              % MATLAB Web browser is NOT available
        msgtext = {['New MATGPR version ' new_version ' has been released.'],...
            '    ',...
            'Please visit http://users.uoa.gr/~atzanis/matgpr/matgpr.html ',...
            'to download and check for details.                           '};
        msg = msgbox(msgtext,'Software Update','warn');
        waitfor(msg);
    end
else
    if response == 1,
        msg = msgbox('Sorry, software updates are not presently available',...
            'matGPR Updates','error');
        waitfor(msg);
    end
end
