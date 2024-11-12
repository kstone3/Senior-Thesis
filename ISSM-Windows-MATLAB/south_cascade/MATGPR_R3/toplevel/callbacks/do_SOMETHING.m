function OPD = do_SOMETHING(IPD, other_argument)
%
% Interface to drive "SOMETHING.m": TEMPLATE FOR THE FAST CREATION OF
% MATGPR CALLBACK FUNCTIONS  
% Author : Andreas Tzanis, 
%          Department of Geophysics, 
%          University of Athens, 
%          atzanis@geol.uoa.gr
%
% Copyright (C) 2005, 2013, Andreas Tzanis. All rights reserved.
%
% This program is copyrighted part of the matGPR software. Do not
% distribute or reverse-engineer and/or modify this program under any
% circumstances without the written permission of the author and
% copyright owner.  
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%

OPD = discardprocdata;                           % discard the current OPD
% Trap common errors ...
% Example 1 : No data in memory ...
if isvoidXPD(IPD),  return,  end;
% Example 2: Unequally spaced data where they shouldn't be! 
if isuneqXPD(IPD),  return,  end;

% Test CONDITIONS for other_argument - delete if not applicable
if  other_argument,                     
    erh = warndlg('ERROR or WARNING MESSAGE', ...
        'MATGPR : ERROR!'); 
    uiwait(erh); 
    return; 
end; 

% Now Proceed ...
OPD = IPD;                                                % Copy IPD to OPD

% Run SOMETHING.m
OPD = SOMETHING(IPD, other_argument);

% Check if Operation has been canceled or aborted and issue a message
if isempty(OPD.d),
    disp('SOMETHING > Operation aborted - No O/P data returned!');
    return
end

% Display data
viewdata(OPD.x,OPD.z,OPD.d,'outdata',OPD.xlab,OPD.zlab);
title('SOMETHING has been done')

% Update processing history
iss = size(OPD.history,1); 
text = 'SOMETHING has indeed been done';
OPD.history(iss+1,1) = cellstr(text); 

return
