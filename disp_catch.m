function [] = disp_catch(CATCHobj,MainFileName,OutVarName)
% This function displays the error messages contained in the object
% returned by catch in a try&catch call.
% The function also takes a string input for diplaying the name
% of the function that initiated the try&catch call. 
% Also, the function can be used to assign the catch object in the
% base-workspace given a varible name string.
%
%Usage:
% disp_catch(CATCHobj,            MainFileName,     OutVarName);  %CATCHobj to display, Name of the file calling disp_catch(), variable name for outputting CATCHobj to base-workspace.
% disp_catch(CATCHobj, 'FileCalling_disp_catch');                 %Display catch messages & stack AND origin of the error
%                                                                 %i.e. the file which called disp_catch, e.g. 'FileCalling_disp_catch'.
%                                                                 %It is best to call it using mfilename from the function that you put it in, i.e. disp_catch(CATCHobj,mfilename)
% disp_catch(CATCHobj, 'FileCalling_disp_catch','CatchVarName');  %Also assign the catch object in base-workspace with variable name 'CatchVarName'
%
%
%How your function and call to disp_catch might look like:
% function Out = myFun(In)
% try
%    something that might go wrong causing an error
% catch CATCHobj
%    disp_catch(CATCHobj,mfilename,'Catch_myFun');
% end
%
%V2.01
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V2.01: (08.09.2015): damn bug in line 48 (cell "diplay stack if not empty") was removed. V2.0: (28.02.2015): initial implementation (changed version from original subfunction of various display tools functions)

disp(' ');
%% depending on inputs output origin of the error output
if(nargin==3)
    assignin('base',OutVarName,CATCHobj); %assign object in base-workspace as variable "varargin{2}" string
    disp(['Error occurred in function "',MainFileName,'"...']);
else
    if(nargin==2)
        disp(['Error occurred in function "',MainFileName,'"...']);
    else
        disp('Error occurred...');
    end
end

%% display identifier & message
disp([CATCHobj.identifier,': ',CATCHobj.message]);

%% diplay stack if not empty
if(~isempty(CATCHobj.stack))
    for IndStack = 1:length(CATCHobj.stack)
        disp(['Stack ',num2str(IndStack),': Error in "',CATCHobj.stack(IndStack).name,'" at line ',num2str(CATCHobj.stack(IndStack).line),'. (',CATCHobj.stack(IndStack).file,')']);
    end
end

%% Done.
disp(' ');

end