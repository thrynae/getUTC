function atomTime=getUTC(override)
%Returns the UTC time in the Matlab datenum format
%
%example syntax:
% disp(datestr(getUTC))
%
% There are several methods implemented in this function:
% - An implementation that requires a C mex function.
%   This method requires write access to a folder and a working C compiler. The compilation result
%   will be stored to a subdirectory of a folder similar to the AddOn path, or the tempdir, or the
%   current folder. Write permission is tested in that order.
% - An implementation using https://www.utctime.net/utc-timestamp.
%   The NIST has a server that returns the time, but it currently blocks API access.
%   This method requires internet access.
% - The local time and timezone offset can be determined with the wmic command or the get-date
%   Powershell function (Windows) or the date command (Linux and Mac).
%   To speed up the usage of this method, you can cache the difference with now() in a persistent
%   variable, that way you avoid the need for a slow system call.
% 
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.1.0                                                         |%
%|  Date:    2022-01-22                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - Some older releases don't support the web implementation.
% - The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10, it seems
%   reasonable to assume that the OS is the cause of the hang. For XP (and older) there is an
%   alternative strategy in place, but this has a higher likelyhood to fail.
% - Similarly, as wmic has been deprecated, a Powershell alternative should be used on newer
%   versions of Windows. The need for this is automatically detected.

if nargin==0
    %Normal flow: first try the cmd method, then the C method, then the web method.
    UTC_epoch_seconds=getUTC_cmd;
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds=getUTC_c;
    end
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds=getUTC_web;
    end
    if isempty(UTC_epoch_seconds)
        error('HJW:getUTC:TimeReadFailed',...
            ['All methods of retrieving the UTC timestamp failed.\nEnsure you ',...
            'have write access to the current folder and check your internet connection.'])
    end
else
    %Override for debug/test, this will not throw an error on fail.
    if override==1
        UTC_epoch_seconds=getUTC_c(false);
    elseif override==2
        UTC_epoch_seconds=getUTC_web;
    elseif override==3
        UTC_epoch_seconds=getUTC_cmd;
    else
        error('non-implemented override')
    end
end
UTC_offset=UTC_epoch_seconds/(24*60*60);
atomTime=UTC_offset+datenum(1970,1,1);
end
function [tf,ME]=CheckMexCompilerExistence
% Returns true if a mex compiler is expected to be installed.
% The method used for R2008a and later is fairly slow, so the flag is stored in a file. Run
% ClearMexCompilerExistenceFlag() to reset this test.
%
% This function may result in false positives (e.g. by detecting an installed compiler that doesn't
% work, or if a compiler is required for a specific language).
% False negatives should be rare.
%
% Based on: http://web.archive.org/web/2/http://www.mathworks.com/matlabcentral/answers/99389
% (this link will redirect to the URL with the full title)
%
% The actual test will be performed in a separate function. That way the same persistent can be
% used for different functions containing this check as a subfunction.

persistent tf_ ME_
if isempty(tf_)
    fun_handle='CheckMexCompilerExistence_persistent';

    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    [p,ME]=CreatePathFolder__CheckMexCompilerExistence_persistent;
    if ~isempty(ME),tf=false;return,end

    fn=fullfile(p,'ClearMexCompilerExistenceFlag.m');
    txt={...
        'function ClearMexCompilerExistenceFlag';
        'fn=create_fn';
        'if exist(fn,''file''),delete(fn),end';
        'end';
        'function fn=create_fn';
        'v=version;v=v(regexp(v,''[a-zA-Z0-9()\.]''));';
        'if ~exist(''OCTAVE_VERSION'', ''builtin'')';
        '    runtime=''MATLAB'';';
        '    type=computer;';
        'else';
        '    runtime=''OCTAVE'';';
        '    arch=computer;arch=arch(1:(min(strfind(arch,''-''))-1));';
        '    if ispc';
        '        if strcmp(arch,''x86_64'')  ,type= ''win_64'';';
        '        elseif strcmp(arch,''i686''),type= ''win_i686'';';
        '        elseif strcmp(arch,''x86'') ,type= ''win_x86'';';
        '        else                      ,type=[''win_'' arch];';
        '        end';
        '    elseif isunix && ~ismac %Essentially this is islinux';
        '        if strcmp(arch,''i686'')      ,type= ''lnx_i686'';';
        '        elseif strcmp(arch,''x86_64''),type= ''lnx_64'';';
        '        else                        ,type=[''lnx_'' arch];';
        '        end';
        '    elseif ismac';
        '        if strcmp(arch,''x86_64''),type= ''mac_64'';';
        '        else                    ,type=[''mac_'' arch];';
        '        end';
        '    end';
        'end';
        'type=strrep(strrep(type,''.'',''''),''-'','''');';
        'flag=[''flag_'' runtime ''_'' v ''_'' type ''.txt''];';
        'fn=fullfile(fileparts(mfilename(''fullpath'')),flag);';
        'end';
        ''};
    fid=fopen(fn,'wt');fprintf(fid,'%s\n',txt{:});fclose(fid);

    fn=fullfile(p,[fun_handle '.m']);
    txt={...
        'function [tf,ME]=CheckMexCompilerExistence_persistent';
        '% Returns true if a mex compiler is expected to be installed.';
        '% The method used for R2008a and later is fairly slow, so the flag is stored in';
        '% a file. Run ClearMexCompilerExistenceFlag() to reset this test.';
        '%';
        '% This function may result in false positives (e.g. by detecting an installed';
        '% compiler that doesn''t work, or if a compiler is required for a specific';
        '% language). False negatives should be rare.';
        '%';
        '% Based on: http://web.archive.org/web/2/';
        '% http://www.mathworks.com/matlabcentral/answers/99389';
        '% (this link will redirect to the URL with the full title)';
        '';
        'persistent tf_ ME_';
        'if isempty(tf_)';
        '    ME_=create_ME;';
        '    fn =create_fn;';
        '    if exist(fn,''file'')';
        '        str=fileread(fn);';
        '        tf_=strcmp(str,''compiler found'');';
        '    else';
        '        % Use evalc to suppress anything printed to the command window.';
        '        [txt,tf_]=evalc(func2str(@get_tf)); %#ok<ASGLU>';
        '        fid=fopen(fn,''w'');';
        '        if tf_,fprintf(fid,''compiler found'');';
        '        else , fprintf(fid,''compiler not found'');end';
        '        fclose(fid);';
        '    end';
        '    ';
        'end';
        'tf=tf_;ME=ME_;';
        'end';
        'function fn=create_fn';
        'v=version;v=v(regexp(v,''[a-zA-Z0-9()\.]''));';
        'if ~exist(''OCTAVE_VERSION'', ''builtin'')';
        '    runtime=''MATLAB'';';
        '    type=computer;';
        'else';
        '    runtime=''OCTAVE'';';
        '    arch=computer;arch=arch(1:(min(strfind(arch,''-''))-1));';
        '    if ispc';
        '        if strcmp(arch,''x86_64'')  ,type= ''win_64'';';
        '        elseif strcmp(arch,''i686''),type= ''win_i686'';';
        '        elseif strcmp(arch,''x86'') ,type= ''win_x86'';';
        '        else                      ,type=[''win_'' arch];';
        '        end';
        '    elseif isunix && ~ismac %Essentially this is islinux';
        '        if strcmp(arch,''i686'')      ,type= ''lnx_i686'';';
        '        elseif strcmp(arch,''x86_64''),type= ''lnx_64'';';
        '        else                        ,type=[''lnx_'' arch];';
        '        end';
        '    elseif ismac';
        '        if strcmp(arch,''x86_64''),type= ''mac_64'';';
        '        else                    ,type=[''mac_'' arch];';
        '        end';
        '    end';
        'end';
        'type=strrep(strrep(type,''.'',''''),''-'','''');';
        'flag=[''flag_'' runtime ''_'' v ''_'' type ''.txt''];';
        'fn=fullfile(fileparts(mfilename(''fullpath'')),flag);';
        'end';
        'function ME_=create_ME';
        'msg={...';
        '    ''No selected compiler was found.'',...';
        '    ''Please make sure a supported compiler is installed and set up.'',...';
        '    ''Run mex(''''-setup'''') for version-specific documentation.'',...';
        '    '''',...';
        '    ''Run ClearMexCompilerExistenceFlag() to reset this test.''};';
        'msg=sprintf(''\n%s'',msg{:});msg=msg(2:end);';
        'ME_=struct(...';
        '    ''identifier'',''HJW:CheckMexCompilerExistence:NoCompiler'',...';
        '    ''message'',msg);';
        'end';
        'function tf=get_tf';
        '[isOctave,v_num]=ver_info;';
        'if isOctave';
        '        % Octave normally comes with a compiler out of the box, but for some';
        '        % methods of installation an additional package may be required.';
        '        tf=~isempty(try_file_compile);';
        '    elseif v_num>=706 % ifversion(''>='',''R2008a'')';
        '        % Just try to compile a MWE. Getting the configuration is very slow.';
        '        [cc,TryNormalCheck]=try_file_compile;';
        '        if TryNormalCheck';
        '            % Something strange happened, so try the normal check anyway.';
        '            try cc=mex.getCompilerConfigurations;catch,cc=[];end';
        '        end';
        '        tf=~isempty(cc);';
        '    else';
        '        if ispc,ext=''.bat'';else,ext=''.sh'';end';
        '        tf=exist(fullfile(prefdir,[''mexopts'' ext]),''file'');';
        'end';
        'end';
        'function [isOctave,v_num]=ver_info';
        '% This is a compact and feature-poor equivalent of ifversion.';
        '% To save space this can be used as an alternative.';
        '% Example: R2018a is 9.4, so v_num will be 904.';
        'isOctave=exist(''OCTAVE_VERSION'', ''builtin'');';
        'v_num=version;';
        'ii=strfind(v_num,''.'');if numel(ii)~=1,v_num(ii(2):end)='''';ii=ii(1);end';
        'v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];';
        'v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);';
        'end';
        'function [cc,TryNormalCheck]=try_file_compile';
        'TryNormalCheck=false;';
        'try';
        '    [p,n]=fileparts(tempname);e=''.c'';';
        '    n=n(regexp(n,''[a-zA-Z0-9_]'')); % Keep only valid characters.';
        '    n=[''test_fun__'' n(1:min(15,end))];';
        '    fid=fopen(fullfile(p,[n e]),''w'');';
        '    fprintf(fid,''%s\n'',...';
        '        ''#include "mex.h"'',...';
        '        ''void mexFunction(int nlhs, mxArray *plhs[],'',...';
        '        ''  int nrhs, const mxArray *prhs[]) {'',...';
        '        ''    plhs[0]=mxCreateString("compiler works");'',...';
        '        ''    return;'',...';
        '        ''}'');';
        '    fclose(fid);';
        'catch';
        '    % If there is a write error in the temp dir, something is wrong.';
        '    % Just try the normal check.';
        '    cc=[];TryNormalCheck=true;return';
        'end';
        'try';
        '    current=cd(p);';
        'catch';
        '    % If the cd fails, something is wrong here. Just try the normal check.';
        '    cc=[];TryNormalCheck=true;return';
        'end';
        'try';
        '    mex([n e]);';
        '    cc=feval(str2func(n));';
        '    clear(n); % Clear to remove file lock.';
        '    cd(current);';
        'catch';
        '    % Either the mex or the feval failed. That means we can safely assume no';
        '    % working compiler is present. The normal check should not be required.';
        '    cd(current);';
        '    cc=[];TryNormalCheck=false;return';
        'end';
        'end';
        ''};
    fid=fopen(fn,'wt');fprintf(fid,'%s\n',txt{:});fclose(fid);
    fun_handle=str2func(fun_handle);

    [tf_,ME_]=feval(fun_handle);
end
tf=tf_;ME=ME_;
end
function [p,ME]=CreatePathFolder__CheckMexCompilerExistence_persistent
%Try creating a folder in either the tempdir or a persistent folder and try adding it to the path
%(if it is not already in there). If the folder is not writable, the current folder will be used.
try
    ME=[];
    p=fullfile(GetWritableFolder,'FileExchange','CheckMexCompilerExistence');
    if isempty(strfind([path ';'],[p ';'])) %#ok<STREMP>
        %This means f is not on the path.
        if ~exist(p,'dir'),mkdir(p);end
        addpath(p,'-end');
    end
catch
    ME=struct('identifier','HJW:CheckMexCompilerExistence:PathFolderFail',...
        'message','Creating a folder on the path to store compiled mex files failed.');
end
end
function UTC_epoch_seconds=getUTC_c(allow_rethrow)
%Use a C implementation, which requires write permission in a folder.
if nargin==0,allow_rethrow=true;end
persistent utc_time_c tempdir_f funname utc_time_fun_handle Compile_attempts_remaining mexfilename
if isempty(utc_time_c)
    %Try creating a folder in either the tempdir or a persistent folder and adding it to the path
    %(if it is not already in there). If the folder is not writable, the current folder will be
    %used.
    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    tempdir_f=fullfile(GetWritableFolder,'FileExchange','getUTC');
    try
        if isempty(strfind([path ';'],[tempdir_f ';'])) %#ok<STREMP>
            %This means f is not on the path.
            if ~exist(tempdir_f,'dir'),mkdir(tempdir_f);end
            addpath(tempdir_f,'-end');
        end
    catch
    end
    
    funname='utc_time';
    [mexfilename,funname]=mexname(funname);
    try utc_time_fun_handle=str2func(funname);catch,end %The try-catch is required for Octave.
    
    %Only allow a few compilation attempts, so this function doesn't cause a lot of disk I/O if
    %there is no working compiler.
    Compile_attempts_remaining=5;
    
    %prepare to write this to a file and compile
    utc_time_c={'#include "mex.h"';
        '#include "time.h"';
        '';
        '/* Abraham Cohn,  3/17/2005 */';
        '/* Philips Medical Systems */';
        '';
        'void mexFunction(int nlhs, mxArray *plhs[], int nrhs,';
        '                 const mxArray *prhs[])';
        '{';
        '  time_t utc;';
        '  ';
        '  if (nlhs > 1) {';
        '    mexErrMsgTxt("Too many output arguments");';
        '  }';
        '  ';
        '  /* Here is a nice ref: www.cplusplus.com/ref/ctime/time.html */';
        '  time(&utc);';
        '  /* mexPrintf("UTC time in local zone: %s",ctime(&utc)); */';
        '  /* mexPrintf("UTC time in GMT: %s",asctime(gmtime(&utc))); */';
        '  ';
        '  /* Create matrix for the return argument. */';
        '  plhs[0] = mxCreateDoubleScalar((double)utc);';
        '   ';
        '}'};
    %(the original had mxCreateScalarDouble)
end

try
    UTC_epoch_seconds=feval(utc_time_fun_handle);
catch
    if exist(mexfilename,'file')
        if allow_rethrow
            ME=lasterror; %#ok<LERR>
            rethrow(ME);
        else
            UTC_epoch_seconds=[];return
        end
    end
    
    %Check if there is compiler in the first place.
    if ~CheckMexCompilerExistence
        Compile_attempts_remaining=0;
        UTC_epoch_seconds=[];return
    end
    
    %Try building missing C file.
    Compile_attempts_remaining=Compile_attempts_remaining-1;
    if Compile_attempts_remaining<0 %Don't endlessly try to compile.
        UTC_epoch_seconds=[];return
    end
    
    if TestFolderWritePermission(tempdir_f)
        f=tempdir_f; %Use the folder in the tempdir to store the mex.
    else
        f=pwd; %Revert to current folder.
    end
    current_folder=cd(f);
    try
        if ~exist(fullfile(f,[funname '.c']),'file')
            fid=fopen(fullfile(f,[funname '.c']),'w');
            for line=1:numel(utc_time_c)
                fprintf(fid,'%s\n',utc_time_c{line});
            end
            fclose(fid);
        end
        
        try
            mex([funname '.c']);
        catch
        end
        %Perform cleanup.
        for ext={'c','o'}
            file=fullfile(f,[funname '.' ext{1}]);if exist(file,'file'),delete(file),end
        end
    catch
    end
    cd(current_folder);
    if exist(mexfilename,'file')
        utc_time_fun_handle=str2func(funname); %Refresh the handle.
        UTC_epoch_seconds=getUTC_c(allow_rethrow); %Use recursion to catch errors.
    else
        %The compiling of the mex function failed
        UTC_epoch_seconds=[];
    end
end
end
function UTC_epoch_seconds=getUTC_cmd
%Use a command line implementation.
%Should return an empty array instead of an error if it fails.
try
    call_type=getUTC_cmd_call_type;
catch
    warning('determination of call type failed')
    UTC_epoch_seconds=[];return
end

try
   switch call_type
       case 'Unix'
           UTC_epoch_seconds=getUTC_cmd_Unix;
       case 'WMIC_sys'
           UTC_epoch_seconds=getUTC_cmd_wmic_sys;
       case 'WMIC_bat'
           UTC_epoch_seconds=getUTC_cmd_wmic_bat;
       case 'PS_get_date'
           UTC_epoch_seconds=getUTC_cmd_PS_get_date;
   end
catch
    UTC_epoch_seconds=[];
end
end
function call_type=getUTC_cmd_call_type
persistent call_type_
if isempty(call_type_)
    if ~ispc
        call_type_='Unix';
    else
        if WinVer<=5
            %The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10, I'm
            %making the assumption that the OS is the cause of the hang.
            call_type_='WMIC_bat';
        elseif WinVer>=10
            %The cmd WMIC interface has been deprecated and seems to have been removed from Windows
            %10 21H1. On older W10 versions we can still use WMIC. Since every version of Windows
            %will either have the WMIC or PowerShell, a PowerShell call seems like a good fallback.
            [wmic_not_available,ignore]=system('wmic /?'); %#ok<ASGLU>
            if wmic_not_available
                call_type_='PS_get_date';
            else
                call_type_='WMIC_sys';
            end
        else
            call_type_='WMIC_sys';
        end
    end
end
call_type=call_type_;
end
function UTC_epoch_seconds=getUTC_cmd_PS_get_date
%Use Powershell to get the UTC time.
[s,str]=system(['powershell $a=get-date;',...
    '$a.ToUniversalTime().ToString(''yyyyMMddHHmmss'')']); %#ok<ASGLU>
str(str<48 | str>57)='';%Remove trailing newline by keeping only digits.
date=mat2cell(str,1,[4 2 2,2 2 2]);date=num2cell(str2double(date));
UTC_epoch_seconds=(datenum(date{:})-datenum(1970,1,1))*24*60*60;
end
function UTC_epoch_seconds=getUTC_cmd_Unix
[status,str]=system('date +%s'); %#ok<ASGLU>
UTC_epoch_seconds=str2double(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_bat
%If a normal system call would hang, use a temporary bat file.
pausetime=1;
fn1=[tempname,'.bat'];fn2=[tempname,'.txt'];
fid=fopen(fn1,'w');
%This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes
fprintf(fid,'%%systemroot%%\\system32\\wbem\\wmic os get LocalDateTime /value > "%s"\r\nexit',fn2);
fclose(fid);
system(['start /min "" cmd /c "' fn1 '"']);
then=now;%Store the current time to adjust for the pausetime and the function time.
pause(pausetime) % Wait for the system call to finish.
str=fileread(fn2);
try delete(fn1);catch,end,try delete(fn2);catch,end % Delete temp files.
UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str)+(now-then)*24*60*60;
end
function UTC_epoch_seconds=getUTC_cmd_wmic_sys
%Use a direct system call (instead of a temporary bat file).
[status,str]=system('wmic os get LocalDateTime /value'); %#ok<ASGLU>
%This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes
UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str)
str=str(str>=43 & str<=57); % Strip irrelevant parts.
date=mat2cell(str(1:21),1,[4 2 2,2 2 2+1+6]);date=str2double(date);
date(5)=date(5)-str2double(str(22:end)); % Add offset.
date=num2cell(date);
UTC_epoch_seconds=(datenum(date{:})-datenum(1970,1,1))*24*60*60;
end
function UTC_epoch_seconds=getUTC_web
%Read the timestamp from a web server.
%For some reason this fails for old Matlab releases (at least ML6.5-R2013b).

persistent UseWebread%Only search the path once per session.
if isempty(UseWebread)
    try UseWebread=~isempty(which(func2str(@webread)));catch,UseWebread=false;end
end

%Skip this function if there is no internet connection (3 timeouts will take a lot of time).
if ~isnetavl,UTC_epoch_seconds=[];return,end

for tries=1:3
    try
        if UseWebread
            data=webread('http://www.utctime.net/utc-timestamp');
        else
            %This probably only works for Octave.
            data=urlread('http://www.utctime.net/utc-timestamp'); %#ok<URLRD>
        end
        break
    catch
    end
end
try
    data(data==' ')='';
    pat='vartimestamp=';
    ind1=strfind(data,pat)+numel(pat);
    ind2=strfind(data,';')-1;
    ind2(ind2<ind1)=[];
    UTC_epoch_seconds=str2double(data(ind1:ind2(1)));
catch
    UTC_epoch_seconds=[];
end
end
function [f,status]=GetWritableFolder(varargin)
%Return a folder with write permission
%
% If the output folder doesn't already exist, this function will attempt to create it. This
% function should provide a reliable and repeatable location to write files.
%
% Syntax:
%   f=GetWritableFolder
%   [f,status]=GetWritableFolder
%   [__]=GetWritableFolder(Name,Value)
%   [__]=GetWritableFolder(options)
%
% Input/output arguments:
% f:
%   Char array with the full path to the writable folder. This does not contain a trailing filesep.
% status:
%   A scalar double ranging from 0 to 3. 0 denotes a failure to find a folder, 1 means the folder
%   is in a folder close to the AddOn folder, 2 that it is a folder in the tempdir, 3 mean that the
%   returned path is a folder in the current directory.
% options:
%   A struct with Name,Value parameters. Missing parameters are filled with the defaults listed
%   below. Using incomplete parameter names or incorrect capitalization is allowed, as long as
%   there is a unique match.
%
% Name,Value parameters:
%   ForceStatus:
%      Retrieve the path corresponding to the status value. Using 0 allows an automatic
%      determination of the location (default=0;).
%    ErrorOnNotFound:
%      Throw an error when failing to find a writeable folder (default=true;).
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.0.0                                                         |%
%|  Date:    2021-02-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - The path returned with status=1 is mostly the same as the addonpath for most releases. Although
%   it is not correct for all release/OS combinations, it should still work. If you have a managed
%   account, this might result in strange behavior.

[success,options,ME]=GetWritableFolder_parse_inputs(varargin{:});
if ~success
    rethrow(ME)
else
    [ForceStatus,ErrorOnNotFound,root_folder_list]=deal(options.ForceStatus,...
        options.ErrorOnNotFound,options.root_folder_list);
end
root_folder_list{end}=pwd;% Set this default here to avoid storing it in a persistent.
if ForceStatus
    status=ForceStatus;f=fullfile(root_folder_list{status},'PersistentFolder');
    try if ~exist(f,'dir'),mkdir(f);end,catch,end
    return
end

% Option 1: use a folder similar to the AddOn Manager.
status=1;f=root_folder_list{status};
try if ~exist(f,'dir'),mkdir(f);end,catch,end
if ~TestFolderWritePermission(f)
    % If the Add-On path is not writable, return the tempdir. It will not be persistent, but it
    % will be writable.
    status=2;f=root_folder_list{status};
    try if ~exist(f,'dir'),mkdir(f);end,catch,end
    if ~TestFolderWritePermission(f)
        % The tempdir should always be writable, but if for some reason it isn't: return the pwd.
        status=3;f=root_folder_list{status};
    end
end

% Add 'PersistentFolder' to whichever path was determined above.
f=fullfile(f,'PersistentFolder');
try if ~exist(f,'dir'),mkdir(f);end,catch,end

if ~TestFolderWritePermission(f)
    % Apparently even the pwd isn't writable, so we will either return an error, or a fail state.
    if ErrorOnNotFound
        error('HJW:GetWritableFolder:NoWritableFolder',...
            'This function was unable to find a folder with write permissions.')
    else
        status=0;f='';
    end
end
end
function [success,options,ME]=GetWritableFolder_parse_inputs(varargin)
%Parse the inputs of the GetWritableFolder function
% This function returns a success flag, the parsed options, and an ME struct.
% As input, the options should either be entered as a struct or as Name,Value pairs. Missing fields
% are filled from the default.

% Pre-assign outputs.
success=false;
ME=struct('identifier','','message','');

persistent default
if isempty(default)
    %Set defaults for options.
    default.ForceStatus=false;
    default.ErrorOnNotFound=false;
    default.root_folder_list={...
        GetPseudoAddonpath;
        fullfile(tempdir,'MATLAB');
        ''};% Overwrite this last element with pwd when called.
end

if nargin==2
    options=default;
    success=true;
    return
end

% Actually parse the Name,Value pairs (or the struct).
[options,replaced]=parse_NameValue(default,varargin{:});

% Test the optional inputs.
for k=1:numel(replaced)
    curr_option=replaced{k};
    item=options.(curr_option);
    ME.identifier=['HJW:GetWritableFolder:incorrect_input_opt_' lower(curr_option)];
    switch curr_option
        case 'ForceStatus'
            try
                if ~isa(default.root_folder_list{item},'char')
                    % This ensures an error for item=[true false true]; as well.
                    error('the indexing must have failed, trigger error')
                end
            catch
                ME.message=sprintf('Invalid input: expected a scalar integer between 1 and %d.',...
                    numel(default.root_folder_list));
                return
            end
        case 'ErrorOnNotFound'
            [passed,options.ErrorOnNotFound]=test_if_scalar_logical(item);
            if ~passed
                ME.message='ErrorOnNotFound should be either true or false.';
                return
            end
        otherwise
            ME.message=sprintf('Name,Value pair not recognized: %s.',curr_option);
            ME.identifier='HJW:GetWritableFolder:incorrect_input_NameValue';
            return
    end
end
success=true;ME=[];
end
function f=GetPseudoAddonpath
% This is mostly the same as the addonpath. Technically this is not correct for all release/OS
% combinations and the code below should be used:
%     addonpath='';
%     try s = Settings;addonpath=get(s.matlab.addons,'InstallationFolder');end %#ok<TRYNC>
%     try s = Settings;addonpath=get(s.matlab.apps,'AppsInstallFolder');end %#ok<TRYNC>
%     try s = settings;addonpath=s.matlab.addons.InstallationFolder.ActiveValue;end %#ok<TRYNC>
%
% However, this returns an inconsistent output:
%     R2011a          <pref doesn't exist>
%     R2015a Ubuntu  $HOME/Documents/MATLAB/Apps
%            Windows %HOMEPATH%\MATLAB\Apps
%     R2018a Ubuntu  $HOME/Documents/MATLAB/Add-Ons
%            Windows %HOMEPATH%\MATLAB\Add-Ons
%     R2020a Windows %APPDATA%\MathWorks\MATLAB Add-Ons
%
% To make the target folder consistent, only one of these options is chosen.
if ispc
    [ignore,appdata]=system('echo %APPDATA%');appdata(appdata<14)=''; %#ok<ASGLU> (remove LF/CRLF)
    f=fullfile(appdata,'MathWorks','MATLAB Add-Ons');
else
    [ignore,home_dir]=system('echo $HOME');home_dir(home_dir<14)=''; %#ok<ASGLU> (remove LF/CRLF)
    f=fullfile(home_dir,'Documents','MATLAB','Add-Ons');
end
end
function [connected,timing]=isnetavl(use_HTML_test_only)
%Check for an internet connection by pinging Google
%
% Syntax:
%   [connected,timing]=isnetavl(use_HTML_test_only)
%
% Input/output arguments:
% connected:
%   A logical denoting the connectivity. Note that this may return unexpected results if Matlab has
%   separate proxy settings and/or if your DNS is having issues.
% timing:
%   This contains the ping time as a double and defaults to 0 if there is no connection.
%   Note that this value will be larger than the true ping time if the HTML fallback is used.
% use_HTML_test_only:
%   If the input is convertible to true the HTML method will be used, which has the benefit of
%   testing if Matlab (or Octave) is able to connect to the internet. This might be relevant if
%   there are proxy settings or firewall rules specific to Matlab/Octave.
%   Note that the ping value will be larger than the true value if the HTML fallback is used.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.0.0                                                         |%
%|  Date:    2021-05-19                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). For the full test matrix, see the HTML doc.
% Compatibility considerations:
% - This is expected to work on all releases.
% - If you have proxy/firewall settings specific to Matlab/Octave, make sure to use the HTML method
%   by providing an input.

if nargin==0,                                             tf=false;
else,[tf1,tf2]=test_if_scalar_logical(use_HTML_test_only);tf=tf1&&tf2;
end
if tf
    %(the timing is not reliable)
    [connected,timing]=isnetavl___ping_via_html;
    return
end

tf=isnetavl__ICMP_is_blocked;
if isempty(tf)
    %Unable to determine if ping is allowed, the connection must be down.
    connected=0;
    timing=0;
else
    if tf
        %Ping is not allowed.
        %(the timing is not reliable)
        [connected,timing]=isnetavl___ping_via_html;
    else
        %Ping is allowed.
        [connected,timing]=isnetavl___ping_via_system;
    end
end
end
function [connected,timing]=isnetavl___ping_via_html
%Ping is blocked by some organizations. As an alternative, the google.com page can be loaded as a
%normal HTML, which should work as well, although it is slower. This also means the ping timing is
%no longer reliable.
persistent UseWebread
if isempty(UseWebread)
    try no_webread=isempty(which(func2str(@webread)));catch,no_webread=true;end
	UseWebread=~no_webread;
end
try
    then=now;
    if UseWebread
        str=webread('http://google.com'); %#ok<NASGU>
    else
        str=urlread('http://google.com'); %#ok<NASGU,URLRD>
    end
    connected=1;
    timing=(now-then)*24*3600*1000;
catch
    connected=0;
    timing=0;
end
end
function [connected,timing]=isnetavl___ping_via_system
if ispc
    try
        %                                   8.8.4.4 will also work
        [ignore_output,b]=system('ping -n 1 8.8.8.8');%#ok<ASGLU> ~
        stats=b(strfind(b,' = ')+3);
        stats=stats(1:3);%[sent received lost]
        if ~strcmp(stats,'110')
            error('trigger error')
        else
            %This branch will error for 'destination host unreachable'.
            connected=1;
            %This assumes there is only one place with '=[digits]ms' in the response, but this code
            %is not language-specific.
            [ind1,ind2]=regexp(b,' [0-9]+ms');
            timing=b((ind1(1)+1):(ind2(1)-2));
            timing=str2double(timing);
        end
    catch
        connected=0;
        timing=0;
    end
elseif isunix
    try
        %                                   8.8.4.4 will also work
        [ignore_output,b]=system('ping -c 1 8.8.8.8');%#ok<ASGLU> ~
        ind=regexp(b,', [01] ');
        if b(ind+2)~='1'
            %This branch includes 'destination host unreachable' errors.
            error('trigger error')
        else
            connected=1;
            %This assumes the first place with '=[digits] ms' in the response contains the ping
            %timing. This code is not language-specific.
            [ind1,ind2]=regexp(b,'=[0-9.]+ ms');
            timing=b((ind1(1)+1):(ind2(1)-2));
            timing=str2double(timing);
        end
    catch
        connected=0;
        timing=0;
    end
else
    error('How did you even get Matlab to work?')
end
end
function [tf,connected,timing]=isnetavl__ICMP_is_blocked
%Check if ICMP 0/8/11 is blocked.
%
%tf is empty if both methods fail.

persistent output
if ~isempty(output)
    tf=output;return
end

%First check if ping works.
[connected,timing]=isnetavl___ping_via_system;
if connected
    %Ping worked and there is an internet connection.
    output=false;
    tf=false;
    return
end

%There are two options: no internet connection, or ping is blocked.
[connected,timing]=isnetavl___ping_via_html;
if connected
    %There is an internet connection, therefore ping must be blocked.
    output=true;
    tf=true;
    return
end

%Both methods failed, internet is down. Leave the value of tf (and the persistent variable) set to
%empty so it is tried next time.
tf=[];
end
function [mex_filename,fun_name]=mexname(fun_name)
%Encode runtime version information in the function name.
% This can be useful if multiple versions of Matlab or Octave need to use the
% same folder to store compiled functions, while not being compatible.
%
% This function replaces a syntax like mex_filename=[fun_name '.' mexext].
%
% Syntax:
%   mex_filename=mexname(fun_name);
%   [mex_filename,updated_fun_name]=mexname(fun_name);
persistent append
if isempty(append)
    v=version;ind=[strfind(v,'.') numel(v)];
    v=sprintf('%02d.%02d',str2double({...
        v(1:(ind(1)-1)           ) ,...
        v(  (ind(1)+1):(ind(2)-1)) }));
    v=['v' strrep(v,'.','_')];
    if ~exist('OCTAVE_VERSION', 'builtin')
        runtime='MATLAB';
        type=computer;
    else
        runtime='OCTAVE';
        arch=computer;arch=arch(1:(min(strfind(arch,'-'))-1));
        if ispc
            if strcmp(arch,'x86_64')  ,type= 'win_64';
            elseif strcmp(arch,'i686'),type= 'win_i686';
            elseif strcmp(arch,'x86') ,type= 'win_x86';
            else                      ,type=['win_' arch];
            end
        elseif isunix && ~ismac %Essentially this is islinux
            if strcmp(arch,'i686')      ,type= 'lnx_i686';
            elseif strcmp(arch,'x86_64'),type= 'lnx_64';
            else                        ,type=['lnx_' arch];
            end
        elseif ismac
            if strcmp(arch,'x86_64'),type= 'mac_64';
            else                    ,type=['mac_' arch];
            end
        end
    end
    type=strrep(strrep(type,'.',''),'-','');
    append=cell(2,1);
    append{1}=['_' runtime '_' v '_' type];
    append{2}=[append{1} '.' mexext];
end

try %Test if fun_name is a valid name.
    if ~isvarname(fun_name),error('trigger catch block'),end
catch
    error('HJW:mexname:InvalidName',...
        'The provided input can''t be a function name')
end

mex_filename=[fun_name append{2}];
fun_name=[fun_name append{1}];
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, ignoring incomplete names and case.
%
%If this fails to find a match, this will throw an error with the offending name as the message.
%
%The second output is a cell containing the field names that have been set.

opts=default;replaced={};
if nargin==1,return,end

%Unwind an input struct to Name,Value pairs.
try
    struct_input=numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input=mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        Names=fieldnames(varargin{1});
        Values=struct2cell(varargin{1});
    else
        %Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names=cellstr(varargin(1:2:end));
        Values=varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    %If this block errors, that is either because a missing Value with the Name,Value syntax, or
    %because the struct input is not a struct, or because an attempt was made to mix the two
    %styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as struct.')
end

%Convert the real fieldnames to a char matrix.
d_Names1=fieldnames(default);d_Names2=lower(d_Names1);
len=cellfun('prodofsize',d_Names2);maxlen=max(len);
for n=find(len<maxlen).' %Pad with spaces where needed
    d_Names2{n}((end+1):maxlen)=' ';
end
d_Names2=vertcat(d_Names2{:});

%Attempt to match the names.
replaced=false(size(d_Names1));
for n=1:numel(Names)
    name=lower(Names{n});
    tmp=d_Names2(:,1:min(end,numel(name)));
    non_matching=numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
    match_idx=find(non_matching==0);
    if numel(match_idx)~=1
        error('parse_NameValue:NonUniqueMatch',Names{n})
    end
    
    %Store the Value in the output struct.
    opts.(d_Names1{match_idx})=Values{n};
    
    replaced(match_idx)=true;
end
replaced=d_Names1(replaced);
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%The char and string test are not case sensitive.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
%- 'enable' or 'disable'
%- 'enabled' or 'disabled'
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    try
        states(end+1,:)=eval('{"on","off"}'); %#ok<EVLCS> 
    catch
    end
end
isLogical=true;
try
    if isa(val,'char') || isa(val,'string')
        try val=lower(val);catch,end
    end
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end
function tf=TestFolderWritePermission(f)
%Returns true if the folder exists and allows Matlab to write files.
%An empty input will generally test the pwd.
%
%examples:
%  fn='foo.txt';if ~TestFolderWritePermission(fileparts(fn)),error('can''t write!'),end

if ~( isempty(f) || exist(f,'dir') )
    tf=false;return
end

fn='';
while isempty(fn) || exist(fn,'file')
    %Generate a random file name, making sure not to overwrite any existing file.
    %This will try to create a file without an extension.
    [ignore,fn]=fileparts(tmpname('write_permission_test_','.txt')); %#ok<ASGLU>
    fn=fullfile(f,fn);
end
try
    %Test write permission.
    fid=fopen(fn,'w');fprintf(fid,'test');fclose(fid);
    delete(fn);
    tf=true;
catch
    %Attempt to clean up.
    if exist(fn,'file'),try delete(fn);catch,end,end
    tf=false;
end
end
function str=tmpname(StartFilenameWith,ext)
%Inject a string in the file name part returned by the tempname function.
if nargin<1,StartFilenameWith='';end
if ~isempty(StartFilenameWith),StartFilenameWith=[StartFilenameWith '_'];end
if nargin<2,ext='';else,if ~strcmp(ext(1),'.'),ext=['.' ext];end,end
str=tempname;
[p,f]=fileparts(str);
str=fullfile(p,[StartFilenameWith f ext]);
end
function out=WinVer
%This returns the main Windows version number (5 for XP, 6 for Vista, etc).
%It will return an empty array for non-Windows or in case of errors.
persistent persistent_val
if ~ispc,out=[];return,end
if isempty(persistent_val)
    try
        [status,str] = system('ver'); %#ok<ASGLU>
        args={'[^0-9]*(\d*).*','$1','tokenize'};
        v=version;v(min(find(v=='.')):end)='';pre_v7=str2double(v)<7; %#ok<MXFND> 
        isOctave=exist('OCTAVE_VERSION','builtin')~=0;
        if ~pre_v7 && ~isOctave
            args(end)=[];%The 'tokenize' option became the default in R14 (v7).
        end
        persistent_val=str2double(regexprep(str,args{:}));
    catch
        try
            [status,str]=system('systeminfo'); %#ok<ASGLU>
            ind1= 1+strfind(str,':');ind1=ind1(3);
            ind2=-1+strfind(str,'.');ind2(ind2<ind1)=[];
            persistent_val=str2double(str(ind1:ind2(1)));
        catch
            persistent_val=[];
        end
    end
end
out=persistent_val;
end