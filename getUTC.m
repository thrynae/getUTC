function atomTime=getUTC(override)
%Returns the UTC time in the Matlab datenum format
%
% example syntax:
%  disp(datestr(getUTC))
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
%   On Windows an NTP query is also implemented (internet connectivity is checked prior to this
%   call).
%
%   To speed up the usage of this method, you can cache the difference with now() in a persistent
%   variable, that way you avoid the need for a slow system call.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.2.0                                                         |%
%|  Date:    2023-11-21                                                    |%
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
%   alternative strategy in place, but this has a higher likelihood to fail.
% - Similarly, as wmic has been deprecated, a Powershell alternative should be used on newer
%   versions of Windows. The need for this is automatically detected.
%
% /=========================================================================================\
% ||                     | Windows             | Linux               | MacOS               ||
% ||---------------------------------------------------------------------------------------||
% || Matlab R2023b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2023a       | W11: Pass           | Ubuntu 22.04: Pass  |                     ||
% || Matlab R2022b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2022a       | W11: Pass           |                     |                     ||
% || Matlab R2021b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2021a       | W11: Pass           |                     |                     ||
% || Matlab R2020b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2020a       | W11: Pass           |                     |                     ||
% || Matlab R2019b       | W11: Pass           | Ubuntu 22.04: Pass  | Monterey: Pass      ||
% || Matlab R2019a       | W11: Fail*          |                     |                     ||
% || Matlab R2018a       | W11: Fail*          | Ubuntu 22.04: Fail* |                     ||
% || Matlab R2017b       | W11: Fail*          | Ubuntu 22.04: Fail* | Monterey: Fail*     ||
% || Matlab R2016b       | W11: Fail*          | Ubuntu 22.04: Fail* | Monterey: Fail*     ||
% || Matlab R2015a       | W11: Pass           | Ubuntu 22.04: Pass  |                     ||
% || Matlab R2013b       | W11: Fail*          |                     |                     ||
% || Matlab R2007b       | W11: Fail*          |                     |                     ||
% || Matlab 6.5 (R13)    | W11: Fail*          |                     |                     ||
% || Octave 8.2.0        | W11: Pass           |                     |                     ||
% || Octave 7.2.0        | W11: Pass           |                     |                     ||
% || Octave 6.2.0        | W11: Pass           | Ubuntu 22.04:       | Catalina:           ||
% || Octave 5.2.0        | W11: Pass           |                     |                     ||
% || Octave 4.4.1        | W11: Pass           |                     | Catalina:           ||
% \=========================================================================================/
%     * See the compatibility considerations. Anything not mentioned there can still be
%       expected to work. Run the included tester function to verify these results for
%       your system and release.
%       Tests involving graphics that were skipped are encoded as fail in this table.

if nargin==0
    % Normal flow: first try the cmd method, then the C method, then the web method.
    UTC_epoch_seconds = getUTC_cmd;
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds = getUTC_c;
    end
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds = getUTC_web;
    end
    if isempty(UTC_epoch_seconds)
        error('HJW:getUTC:TimeReadFailed',...
            ['All methods of retrieving the UTC timestamp failed.\nEnsure you ',...
            'have write access to the current folder and check your internet connection.'])
    end
else
    % Override for debug/test, this will not throw an error on fail.
    switch override
        case 1
            UTC_epoch_seconds = getUTC_c(false);
        case 2
            UTC_epoch_seconds = getUTC_web;
        case 3
            UTC_epoch_seconds = getUTC_cmd;
        otherwise
            if isa(override,'char')
                UTC_epoch_seconds = getUTC_cmd(override);
            else
                error('non-implemented override')
            end
    end
end
UTC_offset = UTC_epoch_seconds/(24*60*60);
atomTime = UTC_offset+datenum(1970,1,1); %#ok<DATNM>
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
    % In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    % putting this code in this block, we are trying to keep these queries to a minimum.
    [p,ME] = CreatePathFolder__CheckMexCompilerExistence_persistent;
    if ~isempty(ME),tf = false;return,end
    
    fn = fullfile(p,'ClearMexCompilerExistenceFlag.m');
    txt = {...
        'function ClearMexCompilerExistenceFlag',...
        'fn = create_fn;',...
        'if exist(fn,''file''),delete(fn),end',...
        'end',...
        'function fn = create_fn',...
        'v = version;v = v(regexp(v,''[a-zA-Z0-9()\.]''));',...
        'if ~exist(''OCTAVE_VERSION'', ''builtin'')',...
        '    runtime = ''MATLAB'';',...
        '    type = computer;',...
        'else',...
        '    runtime = ''OCTAVE'';',...
        '    arch = computer;arch = arch(1:(min(strfind(arch,''-''))-1));',...
        '    if ispc',...
        '        if strcmp(arch,''x86_64'')  ,type =  ''win_64'';',...
        '        elseif strcmp(arch,''i686''),type =  ''win_i686'';',...
        '        elseif strcmp(arch,''x86'') ,type =  ''win_x86'';',...
        '        else                      ,type = [''win_'' arch];',...
        '        end',...
        '    elseif isunix && ~ismac % Essentially this is islinux',...
        '        if strcmp(arch,''i686'')      ,type =  ''lnx_i686'';',...
        '        elseif strcmp(arch,''x86_64''),type =  ''lnx_64'';',...
        '        else                        ,type = [''lnx_'' arch];',...
        '        end',...
        '    elseif ismac',...
        '        if strcmp(arch,''x86_64''),type =  ''mac_64'';',...
        '        else                    ,type = [''mac_'' arch];',...
        '        end',...
        '    end',...
        'end',...
        'type = strrep(strrep(type,''.'',''''),''-'','''');',...
        'flag = [''flag_'' runtime ''_'' v ''_'' type ''.txt''];',...
        'fn = fullfile(fileparts(mfilename(''fullpath'')),flag);',...
        'end',...
        ''};
    fid = fopen(fn,'wt');fprintf(fid,'%s\n',txt{:});fclose(fid);
    
    [tf_,ME_] = CheckMexCompilerExistence_persistent(p);
end
tf = tf_;ME = ME_;
end
function [tf,ME]=CheckMexCompilerExistence_persistent(p)
% Returns true if a mex compiler is expected to be installed.
% The method used for R2008a and later is fairly slow, so the flag is stored in a file. Run
% ClearMexCompilerExistenceFlag() to reset this test.
%
% This function may result in false positives (e.g. by detecting an installed compiler that doesn't
% work, or if a compiler is required for a specific language). False negatives should be rare.
%
% Based on: http://web.archive.org/web/2/http://www.mathworks.com/matlabcentral/answers/99389
% (this link will redirect to the URL with the full title)

persistent tf_ ME_
if isempty(tf_)
    ME_ = create_ME;
    fn  = create_fn(p);
    if exist(fn,'file')
        str = fileread(fn);
        tf_ = strcmp(str,'compiler found');
    else
        % Use evalc to suppress anything printed to the command window.
        [txt,tf_] = evalc(func2str(@get_tf)); %#ok<ASGLU>
        fid = fopen(fn,'w');
        if tf_,fprintf(fid,'compiler found');
        else , fprintf(fid,'compiler not found');end
        fclose(fid);
    end
    
end
tf = tf_;ME = ME_;
end
function fn=create_fn(p)
v = version;v = v(regexp(v,'[a-zA-Z0-9()\.]'));
if ~exist('OCTAVE_VERSION', 'builtin')
    runtime = 'MATLAB';
    type = computer;
else
    runtime = 'OCTAVE';
    arch = computer;arch = arch(1:(min(strfind(arch,'-'))-1));
    if ispc
        if strcmp(arch,'x86_64')  ,type =  'win_64';
        elseif strcmp(arch,'i686'),type =  'win_i686';
        elseif strcmp(arch,'x86') ,type =  'win_x86';
        else                      ,type = ['win_' arch];
        end
    elseif isunix && ~ismac % Essentially this is islinux.
        if strcmp(arch,'i686')      ,type =  'lnx_i686';
        elseif strcmp(arch,'x86_64'),type =  'lnx_64';
        else                        ,type = ['lnx_' arch];
        end
    elseif ismac
        if strcmp(arch,'x86_64'),type =  'mac_64';
        else                    ,type = ['mac_' arch];
        end
    end
end
type = strrep(strrep(type,'.',''),'-','');
flag = ['flag_' runtime '_' v '_' type '.txt'];
fn = fullfile(p,flag);
end
function ME_=create_ME
msg = {...
    'No selected compiler was found.',...
    'Please make sure a supported compiler is installed and set up.',...
    'Run mex(''-setup'') for version-specific documentation.',...
    '',...
    'Run ClearMexCompilerExistenceFlag() to reset this test.'};
msg = sprintf('\n%s',msg{:});msg = msg(2:end);
ME_ = struct(...
    'identifier','HJW:CheckMexCompilerExistence:NoCompiler',...
    'message',msg);
end
function tf=get_tf
[isOctave,v_num] = ver_info;
if isOctave
    % Octave normally comes with a compiler out of the box, but for some methods of installation an
    % additional package may be required.
    tf = ~isempty(try_file_compile);
elseif v_num>=706 % ifversion('>=','R2008a')
    % Just try to compile a MWE. Getting the configuration is very slow. On Windows this is a bad
    % idea, as it starts an interactive prompt. Because this function is called with evalc, that
    % means this function will hang.
    if ispc, TryNormalCheck  = true;
    else,[cc,TryNormalCheck] = try_file_compile;
    end
    if TryNormalCheck
        % Something strange happened, so try the normal check anyway.
        try cc = mex.getCompilerConfigurations;catch,cc=[];end
    end
    tf = ~isempty(cc);
else
    if ispc,ext = '.bat';else,ext = '.sh';end
    tf = exist(fullfile(prefdir,['mexopts' ext]),'file');
end
end
function [isOctave,v_num]=ver_info
% This is a compact and feature-poor equivalent of ifversion.
% To save space this can be used as an alternative.
% Example: R2018a is 9.4, so v_num will be 904.
isOctave = exist('OCTAVE_VERSION', 'builtin');
v_num = version;
ii = strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end) = '';ii = ii(1);end
v_num = [str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
v_num = v_num(1)+v_num(2)/100;v_num = round(100*v_num);
end
function [cc,TryNormalCheck]=try_file_compile
TryNormalCheck = false;
try
    [p,n] = fileparts(tempname);e='.c';
    n = n(regexp(n,'[a-zA-Z0-9_]')); % Keep only valid characters.
    n = ['test_fun__' n(1:min(15,end))];
    fid = fopen(fullfile(p,[n e]),'w');
    fprintf(fid,'%s\n',...
        '#include "mex.h"',...
        'void mexFunction(int nlhs, mxArray *plhs[],',...
        '  int nrhs, const mxArray *prhs[]) {',...
        '    plhs[0]=mxCreateString("compiler works");',...
        '    return;',...
        '}');
    fclose(fid);
catch
    % If there is a write error in the temp dir, something is wrong.
    % Just try the normal check.
    cc = [];TryNormalCheck = true;return
end
try
    current = cd(p);
catch
    % If the cd fails, something is wrong here. Just try the normal check.
    cc = [];TryNormalCheck = true;return
end
try
    mex([n e]);
    cc = feval(str2func(n));
    clear(n); % Clear to remove file lock.
    cd(current);
catch
    % Either the mex or the feval failed. That means we can safely assume no working compiler is
    % present. The normal check should not be required.
    cd(current);
    cc = [];TryNormalCheck = false;return
end
end
function [p,ME]=CreatePathFolder__CheckMexCompilerExistence_persistent
% Try creating a folder in either the tempdir or a persistent folder and try adding it to the path
% (if it is not already in there). If the folder is not writable, the current folder will be used.
try
    ME = [];
    p = fullfile(GetWritableFolder,'FileExchange','CheckMexCompilerExistence');
    if isempty(strfind([path ';'],[p ';'])) %#ok<STREMP>
        % This means f is not on the path.
        if ~exist(p,'dir'),makedir(p);end
        addpath(p,'-end');
    end
catch
    ME = struct('identifier','HJW:CheckMexCompilerExistence:PathFolderFail',...
        'message','Creating a folder on the path to store the compiled function and flag failed.');
end
end
function UTC_epoch_seconds=getUTC_c(allow_rethrow)
% Use a C implementation, which requires write permission in a folder.
if nargin==0,allow_rethrow = true;end
persistent utc_time_c tempdir_f funname utc_time_fun_handle Compile_attempts_remaining mexfilename
if isempty(utc_time_c)
    % Try creating a folder in either the tempdir or a persistent folder and adding it to the path
    % (if it is not already in there). If the folder is not writable, the current folder will be
    % used.
    % In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    % putting this code in this block, we are trying to keep these queries to a minimum.
    tempdir_f = fullfile(GetWritableFolder,'FileExchange','getUTC');
    try
        if isempty(strfind([path ';'],[tempdir_f ';'])) %#ok<STREMP>
            % This means f is not on the path.
            if ~exist(tempdir_f,'dir'),makedir(tempdir_f);end
            addpath(tempdir_f,'-end');
        end
    catch
    end
    
    funname = 'utc_time';
    [mexfilename,funname] = mexname(funname);
    try utc_time_fun_handle = str2func(funname);catch,end % The try-catch is required for Octave.
    
    % Only allow a few compilation attempts, so this function doesn't cause a lot of disk I/O if
    % there is no working compiler.
    Compile_attempts_remaining = 5;
    
    % Prepare to write this to a file and compile.
    utc_time_c = {...
        '#include "mex.h"';
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
    UTC_epoch_seconds = feval(utc_time_fun_handle);
catch
    if exist(mexfilename,'file')
        if allow_rethrow
            ME = lasterror; %#ok<LERR>
            rethrow(ME);
        else
            UTC_epoch_seconds = [];return
        end
    end
    
    % Check if there is compiler in the first place.
    if ~CheckMexCompilerExistence
        Compile_attempts_remaining = 0;
        UTC_epoch_seconds = [];return
    end
    
    % Try building missing C file.
    Compile_attempts_remaining = Compile_attempts_remaining-1;
    if Compile_attempts_remaining<0 % Don't endlessly try to compile.
        UTC_epoch_seconds = [];return
    end
    
    if TestFolderWritePermission(tempdir_f)
        f = tempdir_f; % Use the folder in the tempdir to store the mex.
    else
        f = pwd; % Revert to current folder.
    end
    current_folder = cd(f);
    try
        if ~exist(fullfile(f,[funname '.c']),'file')
            fid = fopen(fullfile(f,[funname '.c']),'w');
            for line=1:numel(utc_time_c)
                fprintf(fid,'%s\n',utc_time_c{line});
            end
            fclose(fid);
        end
        
        try
            % Capture the status message in a variable to suppress it.
            ignore = evalc(['mex([ ''' funname '.c'']);']); %#ok<NASGU>
        catch
        end
        % Perform cleanup.
        for ext={'c','o'}
            file = fullfile(f,[funname '.' ext{1}]);if exist(file,'file'),delete(file),end
        end
    catch
    end
    cd(current_folder);
    if exist(mexfilename,'file')
        utc_time_fun_handle = str2func(funname); % Refresh the handle.
        UTC_epoch_seconds = getUTC_c(allow_rethrow); % Use recursion to catch errors.
    else
        % The compiling of the mex function failed
        UTC_epoch_seconds = [];
    end
end
end
function [UTC_epoch_seconds,call_type]=getUTC_cmd(override_isnetavl,override_call_type)
% Use a command line implementation.
% This should return an empty array instead of an error if it fails.

if nargin==2
    call_type = override_call_type;
else
    if nargin>=1 && ~isempty(override_isnetavl)
        InternetConnected = override_isnetavl;
    else
        InternetConnected = isnetavl;
    end
    try
        call_type = getUTC_cmd_call_type;
        try
            if InternetConnected
                call_type = call_type.online;
            else
                error('trigger offline version')
            end
        catch
            call_type = call_type.offline;
        end
    catch
        warning('determination of call type failed')
        UTC_epoch_seconds = [];return
    end
end
if nargout==2
    try
        call_type = getUTC_cmd_call_type;
    catch
        warning('determination of call type failed')
    end
    UTC_epoch_seconds = [];
    return
end

try
    switch call_type
        case 'Unix'
            UTC_epoch_seconds = getUTC_cmd_Unix;
        case 'NTP_win'
            UTC_epoch_seconds = getUTC_cmd_NTP_win;
        case 'WMIC_sys'
            UTC_epoch_seconds = getUTC_cmd_wmic_sys;
        case 'WMIC_bat'
            UTC_epoch_seconds = getUTC_cmd_wmic_bat;
        case 'PS_get_date'
            UTC_epoch_seconds = getUTC_cmd_PS_get_date;
        otherwise
            error('call type not implemented')
    end
catch
    UTC_epoch_seconds = [];
end
end
function call_type=getUTC_cmd_call_type
persistent call_type_
if isempty(call_type_)
    if ~ispc
        call_type_ = struct('offline','Unix','online','Unix');
    else
        if WinVer<=5
            % The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10,
            % I'm making the assumption that the OS is the cause of the hang.
            call_type_ = 'WMIC_bat';
        elseif WinVer>=10
            % The cmd WMIC interface has been deprecated and seems to have been removed from
            % Windows 10 21H1. On older W10 versions we can still use WMIC. Since every version of
            % Windows will either have the WMIC or PowerShell, a PowerShell call seems like a good
            % fallback.
            [wmic_not_available,ignore] = system('wmic /?'); %#ok<ASGLU>
            if wmic_not_available
                call_type_ = 'PS_get_date';
            else
                call_type_ = 'WMIC_sys';
            end
        else
            call_type_ = 'WMIC_sys';
        end
        call_type_ = struct('offline',call_type_,'online','NTP_win');
    end
end
call_type = call_type_;
end
function UTC_epoch_seconds=getUTC_cmd_PS_get_date
% Use Powershell to get the UTC time.
[s,str] = system(['powershell $a=get-date;',...
    '$a.ToUniversalTime().ToString(''yyyyMMddHHmmss'')']); %#ok<ASGLU>
str(str<48 | str>57) = ''; % Remove trailing newline by keeping only digits.
date = mat2cell(str,1,[4 2 2,2 2 2]);date = num2cell(str2double(date));
UTC_epoch_seconds = (datenum(date{:})-datenum(1970,1,1))*24*60*60; %#ok<DATNM>
end
function UTC_epoch_seconds=getUTC_cmd_Unix
[status,str] = system('date +%s'); %#ok<ASGLU>
UTC_epoch_seconds = str2double(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_bat
% If a normal system call would hang, use a temporary bat file.
pausetime = 1;
fn1 = [tempname,'.bat'];fn2=[tempname,'.txt'];
fid = fopen(fn1,'w');
% This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes.
fprintf(fid,'%%systemroot%%\\system32\\wbem\\wmic os get LocalDateTime /value > "%s"\r\nexit',fn2);
fclose(fid);
system(['start /min "" cmd /c "' fn1 '"']);
then = now; %#ok<TNOW1> Store the current time to adjust for the pausetime and the function time.
pause(pausetime) % Wait for the system call to finish.
str = fileread(fn2);
try delete(fn1);catch,end,try delete(fn2);catch,end % Delete temp files.
UTC_epoch_seconds = getUTC_cmd_wmic_parse_str(str)+(now-then)*24*60*60; %#ok<TNOW1>
end
function UTC_epoch_seconds=getUTC_cmd_wmic_sys
% Use a direct system call (instead of a temporary bat file).
[status,str] = system('wmic os get LocalDateTime /value'); %#ok<ASGLU>
% This command returns YYYYMMDDHHMMSS.milliseconds+UTC_Offset_in_minutes.
UTC_epoch_seconds = getUTC_cmd_wmic_parse_str(str);
end
function UTC_epoch_seconds=getUTC_cmd_wmic_parse_str(str)
str = str(str>=43 & str<=57); % Strip irrelevant parts.
date = mat2cell(str(1:21),1,[4 2 2,2 2 2+1+6]);date = str2double(date);
date(5) = date(5)-str2double(str(22:end)); % Add offset.
date = num2cell(date);
UTC_epoch_seconds = (datenum(date{:})-datenum(1970,1,1))*24*60*60; %#ok<DATNM>
end
function UTC_epoch_seconds=getUTC_cmd_NTP_win
% Use a system call to get the time from the Google NTP server.
persistent LDAP2UTC
if isempty(LDAP2UTC),LDAP2UTC = (datenum(1601,1,1)-datenum(1970,1,1))*(24*60*60);end %#ok<DATNM>
[status,str] = system('w32tm /stripchart /computer:time.google.com /dataonly /samples:1 /rdtsc'); %#ok<ASGLU>
t = regexp_outkeys(str,',\s?([0-9]+)','tokens');
FileTime = str2double(t{end});
UTC_epoch_seconds = LDAP2UTC+FileTime/1e7;
end
function UTC_epoch_seconds=getUTC_web
% Read the timestamp from a web server.
% For some reason this fails for old Matlab releases (at least ML6.5-R2013b).

persistent UseWebread % Only search the path once per session.
if isempty(UseWebread)
    try UseWebread = ~isempty(which(func2str(@webread)));catch,UseWebread = false;end
end

% Skip this function if there is no internet connection (3 timeouts will take a lot of time).
if ~isnetavl,UTC_epoch_seconds = [];return,end

for tries=1:3
    try
        if UseWebread
            data = webread('http://www.utctime.net/utc-timestamp');
        else
            % This probably only works for Octave.
            data = urlread('http://www.utctime.net/utc-timestamp'); %#ok<URLRD>
        end
        break
    catch
    end
end
try
    data(data==' ') = '';
    pat = 'vartimestamp=';
    ind1 = strfind(data,pat)+numel(pat);
    ind2 = strfind(data,';')-1;
    ind2(ind2<ind1) = [];
    UTC_epoch_seconds = str2double(data(ind1:ind2(1)));
catch
    UTC_epoch_seconds = [];
end
end
function [f,status]=GetWritableFolder(varargin)
%Return a folder with write permission
%
% If the output folder doesn't already exist, this function will attempt to create it. This
% function should provide a reliable and repeatable location to write files.
%
% Syntax:
%   f = GetWritableFolder
%   [f,status] = GetWritableFolder
%   [__] = GetWritableFolder(Name,Value)
%   [__] = GetWritableFolder(options)
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

[success,options,ME] = GetWritableFolder_parse_inputs(varargin{:});
if ~success
    rethrow(ME)
else
    [ForceStatus,ErrorOnNotFound,root_folder_list] = deal(options.ForceStatus,...
        options.ErrorOnNotFound,options.root_folder_list);
end
root_folder_list{end} = pwd; % Set this default here to avoid storing it in a persistent.
if ForceStatus
    status = ForceStatus;f=fullfile(root_folder_list{status},'PersistentFolder');
    try if ~exist(f,'dir'),makedir(f);end,catch,end
    return
end

% Option 1: use a folder similar to the AddOn Manager.
status = 1;f = root_folder_list{status};
try if ~exist(f,'dir'),makedir(f);end,catch,end
if ~TestFolderWritePermission(f)
    % If the Add-On path is not writable, return the tempdir. It will not be persistent, but it
    % will be writable.
    status = 2;f = root_folder_list{status};
    try if ~exist(f,'dir'),makedir(f);end,catch,end
    if ~TestFolderWritePermission(f)
        % The tempdir should always be writable, but if for some reason it isn't: return the pwd.
        status = 3;f = root_folder_list{status};
    end
end

% Add 'PersistentFolder' to whichever path was determined above.
f = fullfile(f,'PersistentFolder');
try if ~exist(f,'dir'),makedir(f);end,catch,end

if ~TestFolderWritePermission(f)
    % Apparently even the pwd isn't writable, so we will either return an error, or a fail state.
    if ErrorOnNotFound
        error('HJW:GetWritableFolder:NoWritableFolder',...
            'This function was unable to find a folder with write permissions.')
    else
        status = 0;f = '';
    end
end
end
function [success,options,ME]=GetWritableFolder_parse_inputs(varargin)
%Parse the inputs of the GetWritableFolder function
% This function returns a success flag, the parsed options, and an ME struct.
% As input, the options should either be entered as a struct or as Name,Value pairs. Missing fields
% are filled from the default.

% Pre-assign outputs.
success = false;
ME = struct('identifier','','message','');

persistent default
if isempty(default)
    %Set defaults for options.
    default.ForceStatus = false;
    default.ErrorOnNotFound = false;
    default.root_folder_list = {...
        GetPseudoAddonpath;
        fullfile(tempdir,'MATLAB');
        ''};% Overwrite this last element with pwd when called.
end

if nargin==2
    options = default;
    success = true;
    return
end

% Actually parse the Name,Value pairs (or the struct).
[options,replaced] = parse_NameValue(default,varargin{:});

% Test the optional inputs.
for k=1:numel(replaced)
    curr_option = replaced{k};
    item = options.(curr_option);
    ME.identifier = ['HJW:GetWritableFolder:incorrect_input_opt_' lower(curr_option)];
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
            [passed,options.ErrorOnNotFound] = test_if_scalar_logical(item);
            if ~passed
                ME.message = 'ErrorOnNotFound should be either true or false.';
                return
            end
        otherwise
            ME.message = sprintf('Name,Value pair not recognized: %s.',curr_option);
            ME.identifier = 'HJW:GetWritableFolder:incorrect_input_NameValue';
            return
    end
end
success = true;ME = [];
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
    [ignore,appdata] = system('echo %APPDATA%'); %#ok<ASGLU>
    appdata(appdata<14) = ''; % (remove LF/CRLF)
    f = fullfile(appdata,'MathWorks','MATLAB Add-Ons');
else
    [ignore,home_dir] = system('echo $HOME'); %#ok<ASGLU>
    home_dir(home_dir<14) = ''; % (remove LF/CRLF)
    f = fullfile(home_dir,'Documents','MATLAB','Add-Ons');
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
%   tf = ifversion(test,Rxxxxab)
%   tf = ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Input/output arguments:
% tf:
%   If the current version satisfies the test this returns true. This works similar to verLessThan.
% Rxxxxab:
%   A char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the numeric
%   version (e.g. 6.5, 7, or 9.6). Note that 9.10 is interpreted as 9.1 when using numeric input.
% test:
%   A char array containing a logical test. The interpretation of this is equivalent to
%   eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',23.02) returns true only when run on R2023b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
% ifversion('==',9.10) returns true only when run on R2016b (v9.1) not on R2021a (9.10).
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 1.2.1.1                                                       |%
%|  Date:    2023-10-20                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). You can see the full test matrix below.
% Compatibility considerations:
% - This is expected to work on all releases.

if nargin<2 || nargout>1,error('incorrect number of input/output arguments'),end

% The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
% This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
% remove the potential for float rounding errors.
% Store in persistent for fast recall (don't use getpref, as that is slower than generating the
% variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    % Test if Octave is used instead of Matlab.
    octave = exist('OCTAVE_VERSION', 'builtin');
    
    % Get current version number. This code was suggested by Jan on this thread:
    % https://mathworks.com/matlabcentral/answers/1671199#comment_2040389
    v_num = [100, 1] * sscanf(version, '%d.%d', 2);
    
    % Get dictionary to use for ismember.
    v_dict = {...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;
        'R2023a' 914;'R2023b' 2302};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v = 0.1*Rxxxxab+0.9*fixeps(Rxxxxab);v = round(100*v);
        else
            L = ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf = NaN;return
            else
                v = v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v] = deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v = 0.1*v+0.9*fixeps(v);v = round(100*v);
    else
        [test,v] = deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v = 0.1*v+0.9*fixeps(v);v = round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        % Note that this can't distinguish between 9.1 and 9.10, and will the choose the former.
        v = fixeps(Rxxxxab*100);if mod(v,10)==0,v = fixeps(Rxxxxab)*100+mod(Rxxxxab,1)*10;end
    else
        L = ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf = NaN;return
        else
            v = v_dict{L,2};
        end
    end
end
switch test
    case '==', tf = v_num == v;
    case '<' , tf = v_num <  v;
    case '<=', tf = v_num <= v;
    case '>' , tf = v_num >  v;
    case '>=', tf = v_num >= v;
end
end
function val=fixeps(val)
% Round slightly up to prevent rounding errors using fix().
val = fix(val+eps*1e3);
end
function [connected,timing]=isnetavl(use_HTML_test_only)
%Check for an internet connection by pinging Google
%
% Syntax:
%   [connected,timing] = isnetavl(use_HTML_test_only)
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
%|  Version: 2.0.1                                                         |%
%|  Date:    2022-04-20                                                    |%
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

if nargin==0,                                             tf = false;
else,[tf1,tf2]=test_if_scalar_logical(use_HTML_test_only);tf = tf1&&tf2;
end
if tf
    %(the timing is not reliable)
    [connected,timing] = isnetavl___ping_via_html;
    return
end

tf = isnetavl__ICMP_is_blocked;
if isempty(tf)
    % Unable to determine if ping is allowed, the connection must be down.
    connected = 0;
    timing = 0;
else
    if tf
        % Ping is not allowed.
        % (the timing is not reliable)
        [connected,timing] = isnetavl___ping_via_html;
    else
        % Ping is allowed.
        [connected,timing] = isnetavl___ping_via_system;
    end
end
end
function [connected,timing]=isnetavl___ping_via_html
% Ping is blocked by some organizations. As an alternative, the google.com page can be loaded as a
% normal HTML, which should work as well, although it is slower. This also means the ping timing is
% no longer reliable.
persistent UseWebread
if isempty(UseWebread)
    try no_webread = isempty(which(func2str(@webread)));catch,no_webread = true;end
    UseWebread = ~no_webread;
end
try
    then = now; %#ok<TNOW1>
    if UseWebread
        str = webread('http://google.com'); %#ok<NASGU>
    else
        str = urlread('http://google.com'); %#ok<NASGU,URLRD>
    end
    connected = 1;
    timing = (now-then)*24*3600*1000; %#ok<TNOW1>
catch
    connected = 0;
    timing = 0;
end
end
function [connected,timing]=isnetavl___ping_via_system
if ispc
    try
        %                                     8.8.4.4 will also work
        [ignore_output,b] = system('ping -n 1 8.8.8.8');%#ok<ASGLU> ~
        stats = b(strfind(b,' = ')+3);
        stats = stats(1:3); % [sent received lost]
        if ~strcmp(stats,'110')
            error('trigger error')
        else
            % This branch will error for 'destination host unreachable'.
            connected = 1;
            % This assumes there is only one place with '=[digits]ms' in the response, but this
            % code is not language-specific.
            [ind1,ind2] = regexp(b,' [0-9]+ms');
            timing = b((ind1(1)+1):(ind2(1)-2));
            timing = str2double(timing);
        end
    catch
        connected = 0;
        timing = 0;
    end
elseif isunix
    try
        %                                     8.8.4.4 will also work
        [ignore_output,b] = system('ping -c 1 8.8.8.8');%#ok<ASGLU> ~
        ind = regexp(b,', [01] ');
        if b(ind+2)~='1'
            % This branch includes 'destination host unreachable' errors.
            error('trigger error')
        else
            connected = 1;
            % This assumes the first place with '=[digits] ms' in the response contains the ping
            % timing. This code is not language-specific.
            [ind1,ind2] = regexp(b,'=[0-9.]+ ms');
            timing = b((ind1(1)+1):(ind2(1)-2));
            timing = str2double(timing);
        end
    catch
        connected = 0;
        timing = 0;
    end
else
    error('How did you even get Matlab to work?')
end
end
function [tf,connected,timing]=isnetavl__ICMP_is_blocked
% Check if ICMP 0/8/11 is blocked.
%
% tf is empty if both methods fail.

persistent output
if ~isempty(output)
    tf = output;return
end

% First check if ping works.
[connected,timing] = isnetavl___ping_via_system;
if connected
    % Ping worked and there is an internet connection.
    output = false;
    tf = false;
    return
end

% There are two options: no internet connection, or ping is blocked.
[connected,timing] = isnetavl___ping_via_html;
if connected
    % There is an internet connection, therefore ping must be blocked.
    output = true;
    tf = true;
    return
end

% Both methods failed, internet is down. Leave the value of tf (and the persistent variable) set to
% empty so it is tried next time.
tf = [];
end
function varargout=makedir(d)
% Wrapper function to account for old Matlab releases, where mkdir fails if the parent folder does
% not exist. This function will use the legacy syntax for those releases.
if exist(d,'dir'),return,end % Take a shortcut if the folder already exists.
persistent IsLegacy
if isempty(IsLegacy)
    % The behavior changed after R14SP3 and before R2007b, but since the legacy syntax will still
    % work in later releases there isn't really a reason to pinpoint the exact release.
    IsLegacy = ifversion('<','R2007b','Octave','<',0);
end
varargout = cell(1,nargout);
if IsLegacy
    [d_parent,d_target] = fileparts(d);
    [varargout{:}] = mkdir(d_parent,d_target);
else
    [varargout{:}] = mkdir(d);
end
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
    v = version;ind=[strfind(v,'.') numel(v)];
    v = sprintf('%02d.%02d',str2double({...
        v(1:(ind(1)-1)           ) ,...
        v(  (ind(1)+1):(ind(2)-1)) }));
    v = ['v' strrep(v,'.','_')];
    if ~exist('OCTAVE_VERSION', 'builtin')
        runtime = 'MATLAB';
        type = computer;
    else
        runtime = 'OCTAVE';
        arch = computer;arch = arch(1:(min(strfind(arch,'-'))-1));
        if ispc
            if strcmp(arch,'x86_64')  ,type =  'win_64';
            elseif strcmp(arch,'i686'),type =  'win_i686';
            elseif strcmp(arch,'x86') ,type =  'win_x86';
            else                      ,type = ['win_' arch];
            end
        elseif isunix && ~ismac % Essentially this is islinux
            if strcmp(arch,'i686')      ,type =  'lnx_i686';
            elseif strcmp(arch,'x86_64'),type =  'lnx_64';
            else                        ,type = ['lnx_' arch];
            end
        elseif ismac
            if strcmp(arch,'x86_64'),type =  'mac_64';
            else                    ,type = ['mac_' arch];
            end
        end
    end
    type = strrep(strrep(type,'.',''),'-','');
    append = cell(2,1);
    append{1} = ['_' runtime '_' v '_' type];
    append{2} = [append{1} '.' mexext];
end

try % Test if fun_name is a valid name.
    if ~isvarname(fun_name),error('trigger catch block'),end
catch
    error('HJW:mexname:InvalidName',...
        'The provided input can''t be a function name')
end

mex_filename = [fun_name append{2}];
fun_name = [fun_name append{1}];
end
function [opts,replaced]=parse_NameValue(default,varargin)
%Match the Name,Value pairs to the default option, attempting to autocomplete
%
% The autocomplete ignores incomplete names, case, underscores, and dashes, as long as a unique
% match can be found.
%
% The first output is a struct with the same fields as the first input, with field contents
% replaced according to the supplied options struct or Name,Value pairs.
% The second output is a cellstr containing the field names that have been set.
%
% If this fails to find a match, this will throw an error with the offending name as the message.
%
% If there are multiple occurrences of a Name, only the last Value will be returned. This is the
% same as Matlab internal functions like plot. GNU Octave also has this behavior.
%
% If a struct array is provided, only the first element will be used. An empty struct array will
% trigger an error.

switch numel(default)
    case 0
        error('parse_NameValue:MixedOrBadSyntax',...
            'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
    case 1
        % Do nothing.
    otherwise
        % If this is a struct array, explicitly select the first element.
        default=default(1);
end

% Create default output and return if no other inputs exist.
opts = default;replaced = {};
if nargin==1,return,end

% Unwind an input struct to Name,Value pairs.
try
    struct_input = numel(varargin)==1 && isa(varargin{1},'struct');
    NameValue_input = mod(numel(varargin),2)==0 && all(...
        cellfun('isclass',varargin(1:2:end),'char'  ) | ...
        cellfun('isclass',varargin(1:2:end),'string')   );
    if ~( struct_input || NameValue_input )
        error('trigger')
    end
    if nargin==2
        Names = fieldnames(varargin{1});
        Values = struct2cell(varargin{1});
    else
        % Wrap in cellstr to account for strings (this also deals with the fun(Name=Value) syntax).
        Names = cellstr(varargin(1:2:end));
        Values = varargin(2:2:end);
    end
    if ~iscellstr(Names),error('trigger');end %#ok<ISCLSTR>
catch
    % If this block errors, that is either because a missing Value with the Name,Value syntax, or
    % because the struct input is not a struct, or because an attempt was made to mix the two
    % styles. In future versions of this functions an effort might be made to handle such cases.
    error('parse_NameValue:MixedOrBadSyntax',...
        'Optional inputs must be entered as Name,Value pairs or as a scalar struct.')
end

% The fieldnames will be converted to char matrices in the section below. First an exact match is
% tried, then a case-sensitive (partial) match, then ignoring case, followed by ignoring any
% underscores, and lastly ignoring dashes.
default_Names = fieldnames(default);
Names_char    = cell(1,4);
Names_cell{1} = default_Names;
Names_cell{2} = lower(Names_cell{1});
Names_cell{3} = strrep(Names_cell{2},'_','');
Names_cell{4} = strrep(Names_cell{3},'-','');

% Allow spaces by replacing them with underscores.
Names = strrep(Names,' ','_');

% Attempt to match the names.
replaced = false(size(default_Names));
for n=1:numel(Names)
    name = Names{n};
    
    % Try a case-sensitive match.
    [match_idx,Names_char{1}] = parse_NameValue__find_match(Names_char{1},Names_cell{1},name);
    
    % Try a case-insensitive match.
    if numel(match_idx)~=1
        name = lower(name);
        [match_idx,Names_char{2}] = parse_NameValue__find_match(Names_char{2},Names_cell{2},name);
    end
    
    % Try a case-insensitive match ignoring underscores.
    if numel(match_idx)~=1
        name = strrep(name,'_','');
        [match_idx,Names_char{3}] = parse_NameValue__find_match(Names_char{3},Names_cell{3},name);
    end
    
    % Try a case-insensitive match ignoring underscores and dashes.
    if numel(match_idx)~=1
        name = strrep(name,'-','');
        [match_idx,Names_char{4}] = parse_NameValue__find_match(Names_char{4},Names_cell{4},name);
    end
    
    if numel(match_idx)~=1
        error('parse_NameValue:NonUniqueMatch',Names{n})
    end
    
    % Store the Value in the output struct and mark it as replaced.
    opts.(default_Names{match_idx}) = Values{n};
    replaced(match_idx)=true;
end
replaced = default_Names(replaced);
end
function [match_idx,Names_char]=parse_NameValue__find_match(Names_char,Names_cell,name)
% Try to match the input field to the fields of the struct.

% First attempt an exact match.
match_idx = find(ismember(Names_cell,name));
if numel(match_idx)==1,return,end

% Only spend time building the char array if this point is reached.
if isempty(Names_char),Names_char = parse_NameValue__name2char(Names_cell);end

% Since the exact match did not return a unique match, attempt to match the start of each array.
% Select the first part of the array. Since Names is provided by the user it might be too long.
tmp = Names_char(:,1:min(end,numel(name)));
if size(tmp,2)<numel(name)
    tmp = [tmp repmat(' ', size(tmp,1) , numel(name)-size(tmp,2) )];
end

% Find the number of non-matching characters on every row. The cumprod on the logical array is
% to make sure that only the starting match is considered.
non_matching = numel(name)-sum(cumprod(double(tmp==repmat(name,size(tmp,1),1)),2),2);
match_idx = find(non_matching==0);
end
function Names_char=parse_NameValue__name2char(Names_char)
% Convert a cellstr to a padded char matrix.
len = cellfun('prodofsize',Names_char);maxlen = max(len);
for n=find(len<maxlen).' % Pad with spaces where needed
    Names_char{n}((end+1):maxlen) = ' ';
end
Names_char = vertcat(Names_char{:});
end
function varargout=regexp_outkeys(str,expression,varargin)
%Regexp with outkeys in old releases
%
% On older versions of Matlab the regexp function did not allow you to specify the output keys.
% This function has an implementation of the 'split', 'match', and 'tokens' output keys, so they
% can be used on any version of Matlab or GNU Octave. The 'start' and 'end' output keys were
% already supported as trailing outputs, but are now also explictly supported.
% On releases where this is possible, the builtin is called.
%
% Syntax:
%   out = regexp_outkeys(str,expression,outkey);
%   [out1,...,outN] = regexp_outkeys(str,expression,outkey1,...,outkeyN);
%   [___,startIndex] = regexp_outkeys(___);
%   [___,startIndex,endIndex] = regexp_outkeys(___);
%
% Example:
%  str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
%  words = regexp_outkeys(str,'[ 0-9.]+','split')
%  numbers = regexp_outkeys(str,'[0-9.]*','match')
%  [white,end1,start,end2] = regexp_outkeys(str,' ','match','end')
%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%|                                                                         |%
%|  Version: 2.0.0.1                                                       |%
%|  Date:    2023-09-12                                                    |%
%|  Author:  H.J. Wisselink                                                |%
%|  Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 ) |%
%|  Email = 'h_j_wisselink*alumnus_utwente_nl';                            |%
%|  Real_email = regexprep(Email,{'*','_'},{'@','.'})                      |%
%|                                                                         |%
%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%/%
%
% Tested on several versions of Matlab (ML 6.5 and onward) and Octave (4.4.1 and onward), and on
% multiple operating systems (Windows/Ubuntu/MacOS). You can see the full test matrix below.
% Compatibility considerations:
% - Only the 'match', 'split', 'tokens', 'start', and 'end' options are supported. The additional
%   options provided by regexp are not implemented.
% - Cell array input is not supported.

if nargin<2
    error('HJW:regexp_outkeys:SyntaxError',...
        'No supported syntax used: at least 3 inputs expected.')
    % As an undocumented feature this can also return s1,s2 without any outkeys specified.
end
if ~(ischar(str) && ischar(expression))
    % The extra parameters in varargin are checked inside the loop.
    error('HJW:regexp_outkeys:InputError',...
        'All inputs must be char vectors.')
end
if nargout>nargin
    error('HJW:regexp_outkeys:ArgCount',...
        'Incorrect number of output arguments. Check syntax.')
end

persistent legacy errorstr KeysDone__template
if isempty(legacy)
    % The legacy struct contains the implemented options as field names. It is used in the error
    % message.
    % It is assumed that all Octave versions later than 4.0 support the expanded output, and all
    % earlier versions do not, even if it is very likely most versions will support it.
    
    % Create these as fields so they show up in the error message.
    legacy.start = true;
    legacy.end = true;
    
    % The switch to find matches was introduced in R14 (along with the 'tokenExtents', 'tokens'
    % and 'names' output switches).
    legacy.match = ifversion('<','R14','Octave','<',4);
    legacy.tokens = legacy.match;
    
    % The split option was introduced in R2007b.
    legacy.split = ifversion('<','R2007b','Octave','<',4);
    
    fn = fieldnames(legacy);
    errorstr = ['Extra regexp output type not implemented,',char(10),'only the following',...
        ' types are implemented:',char(10),sprintf('%s, ',fn{:})]; %#ok<CHARTEN>
    errorstr((end-1):end) = ''; % Remove trailing ', '
    
    legacy.any = legacy.match || legacy.split || legacy.tokens;
    
    % Copy all relevant fields and set them to false.
    KeysDone__template = legacy;
    for x=fieldnames(KeysDone__template).',KeysDone__template.(x{1}) = false;end
end

if legacy.any || nargin==2 || any(ismember(lower(varargin),{'start','end'}))
    % Determine s1, s2, and TokenIndices only once for the legacy implementations.
    [s1,s2,TokenIndices] = regexp(str,expression);
end

if nargin==2
    varargout = {s1,s2,TokenIndices};return
end

% Pre-allocate output.
varargout = cell(size(varargin));
done = KeysDone__template; % Keep track of outkey in case of repeats.
% On some releases the Matlab is not convinced split is a variable, so explicitly assign it here.
split = [];
for param=1:(nargin-2)
    if ~ischar(varargin{param})
        error('HJW:regexp_outkeys:InputError',...
            'All inputs must be char vectors.')
    end
    switch lower(varargin{param})
        case 'match'
            if done.match,varargout{param} = match;continue,end
            if legacy.match
                % Legacy implementation.
                match = cell(1,numel(s1));
                for n=1:numel(s1)
                    match{n} = str(s1(n):s2(n));
                end
            else
                [match,s1,s2] = regexp(str,expression,'match');
            end
            varargout{param} = match;done.match = true;
        case 'split'
            if done.split,varargout{param} = split;continue,end
            if legacy.split
                % Legacy implementation.
                split = cell(1,numel(s1)+1);
                start_index = [s1 numel(str)+1];
                stop_index = [0 s2];
                for n=1:numel(start_index)
                    split{n} = str((stop_index(n)+1):(start_index(n)-1));
                    if numel(split{n})==0,split{n} = char(ones(0,0));end
                end
            else
                [split,s1,s2] = regexp(str,expression,'split');
            end
            varargout{param}=split;done.split = true;
        case 'tokens'
            if done.tokens,varargout{param} = tokens;continue,end
            if legacy.tokens
                % Legacy implementation.
                tokens = cell(numel(TokenIndices),0);
                for n=1:numel(TokenIndices)
                    if size(TokenIndices{n},2)~=2
                        % No actual matches for the tokens.
                        tokens{n} = cell(1,0);
                    else
                        for m=1:size(TokenIndices{n},1)
                            tokens{n}{m} = str(TokenIndices{n}(m,1):TokenIndices{n}(m,2));
                        end
                    end
                end
            else
                [tokens,s1,s2] = regexp(str,expression,'tokens');
            end
            varargout{param} = tokens;done.tokens = true;
        case 'start'
            varargout{param} = s1;
        case 'end'
            varargout{param} = s2;
        otherwise
            error('HJW:regexp_outkeys:NotImplemented',errorstr)
    end
end
if nargout>param
    varargout(param+[1 2]) = {s1,s2};
end
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
% The char and string test are not case sensitive.
% (use the first output to trigger an input error, use the second as the parsed input)
%
%  Allowed values:
% - true or false
% - 1 or 0
% - 'on' or 'off'
% - "on" or "off"
% - matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
% - 'enable' or 'disable'
% - 'enabled' or 'disabled'
persistent states
if isempty(states)
    states = {...
        true,false;...
        1,0;...
        'true','false';...
        '1','0';...
        'on','off';...
        'enable','disable';...
        'enabled','disabled'};
    % We don't need string here, as that will be converted to char.
end

% Treat this special case.
if isa(val,'matlab.lang.OnOffSwitchState')
    isLogical = true;val = logical(val);return
end

% Convert a scalar string to char and return an error state for non-scalar strings.
if isa(val,'string')
    if numel(val)~=1,isLogical = false;return
    else            ,val = char(val);
    end
end

% Convert char/string to lower case.
if isa(val,'char'),val = lower(val);end

% Loop through all possible options.
for n=1:size(states,1)
    for m=1:2
        if isequal(val,states{n,m})
            isLogical = true;
            val = states{1,m}; % This selects either true or false.
            return
        end
    end
end

% Apparently there wasn't any match, so return the error state.
isLogical = false;
end
function tf=TestFolderWritePermission(f)
%Returns true if the folder exists and allows Matlab to write files.
% An empty input will generally test the pwd.
%
% examples:
%   fn='foo.txt';if ~TestFolderWritePermission(fileparts(fn)),error('can''t write!'),end

if ~( isempty(f) || exist(f,'dir') )
    tf = false;return
end

fn = '';
while isempty(fn) || exist(fn,'file')
    % Generate a random file name, making sure not to overwrite any existing file.
    % This will try to create a file without an extension.
    [ignore,fn] = fileparts(tmpname('write_permission_test_','.txt')); %#ok<ASGLU>
    fn = fullfile(f,fn);
end
try
    % Test write permission.
    fid = fopen(fn,'w');fprintf(fid,'test');fclose(fid);
    delete(fn);
    tf = true;
catch
    % Attempt to clean up.
    if exist(fn,'file'),try delete(fn);catch,end,end
    tf = false;
end
end
function str=tmpname(StartFilenameWith,ext)
% Inject a string in the file name part returned by the tempname function.
% This is equivalent to the line below:
% str = fullfile(tempdir,[StartFilenameWith '_' strrep(tempname,tempdir,'') ext])
if nargin<1,StartFilenameWith = '';end
if ~isempty(StartFilenameWith),StartFilenameWith = [StartFilenameWith '_'];end
if nargin<2,ext='';else,if ~strcmp(ext(1),'.'),ext = ['.' ext];end,end
str = tempname;
[p,f] = fileparts(str);
str = fullfile(p,[StartFilenameWith f ext]);
end
function out=WinVer
% This returns the main Windows version number (5 for XP, 6 for Vista, etc).
% It will return an empty array for non-Windows or in case of errors.
persistent persistent_val
if ~ispc,out = [];return,end
if isempty(persistent_val)
    try
        [status,str] = system('ver'); %#ok<ASGLU>
        args = {'[^0-9]*(\d*).*','$1','tokenize'};
        if ifversion('>=',7,'Octave','<',0)
            args(end) = []; % The 'tokenize' option became the default in R14 (v7).
        end
        persistent_val = str2double(regexprep(str,args{:}));
    catch
        try
            [status,str] = system('systeminfo'); %#ok<ASGLU>
            ind1 =  1+strfind(str,':');ind1 = ind1(3);
            ind2 = -1+strfind(str,'.');ind2(ind2<ind1) = [];
            persistent_val = str2double(str(ind1:ind2(1)));
        catch
            persistent_val = [];
        end
    end
end
out = persistent_val;
end

