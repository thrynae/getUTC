function pass_part_fail=aaa___getUTC___test(varargin)
%Test function for getUTC (an error for the web is expected for Matlab releases without websave)
%
%Pass:    passes all tests, failures for C method when no compiler is present
%Partial: web method fails (only when expected)
%Fail:    unexpected failures
pass_part_fail = 'pass';
ME = [];

fprintf('Warming up (this may take a long time (>10sec) the first time)...')
getUTC(1); % This will trigger a compile if necessary.
fprintf(' ... ')
getUTC(3); % This will ensure the call type is already determined.
fprintf('Warm-up completed.\n')

t1 = now;                        %#ok<TNOW1>
C = getUTC(1);
t1 = (now-t1)*60*60*24;t2 = now; %#ok<TNOW1>
web =getUTC(2);
t2 = (now-t2)*60*60*24;t3 = now; %#ok<TNOW1>
cmd = getUTC(3);
t3 = (now-t3)*60*60*24;          %#ok<TNOW1>

try
    if isempty(ME)
        ValidateOutput_C(C)
    end
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
end
try
    if isempty(ME)
        pass_part_fail = ValidateOutput_web(web,pass_part_fail);
    end
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
end
try
    if isempty(ME)
        ValidateOutput_cmd(cmd)
    end
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
end
try
    TestEdgeCases_cmd
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
end

SelfTestFailMessage = '';
% Run the self-validator function(s).
SelfTestFailMessage = [SelfTestFailMessage SelfTest__regexp_outkeys];
checkpoint('read');
if ~isempty(SelfTestFailMessage) || ~isempty(ME)
    if nargout==1
        pass_part_fail = 'fail';
    else
        if ~isempty(SelfTestFailMessage)
            error('Self-validator functions returned these error(s):\n%s',SelfTestFailMessage)
        else
            checkpoint('aaa___getUTC___test','char2cellstr','get_trace')
            trace = char2cellstr(get_trace(0,ME.stack));
            error('Test failed with error message:%s',...
                sprintf('\n    %s',ME.message,trace{:}))
        end
    end
else
    fprintf('C function took %.2f seconds\n',t1)
    fprintf('web function took %.2f seconds\n',t2)
    fprintf('system call took %.2f seconds\n',t3)
end
disp(['tester function ' mfilename ' finished '])
if nargout==0,clear,end
end
function ValidateOutput_C(C)
    if isempty(C)
        %Check if there is compiler in the first place.
        checkpoint('aaa___getUTC___test','CheckMexCompilerExistence')
        if CheckMexCompilerExistence
            error('C function failed')
        else
            warning('C function failed [no compiler found]')
        end
    end
end
function pass_part_fail=ValidateOutput_web(web,pass_part_fail)
checkpoint('aaa___getUTC___test','isnetavl')
if ~isnetavl,return;end % The web method can't be tested without a working internet connection.

%Incorrect output (e.g. for v7.10) doesn't matter.
v=version;ind=strfind(v,'.');v=str2double(v(1:(ind(2)-1)));
isOctave=exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isempty(web)
    if ( ~isOctave && ( v<8.4 || ( v>=9.1 && v<=9.6 ) ) ) ...
            || ( isOctave && ismac )
        % Known issue: Matlab 6.5 up to R2013b fail here, let's assume it works again in R2014b, as
        % that was the point where websave was introduced.
        % Later releases (R2016b-R2019a) have also stopped working due to a missing/lapsed CA
        % certificate which prevents SSL connections with this server.
        pass_part_fail='partial: web failed';
    else
        error('web function failed\n')
    end
elseif ~isOctave && v<8
    error('web function did **not** fail as expected\n')
end
end
function ValidateOutput_cmd(cmd)
if isempty(cmd)
    [ignore,call_type] = getUTC_cmd; %#ok<ASGLU>
    error('system call failed with call type ''%s''',call_type)
end
end
function TestEdgeCases_cmd
% Test special input cases. The internal function should never throw an error.

% Force assuming online.
try ME=[]; %#ok<NASGU>
    if isempty(getUTC_cmd(true))
        error('Forced online failed')
    end
catch ME;if isempty(ME),ME = lasterror;end %#ok<LERR>
    checkpoint('aaa___getUTC___test','isnetavl')
    if isnetavl,rethrow(ME),end
end

% Force assuming offline.
if isempty(getUTC_cmd(false))
    error('Forced offline failed')
end

% Use an incorrect call type name.
if ~isempty(getUTC_cmd([],'foobar'))
    error('Incorrect call type not caught.')
end
end
function out=bsxfun_plus(in1,in2)
% Implicit expansion for plus(), but without any input validation.
persistent type
if isempty(type)
    checkpoint('bsxfun_plus','hasFeature')
    type = ...
        double(hasFeature('ImplicitExpansion')) + ...
        double(hasFeature('bsxfun'));
end
if type==2
    % Implicit expansion is available.
    out = in1+in2;
elseif type==1
    % Implicit expansion is only available with bsxfun.
    out = bsxfun(@plus,in1,in2);
else
    % No implicit expansion, expand explicitly.
    % Determine size and find non-singleton dimensions.
    sz1 = ones(1,max(ndims(in1),ndims(in2)));
    sz2 = sz1;
    sz1(1:ndims(in1)) = size(in1);
    sz2(1:ndims(in2)) = size(in2);
    L = sz1~=1 & sz2~=1;
    if ~isequal(sz1(L),sz2(L))
        error('HJW:bsxfun_plus:arrayDimensionsMustMatch',...
            'Non-singleton dimensions of the two input arrays must match each other.')
    end
    if min([sz1 sz2])==0
        % Construct an empty array of the correct size.
        sz1(sz1==0) = inf;sz2(sz2==0) = inf;
        sz = max(sz1,sz2);
        sz(isinf(sz)) = 0;
        % Create an array and cast it to the correct type.
        out = feval(str2func(class(in1)),zeros(sz));
        return
    end
    in1 = repmat(in1,max(1,sz2./sz1));
    in2 = repmat(in2,max(1,sz1./sz2));
    out = in1+in2;
end
end
function c=char2cellstr(str,LineEnding)
% Split char or uint32 vector to cell (1 cell element per line). Default splits are for CRLF/CR/LF.
% The input data type is preserved.
%
% Since the largest valid Unicode codepoint is 0x10FFFF (i.e. 21 bits), all values will fit in an
% int32 as well. This is used internally to deal with different newline conventions.
%
% The second input is a cellstr containing patterns that will be considered as newline encodings.
% This will not be checked for any overlap and will be processed sequentially.

returnChar = isa(str,'char');
str = int32(str); % Convert to signed, this should not crop any valid Unicode codepoints.

if nargin<2
    % Replace CRLF, CR, and LF with -10 (in that order). That makes sure that all valid encodings
    % of newlines are replaced with the same value. This should even handle most cases of files
    % that mix the different styles, even though such mixing should never occur in a properly
    % encoded file. This considers LFCR as two line endings.
    if any(str==13)
        checkpoint('char2cellstr','PatternReplace')
        str = PatternReplace(str,int32([13 10]),int32(-10));
        str(str==13) = -10;
    end
    str(str==10) = -10;
else
    for n=1:numel(LineEnding)
        checkpoint('char2cellstr','PatternReplace')
        str = PatternReplace(str,int32(LineEnding{n}),int32(-10));
    end
end

% Split over newlines.
newlineidx = [0 find(str==-10) numel(str)+1];
c=cell(numel(newlineidx)-1,1);
for n=1:numel(c)
    s1 = (newlineidx(n  )+1);
    s2 = (newlineidx(n+1)-1);
    c{n} = str(s1:s2);
end

% Return to the original data type.
if returnChar
    for n=1:numel(c),c{n} =   char(c{n});end
else
    for n=1:numel(c),c{n} = uint32(c{n});end
end
end
function tf=CharIsUTF8
% This provides a single place to determine if the runtime uses UTF-8 or UTF-16 to encode chars.
% The advantage is that there is only 1 function that needs to change if and when Octave switches
% to UTF-16. This is unlikely, but not impossible.
persistent persistent_tf
if isempty(persistent_tf)
    checkpoint('CharIsUTF8','ifversion')
    if ifversion('<',0,'Octave','>',0)
        % Test if Octave has switched to UTF-16 by looking if the Euro symbol is losslessly encoded
        % with char.
        % Because we will immediately reset it, setting the state for all warnings to off is fine.
        w = struct('w',warning('off','all'));[w.msg,w.ID] = lastwarn;
        persistent_tf = ~isequal(8364,double(char(8364)));
        warning(w.w);lastwarn(w.msg,w.ID); % Reset warning state.
    else
        persistent_tf = false;
    end
end
tf = persistent_tf;
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
    checkpoint('CheckMexCompilerExistence','CreatePathFolder__CheckMexCompilerExistence_persistent')
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
    
    checkpoint('CheckMexCompilerExistence','CheckMexCompilerExistence_persistent')
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
    checkpoint('CreatePathFolder__CheckMexCompilerExistence_persistent','GetWritableFolder')
    p = fullfile(GetWritableFolder,'FileExchange','CheckMexCompilerExistence');
    if isempty(strfind([path ';'],[p ';'])) %#ok<STREMP>
        % This means f is not on the path.
        checkpoint('CreatePathFolder__CheckMexCompilerExistence_persistent','makedir')
        if ~exist(p,'dir'),makedir(p);end
        addpath(p,'-end');
    end
catch
    ME = struct('identifier','HJW:CheckMexCompilerExistence:PathFolderFail',...
        'message','Creating a folder on the path to store the compiled function and flag failed.');
end
end
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers = 1;end
if nargin<2, stack = dbstack;end
stack(1:skip_layers) = [];

% Parse the ML6.5 style of dbstack (the name field includes full file location).
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp = stack(n).name;
        if strcmp(tmp(end),')')
            % Internal function.
            ind = strfind(tmp,'(');
            name = tmp( (ind(end)+1):(end-1) );
            file = tmp(1:(ind(end)-2));
        else
            file = tmp;
            [ignore,name] = fileparts(tmp); %#ok<ASGLU>
        end
        [ignore,stack(n).file] = fileparts(file); %#ok<ASGLU>
        stack(n).name = name;
    end
end

% Parse Octave style of dbstack (the file field includes full file location).
checkpoint('get_trace','ifversion')
persistent isOctave,if isempty(isOctave),isOctave=ifversion('<',0,'Octave','>',0);end
if isOctave
    for n=1:numel(stack)
        [ignore,stack(n).file] = fileparts(stack(n).file); %#ok<ASGLU>
    end
end

% Create the char array with a (potentially) modified stack.
s = stack;
c1 = '>';
str = cell(1,numel(s)-1);
for n=1:numel(s)
    [ignore_path,s(n).file,ignore_ext] = fileparts(s(n).file); %#ok<ASGLU>
    if n==numel(s),s(n).file = '';end
    if strcmp(s(n).file,s(n).name),s(n).file = '';end
    if ~isempty(s(n).file),s(n).file = [s(n).file '>'];end
    str{n} = sprintf('%c In %s%s (line %d)\n',c1,s(n).file,s(n).name,s(n).line);
    c1 = ' ';
end
str = horzcat(str{:});
end
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

if nargin==0
    % Normal flow: first try the cmd method, then the C method, then the web method.
    checkpoint('getUTC','getUTC_cmd')
    UTC_epoch_seconds = getUTC_cmd;
    if isempty(UTC_epoch_seconds)
        checkpoint('getUTC','getUTC_c')
        UTC_epoch_seconds = getUTC_c;
    end
    if isempty(UTC_epoch_seconds)
        checkpoint('getUTC','getUTC_web')
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
            checkpoint('getUTC','getUTC_c')
            UTC_epoch_seconds = getUTC_c(false);
        case 2
            checkpoint('getUTC','getUTC_web')
            UTC_epoch_seconds = getUTC_web;
        case 3
            checkpoint('getUTC','getUTC_cmd')
            UTC_epoch_seconds = getUTC_cmd;
        otherwise
            if isa(override,'char')
                checkpoint('getUTC','getUTC_cmd')
                UTC_epoch_seconds = getUTC_cmd(override);
            else
                error('non-implemented override')
            end
    end
end
UTC_offset = UTC_epoch_seconds/(24*60*60);
atomTime = UTC_offset+datenum(1970,1,1); %#ok<DATNM>
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
    checkpoint('getUTC_c','GetWritableFolder')
    tempdir_f = fullfile(GetWritableFolder,'FileExchange','getUTC');
    try
        if isempty(strfind([path ';'],[tempdir_f ';'])) %#ok<STREMP>
            % This means f is not on the path.
            checkpoint('getUTC_c','makedir')
            if ~exist(tempdir_f,'dir'),makedir(tempdir_f);end
            addpath(tempdir_f,'-end');
        end
    catch
    end
    
    checkpoint('getUTC_c','mexname')
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
    checkpoint('getUTC_c','CheckMexCompilerExistence')
    if ~CheckMexCompilerExistence
        Compile_attempts_remaining = 0;
        UTC_epoch_seconds = [];return
    end
    
    % Try building missing C file.
    Compile_attempts_remaining = Compile_attempts_remaining-1;
    if Compile_attempts_remaining<0 % Don't endlessly try to compile.
        UTC_epoch_seconds = [];return
    end
    
    checkpoint('getUTC_c','TestFolderWritePermission')
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
        checkpoint('getUTC_cmd','isnetavl')
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
            checkpoint('getUTC_cmd______________________________line043','CoverTest')
            UTC_epoch_seconds = getUTC_cmd_Unix;
        case 'NTP_win'
            checkpoint('getUTC_cmd______________________________line046','CoverTest')
            UTC_epoch_seconds = getUTC_cmd_NTP_win;
        case 'WMIC_sys'
            checkpoint('getUTC_cmd______________________________line049','CoverTest')
            UTC_epoch_seconds = getUTC_cmd_wmic_sys;
        case 'WMIC_bat'
            checkpoint('getUTC_cmd______________________________line052','CoverTest')
            UTC_epoch_seconds = getUTC_cmd_wmic_bat;
        case 'PS_get_date'
            checkpoint('getUTC_cmd______________________________line055','CoverTest')
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
        checkpoint('getUTC_cmd','WinVer')
        if WinVer<=5
            % The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10,
            % I'm making the assumption that the OS is the cause of the hang.
            call_type_ = 'WMIC_bat';
            checkpoint('getUTC_cmd','WinVer')
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
checkpoint('getUTC_cmd','regexp_outkeys')
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

checkpoint('getUTC_web','isnetavl')
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
function[v000,v001]=GetWritableFolder(varargin),[v002,v003,v004]=...
GetWritableFolder_f01(varargin{:});if~v002,rethrow(v004),else,[v005,v006,v007]=...
deal(v003.ForceStatus,v003.ErrorOnNotFound,v003.root_folder_list);end,v007{end}=pwd;if v005,...
v001=v005;v000=fullfile(v007{v001},'PersistentFolder');try if~exist(v000,'dir'),...
GetWritableFolder_f02(v000);end,catch,end,return,end,v001=1;v000=v007{v001};try if~exist(v000,...
'dir'),GetWritableFolder_f02(v000);end,catch,end,if~GetWritableFolder_f08(v000),v001=2;v000=...
v007{v001};try if~exist(v000,'dir'),GetWritableFolder_f02(v000);end,catch,end,...
if~GetWritableFolder_f08(v000),v001=3;v000=v007{v001};end,end,v000=fullfile(v000,...
'PersistentFolder');try if~exist(v000,'dir'),GetWritableFolder_f02(v000);end,catch,end,...
if~GetWritableFolder_f08(v000),if v006,error('HJW:GetWritableFolder:NoWritableFolder',...
'This function was unable to find a folder with write permissions.'),else,v001=0;v000='';end,...
end,end
function v002=GetWritableFolder_f00,if ispc,[v000,v001]=system('echo %APPDATA%');v001(v001<14)=...
'';v002=fullfile(v001,'MathWorks','MATLAB Add-Ons');else,[v000,v003]=system('echo $HOME');
v003(v003<14)='';v002=fullfile(v003,'Documents','MATLAB','Add-Ons');end,end
function[v000,v001,v002]=GetWritableFolder_f01(varargin),v000=false;v002=struct('identifier','',...
'message','');persistent v003,if isempty(v003),v003.ForceStatus=false;v003.ErrorOnNotFound=...
false;v003.root_folder_list={GetWritableFolder_f00;fullfile(tempdir,'MATLAB');''};end,if ...
nargin==2,v001=v003;v000=true;return,end,[v001,v004]=GetWritableFolder_f04(v003,varargin{:});
for v005=1:numel(v004),v006=v004{v005};v007=v001.(v006);v002.identifier=...
['HJW:GetWritableFolder:incorrect_input_opt_' lower(v006)];switch v006,case'ForceStatus',try ...
if~isa(v003.root_folder_list{v007},'char'),...
error('the indexing must have failed, trigger error'),end,catch,v002.message=...
sprintf('Invalid input: expected a scalar integer between 1 and %d.',...
numel(v003.root_folder_list));return,end,case'ErrorOnNotFound',[v008,v001.ErrorOnNotFound]=...
GetWritableFolder_f07(v007);if~v008,v002.message=...
'ErrorOnNotFound should be either true or false.';return,end,otherwise,v002.message=...
sprintf('Name,Value pair not recognized: %s.',v006);v002.identifier=...
'HJW:GetWritableFolder:incorrect_input_NameValue';return,end,end,v000=true;v002=[];end
function varargout=GetWritableFolder_f02(v000),if exist(v000,'dir'),return,end,persistent v001,...
if isempty(v001),v001=GetWritableFolder_f09('<','R2007b','Octave','<',0);end,varargout=cell(1,...
nargout);if v001,[v002,v003]=fileparts(v000);[varargout{:}]=mkdir(v002,v003);else,...
[varargout{:}]=mkdir(v000);end,end
function v000=GetWritableFolder_f03(v000),v000=fix(v000+eps*1e3);end
function[v000,v001]=GetWritableFolder_f04(v002,varargin),switch numel(v002),case 0,...
error('parse_NameValue:MixedOrBadSyntax',...
'Optional inputs must be entered as Name,Value pairs or as a scalar struct.'),case 1,otherwise,...
v002=v002(1);end,v000=v002;v001={};if nargin==1,return,end,try v003=numel(varargin)==...
1&&isa(varargin{1},'struct');v004=mod(numel(varargin),2)==0&&all(cellfun('isclass',...
varargin(1:2:end),'char')|cellfun('isclass',varargin(1:2:end),'string'));if~(v003||v004),...
error('trigger'),end,if nargin==2,v005=fieldnames(varargin{1});v006=struct2cell(varargin{1});
else,v005=cellstr(varargin(1:2:end));v006=varargin(2:2:end);end,if~iscellstr(v005),...
error('trigger');end,catch,error('parse_NameValue:MixedOrBadSyntax',...
'Optional inputs must be entered as Name,Value pairs or as a scalar struct.'),end,v007=...
fieldnames(v002);v008=cell(1,4);v009{1}=v007;v009{2}=lower(v009{1});v009{3}=strrep(v009{2},'_',...
'');v009{4}=strrep(v009{3},'-','');v005=strrep(v005,' ','_');v001=false(size(v007));for v010=...
1:numel(v005),v011=v005{v010};[v012,v008{1}]=GetWritableFolder_f05(v008{1},v009{1},v011);if ...
numel(v012)~=1,v011=lower(v011);[v012,v008{2}]=GetWritableFolder_f05(v008{2},v009{2},v011);end,...
if numel(v012)~=1,v011=strrep(v011,'_','');[v012,v008{3}]=GetWritableFolder_f05(v008{3},v009{3},...
v011);end,if numel(v012)~=1,v011=strrep(v011,'-','');[v012,v008{4}]=...
GetWritableFolder_f05(v008{4},v009{4},v011);end,if numel(v012)~=1,...
error('parse_NameValue:NonUniqueMatch',v005{v010}),end,v000.(v007{v012})=v006{v010};v001(v012)=...
true;end,v001=v007(v001);end
function[v000,v001]=GetWritableFolder_f05(v001,v002,v003),v000=find(ismember(v002,v003));if ...
numel(v000)==1,return,end,if isempty(v001),v001=GetWritableFolder_f06(v002);end,v004=v001(:,...
1:min(end,numel(v003)));if size(v004,2)<numel(v003),v004=[v004 repmat(' ',size(v004,1),...
numel(v003)-size(v004,2))];end,v005=numel(v003)-sum(cumprod(double(v004==repmat(v003,size(v004,...
1),1)),2),2);v000=find(v005==0);end
function v000=GetWritableFolder_f06(v000),v001=cellfun('prodofsize',v000);v002=max(v001);for ...
v003=find(v001<v002).',v000{v003}((end+1):v002)=' ';end,v000=vertcat(v000{:});end
function[v000,v001]=GetWritableFolder_f07(v001),persistent v002,if isempty(v002),v002={true,...
false;1,0;'true','false';'1','0';'on','off';'enable','disable';'enabled','disabled'};end,if ...
isa(v001,'matlab.lang.OnOffSwitchState'),v000=true;v001=logical(v001);return,end,if isa(v001,...
'string'),if numel(v001)~=1,v000=false;return,else,v001=char(v001);end,end,if isa(v001,'char'),...
v001=lower(v001);end,for v003=1:size(v002,1),for v004=1:2,if isequal(v001,v002{v003,v004}),v000=...
true;v001=v002{1,v004};return,end,end,end,v000=false;end
function v000=GetWritableFolder_f08(v001),if~(isempty(v001)||exist(v001,'dir')),v000=false;
return,end,v002='';while isempty(v002)||exist(v002,'file'),[v003,v002]=...
fileparts(GetWritableFolder_f10('write_permission_test_','.txt'));v002=fullfile(v001,v002);end,...
try v004=fopen(v002,'w');fprintf(v004,'test');fclose(v004);delete(v002);v000=true;catch,if ...
exist(v002,'file'),try delete(v002);catch,end,end,v000=false;end,end
function v000=GetWritableFolder_f09(v001,v002,v003,v004,v005),if nargin<2||nargout>1,...
error('incorrect number of input/output arguments'),end,persistent v006 v007 v008,if ...
isempty(v006),v008=exist('OCTAVE_VERSION','builtin');v006=[100,1] * sscanf(version,'%d.%d',2);
v007={'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;
'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;'R2009a' 708;
'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;'R2015b' 806;'R2016a' 900;
'R2016b' 901;'R2017a' 902;'R2017b' 903;'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;
'R2020a' 908;'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;'R2023a' 914;
'R2023b' 2302};end,if v008,if nargin==2,warning('HJW:ifversion:NoOctaveTest',...
['No version test for Octave was provided.',char(10),...
'This function might return an unexpected outcome.']),if isnumeric(v002),v009=...
0.1*v002+0.9*GetWritableFolder_f03(v002);v009=round(100*v009);else,v010=ismember(v007(:,1),...
v002);if sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,elseif nargin==4,[v001,v009]=deal(v003,v004);v009=...
0.1*v009+0.9*GetWritableFolder_f03(v009);v009=round(100*v009);else,[v001,v009]=deal(v004,v005);
v009=0.1*v009+0.9*GetWritableFolder_f03(v009);v009=round(100*v009);end,else,if isnumeric(v002),...
v009=GetWritableFolder_f03(v002*100);if mod(v009,10)==0,v009=...
GetWritableFolder_f03(v002)*100+mod(v002,1)*10;end,else,v010=ismember(v007(:,1),v002);if ...
sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,end,switch v001,case'==',v000=v006==v009;case'<',v000=v006 < v009;case'<=',v000=v006 <=...
v009;case'>',v000=v006 > v009;case'>=',v000=v006 >=v009;end,end
function v000=GetWritableFolder_f10(v001,v002),if nargin<1,v001='';end,if~isempty(v001),v001=...
[v001 '_'];end,if nargin<2,v002='';else,if~strcmp(v002(1),'.'),v002=['.' v002];end,end,v000=...
tempname;[v003,v004]=fileparts(v000);v000=fullfile(v003,[v001 v004 v002]);end
function tf=hasFeature(feature)
% Provide a single point to encode whether specific features are available.
persistent FeatureList
if isempty(FeatureList)
    checkpoint('hasFeature','ifversion')
    FeatureList = struct(...
        'HG2'              ,ifversion('>=','R2014b','Octave','<' ,0),...
        'ImplicitExpansion',ifversion('>=','R2016b','Octave','>' ,0),...
        'bsxfun'           ,ifversion('>=','R2007a','Octave','>' ,0),...
        'IntegerArithmetic',ifversion('>=','R2010b','Octave','>' ,0),...
        'String'           ,ifversion('>=','R2016b','Octave','<' ,0),...
        'HTTPS_support'    ,ifversion('>' ,0       ,'Octave','<' ,0),...
        'json'             ,ifversion('>=','R2016b','Octave','>=',7),...
        'strtrim'          ,ifversion('>=',7       ,'Octave','>=',0),...
        'accumarray'       ,ifversion('>=',7       ,'Octave','>=',0));
    checkpoint('hasFeature','CharIsUTF8')
    FeatureList.CharIsUTF8 = CharIsUTF8;
end
tf = FeatureList.(feature);
end
function v000=ifversion(v001,v002,v003,v004,v005),if nargin<2||nargout>1,...
error('incorrect number of input/output arguments'),end,persistent v006 v007 v008,if ...
isempty(v006),v008=exist('OCTAVE_VERSION','builtin');v006=[100,1] * sscanf(version,'%d.%d',2);
v007={'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;
'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;'R2009a' 708;
'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;'R2015b' 806;'R2016a' 900;
'R2016b' 901;'R2017a' 902;'R2017b' 903;'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;
'R2020a' 908;'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;'R2023a' 914;
'R2023b' 2302};end,if v008,if nargin==2,warning('HJW:ifversion:NoOctaveTest',...
['No version test for Octave was provided.',char(10),...
'This function might return an unexpected outcome.']),if isnumeric(v002),v009=...
0.1*v002+0.9*ifversion_f00(v002);v009=round(100*v009);else,v010=ismember(v007(:,1),v002);if ...
sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,elseif nargin==4,[v001,v009]=deal(v003,v004);v009=0.1*v009+0.9*ifversion_f00(v009);v009=...
round(100*v009);else,[v001,v009]=deal(v004,v005);v009=0.1*v009+0.9*ifversion_f00(v009);v009=...
round(100*v009);end,else,if isnumeric(v002),v009=ifversion_f00(v002*100);if mod(v009,10)==0,...
v009=ifversion_f00(v002)*100+mod(v002,1)*10;end,else,v010=ismember(v007(:,1),v002);if ...
sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,end,switch v001,case'==',v000=v006==v009;case'<',v000=v006 < v009;case'<=',v000=v006 <=...
v009;case'>',v000=v006 > v009;case'>=',v000=v006 >=v009;end,end
function v000=ifversion_f00(v000),v000=fix(v000+eps*1e3);end

function[v000,v001]=isnetavl(v002),if nargin==0,v003=false;else,[v004,v005]=isnetavl_f03(v002);
v003=v004&&v005;end,if v003,[v000,v001]=isnetavl_f02;return,end,v003=isnetavl_f01;if ...
isempty(v003),v000=0;v001=0;else,if v003,[v000,v001]=isnetavl_f02;else,[v000,v001]=isnetavl_f00;
end,end,end
function[v003,v006]=isnetavl_f00,if ispc,try [v000,v001]=system('ping -n 1 8.8.8.8');v002=...
v001(strfind(v001,' = ')+3);v002=v002(1:3);if~strcmp(v002,'110'),error('trigger error'),else,...
v003=1;[v004,v005]=regexp(v001,' [0-9]+ms');v006=v001((v004(1)+1):(v005(1)-2));v006=...
str2double(v006);end,catch,v003=0;v006=0;end,elseif isunix,try [v000,v001]=...
system('ping -c 1 8.8.8.8');v007=regexp(v001,', [01] ');if v001(v007+2)~='1',...
error('trigger error'),else,v003=1;[v004,v005]=regexp(v001,'=[0-9.]+ ms');v006=...
v001((v004(1)+1):(v005(1)-2));v006=str2double(v006);end,catch,v003=0;v006=0;end,else,...
error('How did you even get Matlab to work?'),end,end
function[v001,v002,v003]=isnetavl_f01,persistent v000,if~isempty(v000),v001=v000;return,end,...
[v002,v003]=isnetavl_f00;if v002,v000=false;v001=false;return,end,[v002,v003]=isnetavl_f02;if ...
v002,v000=true;v001=true;return,end,v001=[];end
function[v004,v005]=isnetavl_f02,persistent v000,if isempty(v000),try v001=...
isempty(which(func2str(@webread)));catch,v001=true;end,v000=~v001;end,try v002=now;if v000,v003=...
webread('http://google.com');else,v003=urlread('http://google.com');end,v004=1;v005=...
(now-v002)*24*3600*1000;catch,v004=0;v005=0;end,end
function[v000,v001]=isnetavl_f03(v001),persistent v002,if isempty(v002),v002={true,false;1,0;
'true','false';'1','0';'on','off';'enable','disable';'enabled','disabled'};end,if isa(v001,...
'matlab.lang.OnOffSwitchState'),v000=true;v001=logical(v001);return,end,if isa(v001,'string'),...
if numel(v001)~=1,v000=false;return,else,v001=char(v001);end,end,if isa(v001,'char'),v001=...
lower(v001);end,for v003=1:size(v002,1),for v004=1:2,if isequal(v001,v002{v003,v004}),v000=true;
v001=v002{1,v004};return,end,end,end,v000=false;end
function varargout=makedir(d)
% Wrapper function to account for old Matlab releases, where mkdir fails if the parent folder does
% not exist. This function will use the legacy syntax for those releases.
if exist(d,'dir'),return,end % Take a shortcut if the folder already exists.
persistent IsLegacy
if isempty(IsLegacy)
    % The behavior changed after R14SP3 and before R2007b, but since the legacy syntax will still
    % work in later releases there isn't really a reason to pinpoint the exact release.
    checkpoint('makedir','ifversion')
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
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
% Any input is converted to a row vector.

in = reshape(in,1,[]);
out = in;
if numel(pattern)==0 || numel(pattern)>numel(in)
    % Return input unchanged (apart from the reshape), as strrep does as well.
    return
end

L = true(size(in));
L((end-numel(pattern)+2):end) = false; % Avoid partial matches
for n=1:numel(pattern)
    % For every element of the pattern, look for matches in the data. Keep track of all possible
    % locations of a match by shifting the logical vector.
    % The last n elements should be left unchanged, to avoid false positives with a wrap-around.
    L_n = in==pattern(n);
    L_n = circshift(L_n,[0 1-n]);
    L_n(1:(n-1)) = L(1:(n-1));
    L = L & L_n;
    
    % If there are no matches left (even if n<numel(pat)), the process can be aborted.
    if ~any(L),return,end
end

if numel(rep)==0
    out(L)=[];
    return
end

% For the replacement, we will create a shadow copy with a coded char array. Non-matching values
% will be coded with a space, the first character of a match will be encoded with an asterisk, and
% trailing characters will be encoded with an underscore.
% In the next step, regexprep will be used to perform the replacement, after which indexing can be
% used to compose the final array.
if numel(pattern)>1
    checkpoint('PatternReplace','bsxfun_plus')
    idx = bsxfun_plus(find(L),reshape(1:(numel(pattern)-1),[],1));
else
    idx = find(L);
end
idx = reshape(idx,1,[]);
str = repmat(' ',1,numel(in));
str(idx) = '_';
str( L ) = '*';
NonMatchL = str==' ';

% The regular expression will take care of the lazy pattern matching. This also shifts the number
% of underscores to the length of the replacement array.
str = regexprep(str,'\*_*',['*' repmat('_',1,numel(rep)-1)]);

% We can paste in the non-matching positions. Positions where the replacement should be inserted
% may or may not be correct.
out(str==' ') = in(NonMatchL);

% Now we can paste in all the replacements.
x = strfind(str,'*');
checkpoint('PatternReplace','bsxfun_plus')
idx = bsxfun_plus(x,reshape(0:(numel(rep)-1),[],1));
idx = reshape(idx,1,[]);
out(idx) = repmat(rep,1,numel(x));

% Remove the elements beyond the range of what the resultant array should be.
out((numel(str)+1):end) = [];
end
function varargout=regexp_outkeys(v000,v001,varargin),if nargin<2,...
error('HJW:regexp_outkeys:SyntaxError','No supported syntax used: at least 3 inputs expected.'),...
end,if~(ischar(v000)&&ischar(v001)),error('HJW:regexp_outkeys:InputError',...
'All inputs must be char vectors.'),end,if nargout>nargin,error('HJW:regexp_outkeys:ArgCount',...
'Incorrect number of output arguments. Check syntax.'),end,persistent v002 v003 v004,if ...
isempty(v002),v002.start=true;v002.end=true;v002.match=regexp_outkeys_f01('<','R14','Octave',...
'<',4);v002.tokens=v002.match;v002.split=regexp_outkeys_f01('<','R2007b','Octave','<',4);v005=...
fieldnames(v002);v003=['Extra regexp output type not implemented,',char(10),...
'only the following',' types are implemented:',char(10),sprintf('%s, ',v005{:})];
v003((end-1):end)='';v002.any=v002.match||v002.split||v002.tokens;v004=v002;for v006=...
fieldnames(v004).',v004.(v006{1})=false;end,end,if v002.any||nargin==...
2||any(ismember(lower(varargin),{'start','end'})),[v007,v008,v009]=regexp(v000,v001);end,if ...
nargin==2,varargout={v007,v008,v009};return,end,varargout=cell(size(varargin));v010=v004;v011=...
[];for v012=1:(nargin-2),if~ischar(varargin{v012}),error('HJW:regexp_outkeys:InputError',...
'All inputs must be char vectors.'),end,switch lower(varargin{v012}),case'match',if v010.match,...
varargout{v012}=v013;continue,end,if v002.match,v013=cell(1,numel(v007));for v014=1:numel(v007),...
v013{v014}=v000(v007(v014):v008(v014));end,else,[v013,v007,v008]=regexp(v000,v001,'match');end,...
varargout{v012}=v013;v010.match=true;case'split',if v010.split,varargout{v012}=v011;continue,...
end,if v002.split,v011=cell(1,numel(v007)+1);v015=[v007 numel(v000)+1];v016=[0 v008];for v014=...
1:numel(v015),v011{v014}=v000((v016(v014)+1):(v015(v014)-1));if numel(v011{v014})==0,v011{v014}=...
char(ones(0,0));end,end,else,[v011,v007,v008]=regexp(v000,v001,'split');end,varargout{v012}=...
v011;v010.split=true;case'tokens',if v010.tokens,varargout{v012}=v017;continue,end,if ...
v002.tokens,v017=cell(numel(v009),0);for v014=1:numel(v009),if size(v009{v014},2)~=2,v017{v014}=...
cell(1,0);else,for v018=1:size(v009{v014},1),v017{v014}{v018}=v000(v009{v014}(v018,...
1):v009{v014}(v018,2));end,end,end,else,[v017,v007,v008]=regexp(v000,v001,'tokens');end,...
varargout{v012}=v017;v010.tokens=true;case'start',varargout{v012}=v007;case'end',...
varargout{v012}=v008;otherwise,error('HJW:regexp_outkeys:NotImplemented',v003),end,end,if ...
nargout>v012,varargout(v012+[1 2])={v007,v008};end,end
function v000=regexp_outkeys_f00(v000),v000=fix(v000+eps*1e3);end
function v000=regexp_outkeys_f01(v001,v002,v003,v004,v005),if nargin<2||nargout>1,...
error('incorrect number of input/output arguments'),end,persistent v006 v007 v008,if ...
isempty(v006),v008=exist('OCTAVE_VERSION','builtin');v006=[100,1] * sscanf(version,'%d.%d',2);
v007={'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;'R14SP3' 701;
'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;'R2008a' 706;'R2008b' 707;'R2009a' 708;
'R2009b' 709;'R2010a' 710;'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;'R2015b' 806;'R2016a' 900;
'R2016b' 901;'R2017a' 902;'R2017b' 903;'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;
'R2020a' 908;'R2020b' 909;'R2021a' 910;'R2021b' 911;'R2022a' 912;'R2022b' 913;'R2023a' 914};end,...
if v008,if nargin==2,warning('HJW:ifversion:NoOctaveTest',...
['No version test for Octave was provided.',char(10),...
'This function might return an unexpected outcome.']),if isnumeric(v002),v009=...
0.1*v002+0.9*regexp_outkeys_f00(v002);v009=round(100*v009);else,v010=ismember(v007(:,1),v002);
if sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,elseif nargin==4,[v001,v009]=deal(v003,v004);v009=0.1*v009+0.9*regexp_outkeys_f00(v009);
v009=round(100*v009);else,[v001,v009]=deal(v004,v005);v009=...
0.1*v009+0.9*regexp_outkeys_f00(v009);v009=round(100*v009);end,else,if isnumeric(v002),v009=...
regexp_outkeys_f00(v002*100);if mod(v009,10)==0,v009=regexp_outkeys_f00(v002)*100+mod(v002,...
1)*10;end,else,v010=ismember(v007(:,1),v002);if sum(v010)~=1,warning('HJW:ifversion:NotInDict',...
'The requested version is not in the hard-coded list.'),v000=NaN;return,else,v009=v007{v010,2};
end,end,end,switch v001,case'==',v000=v006==v009;case'<',v000=v006 < v009;case'<=',v000=v006 <=...
v009;case'>',v000=v006 > v009;case'>=',v000=v006 >=v009;end,end

function SelfTestFailMessage=SelfTest__regexp_outkeys
% Run a self-test to ensure the function works as intended.
% This is intended to test internal function that do not have stand-alone testers, or are included
% in many different functions as subfunction, which would make bug regression a larger issue.

checkpoint('SelfTest__regexp_outkeys','regexp_outkeys')
ParentFunction = 'regexp_outkeys';
% This flag will be reset if an error occurs, but otherwise should ensure this test function
% immediately exits in order to minimize the impact on runtime.
if nargout==1,SelfTestFailMessage='';end
persistent SelfTestFlag,if ~isempty(SelfTestFlag),return,end
SelfTestFlag = true; % Prevent infinite recursion.

test_number = 0;ErrorFlag = false;
while true,test_number=test_number+1;
    switch test_number
        case 0 % (test template)
            try ME=[];
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 1
            % Test if all implemented output keys will return a value.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
                [val1,val2,val3] = regexp_outkeys(str,'( )','split','match','tokens');
                if isempty(val1) || isempty(val2) || isempty(val3)
                    error('one of the implemented outkeys is empty')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 2
            % Test if adding the start and end indices as outputs does not alter the others.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
                [a1,a2,a3,start_,end_] = regexp_outkeys(str,'( )','split','match','tokens');
                [b1,b2,b3] = regexp_outkeys(str,'( )','split','match','tokens');
                if isempty(start_) || isempty(end_) || ...
                        ~isequal(a1,b1) || ~isequal(a2,b2) || ~isequal(a3,b3)
                    error('one of the implemented outkeys is empty')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 3
            % Confirm a regex without tokens will have an empty tokens output.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
                NoTokenMatch = regexp_outkeys(str,' ','tokens');
                expected = repmat({cell(1,0)},1,6);
                if ~isequal(NoTokenMatch,expected)
                    error('no tokens in regex did not return empty result')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 4
            % Check the split option, including trailing empty.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
                SpaceDelimitedElements = regexp_outkeys(str,' ','split');
                expected = {'lorem1','ipsum1.2','dolor3','sit','amet','99',char(ones(0,0))};
                if ~isequal(SpaceDelimitedElements,expected)
                    error(['space delimited elements did not match expected result' char(10) ...
                        '(perhaps the trailing empty is 1x0 instead of 0x0)']) %#ok<CHARTEN>
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 5
            % Check the split option, including trailing empty.
            try ME=[];
                SpaceDelimitedElements = regexp_outkeys('',' ','split');
                expected = {char(ones(0,0))};
                if ~isequal(SpaceDelimitedElements,expected)
                    size(SpaceDelimitedElements{end}),size(expected{end}),keyboard
                    error('split on empty str did not return 0x0 empty')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 6
            % Check the extraction of a matched pattern.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 99 ';
                RawTokens = regexp_outkeys(str,'([a-z]+)[0-9]','tokens');
                words_with_number = horzcat(RawTokens{:});
                expected = {'lorem','ipsum','dolor'};
                if ~isequal(words_with_number,expected)
                    error('actual results did not match expected result')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 7
            % Check the extraction of a matched pattern.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 9x9 ';
                numbers = regexp_outkeys(str,'[0-9.]*','match');
                expected = {'1','1.2','3','9','9'};
                if ~isequal(numbers,expected)
                    error('actual results did not match expected result')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 8
            % Check the addition of start and end as tokens.
            try ME=[];
                str = 'lorem1 ipsum1.2 dolor3 sit amet 9x9 ';
                [ignore,end1,start1,end2] = regexp_outkeys(str,' .','match','end'); %#ok<ASGLU>
                [start2,end3] = regexp_outkeys(str,' .');
                expected = [7 16 23 27 32];
                if ~isequal(start1,expected) || ~isequal(start2,expected)
                    error('start indices did not match expectation')
                end
                expected = expected+1;
                if ~isequal(end1,expected) || ~isequal(end2,expected) || ~isequal(end3,expected)
                    error('end indices did not match expectation')
                end
            catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
                ErrorFlag = true;break
            end
        case 9
            % Check multi-element tokens.
            [t,s] = regexp_outkeys('2000/12/31','(\d\d\d\d)/(\d\d)/(\d\d)','tokens','split');
            expected1 = {{'2000','12','31'}};
            expected2 = {char(ones(0,0)),char(ones(0,0))};
            if ~isequal(t,expected1) || ~isequal(s,expected2)
                error('result did not match expectation for multiple tokens')
            end
        case 10
            % Check multi-element tokens.
            t = regexp_outkeys('12/34 56/78','(\d\d)/(\d\d)','tokens');
            expected = {{'12','34'},{'56','78'}};
            if ~isequal(t,expected)
                error('result did not match expectation for multiple tokens')
            end
        otherwise % No more tests.
            break
    end
end
if ErrorFlag
    SelfTestFlag = [];
    if isempty(ME)
        if nargout==1
            SelfTestFailMessage=sprintf('Self-validator %s failed on test %d.\n',...
                ParentFunction,test_number);
        else
            error('self-test %d failed',test_number)
        end
    else
        if nargout==1
            SelfTestFailMessage=sprintf(...
                'Self-validator %s failed on test %d.\n   ID: %s\n   msg: %s\n',...
                ParentFunction,test_number,ME.identifier,ME.message);
        else
            error('self-test %d failed\n   ID: %s\n   msg: %s',...
                test_number,ME.identifier,ME.message)
        end
    end
end
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
    checkpoint('TestFolderWritePermission','tmpname')
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
        checkpoint('WinVer','ifversion')
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

function out=checkpoint(caller,varargin)
% This function has limited functionality compared to the debugging version.
% (one of the differences being that this doesn't read/write to a file)
% Syntax:
%   checkpoint(caller,dependency)
%   checkpoint(caller,dependency_1,...,dependency_n)
%   checkpoint(caller,checkpoint_flag)
%   checkpoint('reset')
%   checkpoint('read')
%   checkpoint('write_only_to_file_on_read')
%   checkpoint('write_to_file_every_call')

persistent data
if isempty(data)||strcmp(caller,'reset')
    data = struct('total',0,'time',0,'callers',{{}});
end
if strcmp(caller,"read")
    out = data.time;return
end
if nargin==1,return,end
then = now;
for n=1:numel(varargin)
    data.total = data.total+1;
    data.callers = sort(unique([data.callers {caller}]));
    if ~isfield(data,varargin{n}),data.(varargin{n})=0;end
    data.(varargin{n}) = data.(varargin{n})+1;
end
data.time = data.time+ (now-then)*( 24*60*60*1e3 );
data.time = round(data.time);
end

