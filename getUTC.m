function atomTime=getUTC(debug_test)
%Returns the UTC time. The value is in Matlab datenum format.
%
%example syntax:
% disp(datestr(getUTC))
%
% There are two methods implemented in this function:
% - An implementation that requires a C mex function.
%   This method requires write access to a subdirectory of the tempdir (or the current folder) and
%   a working C compiler.
%   (you may want to compile the mex function to the same folder you store this m-file)
% - An implementation using https://www.utctime.net/utc-timestamp. The NIST has a server that
%   returns the time, but it currently blocks API access.
%   This method requires internet access.
% 
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020b     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  partial #1 |  partial #1      |  not tested          |
% | ML 6.5 (R13)  |  partial #1 |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% note #1: web method doesn't work
%
% Version: 1.1.0
% Date:    2020-10-29
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin==0
    %normal flow: first try the C method, then the web
    UTC_epoch_seconds=getUTC_c;
    if isempty(UTC_epoch_seconds)
        UTC_epoch_seconds=getUTC_web;
        if isempty(UTC_epoch_seconds)
            error('HJW:getUTC:TimeReadFailed',...
                ['Both methods of retrieving the UTC timestamp failed.\nEnsure you ',...
                'have write access to the current folder and check your internet connection.'])
        end
    end
else
    %debug/test, will not throw an error on fail
    if debug_test==1
        UTC_epoch_seconds=getUTC_c;
    else
        UTC_epoch_seconds=getUTC_web;
    end
end
UTC_offset=UTC_epoch_seconds/(24*60*60);
atomTime=UTC_offset+datenum(1970,1,1);
end
function UTC_epoch_seconds=getUTC_c
%Use a C implementation, which requires write permission in the tempdir or the current directory.
%Should return an empty array instead of an error if it fails.

persistent utc_time_c tempdir_f funname utc_time_fun_handle Compile_attempts_remaining
if isempty(utc_time_c)
    funname=['utc_time_' mymexext];
    funname=strrep(funname,['.' mexext],'');
    try utc_time_fun_handle=str2func(funname);catch,end %The try-catch required in Octave.
    
    %Try creating a folder in the tempdir and adding it to the path (if it is not already in
    %there). If the folder is not writeable, the current folder will be used.
    %In some release-runtime combinations addpath has a permanent effect, in others it doesn't. By
    %putting this code in this block, we are trying to keep these queries to a minimum.
    tempdir_f=fullfile(tempdir,'MATLAB','FileExchange','getUTC');
    try
        if isempty(strfind([path ';'],[tempdir_f ';'])) %#ok<STREMP> f is not on the path
            if ~exist(tempdir_f,'dir'),mkdir(tempdir_f),end
            addpath(tempdir_f,'-end')
        end
    catch
    end
    
    %Only allow a few compilation attempts, so this function doesn't cause a lot of disk I/O if
    %there is no working compiler.
    Compile_attempts_remaining=5;
    
    %prepare to write this to a file and compile
    utc_time_c={'#include "mex.h"'
        '#include "time.h"'
        ''
        '/* Abraham Cohn,  3/17/2005 */'
        '/* Philips Medical Systems */'
        ''
        'void mexFunction(int nlhs, mxArray *plhs[], int nrhs,'
        '                 const mxArray *prhs[])'
        '{'
        '  time_t utc;'
        '  '
        '  if (nlhs > 1) {'
        '    mexErrMsgTxt("Too many output arguments");'
        '  }'
        '  '
        '  /* Here is a nice ref: www.cplusplus.com/ref/ctime/time.html */'
        '  time(&utc);'
        '  /* mexPrintf("UTC time in local zone: %s",ctime(&utc)); */'
        '  /* mexPrintf("UTC time in GMT: %s",asctime(gmtime(&utc))); */'
        '  '
        '  /* Create matrix for the return argument. */'
        '  plhs[0] = mxCreateDoubleScalar((double)utc);'
        '   '
        '}'};
    %the original had mxCreateScalarDouble
end

try
    UTC_epoch_seconds=feval(utc_time_fun_handle);
catch
    if exist(['utc_time.' mymexext],'file')
        ME=lasterror; %#ok<LERR>
        rethrow(ME);
    end
    
    %Try building missing C file
    Compile_attempts_remaining=Compile_attempts_remaining-1;
    if Compile_attempts_remaining<0 % Don't endlessly try to compile.
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
        %cleanup
        for ext={'c','o'}
            file=fullfile(f,[funname '.' ext{1}]);if exist(file,'file'),delete(file),end
        end
    catch
    end
    cd(current_folder);
    if exist([funname '.' mexext],'file')
        utc_time_fun_handle=str2func(funname);%refresh handle
        UTC_epoch_seconds=feval(utc_time_fun_handle);
    else
        %the compiling of the mex function failed
        UTC_epoch_seconds=[];
    end
end
end
function UTC_epoch_seconds=getUTC_web
%read the timestamp from a web server
%this fails for ML6.5 for some reason

%skip this function if there is no internet connection (3 timeouts will take a lot of time)
if ~isnetavl,UTC_epoch_seconds=[];return,end
for tries=1:3
    try
        if exist('webread','file')
            data=webread('http://www.utctime.net/utc-timestamp');
        else
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
function [connected,timing]=isnetavl
% Ping to one of Google's DNSes.
% Optional second output is the ping time (0 if not connected).
%
% Includes a fallback to HTML if usage of ping is not allowed. This increases the measured ping.
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020b     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.2.2
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

tf=isnetavl__ICMP_is_blocked;
if isempty(tf)
    %Unable to determine if ping is allowed, the connection must be down
    connected=0;
    timing=0;
else
    if tf
        %ping not allowed
        %(the timing is not reliable)
        [connected,timing]=isnetavl___ping_via_html;
    else
        %ping allowed
        [connected,timing]=isnetavl___ping_via_system;
    end
end
end
function [connected,timing]=isnetavl___ping_via_html
%Ping is blocked by some organizations. As an alternative, the google.com page can be loaded as a
%normal HTML, which should work as well, although it is slower. This also means the ping timing is
%no longer reliable.
try
    then=now;
    if exist('webread','file')
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
            %This branch will error for 'destination host unreachable'
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
            %This branch includes 'destination host unreachable' errors
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
%Check if ICMP 0/8/11 is blocked
%
%tf is empty if both methods fail

persistent output
if ~isempty(output)
    tf=output;return
end

%First check if ping works
[connected,timing]=isnetavl___ping_via_system;
if connected
    %Ping worked and there is an internet connection
    output=false;
    tf=false;
    return
end

%There are two options: no internet connection, or ping is blocked
[connected,timing]=isnetavl___ping_via_html;
if connected
    %There is an internet connection, therefore ping must be blocked
    output=true;
    tf=true;
    return
end

%Both methods failed, internet is down. Leave the value of tf (and the persistent variable) set to
%empty so it is tried next time.
tf=[];
end
function ext=mymexext
v=version;ind=strfind(v,'.');v(ind(2):end)='';v=['v' strrep(v,'.','_')];
if ~exist('OCTAVE_VERSION', 'builtin')
    type=computer;
else
    arch=computer;arch=arch(1:(strfind(arch,'-')-1));
    if ispc
        if strcmp(arch,'x86_64')  ,type= 'win_64';
        elseif strcmp(arch,'i686'),type= 'win_i686';
        elseif strcmp(arch,'x86') ,type= 'win_x86';
        else                      ,type=['win_' arch];
        end
    elseif isunix && ~ismac %essentially this is islinux
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
ext=[v '_' type '.' mexext];
end
function tf=TestFolderWritePermission(f)
%returns true if the folder exists and allows Matlab to write files
%an empty input will generally test the pwd
%
%examples:
%  fn='foo.txt';if ~TestFolderWritePermission(fileparts(fn)),error('can''t write!'),end

%test existance
if ~( isempty(f) || exist(f,'dir') )
    tf=false;return
end

%test write permission
fn='';
while isempty(fn) || exist(fn,'file')
    %generate a random file name, making sure not to overwrite an exisiting file
    [ignore,fn]=fileparts(tmpname('write_permission_test_','.txt')); %#ok<ASGLU>
    fn=fullfile(f,fn);
end
try
    fid=fopen(fn,'w');fprintf(fid,'test');fclose(fid);
    delete(fn);
    tf=true;
catch
    try delete(fn);catch,end
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