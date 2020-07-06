function aaa___test___getUTC(varargin)
%(error for ML6.5 and R2011a web expected)

% suppress v6.5 warning
v=version;v(strfind(v,'.'):end)='';v=str2double(v);
if v<=7
    warning off MATLAB:m_warning_end_without_block;clc
end

%if nargin==0,RunTestHeadless=false;else,RunTestHeadless=true;end
getUTC(1);%triggers compile if needed

t1=now;                      %tic
C=getUTC(1);
t1=(now-t1)*60*60*24;t2=now; %toc,tic
web=getUTC(2);
t2=(now-t2)*60*60*24;%       %toc

isOctave=exist('OCTAVE_VERSION', 'builtin') ~= 0;
v=version;ind=strfind(v,'.');v=str2double(v(1:(ind(1)-1)));%only the value before the first period
fail=false;
if isempty(C)
    fail=true;fprintf('C function failed\n')
end
if isempty(web)
    if ~isOctave && v<8
        %Known issue: Matlab 6.5 and R2011a fail here, let's assume it works again in R2012b
        fprintf('web function failed as expected\n')
    else
        fail=true;fprintf('web function failed\n')
    end
elseif ~isOctave && v<8
    fail=true;fprintf('web function did not fail as expected, update compatibility matrix\n')
end
if fail
    error('test failed')
else
    fprintf('C function took %.2f seconds\n',t1)
    fprintf('web function took %.2f seconds\n',t2)
    fprintf('test completed\n')
end
end
function atomTime=getUTC(debug_test)
%Returns the UTC time. The value is in Matlab datenum format.
%
%example syntax:
% disp(datestr(getUTC))
%
% There are two methods implemented in this function:
% - An implementation that requires a C mex function.
%   This method requires write access to the current folder and a working C compiler.
%   (you may want to compile the mex function to the same folder you store this m-file)
% - An implementation using https://www.utctime.net/utc-timestamp. The NIST has a server that
%   returns the time, but it currently blocks API access.
%   This method requires internet access.
% 
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020a     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  partial #1 |  partial #1      |  not tested          |
% | ML 6.5 (R13)  |  partial #1 |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% note #1: web method doesn't work
%
% Version: 1.0
% Date:    2020-05-20
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
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
%Use a C implementation, which requires write permission in the current directory.
%Should return an empty array instead of an error if it fails.

persistent utc_time_c
if isempty(utc_time_c)
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
    UTC_epoch_seconds=utc_time;
catch
    %build missing C file
    if exist(['utc_time.' mexext],'file')
        ME=lasterror; %#ok<LERR>
        rethrow(ME);
    end
    if ~exist('utc_time.c','file')
        fid=fopen('utc_time.c','w');
        for line=1:numel(utc_time_c)
            fprintf(fid,'%s\n',utc_time_c{line});
        end
        fclose(fid);
    end
    mex('utc_time.c');
    delete('utc_time.c')%cleanup
    if exist('utc_time.o','file'),delete('utc_time.o'),end%cleanup on Octave
    if exist(['utc_time.' mexext],'file')
        UTC_epoch_seconds=utc_time;
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
% | ML R2020a     |  works      |  not tested      |  not tested          |
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
% Licence: CC by-nc-sa 4.0 ( http://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email=  'h_j_wisselink*alumnus_utwente_nl';
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