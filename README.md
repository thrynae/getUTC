# getUTC documentation
[![View ifversion on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/75688-getutc)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=thrynae/getUTC)

Table of contents

- Description section:
- - [Description](#description)
- Matlab/Octave section:
- - [Syntax](#syntax)
- - [Output arguments](#output-arguments)
- - [Compatibility, version info, and licence](#compatibility-version-info-and-licence)

## Description
For some applications it is necessary to know the current UTC time. This function has three methods to determine the UTC timestamp: an implementation written in C (which will be compiled), a call to a website, and a system call.

There are several methods implemented in this function:
- An implementation that requires a C mex function (which will be compiled).<br>This method requires write access to a folder and a working C compiler. The compilation result will be stored to a subdirectory of a folder similar to the AddOn path, or the tempdir, or the current folder. Write permission is tested in that order.
- An implementation using https://www.utctime.net/utc-timestamp.<br>The NIST has a server that returns the time, but it currently blocks API access, hence this replacement.<br>This method requires internet access.
- The local time and timezone offset can be determined with the `wmic` command or the `get-date` Powershell function (Windows) or the `date` command (Linux and Mac).<br>On Windows an NTP query is also implemented (internet connectivity is checked prior to this call).

To speed up the usage of this method, you can cache the difference with `now()` in a persistent variable, that way you avoid the need for a slow system call.

## Matlab/Octave

### Syntax

    atomTime = getUTC

### Output arguments

|Argument|Description|
|---|---|
|atomTime|The current UTC in the numeric format produced by `datenum`. Depending on which method is used, the accuracy may be reduced by longer execution times (up to a second).|

### Compatibility, version info, and licence
Compatibility considerations:
- Some older releases don't support the web implementation.
- The normal system call hangs on ML7.1 on XP. Since ML6.5 works fine on Windows 10, it seems reasonable to assume that the OS is the cause of the hang. For XP (and older) there is an alternative strategy in place, but this has a higher likelihood to fail.
- Similarly, as `wmic` has been deprecated, a Powershell alternative should be used on newer versions of Windows. The need for this is automatically detected.

|Test suite result|Windows|Linux|MacOS|
|---|---|---|---|
|Matlab R2023b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2023a|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it></it>|
|Matlab R2022b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2022a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2021b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2021a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2020b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2020a|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Matlab R2019b|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it>Monterey : Pass</it>|
|Matlab R2019a|<it>W11 : Partial</it>|<it></it>|<it></it>|
|Matlab R2018a|<it>W11 : Partial</it>|<it>ubuntu_22.04 : Partial</it>|<it></it>|
|Matlab R2017b|<it>W11 : Partial</it>|<it>ubuntu_22.04 : Partial</it>|<it>Monterey : Partial</it>|
|Matlab R2016b|<it>W11 : Partial</it>|<it>ubuntu_22.04 : Partial</it>|<it>Monterey : Partial</it>|
|Matlab R2015a|<it>W11 : Pass</it>|<it>ubuntu_22.04 : Pass</it>|<it></it>|
|Matlab R2013b|<it>W11 : Partial</it>|<it></it>|<it></it>|
|Matlab R2007b|<it>W11 : Partial</it>|<it></it>|<it></it>|
|Matlab 6.5 (R13)|<it>W11 : Partial</it>|<it></it>|<it></it>|
|Octave 8.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 7.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 6.2.0|<it>W11 : Pass</it>|<it>ubuntu_22.04 : </it>|<it>Catalina : </it>|
|Octave 5.2.0|<it>W11 : Pass</it>|<it></it>|<it></it>|
|Octave 4.4.1|<it>W11 : Pass</it>|<it></it>|<it>Catalina : </it>|

    Version: 2.2.0
    Date:    2023-11-21
    Author:  H.J. Wisselink
    Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
    Email = 'h_j_wisselink*alumnus_utwente_nl';
    Real_email = regexprep(Email,{'*','_'},{'@','.'})

### Test suite

The tester is included so you can test if your own modifications would introduce any bugs. These tests form the basis for the compatibility table above. Note that functions may be different between the tester version and the normal function. Make sure to apply any modifications to both.
