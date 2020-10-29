[![View getUTC on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/75688-getutc)

For some applications it is necessary to know the current UTC time. This function has an implementation written in C that it will automatically compile.  
If a subdirectory of the `tempdir` does not have writing permissions, the current folder will be used as the build target. If the building process fails, the function will attempt to use a website that provides the UTC timestamp. This used to be possible with a NIST server, but that no longer works. So to use this function you will either need write access to a folder and a working C compiler, or internet connectivity.

Licence: CC by-nc-sa 4.0