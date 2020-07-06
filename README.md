For some applications it is necessary to know the current UTC time. This function has an implementation written in C that it will automatically compile.  
If the current folder does not have writing permissions or the building process fails, it will attempt to use a website that provides the UTC timestamp. This used to be possible with a NIST server, but that no longer works. So to use this function you will either need write access to a folder (you could `cd` to the `tempdir` and back), or internet connectivity.

Licence: CC by-nc-sa 4.0