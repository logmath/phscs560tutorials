function x = binfileload( path, IDname,IDnum,CHnum,N,Nstart )
% x = binfileload( path, IDname,IDnum,CHnum,N,Nstart )
% This function loads single-precision, little-endian binary files without header information.
% The file name has the format: IDnameIDnum_CHnum.bin, where IDnum and CHnum have %03.0f format. 
% Example: % ID001_004.bin
% path: file path, e.g., 'C:\Data'
% IDname: Root test name, e.g., 'ID'
% IDnum: Test number, e.g., 4
% CHnum: Channel number, e.g., 12
% Nstart: number of samples to offset from beginning of file.  Default is beginning of file
% N: Number of samples to read.  Default is the entire file
%
% Note - the path and IDname variables must be characters (' '), not strings (" ")
%
% Author: Kent Gee, 11/14/13
% Update: Includes provisions for four-digit ID numbers native to AFR
% version 12.16.2.76 and newer. LM 7 JUL 2020

if nargin<6
    Nstart=0;
end

if nargin<5
    N=inf;
end


filename=[path,filesep,IDname,sprintf('%03.0f',IDnum),'_',sprintf('%03.0f',CHnum),'.bin'];

% Error Checking: See if file exists
if ~isfile(filename) % Check to see if the three-digid ID number does not work, then try the four digit ID number.
    
    fourDigitFilename=[path,filesep,IDname,sprintf('%04.0f',IDnum),'_',sprintf('%03.0f',CHnum),'.bin'];
    
    if ~isfile(fourDigitFilename)
        disp('binfileload: The following file name is not found (with 3 or 4-digit ID number):')
        disp(filename) % display the 3-digit filename (for clarity)
    else
        disp('binfileload: 4-digit ID detected in file, using 4-digit ID number')
        filename = fourDigitFilename;
    end
    
end


% Read in data

fid=fopen(filename,'r');
Nstart=Nstart*4;   % Convert from samples to bytes
fseek(fid,Nstart,'bof');
x=fread(fid,N,'single');
fclose(fid);


end

