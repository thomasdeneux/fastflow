function  loadsum2_global(filename, data_type, xsize, ysize, nframes, header_size)
%function  loadsum2_global(filename, data_type, xsize, ysize, nframes, header_size)

%LOADSUM Loads frames from sum-files into matlab.
%	Y = LOADSUM2_REF('Y',FILENAME, FRAMES, SUMFRAMES, DATA_TYPE, XSIZE, YSIZE, HEADER_SIZE, XINDICES, YINDICES)
%   'Y': name of Y in caller workspace
%	FILENAME : name of file (many formats supported)
%	FRAMES 	 : a vector containing locations of the desired frames in the
%		   file, starting at 0.
%	SUMFRAMES: an optional flag (default 1) which determines whether the
%		   routine will sum all the FRAMES into a single frame or 
%		   return an array with all the frames.
%	DATA_TYPE: 0 for ushort  (16bit) (CCD, raw VDAQ, ORA) -- default
%		   1 for ulong   (32bit) (sum16)
%		   2 for float   (32bit) (sumsum, convert)
%		   3 for double  (64bit) (matlab)
%		   4 for signed char     (Fuji)
%		   5 for unsigned char   (---)
%		   6 for signed short	 (DYEDAQ)
%	XSIZE & YSIZE: dimensions of the frames. Should be provided either as
%		       parameters or as global variables.
%	HEADER_SIZE:            Optional size (in bytes) of header (default 0).
%   XINDICES, YINDICES:     Optional subselection of voxels 
%    [=ajout T. Deneux]     global indices if XINDICES only (default: 1:xsize and 1:ysize) 
%
%	The frames are returned as XSIZE*YSIZE-long column vectors.
%
%	Doron. 19 Aug 1993.
%	Last revision: 10 Dec 1997.
%
%---
% ne doit pas recr�er de l'espace m�moire - trop chiant ce faux passage par
% r�f�rence - Thomas 7 Dec 2005

global Y

if isempty(Y) || any(size(Y)~=[xsize ysize nframes])
    error('global variable ''Y'' does not have the good dimensions')
end

if (data_type == 0);
	type = 'ushort';
	bytes_per_pixel = 2;
elseif (data_type == 1),
	type = 'ulong';
	bytes_per_pixel = 4;
elseif (data_type == 2),
	type = 'float';
	bytes_per_pixel = 4;
elseif (data_type == 3),
        type = 'double';
        bytes_per_pixel = 8;
elseif (data_type == 4),
	type = 'schar';
	bytes_per_pixel = 1;
elseif (data_type == 5),
	type = 'uchar';
	bytes_per_pixel = 1;
elseif (data_type == 6);
	type = 'short';
	bytes_per_pixel = 2;
else
	error('Wrong data_type!');
end;
frame_size = xsize * ysize;

Y = reshape(Y,frame_size, nframes);

[file, msg] = fopen(filename, 'r', 'l'); %assume files are little-endian (PC or VAX)
if (file == -1), error(msg); end;

frames = 0:nframes-1;
for i=1:nframes,
	offset = frames(i) * frame_size * bytes_per_pixel + header_size;
	status = fseek(file, offset, 'bof');
	if status, error(['Could not seek to frame ' num2str(frames(i))]); end;
	[frame, count] = fread(file, frame_size, type);
	if (count ~= frame_size), error(['Could not read frame ' num2str(frames(i))]); end;
    Y(:, i) = frame;
end

Y = reshape(Y,[xsize ysize nframes]);

fclose(file);
