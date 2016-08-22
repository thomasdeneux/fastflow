function b = lut_design(A,option)
% function b = lut_design(A,option)
%---
% mask for performing forks detection in binary image
% Input:
% - A is a 3x3 matrix
% - option: condition over neighbourghs number -> 'n', '+n' or '-n'
% Output:
% - b is true or false

% Input
if nargin<2, option='+3'; end
n = abs(str2num(option));
pflag = findstr(option,'+');
nflag = findstr(option,'-');

b = false;

% center must be on
if ~A(5), return, end

% counting sign changes around the border
count = 0;
for i=[2 3 6 9 8 7 4 1];
    if A(i)
        count = count+1;
    end
end

% Output
if pflag
    b = (count>=n);
elseif nflag
    b = (count<=n);
else
    b = (count==n);
end
