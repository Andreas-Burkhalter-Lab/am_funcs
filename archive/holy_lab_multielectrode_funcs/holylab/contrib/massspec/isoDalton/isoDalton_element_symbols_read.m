function names = isoDalton_element_symbols_read()
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_element_symbols_read.m
% Description:  	This file parses the atomic name and symbol file isoDalton_element_symbols.txt
%                 The data file was created from http://www.chem.qmul.ac.uk/iupac/AtWt/
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  1:31:42 PM			
%
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  None (reads a file that should be present)
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	The cell array element with the following fields:
%
%        element{atomic_number}.symbol    % atomic_symbol;
%        element{atomic_number}.name      % element name
%
%---------------------------------------------------------------------------------------------------
%
% Modifications (give date, author, and description)
%
% None
%
%
% Please send bug reports and enhancement requests to isoDalton@snidertech.com
%
%---------------------------------------------------------------------------------------------------
%            
%    Copyright (C) 2007  Ross K. Snider
%
%    This software is associated with the following paper:
%    Snider, R.K. Efficient Calculation of Exact Mass Isotopic Distributions
%    J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.
%    The digital object identifier (DOI) link to paper:  http://dx.doi.org/10.1016/j.jasms.2007.05.016
%
%    This library is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 2.1 of the License, or (at your option) any later version.
%
%    This library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with this library; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%    Ross Snider
%    Snider Technology, Inc.
%    40 Arrowhead Trail
%    Bozeman, MT  59718
%    ross@snidertech.com
%
%---------------------------------------------------------------------------------------------------

fid=fopen('isoDalton_element_symbols.txt');
if fid < 0
    error('Cannot open isoDalton_element_symbols.txt');
end

line = fgetl(fid);

skip_header = 1;    % skip header  (lines starting with %)
while skip_header == 1 & length(line) > 0
    if line(1) == '%';
        line = fgetl(fid);
    else
        skip_header = 0;
    end
end

while length(line) == 0   % skip blank lines
    line = fgetl(fid);
end


skip_flag = 1;
while 1
    if skip_flag == 0
        line = fgetl(fid);
        if ~ischar(line), break, end
    end

    [w, line] = strtok(line);
    [s, line] = strtok(line);
    [n, line] = strtok(line);
   
   i = str2num(w);
   names{i}.symbol = s;
   names{i}.name = n;

   skip_flag = 0;    
end

fclose(fid);



