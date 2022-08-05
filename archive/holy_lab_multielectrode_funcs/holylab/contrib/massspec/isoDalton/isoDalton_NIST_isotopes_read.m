function elements = isoDalton_NIST_isotopes_read()
%--------------------------------------------------------------------------------------------------
%
% Filename:     	NIST_isotope_read.m
% Description:  	This file parses the NIST isotope data file isoDalton_NIST_isotopes.txt
%                 The isoDalton_NIST_isotopes.txt data file was created from http://physics.nist.gov/PhysRefData/Compositions/index.html
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  1:24:52 PM			
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
% Output:  	The cell array with the following fields:
%
%        element{atomic_number}.standard_atomic_weight                       % standard_atomic_weight;
%        element{atomic_number}.isotope{mass_number}.symbol                  % atomic_symbol;
%        element{atomic_number}.isotope{mass_number}.relative_atomic_mass    % relative_atomic_mass;
%        element{atomic_number}.isotope{mass_number}.isotopic_composition    % isotopic_composition in fractional format (not percent)
%        element{atomic_number}.isotope{mass_number}.notes                   %  notes;
%
%        (and yes, this is inefficient storage...)
%
%---------------------------------------------------------------------------------------------------
%
% Modifications (give date, author, and description)
%
% None
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


element = cell(116,1);


fid=fopen('isoDalton_NIST_all_isotopes.txt');
if fid < 0
    error('Cannot open isoDalton_NIST_all_isotopes.txt');
end

line = fgetl(fid);
line_count = 1;

skip_header = 1;    % skip header  (lines starting with %)
while skip_header == 1 && length(line) > 0
    if line(1) == '%';
        line = fgetl(fid); line_count = line_count + 1;
    else
        skip_header = 0;
    end
end

while length(line) == 0   % skip blank lines
    line = fgetl(fid); line_count = line_count + 1;
end

skip_flag = 1;
while 1
    if skip_flag == 0
        line = fgetl(fid); line_count = line_count + 1;
        if ~ischar(line), break, end
    end
    
    stringv = 'Atomic Number =';
    k1 = findstr(line,stringv);      
    if length(k1) > 0
        atomic_number = str2double(line(k1+length(stringv):end));
    end
    
    line = fgetl(fid); line_count = line_count + 1;
    stringv = 'Atomic Symbol =';
    k1 = findstr(line,stringv);      
    if length(k1) > 0
        atomic_symbol = line(k1+length(stringv):end);
%        line
%        line_count
    end
    
    line = fgetl(fid); line_count = line_count + 1;
    stringv = 'Mass Number =';
    k1 = findstr(line,stringv);      
    if length(k1) > 0
        mass_number = str2double(line(k1+length(stringv):end));
    else
%        line
%        line_count
    end
    
    line = fgetl(fid); line_count = line_count + 1;
    stringv1 = 'Relative Atomic Mass =';
    k1 = findstr(line,stringv1);
    stringv2 = '(';
    k2 = findstr(line,stringv2);      
    if length(k2) > 0
        relative_atomic_mass = str2double(line(k1+length(stringv1):k2-1));
    else
        stringv2 = line(k1+length(stringv1)+1:end);
        if sum(isletter(stringv2)) > 0
            relative_atomic_mass = str2double(stringv2);      
        else
            relative_atomic_mass = -realmax;  % not given, so give an unrealistic value that should be apperent in any result that is using this isotope
        end
    end
    
    line = fgetl(fid); line_count = line_count + 1;
    stringv1 = 'Isotopic Composition =';
    k1 = findstr(line,stringv1);
    stringv2 = '(';
    k2 = findstr(line,stringv2);      
    if length(k2) > 0
        isotopic_composition = str2double(line(k1+length(stringv1):k2-1));
    else
        stringv2 = line(k1+length(stringv1)+1:end);
        if length(stringv2) > 0
            isotopic_composition = str2double(stringv2);      
        else
            isotopic_composition = 0;
        end
    end

    line = fgetl(fid); line_count = line_count + 1;
    stringv1 = 'Standard Atomic Weight =';
    k1 = findstr(line,stringv1);
    stringv2 = '(';
    k2 = findstr(line,stringv2);      
    if length(k2) > 0
        standard_atomic_weight = str2double(line(k1+length(stringv1):k2-1));
    else
        stringv2 = '[';
        k2 = findstr(line,stringv2);      
        stringv3 = ']';
        k3 = findstr(line,stringv3);     
        if length(k2) > 0 && length(k3) > 0
             standard_atomic_weight = str2double(line(k2+1:k3-1));
        else
             standard_atomic_weight = [];
        end
    end
    
    line = fgetl(fid); line_count = line_count + 1;
    stringv = 'Notes =';
    k1 = findstr(line,stringv);      
    notes = line(k1+length(stringv):end);
    
    line = fgetl(fid); line_count = line_count + 1;  % skip blank line
   
    if length(atomic_number) > 0 && length(mass_number) > 0
        element{atomic_number}.standard_atomic_weight = standard_atomic_weight;
        element{atomic_number}.isotope{mass_number}.symbol = atomic_symbol;
        element{atomic_number}.isotope{mass_number}.relative_atomic_mass =  relative_atomic_mass;
        element{atomic_number}.isotope{mass_number}.isotopic_composition =  isotopic_composition/100;   % put in fractional format
        element{atomic_number}.isotope{mass_number}.notes =  notes;
    else
        %line
        %line_count
    end

    skip_flag = 0;
end
fclose(fid);
elements = element;


