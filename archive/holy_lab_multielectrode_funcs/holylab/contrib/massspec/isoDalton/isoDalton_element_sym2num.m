function atomic_numbers = isoDalton_element_sym2num(element_symbols,names)
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_element_sym2num.m
% Description:  	The file returns atomic numbers of the element symbols
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  4:58:15 PM		
%
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  String of elements, example: element_symbols = 'H He Li Be B C' 
%         and array names = isoDalton_element_symbols_read();
%
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	returns atomic numbers of the symbols
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


symbol_index = 1;
[s, element_symbols] = strtok(element_symbols);
symbols{symbol_index} = s;
symbol_index = symbol_index + 1;
while length(symbols) > 0
    [s, element_symbols] = strtok(element_symbols);
    if length(s) > 0
        symbols{symbol_index} = s;
        symbol_index = symbol_index + 1;
    else
        break;
    end
end
symbol_count = symbol_index - 1;
name_count = 112;

atomic_numbers = zeros(1,symbol_count);
for i=1:symbol_count
    found_flag = 0;
    for j=1:name_count
        if strcmp(lower(symbols{i}),lower(names{j}.symbol)) == 1
            atomic_numbers(i) = j;
            found_flag = 1;
            break;
        end
    end
    if found_flag == 0
        error([symbols{i} ' is not an element'])
    end
end    
   
