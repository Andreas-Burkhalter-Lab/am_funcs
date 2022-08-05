function molecule = isoDalton_get_molecule_isotopes(molecule_string);
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_get_molecule_isotopes.m
% Description:  	The file returns isotope information for elements in the molecule string
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  7:18:43 PM			
%
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  String of elements <Ex> (any of the first 112 on the periodic chart of the elements)
%         paired with the number of atoms for each element <Cx>, i.e. ['E1 C1 E2 C2 E3 C3 etc.'].  
%         Spaces are not required, i.e. all the following cases are allowed, although the top case is 
%         the most readable.:
%
%         molecule_string = 'C254 H377 N65 O75 S6 Fe2';
%         molecule_string = 'C254H377N65O75S6Fe2';
%         molecule_string = 'C 254 H 377 N 65 O 75 S 6 Fe 2';
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	The cell array molecule, with the following fields:
%
%        molecule.element_count     % number of elements in molecule
%                                     where the element field is a cell array with fields shown below. 
%                                     An example is given for the first element
%
%        molecule.element{1} = 
%                                symbol: 'C'
%                                  name: 'Carbon'
%                         atomic_number: 6
%                       number_of_atoms: 254
%                         isotope_count: 15
%       isotopic_relative_atomic_masses: [1x15 double]
%                  isotopic_composition: [0 0 0 0 0.9893 0.0107 0 0 0 0 0 0 0 0 0]
%                 isotopic_mass_numbers: [8 9 10 11 12 13 14 15 16 17 18 19 20 21 22]
%          dominant_isotope_mass_number: 12
%          dominant_isotope_composition: 0.9893
%                dominant_isotope_index: 5
%                non_zero_isotope_count: 2
%    non_zero_isotope_composition_index: [5 6]
%    
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

elements = isoDalton_get_isotope_info();  % get the isotopes
names = isoDalton_element_symbols_read();


Nms = length(molecule_string);
s_letter = isletter(molecule_string);
s_space   = isspace(molecule_string);

symbol_start = 0;
symbol_end   = 0;
symbol_count = 1;
for i=1:Nms
    if i == 1 && s_letter(1) == 1
        symbol_start = 1;
    end
    if i == Nms && s_letter(i) == 1 
        symbol_end = i;
    end
    if i == Nms && s_letter(i-1) == 1 && s_letter(i) == 0
        symbol_end = i-1;
    end
    if i > 1 && i < Nms
        if s_letter(i) == 1 &&  s_letter(i-1) == 0
            symbol_start = i;
        end
        if s_letter(i) == 0 &&  s_letter(i-1) == 1
            symbol_end = i-1;
        end
    end
    if symbol_start > 0 && symbol_end > 0
        symbols{symbol_count}.symbol_string = molecule_string(symbol_start:symbol_end);
        symbols{symbol_count}.symbol_start  = symbol_start;
        symbols{symbol_count}.symbol_end    = symbol_end;
        symbol_count = symbol_count + 1;
        symbol_start = 0;
        symbol_end   = 0;
    end    
end


Nsym = length(symbols);
for i=1:Nsym-1
    number_string = molecule_string(symbols{i}.symbol_end+1 : symbols{i+1}.symbol_start-1);
    number_space = isspace(number_string);
    N2 = length(number_space);
    k1 = 1;
    while number_space(k1) == 1
        k1 = k1 + 1;
    end
    k2 = N2;
    while number_space(k2) == 1
        k2 = k2 - 1;
    end
    element_count{i} = str2num(number_string(k1:k2));
end
number_string = molecule_string(symbols{Nsym}.symbol_end+1 : end);
number_space = isspace(number_string);
N2 = length(number_space);
k1 = 1;
while number_space(k1) == 1
    k1 = k1 + 1;
end
k2 = N2;
while number_space(k2) == 1
    k2 = k2 - 1;
end
element_count{Nsym} = str2num(number_string(k1:k2));


atomic_numbers = zeros(1,Nsym);
for i=1:Nsym
    atomic_numbers(i) = isoDalton_element_sym2num(symbols{i}.symbol_string,names);
end
[y,an_index]=sort(atomic_numbers);  % get index so the elements will be listed in atomic order weight

molecule.element_count = Nsym;
for i=1:Nsym
    an = atomic_numbers(i);
    Nic = elements{an}.isotope_count;
    isotope_weights = zeros(1,Nic);
    isotope_comp    = zeros(1,Nic);
    mass_numbers    = zeros(1,Nic);
    for j = 1:Nic
        isotope_weights(j) = elements{an}.isotopes{j}.relative_atomic_mass;
        isotope_comp(j)    = elements{an}.isotopes{j}.isotopic_composition;
        mass_numbers(j)    = elements{an}.isotopes{j}.mass_number;
    end
    molecule.element{an_index(i)}.symbol                              = symbols{i}.symbol_string;
    molecule.element{an_index(i)}.name                                = elements{an}.name;
    molecule.element{an_index(i)}.atomic_number                       = atomic_numbers(i);
    molecule.element{an_index(i)}.number_of_atoms                     = element_count{i};
    molecule.element{an_index(i)}.isotope_count                       = Nic;
    molecule.element{an_index(i)}.isotopic_relative_atomic_masses     = isotope_weights;
    molecule.element{an_index(i)}.isotopic_composition                = isotope_comp;
    molecule.element{an_index(i)}.isotopic_mass_numbers               = mass_numbers;
    molecule.element{an_index(i)}.dominant_isotope_mass_number        =  elements{an}.dominant_isotope_mass_number;
    molecule.element{an_index(i)}.dominant_isotope_composition        =  elements{an}.dominant_isotope_composition;
    molecule.element{an_index(i)}.dominant_isotope_index              =  elements{an}.dominant_isotope_index;        
    molecule.element{an_index(i)}.non_zero_isotope_count              =  elements{an}.non_zero_isotope_count;
    molecule.element{an_index(i)}.non_zero_isotope_composition_index  =  elements{an}.non_zero_isotope_composition_index;
end








