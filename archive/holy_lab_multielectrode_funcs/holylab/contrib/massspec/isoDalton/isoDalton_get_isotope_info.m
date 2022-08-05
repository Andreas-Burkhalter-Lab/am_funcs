function elements = isoDalton_get_isotope_info();
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_get_isotope_info.m
% Description:  	The file returns isotope information for elements up to and including atomic number 112
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  3:31:32 PM			
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  None (reads files that should be present)
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	The cell array elements, indexed by atomic number, with the following fields:
%           (copper, atomic number 29, is taken for example):
%
%        elements{atomic_number}
%        elements{29} = 
%                                      name: 'Copper'
%                                    symbol: 'Cu'
%                    standard_atomic_weight: 63.5460
%                             isotope_count: 29
%                      isotope_mass_numbers: [52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80]
%                                  isotopes: {1x29 cell}
%              dominant_isotope_mass_number: 63
%              dominant_isotope_composition: 0.69170000000000
%                    dominant_isotope_index: 12
%                    non_zero_isotope_count: 2
%        non_zero_isotope_composition_index: [12 14]                        
%                           composition_sum: 1
%
%        elements{29}.isotopes{12} = 
%                                    symbol: ' Cu'
%                               mass_number: 63
%                      relative_atomic_mass: 62.92960110000000
%                      isotopic_composition: 0.69170000000000
%
%        elements{29}.isotopes{14} = 
%                                    symbol: ' Cu'
%                               mass_number: 65
%                      relative_atomic_mass: 64.92779370000000
%                      isotopic_composition: 0.30830000000000
%    
%
%---------------------------------------------------------------------------------------------------
%
% Modifications (give date, author, and description)
%
% 1.  June 26, 2007,  Ross K. Snider
%     Added ability to load in precomputed element file
%     NOTE: path to file is hard coded......
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

data_dir = 'D:\users\ross\Matlab\isoDalton\';
data_file = 'isoDalton_element_info.mat';

if exist([data_dir data_file],'file') == 2 % load precomputed mappings
    display('Loading NIST isotope information')
    display(['To add custom isotopic compositions, Modify: isoDalton_modify_isotope_composition.m, and Delete file: ' data_dir data_file])
    eval(['load ' data_dir data_file ';']);       % ---- load mappings  ------
elseif exist([data_dir data_file],'file') == 0
    display('Creating NIST isotope information file')
    tic

    element = isoDalton_NIST_isotopes_read();    % get the NIST isotope information
    names   = isoDalton_element_symbols_read();  % get element symbols

    %check = [];
    N_elements = 112;  % number of elements
    elements = cell(N_elements,1);
    for i=1:N_elements
        elements{i}.name   = names{i}.name;
        elements{i}.symbol = names{i}.symbol;
        elements{i}.standard_atomic_weight = element{i}.standard_atomic_weight;

        %------------------------------------------------
        % put isotope information in a more useable form
        %------------------------------------------------
        Ni = length(element{i}.isotope);
        icount = 0;
        mass = [];
        for j=1:Ni
            if length(element{i}.isotope{j}) > 0 && length(element{i}.isotope{j}.isotopic_composition) > 0
                mass = [mass j];
                icount = icount + 1;
            end
        end
        if icount == 0;
            icount = 1;
            mass = Ni;
        end
        elements{i}.isotope_count = icount;
        elements{i}.isotope_mass_numbers = mass;

        if elements{i}.isotope_count > 1
            Ni = length(element{i}.isotope);
            icount = 1;
            max_val = -realmax;
            max_num = 0;
            max_index = 0;
            comp_sum = 0;
            nzmass =[];
            min_dist = realmax;
            for j=1:Ni
                if length(element{i}.isotope{j}) > 0 && length(element{i}.isotope{j}.isotopic_composition) > 0
                    elements{i}.isotopes{icount}.symbol = element{i}.isotope{j}.symbol;
                    elements{i}.isotopes{icount}.mass_number = j;
                    elements{i}.isotopes{icount}.relative_atomic_mass = element{i}.isotope{j}.relative_atomic_mass;
                    elements{i}.isotopes{icount}.isotopic_composition = element{i}.isotope{j}.isotopic_composition;
                    comp_sum = comp_sum + elements{i}.isotopes{icount}.isotopic_composition;
                    if max_val < elements{i}.isotopes{icount}.isotopic_composition
                        max_val = elements{i}.isotopes{icount}.isotopic_composition;
                        max_num = elements{i}.isotopes{icount}.mass_number;  % find dominant isotope
                        max_index = icount;
                    end
                    if element{i}.isotope{j}.isotopic_composition > 0
                        nzmass =[nzmass icount];
                    end
                    if min_dist > abs(elements{i}.standard_atomic_weight - element{i}.isotope{j}.relative_atomic_mass)
                        min_dist = abs(elements{i}.standard_atomic_weight - element{i}.isotope{j}.relative_atomic_mass);
                        min_dist_index = icount;
                        min_dist_mass_num = j;
                    end
                    icount = icount + 1;
                end
            end
            if max_val > 0
                elements{i}.dominant_isotope_mass_number = max_num;
                elements{i}.dominant_isotope_composition = max_val;
                elements{i}.dominant_isotope_index = max_index;
                elements{i}.non_zero_isotope_count = length(nzmass);
                elements{i}.non_zero_isotope_composition_index = nzmass;
                elements{i}.composition_sum = comp_sum;
            else
                elements{i}.dominant_isotope_mass_number = min_dist_mass_num;   % unstable element case
                elements{i}.dominant_isotope_composition = 1.0;
                elements{i}.dominant_isotope_index = min_dist_index;
                elements{i}.non_zero_isotope_count = 1;
                elements{i}.non_zero_isotope_composition_index = min_dist_index;
                elements{i}.composition_sum = 1.0;
            end
        else
            icount = 1;
            comp_sum = 0;
            elements{i}.isotopes{icount}.symbol = element{i}.isotope{Ni}.symbol;
            elements{i}.isotopes{icount}.mass_number = Ni;
            elements{i}.isotopes{icount}.relative_atomic_mass = element{i}.isotope{Ni}.relative_atomic_mass;
            elements{i}.isotopes{icount}.isotopic_composition = 1.0;
            elements{i}.dominant_isotope_mass_number = Ni;
            elements{i}.dominant_isotope_composition = 1.0;
            elements{i}.dominant_isotope_index = 1;
            elements{i}.non_zero_isotope_count = 1;
            elements{i}.non_zero_isotope_composition_index = 1;
        end
    end
    %-----------------------------------------------
    % Implement user defined isotopic compositions
    %-----------------------------------------------
    elements = isoDalton_modify_isotope_composition(elements);

    eval(['save ' data_dir data_file ' elements']);       % ---- load mappings  ------
    t=toc
end

