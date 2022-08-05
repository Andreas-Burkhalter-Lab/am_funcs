function elements = isoDalton_modify_isotope_composition(elements)
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_modify_isotope_composition.m
% Description:  	The file modifies isotopic composition based on user input in this file
% Author:					Ross K. Snider
% Creation Date:	Wednesday  -  March 07, 2007  -  3:29:27 PM		
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  elements{} structure from isoDalton_get_isotope_info()
%         User entered modify{} cell array structure for each modified element (added below the **User Input Here** line in this file, line 72)
%         the modify{} cell array structure needs 3 fields: 
%             modify{k}.atomic_number            - Atomic Number of the kth element being modified
%             modify{k}.isotopic_mass_numbers    - Isotopic mass numbers whose composition is being modified
%             modify{k}.isotopic_composition     - The desired composition for each mass number, values should be normalized to 1.0
%
%             Example:
%             %--------------------------------------
%             % Modify Carbon to add C14 
%             %--------------------------------------
%             % Note: Carbon (atomic mass = 6) has 15 isotopes and all isotopic mass numbers are listed below
%             %       Isotopic composition given is arbitrary, the user needs to put in the appropriate values based on their own situation
%             modify{1}.atomic_number = 6;
%             modify{1}.isotopic_mass_numbers = [8  9  10  11  12       13       14       15   16   17  18  19  20  21  22];  % if mass number is not entered, the isotopic composition will be assumed to be zero
%             modify{1}.isotopic_composition  = [0  0  0   0   0.9844   0.0106   0.005    0    0    0   0   0   0   0   0 ];   % must be paired with mass numbers, i.e. same indexes
%
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	modified elements{} structure with isotopic composition changed according to the entered modify{} structure
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
modify = [];
%**************************************************************
%*----------     User Input Below this Line -------------------
%**************************************************************
% Create the appropriate number of modify{} cell arrays for each desired element as given by the examples below.
% The isotopic compositions should be normalized to 1.0   Normalization will be performed regardless of the given composition values.
% Make sure that the modify{} cell arrays are uncommented (remove a single '%' at the start of the line) and that they are indexed in integer order starting at 1
%********************************************************************


%% Example:
%%--------------------------------------
%% Modify Carbon to add C14 
%%--------------------------------------
%% Note: Carbon (atomic mass = 6) has 15 isotopes and all isotopic mass numbers are listed below
%%       Isotopic composition given is arbitrary, the user needs to put in the appropriate values based on their own situation
%modify{1}.atomic_number = 6;
%modify{1}.isotopic_mass_numbers = [8  9  10  11  12       13       14       15   16   17  18  19  20  21  22];  % if mass number is not entered, the isotopic composition will be assumed to be zero
%modify{1}.isotopic_composition  = [0  0  0   0   0.9844   0.0106   0.005    0    0    0   0   0   0   0   0 ];   % must be paired with mass numbers, i.e. same indexes
%%-------------------------------------------------------------------------
%% Also modify Silver (atomic mass = 47) to add the unstable isotope Ag108
%%-------------------------------------------------------------------------
%% Note: Silver has 34 isotopes and just the desired non-zero isotopic compositions are given below
%%       Isotopic composition given is arbitrary, the user needs to put in the appropriate values based on their own situation
%modify{2}.atomic_number = 47;
%modify{2}.isotopic_mass_numbers = [107     108     109   ]; % if mass number is not entered, the isotopic composition will be assumed to be zero
%modify{2}.isotopic_composition  = [0.5133  0.0099  0.4768]; % must be paired with mass numbers, i.e. same indexes
%disp(['Isotope composition has been modified.]);






%**************************************************************
%*---------- Modifications Implemented Here -------------------
%*----------    Do Not Change Code Below    -------------------
%**************************************************************
Nmods = length(modify);  % number of modified elements - this will be zero if the modify{} cell arrays are commented out or none exist
if Nmods > 0
    for k=1:Nmods
        atomic_number = modify{k}.atomic_number;
        isotope_count = elements{atomic_number}.isotope_count;
        %-------- Normalize given isotopic composition -----------
        isotopic_composition = modify{k}.isotopic_composition;
        modify{k}.isotopic_composition = isotopic_composition/sum(isotopic_composition);
        %---- zero isotopic composition ------
        for i=1:isotope_count
            elements{atomic_number}.isotopes{i}.isotopic_composition = 0;
        end
        %--- Put in desired isotopic composition values ------
        Nmass = length(modify{k}.isotopic_mass_numbers);
        for j=1:Nmass
            mass_number = modify{k}.isotopic_mass_numbers(j);
            composition = modify{k}.isotopic_composition(j);
            for i=1:isotope_count
                if mass_number == elements{atomic_number}.isotopes{i}.mass_number
                    elements{atomic_number}.isotopes{i}.isotopic_composition = composition;
                    break;
                end
            end
        end
        %--- Update elements{} structure fields  ------
        max_val = -realmax;
        max_num = 0;
        max_index = 0;     
        nzmass =[];   
        for i=1:isotope_count
            if max_val < elements{atomic_number}.isotopes{i}.isotopic_composition
                max_val = elements{atomic_number}.isotopes{i}.isotopic_composition;
                max_num = elements{atomic_number}.isotopes{i}.mass_number;  % find dominant isotope
                max_index = i;
            end
            if elements{atomic_number}.isotopes{i}.isotopic_composition > 0
                nzmass =[nzmass i];
            end
        end
        elements{atomic_number}.dominant_isotope_mass_number = max_num;
        elements{atomic_number}.dominant_isotope_composition = max_val;
        elements{atomic_number}.dominant_isotope_index = max_index;        
        elements{atomic_number}.non_zero_isotope_count = length(nzmass);
        elements{atomic_number}.non_zero_isotope_composition_index = nzmass;
        elements{atomic_number}.composition_sum = sum(modify{k}.isotopic_composition);
    end
end

