function states = isoDalton_exact_mass(molecule_string,maxstates)
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_exact_mass.m
% Description:  	The file returns the exact mass isotopic distributions of a molecule 
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  May 5, 2006  -  1:55:47 PM		
%
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Input:  molecule_string :
%         String of elements <Ex> (any of the first 112 on the periodic chart of the elements)
%         paired with the number of atoms for each element <Cx>, i.e. ['E1 C1 E2 C2 E3 C3 etc.'].  
%         Spaces are not required, i.e. all the following cases are allowed, although the top case is 
%         the most readable.:
%         molecule_string = 'C254 H377 N65 O75 S6 Fe2';
%         molecule_string = 'C254H377N65O75S6Fe2';
%         molecule_string = 'C 254 H 377 N 65 O 75 S 6 Fe 2';
%
%         maxstates : maximum number of mass states
%                     set maxstates = realmax; if one wishes to keep *ALL* the states.  Note that this
%                     is practical only for small molecules, since for large molecules this will slow
%                     to a catatonic crawl and run out of memory.
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	returns a two column matrix where the first column contains the exact masses
%           (most probable ones) and the second column are the probabilities (pruned).
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


%---------------------------
% Get isotopes of molecule
%---------------------------
molecule = isoDalton_get_molecule_isotopes(molecule_string);


Nel = molecule.element_count;  % number of different elements in molecule
isotope_count = zeros(Nel,1);
for i=1:Nel
    isotope_count(i) = molecule.element{i}.non_zero_isotope_count;   % number of isotopes for each element
end
[v,element_index] = sort(isotope_count);    % element_index is used to order the elements in increasing number of isotopes
                                            % we want to start with the elements with fewest isotopes to minimize growth of the trellis

%-----------------------------------------
% Get mass and probability spanning info
%-----------------------------------------
term_lightest       = 0;
term_heaviest       = 0;
term_most_probable  = 0;
term_least_probable = 0;
for i=1:Nel  % elements in molecule
    Noa = molecule.element{i}.number_of_atoms;
    Nic = molecule.element{i}.non_zero_isotope_count;   % isotope count
    mass_min = realmax;
    mass_max = -realmax;
    prob_min = realmax;
    prob_max = -realmax;
    for j=1:Nic
        isotope_index        = molecule.element{i}.non_zero_isotope_composition_index(j);
        isotope_mass         = molecule.element{i}.isotopic_relative_atomic_masses(isotope_index);
        isotope_probability  = molecule.element{i}.isotopic_composition(isotope_index);
        if mass_min > isotope_mass
            mass_min = isotope_mass;
        end
        if mass_max < isotope_mass
            mass_max = isotope_mass;
        end
        if prob_min > isotope_probability
            prob_min = isotope_probability;
        end
        if prob_max < isotope_probability
            prob_max = isotope_probability;
        end
    end
%    mass_min
%    mass_max
%    prob_min
%    prob_max
    term_lightest       = term_lightest + Noa*mass_min;
    term_heaviest       = term_heaviest + Noa*mass_max;
    term_most_probable  = term_most_probable  + Noa*log10(prob_max);
    term_least_probable = term_least_probable + Noa*log10(prob_min);
end
disp(['*************************************************************'])
disp(['Information regarding molecule ' molecule_string])
disp(['The lightest mass term = ' num2str(term_lightest,'%10.2f') ' daltons'])
disp(['The heaviest mass term = ' num2str(term_heaviest,'%10.2f') ' daltons'])
distribution_span = term_heaviest - term_lightest;
disp(['The isotopic distribution spans = ' num2str(distribution_span,'%10.2f') ' daltons'])
term_most_probable = 10^term_most_probable;
disp(['The most probable term = ' num2str(term_most_probable,'%10.5f')])
disp(['The least probable term (log10) = ' num2str(term_least_probable,'%10.2f')])
disp(['*************************************************************'])

                                           
tic
%---------------------------------
% Setup initial states
%---------------------------------
Ni = molecule.element{element_index(1)}.non_zero_isotope_count;
states1=zeros(Ni,2);
for i=1:Ni
    isotope_index = molecule.element{element_index(1)}.non_zero_isotope_composition_index(i);
    states1(i,1) = molecule.element{element_index(1)}.isotopic_relative_atomic_masses(isotope_index);
    states1(i,2) = molecule.element{element_index(1)}.isotopic_composition(isotope_index);
end

%---------------------------------
% Trellis
%---------------------------------
display(['Element ' molecule.element{element_index(1)}.symbol ' (' num2str(1) ' of ' num2str(Nel) ')  Count ' num2str(1) ' of ' num2str(molecule.element{element_index(1)}.number_of_atoms) '  Current states: ' num2str(Ni) ])
for iel=1:Nel  % elements in molecule
    Noa = molecule.element{element_index(iel)}.number_of_atoms;  % number of atoms of element
    for ioa=1:Noa % atoms of each element
        if ~(iel==1 && ioa==1)   % we already have these initial states
            Nst1 = length(states1(:,1));  % number of states in previous step
            Nic = molecule.element{element_index(iel)}.non_zero_isotope_count;   % isotope count
            states2=zeros(Nst1*Nic,2);  % number of new states is the product of the old states and number of isotopes of atom being added
            states2_index = 1;
            for iic=1:Nic
                for ist=1:Nst1
                    isotope_index = molecule.element{element_index(iel)}.non_zero_isotope_composition_index(iic);
                    states2(states2_index,1) = states1(ist,1) + molecule.element{element_index(iel)}.isotopic_relative_atomic_masses(isotope_index); % mass of new states
                    states2(states2_index,2) = states1(ist,2) * molecule.element{element_index(iel)}.isotopic_composition(isotope_index);  % probability of new states
                    states2_index = states2_index + 1;
                end
            end
            Nst2 = states2_index-1;
            states2 = sortrows(states2,1);  % sort new states by mass

            %-------------------------------------
            % combine states with identical masses
            %-------------------------------------
            [umass, index1, index2] = unique(states2(:,1));  % unique also sorts
            Num = length(umass);
            states3 = zeros(Num,2);
            states3(:,1) = umass;  % put unique masses in first column
            for km = 1:Num
                if km == 1
                    index = 1:index1(km);
                else
                    index = index1(km-1)+1:index1(km);
                end
                states3(km,2) = sum(states2(index,2));
            end
            %------------------------------------------------------------------------------------------
            % combine states with mass differences equal to or less than eps() (i.e. Matlab precision)
            %------------------------------------------------------------------------------------------
            eps_threshold = eps(max(states3(:,1)));
            state_diff = diff(states3(:,1));  % mass differences
            mindif=min(state_diff);
            if mindif <= eps_threshold   % combine states that are equal to eps threshold.
                states4 = zeros(size(states3));
                state_count = 1;
                k=1;
                Ns = length(state_diff);
                while k <= Ns
                    if state_diff(k) <= eps_threshold
                        m1 = states3(k,1);
                        m2 = states3(k+1,1);
                        p1 = states3(k,2);
                        p2 = states3(k+1,2);
                        states4(state_count,1) = m1;
                        states4(state_count,2) = p1+p2;
                        state_count = state_count + 1;
                        k=k+1;
                    else
                        states4(state_count,1) = states3(k,1);
                        states4(state_count,2) = states3(k,2);
                        state_count = state_count + 1;
                    end
                    k = k + 1;
                end
                if state_diff(Ns) > eps_threshold
                    states4(state_count,1) = states3(Ns+1,1);
                    states4(state_count,2) = states3(Ns+1,2);
                    state_count = state_count + 1;
                end        
                states3 = states4(1:state_count-1,:);
            end
            %----------------
            % prune states
            %----------------
            Num = length(states3);
            if maxstates < Num
                states3 = sortrows(states3,-2);  % sort new states by probability
                states3 = states3(1:maxstates,:);   % keep only the most probably number of states (exact mass)
            end            
            Num = length(states3);
            display(['Element ' molecule.element{element_index(iel)}.symbol ' (' num2str(iel) ' of ' num2str(Nel) ')  Count ' num2str(ioa) ' of ' num2str(Noa) '  Current states: ' num2str(Num) ])
%            [iel Nel ioa Noa Num]
            %plot_states(states3)
            %pause
            %----------------
            % update states
            %----------------
            states1 = states3;  % the new states become the old states            
        end
    end 
end 
timev = toc;  % should make sure there are no printing inside loops
if maxstates == realmax
   display(['Time to process REALMAX maximum states took ' num2str(timev) ' seconds.'])
else
   display(['Time to process ' num2str(maxstates) ' maximum states took ' num2str(timev) ' seconds.'])
end

%states = sortrows(states3,1);  % sort by ascending mass
states = sortrows(states3,-2);  % sort by descending probability 
%timev = toc


