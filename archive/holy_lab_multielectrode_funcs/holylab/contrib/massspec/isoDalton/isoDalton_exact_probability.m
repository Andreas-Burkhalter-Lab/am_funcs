function states =  isoDalton_exact_probability(molecule_string,maxstates,massdiff_threshold,minprob_threshold)
%--------------------------------------------------------------------------------------------------
%
% Filename:     	isoDalton_exact_probability.m
% Description:  	The file returns the exact probability isotopic distributions of a molecule 
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
%                     to a crawl and run out of memory.
%
%         massdiff_threshold : mass difference threshold
%                     If the number of required states is greater than maxstates, then all states with
%                     a mass difference less than this threshold are merged.
%
%         minprob_threshold : minimum probability threshold
%                     If the number of required states is greater than maxstates, and the mass difference
%                     threshold case is not satisfied, then neighboring states with probabilities less than
%                     this threshold are merged.
%
%---------------------------------------------------------------------------------------------------
%
% Output:  	returns a two column matrix where the first column contains the masses
%           and the second column are the probabilities (merged, i.e. added).
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

if nargin == 2
    massprob_product = 1;
else
    massprob_product = 0;    
end

%---------------------------
% Get isotopes of molecule
%---------------------------
molecule = isoDalton_get_molecule_isotopes(molecule_string);


Nel = molecule.element_count;  % number of different elements in molecule
isotope_count = zeros(Nel,1);
for i=1:Nel
    isotope_count(i) = molecule.element{i}.non_zero_isotope_count;   % number of isotopes for each element with non zero composition values
end
[v,element_index] = sort(isotope_count);    % element_index is used to order the elements in increasing number of isotopes
                                            % we want to start with the elements with fewest isotopes to minimize growth of the trellis
                                            
Natoms = 0;
for iel=1:Nel  % elements (atoms) in molecule
    Natoms = Natoms + molecule.element{element_index(iel)}.number_of_atoms;
end

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

%pause
tic
%---------------------------------
% Setup initial states
%---------------------------------
Ni = molecule.element{element_index(1)}.non_zero_isotope_count;
states1=zeros(Ni,2);
for i=1:Ni
    isotope_index = molecule.element{element_index(1)}.non_zero_isotope_composition_index(i);
    states1(i,1)  = molecule.element{element_index(1)}.isotopic_relative_atomic_masses(isotope_index);
    states1(i,2)  = molecule.element{element_index(1)}.isotopic_composition(isotope_index);
end

max_prob = zeros(Natoms,1);
time_step = 1;
mp = max(states1(:,2));
max_prob(time_step) = mp;
time_step = time_step + 1;
states1(:,2) = states1(:,2)/mp;   % normalize state probabilities

%---------------------------------
% Trellis
%---------------------------------
for iel=1:Nel  % elements in molecule
    Noa = molecule.element{element_index(iel)}.number_of_atoms;
    for ioa=1:Noa % atoms of each element
        if ~(iel==1 && ioa==1)   % we already have these initial states
            Nst1 = length(states1(:,1));  % number of previous states
            Nic = molecule.element{element_index(iel)}.non_zero_isotope_count;   % isotope count
            states2=zeros(Nst1*Nic,2);
            states2_index = 1;
            for iic=1:Nic
                for ist=1:Nst1
                    isotope_index = molecule.element{element_index(iel)}.non_zero_isotope_composition_index(iic);
                    states2(states2_index,1) = states1(ist,1) + molecule.element{element_index(iel)}.isotopic_relative_atomic_masses(isotope_index); % mass
                    states2(states2_index,2) = states1(ist,2) * molecule.element{element_index(iel)}.isotopic_composition(isotope_index);  % probability
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
            if massprob_product == 1
                %---------------------------------------------------------------------------------------------------------------
                % absorb states with smallest massdiff probability product within an isotope peak cluster, i.e. massdiff < 0.5
                %
                % merging strategies implemented here........
                %
                %---------------------------------------------------------------------------------------------------------------
                Nst3 = length(states3(:,1));
                while Nst3 > maxstates 
                    states3 = sortrows(states3,1);  % sort based on acending masses
                    state_diff = diff(states3(:,1));  % mass differences
                    index = find(state_diff > 0.5);
                    mdp_product = state_diff .* states3(1:end-1,2) .* states3(2:end,2);  % massdiff probability product
                    mdp_product(index) = ones(length(index),1)*realmax;  % make the product huge for mass differences > 0.5
                    [mdp_min,min_index]=min(mdp_product);  %  min massdiff probability product (within clusters)
                    m1 = states3(min_index,1);
                    m2 = states3(min_index+1,1);
                    p1 = states3(min_index,2);
                    p2 = states3(min_index+1,2);
                    if p1 > 0 && p2 > 0  % worry about small numbers
                        psum = 1/(p1+p2);
                        if isfinite(psum)
                            p1n = p1*psum;
                            p2n = p2*psum;
                        else
                            lp1 = log10(p1);
                            lp2 = log10(p2);
                            lmp = min([lp1 lp2]);
                            lp1 = lp1 - lmp;
                            lp2 = lp2 - lmp;
                            p11 = 10^lp1;
                            p22 = 10^lp2;
                            psum = 1/(p11+p22);
                            p1n = p11*psum;
                            p2n = p22*psum;
                        end
                        states3(min_index,1)  = states3(min_index,1)*p1n + states3(min_index+1,1)*p2n;  % centroid mass
                        states3(min_index,2)  = p1 + p2;
                        states3(min_index+1,:)=[];
                    elseif p1 > 0 
                        states3(min_index+1,:)=[];    % just keep the p1 state and delete p2 since it is zero                
                    elseif p2 > 0 
                        states3(min_index,:)  =[];    % just keep the p2 state and delete p1 since it is zero                
                    else  % both zero
                        states3(min_index,:)  =[];    % eliminate both
                        states3(min_index+1,:)=[]; 
                        Nst3 = Nst3 - 1;   % we are eliminating two states so subtract one here                     
                    end
                    Nst3 = Nst3 - 1;
    %                pause
                end   
            else
                %---------------------------------------------------------------------------------------------------------------
                % absorb states with mass differences less than threshold or probabilities less than threshold
                %
                % merging strategies implemented here........  Note: code needs to be cleaned up.....
                %
                %---------------------------------------------------------------------------------------------------------------
                states1 = states3;  % the new states become the old states
                Nst1 = length(states1(:,1)); 
                combine_flag = 1;
                while combine_flag == 1
                    if maxstates < Nst1
                        states1 = sortrows(states1,1);  % sort based on acending masses
                        state_diff = diff(states1(:,1));  % mass differences
                        [mindif,index]=min(state_diff);
                        if mindif < massdiff_threshold   % combine states that are closer together than massdiff threshold.
                            m1 = states1(index,1);
                            m2 = states1(index+1,1);
                            p1 = states1(index,2);
                            p2 = states1(index+1,2);
                            psum = 1/(p1+p2);
                            p1n = p1*psum;
                            p2n = p2*psum;
                            states1(index,1) = states1(index,1)*p1n + states1(index+1,1)*p2n;
                            states1(index,2) = p1 + p2;
                            states1(index+1,:)=[];
                            Nst1 = Nst1 - 1;
                        else  % combine neighboring states that have combined probabilities less the minprob threshold
                            [minprob,pindex]=min(states1(:,2));
                            if pindex > 1
                                minprob1 = sum(states1(pindex-1:pindex,2));
                            else
                                minprob1 = 1;
                            end
                            if pindex < Nst1
                                minprob2 = sum(states1(pindex:pindex+1,2));
                            else
                                minprob2 = 1;
                            end
                            if minprob1 < minprob_threshold || minprob2 < minprob_threshold
                                if minprob1 <= minprob2
                                    pindex1 = pindex-1;
                                    pindex2 = pindex;
                                else
                                    pindex1 = pindex;
                                    pindex2 = pindex+1;
                                end
                                m1 = states1(pindex1,1);
                                m2 = states1(pindex2,1);
                                p1 = states1(pindex1,2);
                                p2 = states1(pindex2,2);
                                psum = p1+p2;
                                if psum > 0
                                    p3 = 1/psum;
                                    if isfinite(p3)                                    
                                        p1n = p1*p3;
                                        p2n = p2*p3;
                                        m = states1(pindex1,1)*p1n + states1(pindex2,1)*p2n;                                    
                                        states1(pindex1,1) = m;
                                        states1(pindex1,2) = psum;                                
                                        states1(pindex2,:)=[];
                                        Nst1 = Nst1 - 1;                                
                                    else
                                        states1(pindex1:pindex2,:)=[];  % they are both too small so toss them out
                                        Nst1 = Nst1 - 2;
                                    end
                                else
                                    states1(pindex1:pindex2,:)=[];  % they are both too small so toss them out
                                    Nst1 = Nst1 - 2;
                                end                            
                            else                        
                                display('Adding a state since the combine state criteria is not met');
                                states1
                                maxstates = maxstates + 1
                                mindif
                                Nst1
                                min(states1(:,1))
                                max(states1(:,1))
                                imin
                                imax
                                display('Paused here'); pause
                            end
                        end
    %                    pause
                    else
                        %states1 = states1(1:Nst1,:);  
                        break;
                    end
                end
                states3 = states1;
            end
            %[iel Nel ioa Noa]
            %plot_states(states3)
            %pause
%            minprob = min(states3(:,2))
%            pause
            %----------------
            % update states
            %----------------
            states1 = states3;  % the new states become the old states            
        end
    end 
end 
timev = toc;
display(['Time to process ' num2str(maxstates) ' maximum states took ' num2str(timev) ' seconds.'])


%states = sortrows(states3,1);  % sort by ascending mass
states = sortrows(states3,-2);  % sort by descending probability 


