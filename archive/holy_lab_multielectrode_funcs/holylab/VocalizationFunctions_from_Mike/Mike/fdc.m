function [h,bins,binwidth,g] = fdc(x,varargin)
%FDC Freedman-Diaconis Choice Histogram
%   Function fdc calculates a histogram based upon the Freedman-Diaconis
%   choice for optimal bin width. This is calculated as follows:

%   optimal_bin_width = (2*interquartile_range)/(N^(1/3))
%   where N = length(x);
%   and interquartile_range = iqr(x); 
%   The Freedman-Diaconis uses the interquartile range as the estimator for
%   optimal bin rather than the sample standard deviation.
%   The number of bins is then computed as:
%   nbins = ceil((max(x) - min(x))/optimal_bin_width)
%   where ceil() is the ceiling function.

%   Use:
%   fdc(x): Calculate a histogram with Freedman-Diaconis bins/bin width, and
%   plot as counts x bin. Input x is a unidimensional vector of length n.
%   Optional arguments:
%   fdc(x,'ystyle',ystyle)
%   Where ystyle is: 
%       'fraction' - fraction (0 to 1) of all observations
%       'count'   - # of observations (default)
%       'percent' - fraction * 100
%   fdc(x,'xstyle',xstyle)
%   Where xstyle is:
%       'semilog' - plots log10(nbins) with tickmark values from bins
%                                       (untransformed bin values)
%       'log' - plots log10(bins)
%       'full' - plots untransformed bins as x axis
%   fdc(x,'noplot')
%       calculates outputs without plotting graph
%   fdc(x, 'titletext',titletext)
%       Prints custom title.
%   fdc(x, 'measure',measure)
%       Name for variable being counted: counts, whistles, heights, weights,
%       etc.
%   fdc(x, 'unit', unit)
%       Name for the unit being used: kHz, seconds, volts, etc.
    
%  (2013) Michael A. Rieger

    if (isnumeric(x) && (size(x,1) == 1)) % Checks that x is a vector
       n = length(x); % Calculates N
       binwidth = (2*iqr(x))/(n^(1/3)); % Calculates FDC
       nbins = ceil((max(x) - min(x))/binwidth); % Calculates # of bins
       bins = (min(x)):binwidth:(max(x));
       if ((nargin > 1) && isempty(intersect(varargin,...
               {'ystyle','xstyle','noplot','titletext','measure','unit'})))
           % Checks that only proper arguments are specified
           error('Argument improperly specified.'); 
       elseif ((nargin > 1) && (sum(strcmp('noplot',varargin)) == 1))
           % Set dummy variable 'noplot' for later.
           h = hist(x,bins);
           noplot = 1;
           
       elseif ((nargin > 1) && (sum(strcmp('noplot',varargin)) == 0 ))
           
           h = hist(x,bins); % Calculates the histogram h and bin centers
           
           for i = 1:length(varargin)
           
               if strcmp('ystyle',varargin{i})
                   if strcmp('fraction',varargin{i+1})
                       ydata = h./sum(h);
                       ytext = 'Fraction of all observations';
                   elseif strcmp('count',varargin{i+1})
                       ydata = h;
                       ytext = 'Count';
                   elseif strcmp('percent',varargin{i+1})
                       ydata = (h./sum(h)).*100;
                       ytext = 'Percent of all observations';
                   else
                       error('ystyle improperly specified.');
                   end
               elseif strcmp('xstyle',varargin{i})
                   if strcmp('semilog',varargin{i+1})
                       xdata = log10(bins);
                       xticks = []; % Makes a variable xticks to look for later
                   elseif strcmp('log',varargin{i+1})
                       xdata = log10(bins);                       
                   elseif strcmp('full',varargin{i+1})
                       xdata = bins;
                   else
                       error('xstyle improperly specified.');
                   end
               elseif strcmp('titletext',varargin{i})
                   titletext = varargin{i+1};
               elseif strcmp('measure',varargin{i})
                   measure = varargin{i+1};
               elseif strcmp('unit',varargin{i})
                   unit = varargin{i+1};
               end
           
           end
           
       elseif (nargin == 1)
           
           h = hist(x,bins);
           ydata = h;
           xdata = bins;
       
       end
    else
       error('Sorry, input must be vector.');
    end

    
    % Make Plot
    if ~exist('noplot','var') % if noplot is a variable, then don't plot
        if ~exist('xdata','var')
           xdata = bins;
        end
    
        if ~exist('ydata','var')
            ydata = h;
        end
        g = figure;
        bar(xdata,ydata); 
        
        if exist('titletext','var')
         title({titletext; ['N = ', num2str(sum(h))];...
            ['Bin number = ',num2str(nbins)]});
        else
        title({['N = ',num2str(sum(h))];...
            ['Bin number  = ',num2str(nbins)]});
        end
    
         if exist('xticks','var')
          xticks = get(gca,'XTick');
          xticks = 10.^(xticks);
          xticks = round(xticks.*1000)/1000;
          set(gca,'XTickLabel',xticks);
         end
    
        if exist('unit','var') && exist('measure','var')
            xlabel({measure;...
            ['Bin width = ',num2str(round(binwidth*1000)/1000),' ',unit]});
        elseif exist('unit','var') && ~exist('measure','var')
            xlabel(['Bin width = ',num2str(round(binwidth*1000)/1000),' ',unit]);
        else
            xlabel(['Bin width = ',num2str(round(binwidth*1000)/1000)]);
        end
    
        if exist('ytext','var')
            ylabel(ytext);
        else
            ylabel('Count');
        end
    end
end
