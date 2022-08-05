% Plotting and graphical utility functions
%
% Zooming:
%   sliderwindow  - GUI for zooming in on graph.
%   slidwincmenu  - Put a context menu on a graph to enable sliderwindow.
%
% Axis manipulation:
%   SplitVert     - Split current axis vertically into several axes.
%   SplitHoriz    - Split current axis horizontally into several axes.
%   ghostaxis     - Create a convenient axis for annotation.
%   data2norm     - Convert data units into figure-normal units (annotation).
%   LetterFigures - Put parts labels ("A", "B", etc.) into figure.
%   CalcSubplotDims - Determine grid of subplots, optional overall labels.
%   raisexlabel   - Shift all x labels by given # of pixels.
%   vertsubplot   - Like subplot, except axes indexed vertically.
%   tickleneq     - Equal-length tick marks across axes in a figure.
%
% Color manipulation:
%   lightencolors - Generate pale versions of supplied colors.
%
% Graphical export:
%   exportfig     - Export a figure (PS, TIFF).
%   applytofig    - Apply export options to figure on screen.
%   previewfig    - Preview figure with export options.
%   restorefig    - Restore a figure's properties.
%
% Plotting functions and utilities:
%   errorbarlog   - Error bar plots for logarithmic axes.
%   xerrorbar     - Put error bars on x coordinate.
%   fillmm2       - Fill region between min and max.
%   PlotDots      - A more flexible version of MATLAB's dotted lines
%   points        - Scatter plot with intensity code.
%   addbreaks     - Break lines at gaps in data.
%   image_textmask  - Write text into images. 
%
% Graphical object selection:
%   SelectPoints  - Get handles of points in scatterplot.
%   GetSelRect    - Get a selection rectangle in axis coordinates.
%   moveline      - Mouse dragging of a vertical line (e.g., set threshold)
%   selectcircle  - Create a selection circle
%
% GUIs:
%   appdatadefault - supply default values for appdata
