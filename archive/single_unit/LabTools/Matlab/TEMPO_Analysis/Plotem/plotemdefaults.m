% plotemdefaults.m
% Script that specifies default paths and parameters for the
%   function plotem.m.
% HOW TO USE:  Under "Setting Default Parameters", under the 
%   condition that specifies the data you want to analyze (e.g. 
%   if filename starts with 'ds'), specify each variable; if 
%   you don't want to specify it, set it to the default, as 
%   is shown in the last ('else') condition.  Nothing else usually
%   needs to be changed, unless adding a new paradigm, in which
%   case add the strings that precede the x and y positions of
%   the fixation point in the .s file under "Setting Paradigm
%   Specific Variables".
% VARIABLES:
%   datapath		cell array of paths to data files (*.s, 
%			*.ras, *.e1h, etc); format:  {'path1';'path2'};
%			can be any length
%   indpath		path to indices file for paradigm used
%   regflg		flag for whether to attempt to regularize the
%   			data; can have a value of 0, 1 or 2.
%				0 = don't regularize the data
%				1 = find regular samples for whole trial
%				2 = find longest stretch of regularly
%				spaced samples (used if regularizing
%				whole trial is impossible)
%   truncatestring	event to truncate data at; should be one of
%				the events present in eventnames below.
%   tzerostring		event to set time=0 to; should be one of the
%				events present in eventnames below.
%   algsetname		name of variable containing a set of parameters
%				for the saccade-picking algorithm (msacalg.m);
%				program looks for this variable in the file
%				algparamsets.mat in /HOME/crista/matlab/eyes.
%   defaulttrlsperax	number of trials to plot in each axes on the screen;
%				range 1 to 20.
%   defaultaxesmode	block structure to use in plotting data; options include:
%				1 = 1 trial at a time
%				2 = 1 block of n trials at a time
%				3 = 2 blocks of n trials each
%				4 = 3 blocks of n trials each
%   eventnames={'None';
%	    'SIN_FPON';
%	    'SIN_FIX';
%	    'SIN_STON';
%	    'SIN_STOFF';
%	    'SIN_TRGON';
%	    'SIN_TRGOFF';
%	    'SIN_FPOFF';
%	    'SIN_SAC';
%	    'SIN_TRGAC';
%	    'SIN_DTON';
%	    'SIN_DTOFF';
%	    'SIN_SACST';
%	    'SIN_SACEN'};
% Written by Crista, 9/25/98.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Getting Variable to Test Against %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Getting first two letters of file name, i.e. monkey name.
if length(filename)>=2
	monkey=filename(1:2);
else
	monkey='xx';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Default Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% MONKEY:  ben 
if strcmp(monkey,'be')
	datapath={'Z:\Data\Rex\Ben\spikes'};
	indpath='Z:\LabTools\Matlab\PlotEM';
	regflg=2;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
   tempo=0;

%%% MONKEY:  Radar
elseif strcmp(monkey,'rd')
	datapath={'Z:\Data\Rex\Radar\spikes'};
	indpath='Z:\LabTools\Matlab\PlotEM';
	regflg=2;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
	defaultaxesmode=3;
   tempo=0;
   
%%% MONKEY:  Wally
elseif strcmp(monkey,'w0')
	datapath={'Z:\Data\Rex\Wally\spikes'};
	indpath='Z:\LabTools\Matlab\PlotEM';
	regflg=2;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
   tempo=0;
   
%%% MONKEY:  TEMPO TEST MONKEY
elseif strcmp(monkey,'m0')
	datapath={'Z:\Data\Tempo\Test\'};			%need \ at end of paths for tempo data
	indpath='Z:\LabTools\Matlab\PlotEM\';
	regflg=0;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
	tempo=1;
   
%%% MONKEY:  TEMPO Ben the Monkey
elseif strcmp(monkey,'m1')
   %datapath={'F:\'};
  	datapath={'Z:\Data\Tempo\Ben\Raw\'};			%need \ at end of paths for tempo data
	indpath='Z:\LabTools\Matlab\PlotEM\';
	regflg=0;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
   tempo=1;
   
%%% MONKEY:  TEMPO Jerry the Monkey
elseif strcmp(monkey,'m2')
 %  datapath={'E:\DATA\'};	
	datapath={'Z:\Data\Tempo\Jerry\Raw\'};			%need \ at end of paths for tempo data
	indpath='\\RHONE\LabTools\Matlab\PlotEM\';
	regflg=0;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
	tempo=1;
  %%% MONKEY:  TEMPO Baskin the Monkey
elseif strcmp(monkey,'m4')
 %  datapath={'E:\DATA\'};	
	datapath={'Z:\Data\Tempo\Baskin\Raw\'};			%need \ at end of paths for tempo data
	indpath='\\RHONE\LabTools\Matlab\PlotEM\';
	regflg=0;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
	tempo=1;
  %%% MONKEY:  TEMPO Robbins the Monkey
elseif strcmp(monkey,'m3')
 %  datapath={'E:\DATA\'};	
	datapath={'Z:\Data\Tempo\Robbins\Raw\'};			%need \ at end of paths for tempo data
	indpath='\\RHONE\LabTools\Matlab\PlotEM\';
	regflg=0;
	truncatestring='SIN_STOFF';
	truncatestring2='SIN_STON';
	tzerostring='SIN_STON';
	algsetname='algparams3';
	defaulttrlsperax=5;
   defaultaxesmode=3;
	tempo=1;  
%%% NO MONKEY SPECIFIED
else
	datapath={};
	indpath=[];
	regflg=0;
	truncatestring='None';
	truncatestring2='None';
	tzerostring='None';
	algsetname=[];
	defaulttrlsperax=5;
	defaultaxesmode=3;
   tempo = 0;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Paradigm Specific Variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tempo == 0
	fpstrings363={'FP X-Coord:';'FP Y-Coord:'};
	fpstrings364={'FP X-Coord:';'FP Y-Coord:'};
	fpstrings502={'Fx';'Fy'};
	fpstrings667={'Fx';'Fy'};
	fpstrings802={'FPXCTR:';'FPYCTR:'};
	fpstrings803={'FPXCTR:';'FPYCTR:'};
	fpstrings805={'FPXCTR:';'FPYCTR:'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Setting Other Paths %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

algsetpath='Z:\LabTools\Matlab\PlotEM\algparamsets';

