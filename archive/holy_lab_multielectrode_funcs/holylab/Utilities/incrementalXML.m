classdef incrementalXML
  % incrementalXML: parse XML files line-by-line, to avoid memory problems in large files
  properties (SetAccess = private,GetAccess = public)
    filename = '';
    open_node = {};
    index = 0;
  end
  properties (SetAccess = private,GetAccess = private)
    fid = -1;
  end
  methods
    %% Constructor
    function ixml = incrementalXML(fn)
      ixml.filename = fn;
      [ixml.fid,msg] = fopen(fn,'r');
      if (ixml.fid < 0)
        error(msg)
      end
      % Check that this is an xml file
      l = fgetl(ixml.fid);
      if (l == -1)
        error('Empty file');
      end
      if ~strncmp(l,'<?xml',5)
        error('Not an xml file');
      end
    end % constructor
    %% Return the currently-building node
    function tmp = current(ixml)
      tmp = ixml.open_node{ixml.index};
    end
    %% Read another line
    function [ixml,cnode] = next(ixml)
      l = fgetl(ixml.fid);
      if (l == -1)
        % end of file
        ixml = -1;
        cnode = [];
        return
      end
      l = strtrim(l);  % remove whitespace
      if l(1) ~= '<'
        error('Malformed XML file')
      end
      matchstart = 2;
      % Read the node name (the first space- or >-delimited word)
      matchend = regexp(l,'[\s>]','once');
      if isempty(matchend)
        error('Empty delimiter match1');
      end
      name = l(matchstart:matchend-1);
      % If the name closes a node, return the node and move the pointer
      % back
      if name(1) == '/'
        cnode = ixml.open_node{ixml.index};
        if ~strcmp(cnode.name,name(2:end))
          error('Names don''t match')
        end
        ixml.index = ixml.index-1;
        return
      end
      % Start a new node
      tmp = xmlNode(name);
      l = l(matchend:end);
      % Parse attributes
      [attrstart,attrend] = regexp(l,'\s[^>]*?="[^>]*?"');
      if ~isempty(attrstart)
        tmp.attributes = repmat(struct('name','','value',''),1,length(attrstart));
        for i = 1:length(attrstart)
          thisattr = l(attrstart(i)+1:attrend(i));
          indx = regexp(thisattr,'=','once');
          tmp.attributes(i).name = thisattr(1:indx-1);
          tmp.attributes(i).value = thisattr(indx+2:end-1);
        end
      end
      if ~isempty(attrend)
        l = l(attrend(end)+1:end);
      end
      % Look for end of node start
      if strcmp(l,'>')
        % Append the new node
        ixml.index = ixml.index+1;
        ixml.open_node{ixml.index} = tmp;
        cnode = [];
        return
      end
      % Look for closure
      if strcmp(l,' />')
        cnode = tmp;
        return
      end
      % Look for content
      if l(1) == '>'
        matchstart = regexp(l,'<','once');
        tmp.content = l(2:matchstart);
        cnode = tmp;
        % Trust that the closure is all that follows
        return;
      end
      error('Something unexpected happened')
    end
    %% Advance to the opening of a particular block
    function [ixml,cnode] = advance(ixml,str)
      cnode = [];
      while true
        if isequal(ixml,-1)
          error(['End of file reached without finding ' str])
        end
        if ixml.index ~= 0 && strcmp(ixml.open_node{ixml.index}.name,str)
          break
        end
        if ~isempty(cnode) && strcmp(cnode.name,str)
          break
        end
        [ixml,cnode] = ixml.next;
      end
    end
  end
end
