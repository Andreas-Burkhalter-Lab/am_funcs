classdef xmlNode
  % xmlNode: represents a node in an XML document
  properties (SetAccess = public, GetAccess = public)
    name = '';
    attributes = [];
    content = '';
  end
  methods
    %% Constructor
    function n = xmlNode(str)
      n.name = str;
    end
    %% Has attribute
    function tf = has_attribute(node,str)
      if isempty(node.attributes)
        tf = false;
      else
        flag = strcmp(str,{node.attributes.name});
        tf = sum(flag) > 0;
      end
    end
    %% Attribute value
    function v = value(node,str)
      if ~isempty(node.attributes)
        flag = strcmp(str,{node.attributes.name});
        if sum(flag) ~= 1
          error(['Did not find unique attribute with name ' str]);
        end
        v = node.attributes(flag).value;
      end
    end
  end
end
