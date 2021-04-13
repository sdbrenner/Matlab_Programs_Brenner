function SA = struct2structarray(S)
% STRUCT2STRUCTARRAY re-organizes a regular structure into a
% structure-array format
%
%   SA = struct2structarray(S) converts from a scalar structure S,
%   containing a number of fields of size 1xN, or 1x1 to a
%   structure-array SA of size 1xN containing the fields of S.
%
%   There are likely restrictions on field sizes and data types in S. This
%   function has not been extensively tested.
%
%   S.D.Brenner, 2020

% Based on answers from MatlabCentral:
% https://www.mathworks.com/matlabcentral/answers/46236-efficient-way-to-change-a-struct-of-mixed-type-arrays-into-a-structure-array

B = [fieldnames(S).';cellfun(@num2cell,struct2cell(S).','un',0)];
SA = struct(B{:});