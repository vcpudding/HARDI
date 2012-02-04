function div=vectorDivergence(varargin)
%DIVERGENCE  Divergence of a vector field.
%   DIV = DIVERGENCE(X,Y,Z,U,V,W) computes the divergence of a 3-D
%   vector field U,V,W. The arrays X,Y,Z define the coordinates for
%   U,V,W and must be monotonic and 3-D plaid (as if produced by
%   MESHGRID).
%   
%   DIV = DIVERGENCE(U,V,W) assumes 
%         [X Y Z] = meshgrid(1:N, 1:M, 1:P) where [M,N,P]=SIZE(U). 
%
%   DIV = DIVERGENCE(X,Y,U,V) computes the divergence of a 2-D
%   vector field U,V. The arrays X,Y define the coordinates for U,V
%   and must be monotonic and 2-D plaid (as if produced by
%   MESHGRID). 
%   
%   DIV = DIVERGENCE(U,V) assumes 
%         [X Y] = meshgrid(1:N, 1:M) where [M,N]=SIZE(U). 
%   
%   Example:
%      load wind
%      div = divergence(x,y,z,u,v,w);
%      slice(x,y,z,div,[90 134],[59],[0]); shading interp
%      daspect([1 1 1])
%      camlight
%
%   See also STREAMTUBE, CURL, ISOSURFACE.

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 1.4.4.1 $  $Date: 2005/04/28 19:56:16 $

error(nargchk(2,3,nargin,'struct'));
[u v w] = parseargs(nargin,varargin);

% Take this out when other data types are handled
u = double(u);
v = double(v);
w = double(w);

if isempty(w)  % 2D

  [px junk] = vectorGradient(u); %#ok
  [junk qy] = vectorGradient(v); %#ok
  div = px+qy;
  
else   %3D
    
  [px junk junk] = vectorGradient(u); %#ok
  [junk qy junk] = vectorGradient(v); %#ok
  [junk junk rz] = vectorGradient(w); %#ok
  
  div = px+qy+rz;
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, v, w] = parseargs(nin, vargin)

w = [];

if nin==2         % divergence(u,v)
  u = vargin{1};
  v = vargin{2};
elseif nin==3     % divergence(u,v,w)
  u = vargin{1};
  v = vargin{2};
  w = vargin{3};
else
  error('MATLAB:divergence:WrongNumberOfInputs',...
        'Wrong number of input arguments.'); 
end
