function [dist,PP0] = pointTriangleDistance(TRI,P)
% calculate distance between a point and a triangle in 3D
%
% Copyright (c) 2009, Gwendolyn Fischer
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% SYNTAX
%   dist = pointTriangleDistance(TRI,P)
%   [dist,PP0] = pointTriangleDistance(TRI,P)
%
% DESCRIPTION
%   Calculate the distance of a given point P from a triangle TRI.
%   Point P is a row vector of the form 1x3. The triangle is a matrix
%   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
%   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
%   to the triangle TRI.
%   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
%   closest point PP0 to P on the triangle TRI.
%
% Author: Gwendolyn Fischer
% Release: 1.0
% Release date: 09/02/02
% Release: 1.1 Fixed Bug because of normalization

% Possible extention could be a version tailored not to return the distance
% and additionally the closest point, but instead return only the closest
% point. Could lead to a small speed gain.

% Example:
% %% The Problem
% P0 = [0.5 -0.3 0.5];
% 
% P1 = [0 -1 0];
% P2 = [1  0 0];
% P3 = [0  0 0];
% 
% vertices = [P1; P2; P3];
% faces = [1 2 3];
% 
% %% The Engine
% [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0);
%
% %% Visualization
% [x,y,z] = sphere(20);
% x = dist*x+P0(1);
% y = dist*y+P0(2);
% z = dist*z+P0(3);
% 
% figure
% hold all
% patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8);
% plot3(P0(1),P0(2),P0(3),'b*');
% plot3(PP0(1),PP0(2),PP0(3),'*g')
% surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
% view(3)

% The algorithm is based on 
% "David Eberly, 'Distance Between Point and Triangle in 3D',
% Geometric Tools, LLC, (1999)"
% http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
%
%        ^t
%  \     |
%   \reg2|
%    \   |
%     \  |
%      \ |
%       \|
%        *P2
%        |\
%        | \
%  reg3  |  \ reg1
%        |   \
%        |reg0\ 
%        |     \ 
%        |      \ P1
% -------*-------*------->s
%        |P0      \ 
%  reg4  | reg5    \ reg6



% rewrite triangle in normal form
B = TRI(1,:);
E0 = TRI(2,:)-B;
E1 = TRI(3,:)-B;


D = B - P;
a = dot(E0,E0);
b = dot(E0,E1);
c = dot(E1,E1);
d = dot(E0,D);
e = dot(E1,D);
f = dot(D,D);

det = a*c - b*b;
s   = b*e - c*d;
t   = b*d - a*e;

% Terible tree of conditionals to determine in which region of the diagram
% shown above the projection of the point into the triangle-plane lies.
if (s+t) <= det
  if s < 0
    if t < 0
      %region4
      if (d < 0)
        t = 0;
        if (-d >= a)
          s = 1;
        else
          s = -d/a;
        end
      else
        s = 0;
        if (e >= 0)
          t = 0;
        else
          if (-e >= c)
            t = 1;
          else
            t = -e/c;
          end
        end
      end %of region 4
    else
      % region 3
      s = 0;
      if e >= 0
        t = 0;
      else
        if -e >= c
          t = 1;
        else
          t = -e/c;
        end
      end
    end %of region 3 
  else
    if t < 0
      % region 5
      t = 0;
      if d >= 0
        s = 0;
      else
        if -d >= a
          s = 1;
        else
          s = -d/a;
        end
      end
    else
      % region 0
      invDet = 1/det;
      s = s*invDet;
      t = t*invDet;
    end
  end
else
  if s < 0
    % region 2
    tmp0 = b + d;
    tmp1 = c + e;
    if tmp1 > tmp0 % minimum on edge s+t=1
      numer = tmp1 - tmp0;
      denom = a - 2*b + c;
      if numer >= denom
        s = 1;
        t = 0;
      else
        s = numer/denom;
        t = 1-s;
      end
    else          % minimum on edge s=0
      s = 0;
      if tmp1 <= 0
        t = 1;
      else
        if e >= 0
          t = 0;
        else
          t = -e/c;
        end
      end
    end %of region 2
  else
    if t < 0
      %region6 
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
        numer = tmp1 - tmp0;
        denom = a-2*b+c;
        if (numer >= denom)
          t = 1;
          s = 0;
        else
          t = numer/denom;
          s = 1 - t;
        end
      else  
        t = 0;
        if (tmp1 <= 0)
            s = 1;
        else
          if (d >= 0)
              s = 0;
          else
              s = -d/a;
          end
        end
      end
      %end region 6
    else
      % region 1
      numer = c + e - b - d;
      if numer <= 0
        s = 0;
        t = 1;
      else
        denom = a - 2*b + c;
        if numer >= denom
          s = 1;
          t = 0;
        else
          s = numer/denom;
          t = 1-s;
        end
      end %of region 1
    end
  end
end

PP0 = B + s*E0 + t*E1;

dist = norm(P - PP0);