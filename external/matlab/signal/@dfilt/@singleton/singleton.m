function varargout = singleton(varargin)
%SINGLETON Singleton filter virtual class.
%   SINGLETON is a virtual class---it is never intended to be instantiated.
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.

error(message('signal:dfilt:singleton:singleton:SigErr', 'SINGLETON')); 
