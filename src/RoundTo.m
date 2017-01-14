function [ a ] = RoundTo( a, tol )
%ROUNDTO Rounds the provided matrix to the desired tolerance
%
% Usage: [ a ] = RoundTo( a, tol )
%
% example: RoundTo(a, 0.001)
%
% see also: round
%
% by Bryan Palmintier, 2010

a = round(a./tol).*tol;

end
