%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	This file is part of the Snoopy code.

%   Snoopy code is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    Snoopy code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the data 

clear all;
%full_file=importdata('rates_16_6/rates_r00_copy.txt');
full_file=importdata('../../rates_r00.txt');

timevar=full_file.data;
nvar=size(timevar,2);
%var_name=strread(full_file.textdata{2},'%s',nvar);

% Assign variables to columns of data 
%for i=1:nvar
   %assignin('base',var_name{i},timevar(:,i)); 
%end

% Create a table with all variables
%tblA = table(KX,KY,KZ,TIME,EM,GR);

tblA = table(timevar(:,1),timevar(:,2),timevar(:,3),timevar(:,4),timevar(:,5),timevar(:,6));
KX=timevar(:,1);
% Sort the rows of the table
tblB = sortrows(tblA);