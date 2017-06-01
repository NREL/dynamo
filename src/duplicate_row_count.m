function [d_points, d_values]=duplicate_row_count(points, values)
   %duplicate_row_count takes as inputs:
   %            - points: array of points
   %            - values: array of values       
   %The function looks for duplicate points in the list and count the
   %number of occurence of these points. The function also needs the
   %list of values such that it can split the value list in a consistent way.    
   % RETURNS:
   %         - d_points: 2D array where each line is a point and there is no
   %                     duplicate points. For each line there is an extra column
   %                     at the end that gives the number of occurence of the point
   %                     in the initial array.
   %         - d_values: Cell array that contains as many cells as unique points.
   %                     d_values{k} is an array storing all the values of
   %                     the unique point in row k of d_points.
   %                     
   %     See also update_scheduler
   %     
   %     Code by Nicolas Gensollen, 06/01/2017
   %
   % HISTORY
   % ver     date    time       who     changes made
   % ---  ---------- -----  ----------- ---------------------------------------
   % 0.1  2017-06-01 11:00  NicolasG     Initial Version adapted from Dynamo Python
   
   [unique_points, ~, J]=unique(points, 'rows');
   Count=accumarray(J,1);
   d_points=[unique_points Count];
   d_values={};
   for k=1:size(unique_points,1)
       d_values{k}=[];
       for idx=1:size(points, 1)
           if points(idx,:)==unique_points(k,:)
               d_values{k}=[d_values{k} values(idx)];
           end
       end
   end
           
end