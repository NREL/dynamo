function [queries_points, queries_values]=update_scheduler(points, values)
       % update_scheduler takes as INPUTS:
       %      - points: array of points
       %      - values: array of values
       %     
       % If there is no duplicates in points, it basically does nothing at all
       % because the update can be done in a single query.
       % If there are duplicates, the function gives a warning, but split the
       % update in multiple update queries. Each sub-query will not contain
       % duplicates. The purpose of this function is to make sure that, if
       % some points appear multiple times, they will be considered with appropriate
       % step size values.
       % 
       % RETURNS:
       %     - queries_points: Cell array of arrays holding the points
       %     - queries_values: Cell array of arrays holding the values
       %     
       % USAGE: >>>[p_list,v_list]=update_scheduler(points, values)
       %        >>>for idx=1:size(p_list,2) 
       %        ...     my_func.update(p_list{idx}, v_list{idx})
       %        
       % See also duplicate_row_count
       % 
       % Code by Nicolas Gensollen, 06/01/2017
       % HISTORY
       % ver     date    time       who     changes made
       % ---  ---------- -----  ----------- ---------------------------------------
       % 0.1  2017-06-01 13:00  NicolasG    Initial Version adapted from Dynamo Python
       
       [counts, vals]=duplicate_row_count(points, values);
       
       %If we have no duplicate, then simply return
       %the arguments
       if all(counts(:,end)==1)
           queries_points=points;
           queries_values=values;
           return
       end
       
       %Otherwise, there are duplicates
       %Give a warning first
       warning('ADP:FaDiscrete:Duplicate points in update')
       
       queries_points={};
       queries_values={};
       query_idx=0;
       
       new_query=true;
       
       while new_query==true
           query_idx=query_idx+1;
           queries_points{query_idx}=zeros(0,size(points,2));
           queries_values{query_idx}=[];
           new_query=false;
           for idx=1:size(counts,1)
               if counts(idx,end)>0
                   new_query=true;
                   queries_points{query_idx}=[queries_points{query_idx}; counts(idx,1:end-1)];   
                   queries_values{query_idx}=[queries_values{query_idx}; vals{idx}(end)];
                   counts(idx,end)=counts(idx,end)-1;
                   vals{idx}=vals{idx}(1:end-1);
               end
           end
       end
       queries_points=queries_points(1:end-1);
       queries_values=queries_values(1:end-1);
                       
end