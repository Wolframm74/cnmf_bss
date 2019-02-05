function [Y_lk]=projfunc_wrapper_matlab_Ylk(Y_lk)

[L_static, K_static]=size(Y_lk);

% sparseness=0.6;

compute_flag=0;

before_sparseness_sum=0;
after_sparseness_sum=0;

for l=1:L_static

    %z_col_l=squeeze(Z_ol(:,l));
    %y_col_k=squeeze(Y_lk(:,k));
    y_row_l=squeeze(Y_lk(l, :));
    k2=norm(y_row_l,2);
    
    k1_current=norm(y_row_l, 1);
    
    current_sparseness=(sqrt(K_static)-(k1_current/k2))/(sqrt(K_static)-1);
    
    %if the current sparseness is less than 0.75, then compute a sparseness
    %that is 10% sparser than its current sparseness. 
    
    target_sparseness=0.5;
    sparseness_diff=0.01;    
    
    if (current_sparseness<=(target_sparseness-sparseness_diff))
        
            sparseness_to_pass=current_sparseness+sparseness_diff;
            compute_flag=1;
            
    elseif (current_sparseness <target_sparseness)
            
            sparseness_to_pass=target_sparseness;
            compute_flag=1;
            
    elseif (current_sparseness==target_sparseness)
        
            compute_flag=0;
        
    elseif (current_sparseness>target_sparseness)
            
            sparseness_to_pass=target_sparseness;
            compute_flag=1; 
            
    end    
            
    %compute the l1 norm that is required to achieve the desired sparseness
    %value
    
    if (compute_flag)

    before_sparseness_sum=before_sparseness_sum+current_sparseness;
    after_sparseness_sum=after_sparseness_sum+sparseness_to_pass;        
        
    k1=k2*(sqrt(K_static)-sparseness_to_pass*(sqrt(K_static)-1));
    
    y_col_l=transpose(y_row_l);
    
    [v, usediters]=projfunc_local(y_col_l, k1, k2, 1);
    
    %Z_ol(:,l)=v;
    Y_lk(l,:)=transpose(v);
    
    end
    
end     

input_sparseness_avg=before_sparseness_sum/L_static
output_sparseness_avg=after_sparseness_sum/L_static

%Y_lk

end

function [v,usediters] = projfunc_local( s, k1, k2, nn )

% Solves the following problem:
% Given a vector s, find the vector v having sum(abs(v))=k1 
% and sum(v.^2)=k2 which is closest to s in the euclidian sense.
% If the binary flag nn is set, the vector v is additionally
% restricted to being non-negative (v>=0).
%    
% Written 2.7.2004 by Patrik O. Hoyer
%
    
% Problem dimension
N = length(s);

% If non-negativity flag not set, record signs and take abs
if ~nn,
    isneg = s<0;
    s = abs(s);
end

% Start by projecting the point to the sum constraint hyperplane
v = s + (k1-sum(s))/N;

% Initialize zerocoeff (initially, no elements are assumed zero)
zerocoeff = [];

j = 0;
while 1,

    % This does the proposed projection operator
    midpoint = ones(N,1)*k1/(N-length(zerocoeff));
    midpoint(zerocoeff) = 0;
    w = v-midpoint;
    a = sum(w.^2);
    b = 2*w'*v;
    c = sum(v.^2)-k2;
    alphap = (-b+real(sqrt(b^2-4*a*c)))/(2*a);
    v = alphap*w + v;
    
    if all(v>=0),
	% We've found our solution
	usediters = j+1;
	break;
    end
        
    j = j+1;
        
    % Set negs to zero, subtract appropriate amount from rest
    zerocoeff = find(v<=0);
    v(zerocoeff) = 0;
    tempsum = sum(v);
    v = v + (k1-tempsum)/(N-length(zerocoeff));
    v(zerocoeff) = 0;
            
end

% If non-negativity flag not set, return signs to solution
if ~nn,
    v = (-2*isneg + 1).*v;
end

% Check for problems
if max(max(abs(imag(v))))>1e-10,
    error('Somehow got imaginary values!');
end

end

