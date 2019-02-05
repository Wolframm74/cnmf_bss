function [T_fk]=projfunc_wrapper_matlab_Tfk(T_fk)

display('INSIDE PROJFUNC_WRAPPER_MATLAB_T_fk!!!!!!!!!!!!!!!1');

[F_static, K_static]=size(T_fk)

T_fk=rand(F_static,K_static);

% sparseness=0.6;

compute_flag=0;

before_sparseness_sum=0;
after_sparseness_sum=0;

for k=1:K_static

    display([ 'Line 16: The current value of k= ', num2str(k) ]); 
    
    t_col_k=squeeze(T_fk(:,k))

    figure, plot(t_col_k)
    title(['t_col_k: k=1', num2str(k)]);
    
    k2=norm(t_col_k,2);
    
    k1_current=norm(t_col_k, 1);
    
    current_sparseness=(sqrt(F_static)-(k1_current/k2))/(sqrt(F_static)-1)
    
    %if the current sparseness is less than 0.6, then compute a sparseness
    %that is 10% sparser than its current sparseness. 

    display([ 'Line 28: The current value of k= ', num2str(k) ]);     
    
    if (current_sparseness<=0.5)
        
            target_sparseness=current_sparseness+0.1;
            compute_flag=1;
            
    elseif (current_sparseness <0.6)
            
            target_sparseness=0.6;
            compute_flag=1;
            
    elseif (current_sparseness==0.6)
        
            compute_flag=0;
        
    elseif (current_sparseness>0.6)
            
            target_sparseness=0.6;
            compute_flag=1; 
            
    end    
    
    display([ 'Line 51: The current value of k= ', num2str(k) ]); 
        
    %compute the l1 norm that is required to achieve the desired sparseness
    %value
    
    if (compute_flag)

    display([ 'Line 58: The current value of k= ', num2str(k) ]);         
        
    before_sparseness_sum=before_sparseness_sum+current_sparseness;
    after_sparseness_sum=after_sparseness_sum+target_sparseness;        
        
    k1=k2*(sqrt(F_static)-target_sparseness*(sqrt(F_static)-1));

    display([ 'Line 65: The current value of k= ', num2str(k) ]);     
    
    [v, usediters]=projfunc_local(t_col_k, k1, k2, 1);

    display([ 'Line 69: The current value of k= ', num2str(k) ]);     
    
    figure, plot(v)
    title(['v_vector: k=1', num2str(k)]);
    
    T_fk(:,k)=v;
    
    end

    display([ 'Line 75: The current value of k= ', num2str(k) ]); 
    
end     

input_sparseness_avg=before_sparseness_sum/K_static
output_sparseness_avg=after_sparseness_sum/K_static

%T_fk

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

