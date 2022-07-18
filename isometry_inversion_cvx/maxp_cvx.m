function [maxp S F]=maxp_cvx(din, dout, d, D, Jin, Jout, k, choose_protocol, isComplex)
% Code that evaluates the maximal probability of transforming the set of
% k isometries encoded in Choi operators Jin and obtain the isometry
% operators Jout in parallel/sequential/general protocol.
% 'din' corresponds to the input dimension of the output operation
% 'dout' corresponds to the input dimension of the output operation
% Jin and Jout are tenors with 'n' Choi operators that should be written as
% Jin(:,:,i), where the variable 'i' ranges from 1 to 'n'
% 'k' corresponds to the number of uses of the input-operation
% if isComplex is '0', the code only consider supermaps with real numbers
% is isComplex is '1', the code considers general supermaps with complex
% numbers

n_in=size(Jin,3);   %Count the number of choi operators
     cvx_begin SDP quiet
    %cvx_solver MOSEK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 1 copy
if k==1
    %Declare the SDP variable related to the probability p
    variable p
    % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*dout,din*d*D*dout) complex semidefinite
        variable F(din*d*D*dout,din*d*D*dout) complex semidefinite
    else
        variable S(din*d*D*dout,din*d*D*dout) semidefinite
        variable F(din*d*D*dout,din*d*D*dout) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        PartialTrace(S * Tensor(eye(din) , transpose(Jin(:,:,i)) , eye(dout)) ,[2 3],[din d D dout]) == p*Jout(:,:,i);
    end
end %end if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 2 copies
if k==2
    %Declare the SDP variable related to the probability p
    variable p
   % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*d*D*dout,din*d*D*d*D*dout) complex semidefinite
        variable F(din*d*D*d*D*dout,din*d*D*d*D*dout) complex semidefinite
    else
        variable S(din*d*D*d*D*dout,din*d*D*d*D*dout) semidefinite
        variable F(din*d*D*d*D*dout,din*d*D*d*D*dout) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout)) ,[2 3 4 5],[din d D d D dout]) == p*Jout(:,:,i);
    end
end %end if k==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 3 copies
if k==3
    %Declare the SDP variable related to the probability p
    variable p
    % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout) complex semidefinite
        variable F(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout) complex semidefinite
    else
        variable S(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout) semidefinite
        variable F(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout)) ,[2 3 4 5 6 7],[din d D d D d D dout]) == p*Jout(:,:,i);
    end
end  % End if k==3
%%%%%%%%%%%%%%%%%%%%%%%%
%The elements S and F should sum to a parallel/sequential/general superchannel
C=S+F;
if choose_protocol==1
    is_parallel_protocol_cvx(C,din,dout,d,D,k);
end
if choose_protocol==2
    is_sequential_protocol_cvx(C,din,dout,d,D,k);
end
if choose_protocol==3
    is_general_protocol_cvx(C,din,dout,d,D,k);
end

maximise p
%finish SDP
    cvx_end
    maxp=p;
end
