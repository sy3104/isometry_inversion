function [maxp S F]=maxp_isometry_inversion_sod_cvx(d, D, Jin,k,choose_protocol,isComplex)
% Code that evaluates the maximal probability of transforming the set of
% k isometries encoded in Choi operators Jin and obtain the inverse maps
% in "success-or-draw" manner.
% Jin is a tensor with 'n' Choi operators that should be written as
% Jin(:,:,i), where the variable 'i' ranges from 1 to 'n'
% 'k' corresponds to the number of uses of the input-operation
% Choose the protocol that will be used
% 1 corresponds to parallel
% 2 corresponds to sequential
% 3 corresponds to general (potentially with indefinite causal order)
% if isComplex is '0', the code only consider supermaps with real numbers
% is isComplex is '1', the code considers general supermaps with complex
% numbers

n_in=size(Jin,3);  %Count the number of choi operators
din=D;             %dimension of the input space of the output operation
dout_S=d;          %dimension of the output space of the successful output operation
dout_F=D;          %dimension of the output space of the failure output operation
max_ent_times_d=IsotropicState(d,1)*d;
    cvx_begin SDP quiet
    %cvx_solver sedumi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 1 copy
if k==1
    %Declare the SDP variable related to the probability p
    variable p
    % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*dout_S,din*d*D*dout_S) complex semidefinite
        variable F(din*d*D*dout_F,din*d*D*dout_F) complex semidefinite
    else
        variable S(din*d*D*dout_S,din*d*D*dout_S) semidefinite
        variable F(din*d*D*dout_F,din*d*D*dout_F) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout_S is the Choi operator of the successful output operation.
        %Jout_F is the Choi operator of the failure output operation.
        %Impose the condition that
        %the link product of Jout_S and Jin is proportional to the maximally entangled state.
        %the link product of Jout_F and Jin is proportional to Jin.
        Jout_S = PartialTrace(S * Tensor(eye(din) , transpose(Jin(:,:,i)) , eye(dout_S)) ,[2 3],[din d D dout_S]);
        Jout_F = PartialTrace(F * Tensor(eye(din) , transpose(Jin(:,:,i)) , eye(dout_F)) ,[2 3],[din d D dout_F]);
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_S)), [2], [d D dout_S])*Tensor(eye(d), Jout_S), [2], [d D dout_S])== p*max_ent_times_d;
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_F)), [2], [d D dout_F])*Tensor(eye(d), Jout_F), [2], [d D dout_F])== (1-p)*Jin(:,:,i);
    end
end %end if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 2 copies
if k==2
    %Declare the SDP variable related to the probability p
    variable p
   % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*d*D*dout_S,din*d*D*d*D*dout_S) complex semidefinite
        variable F(din*d*D*d*D*dout_F,din*d*D*d*D*dout_F) complex semidefinite
    else
        variable S(din*d*D*d*D*dout_S,din*d*D*d*D*dout_S) semidefinite
        variable F(din*d*D*d*D*dout_F,din*d*D*d*D*dout_F) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout_S is the Choi operator of the successful output operation.
        %Jout_F is the Choi operator of the failure output operation.
        %Impose the condition that
        %the link product of Jout_S and Jin is proportional to the maximally entangled state.
        %the link product of Jout_F and Jin is proportional to Jin.
        Jout_S = PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout_S)) ,[2 3 4 5],[din d D d D dout_S]);
        Jout_F = PartialTrace(F * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout_F)) ,[2 3 4 5],[din d D d D dout_F]);
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_S)), [2], [d D dout_S])*Tensor(eye(d), Jout_S), [2], [d D dout_S])== p*max_ent_times_d;
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_F)), [2], [d D dout_F])*Tensor(eye(d), Jout_F), [2], [d D dout_F])== (1-p)*Jin(:,:,i);
    end
end %end if k==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 3 copies
if k==3
    %Declare the SDP variable related to the probability p
    variable p
    % Declare the SDP variables related to S and F
    if isComplex==1
        variable S(din*d*D*d*D*d*D*dout_S,din*d*D*d*D*d*D*dout_S) complex semidefinite
        variable F(din*d*D*d*D*d*D*dout_F,din*d*D*d*D*d*D*dout_F) complex semidefinite
    else
        variable S(din*d*D*d*D*d*D*dout_S,din*d*D*d*D*d*D*dout_S) semidefinite
        variable F(din*d*D*d*D*d*D*dout_F,din*d*D*d*D*d*D*dout_F) semidefinite
    end
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout_S is the Choi operator of the successful output operation.
        %Jout_F is the Choi operator of the failure output operation.
        %Impose the condition that
        %the link product of Jout_S and Jin is proportional to the maximally entangled state.
        %the link product of Jout_F and Jin is proportional to Jin.
        Jout_S = PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout_S)) ,[2 3 4 5 6 7],[din d D d D d D dout_S]);
        Jout_F = PartialTrace(F * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout_F)) ,[2 3 4 5 6 7],[din d D d D d D dout_F]);
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_S)), [2], [d D dout_S])*Tensor(eye(d), Jout_S), [2], [d D dout_S])== p*max_ent_times_d;
        PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(dout_F)), [2], [d D dout_F])*Tensor(eye(d), Jout_F), [2], [d D dout_F])== (1-p)*Jin(:,:,i);
    end
end  % End if k==3
%%%%%%%%%%%%%%%%%%%%%%%%
%The elements S and F should sum to a parallel/sequential/general superchannel
C=PartialTrace(S, [4], [din d^k D^k dout_S]) +PartialTrace(F, [4], [din d^k D^k dout_F]);
if choose_protocol==1
    is_parallel_protocol_cvx(C,din,1,d,D,k);
end
if choose_protocol==2
    is_sequential_protocol_cvx(C,din,1,d,D,k);
end
if choose_protocol==3
    is_general_protocol_cvx(C,din,1,d,D,k);
end

maximise p
%finish SDP
    cvx_end
    maxp=p;
end
