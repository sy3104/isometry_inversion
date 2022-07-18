function [maxp S F]=maxp_isometry_inversion_yalmip(d, D, Jin,k,choose_protocol,isComplex)
% Code that evaluates the maximal probability of transforming the set of
% k isometries encoded in Choi operators Jin and obtain the inverse maps.
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
dout=d;            %dimension of the output space of the output operation
max_ent_times_d=IsotropicState(d,1)*d;
    yalmip('clear');
    ops = sdpsettings('solver','scs','scs.eps',1e-5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 1 copy
if k==1
    %Declare the SDP variable related to the probability p
    sdpvar p;
    % Declare the SDP variables related to S and F
    if isComplex==1
        S=sdpvar(din*d*D*dout,din*d*D*dout,'hermitian','complex');
        F=sdpvar(din*d*D*dout,din*d*D*dout,'hermitian','complex');
    else
        S=sdpvar(din*d*D*dout,din*d*D*dout,'symmetric','real');
        F=sdpvar(din*d*D*dout,din*d*D*dout,'symmetric','real');
    end

    Constraints=[S>=0, F>=0];
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout is the Choi operator of the output operation.
        %Impose the condition that the link product of Jout and Jin is
        %proportional to the maximally entangled state.
        Jout = PartialTrace(S * Tensor(eye(din) , transpose(Jin(:,:,i)) , eye(dout)) ,[2 3],[din d D dout]);
        Constraints = [Constraints, PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(d)), [2], [d D d])*Tensor(eye(d), Jout), [2], [d D d])== p*max_ent_times_d];
    end
end %end if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 2 copies
if k==2
    %Declare the SDP variable related to the probability p
    sdpvar p;
    % Declare the SDP variables related to S and F
    if isComplex==1
        S=sdpvar(din*d*D*d*D*dout,din*d*D*d*D*dout,'hermitian','complex');
        F=sdpvar(din*d*D*d*D*dout,din*d*D*d*D*dout,'hermitian','complex');
    else
        S=sdpvar(din*d*D*d*D*dout,din*d*D*d*D*dout,'symmetric','real');
        F=sdpvar(din*d*D*d*D*dout,din*d*D*d*D*dout,'symmetric','real');
    end

    Constraints=[S>=0, F>=0];
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout is the Choi operator of the output operation.
        %Impose the condition that the link product of Jout and Jin is
        %proportional to the maximally entangled state.
        Jout = PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout)) ,[2 3 4 5],[din d D d D dout]);
        Constraints = [Constraints, PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(d)), [2], [d D d])*Tensor(eye(d), Jout), [2], [d D d])== p*max_ent_times_d];
    end
end %end if k==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 3 copies
if k==3
    %Declare the SDP variable related to the probability p
    sdpvar p;
    % Declare the SDP variables related to S and F
    if isComplex==1
        S=sdpvar(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout,'hermitian','complex');
        F=sdpvar(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout,'hermitian','complex');
    else
        S=sdpvar(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout,'symmetric','real');
        F=sdpvar(din*d*D*d*D*d*D*dout,din*d*D*d*D*d*D*dout,'symmetric','real');
    end

    Constraints=[S>=0, F>=0];
% Impose that the relation between input and output must be satisfied
% This is made exploiting the Choi isomorphism
    for i=1:n_in
        %Jout is the Choi operator of the output operation.
        %Impose the condition that the link product of Jout and Jin is
        %proportional to the maximally entangled state.
        Jout = PartialTrace(S * Tensor(eye(din) , transpose(Tensor(Jin(:,:,i),k)) , eye(dout)) ,[2 3 4 5 6 7],[din d D d D d D dout]);
        Constraints = [Constraints, PartialTrace(PartialTranspose(Tensor(Jin(:,:,i), eye(d)), [2], [d D d])*Tensor(eye(d), Jout), [2], [d D d])== p*max_ent_times_d];
    end
end  % End if k==3
%%%%%%%%%%%%%%%%%%%%%%%%
%The elements S and F should sum to a parallel/sequential/general superchannel
C=S+F;
if choose_protocol==1
    Constraints = is_parallel_protocol_yalmip(C,din,dout,d,D,k,Constraints);
end
if choose_protocol==2
    Constraints = is_sequential_protocol_yalmip(C,din,dout,d,D,k,Constraints);
end
if choose_protocol==3
    Constraints = is_general_protocol_yalmip(C,din,dout,d,D,k,Constraints);
end

    sol=optimize(Constraints, -p, ops);
    maxp=p;
end
