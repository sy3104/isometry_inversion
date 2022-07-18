function is_parallel_protocol_cvx(C,din,dout,d,D,k)
% Code that imposes SDP constrains on a variable C to enforce that C is a valid
% sequential superchannel
% 'din' corresponds to the dimension of the input space of the output operation
% 'dout' corresponds to the dimension of the output space of the output operation
% 'd' corresponds to the dimension of the input space of the input operation
% 'D' corresponds to the dimension of the output space of the input operation
% 'k' corresponds to the number of uses of the input-operation

d1=din;
d2=d^k;
d3=D^k;
d4=dout;

    cvx_begin SDP
% Regroup the input spaces and output spaces of the input operation.
if k==2
    C=PermuteSystems(C, [1 2 4 3 5 6], [din d D d D dout]);
end
if k==3  
    C=PermuteSystems(C, [1 2 4 6 3 5 7 8], [din d D d D d D dout]);
end

%Impose the condition that C is a superchannel that makes a single use of the input-operation
    PartialTrace(C,4,[d1 d2 d3 d4])==kron( PartialTrace(C,[3 4],[d1 d2 d3 d4]),eye(d3)/d3);
    PartialTrace(C,[2 3 4],[d1 d2 d3 d4])==eye(d1)*d3;
    cvx_end
end
