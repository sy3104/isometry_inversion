function is_sequential_protocol_cvx(C,din,dout,d,D,k)
% Code that imposes SDP constrains on a variable C to enforce that C is a valid
% sequential superchannel
% 'din' corresponds to the dimension of the input space of the output operation
% 'dout' corresponds to the dimension of the output space of the output operation
% 'd' corresponds to the dimension of the input space of the input operation
% 'D' corresponds to the dimension of the output space of the input operation
% 'k' corresponds to the number of uses of the input-operation

    cvx_begin SDP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 1 copy
if k==1
    PartialTrace(C,4,[din d D dout])==kron(PartialTrace(C,[3 4],[din d D dout]),eye(D)/D);
    PartialTrace(C,[2 3 4],[din d D dout])==eye(din)*D;
end %end if k==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 2 copies
if k==2
     PartialTrace(C,6,[din d D d D dout]) == kron(PartialTrace(C,[5 6],[din d D d D dout]),eye(D)/D);
     PartialTrace(C,[4 5 6],[din d D d D dout]) == kron(PartialTrace(C,[3 4 5 6],[din d D d D dout]),eye(D)/D);
     PartialTrace(C,[2 3 4 5 6],[din d D d D dout])==eye(din)*D^2;
end %end if k==2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% case of 3 copies
if k==3  
    PartialTrace(C,8, [din d D d D d D dout])==kron(PartialTrace(C,[7 8], [din d D d D d D dout]),eye(D)/D);
    PartialTrace(C,[6 7 8], [din d D d D d D dout])==kron(PartialTrace(C,[5 6 7 8], [din d D d D d D dout]),eye(D)/D);
    PartialTrace(C,[4 5 6 7 8], [din d D d D d D dout])==kron(PartialTrace(C,[3 4 5 6 7 8], [din d D d D d D dout]),eye(D)/D);
    PartialTrace(C,[2 3 4 5 6 7 8], [din d D d D d D dout])==eye(din)*D^3;     
end %end if k==3
    cvx_end
end
