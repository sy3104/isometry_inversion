clear
%Code that obtains the maximal success probability of transforming 'k' uses
%of a set of 'n' random isometry operations from 'd'-dimensional space to
%'D'-dimensional space into its complex conjugate map.
%Please set the parameters 'k','d','D','n','choose protocol', and 'isComplex'
%The interpreter cvx is used to conduct SDP.


%Set the dimension d and D of the isometry

%Set the number k of uses of the input-isometry

%Choose the protocol that will be used
%1 corresponds to parallel
%2 corresponds to sequential
%3 corresponds to general (potentially with indefinite causal order)
     
%Set the number of random isometries that will be used.
%If n==0, we set the dimension of the space spanned by k+1 isometry
%operations from d-dimensional space to D-dimensional space.

%Decide if supermaps can use complex numbers of if we restrict only to real
%numbers.
%Set isComplex=0 to consider supermaps only with complex numbers
%From numerical evidence, we see that all optimal protocols used
%supermaps with only real numbers. When the SDP variable associated to the
%superinstruments is set to be real, the SDPs are faster and require less
%RAM memory
%Set isComplex==1 to consider general complex valued supermaps

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Start of Adjustable Parameters  %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        d=2                %Dimension 'd'
        D=3                %Dimension 'D'
        k=1                %Number of uses 'k'
        n=0;               %Number of linear independent V^k, n==0 implies full basis
        choose_protocol=1; %1 for parallel, 2 for sequential, 3 for general
        isComplex=0;       %Set isComplex=0 to restrict to supermaps with real coefficients

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%   End of Adjustable Parameters   %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if n==0
% Set 'n' as the dimension of the space spanned by k isometry operations
% from d-dimensional space to D-dimensional space.
    n=find_dimension_isometry(d,D,k)
else
    n=n
        if find_dimension_isometry(d,D,k)>n
            'The dimension of the space spanned by V^k is smaller than the number of isometries n'
            'The value pmax that will be displayed is just an upper bound'
        end
end

%Create a maximally entangled state and multiply by the dimension 'd'
    max_ent_times_d=IsotropicState(d,1)*d;
%Declare variables that contain the Choi input and Choi output operations
    Jin=zeros(d*D,d*D,n);  
    Jout=zeros(d*D,d*D,n);
    list_V=zeros(D,d,n);
%Asign values for Choi inputs and Choi outputs
for i=1:n 
    %Creates a random isometry from d-dimensional space to D-dimensional
    %space by truncating a random D-dimensional unitary sampled by the Haar measure
        V=RandomUnitary(D);
        for a=d+1:D
            V(:,d+1)=[];
        end
    %Creates the Choi operator for the unitary U and stock at a tensor for
    %the input
        Jin(:,:,i)=kron(eye(d),V) * max_ent_times_d * kron(eye(d),V');
    %Define the Target isometry as the complex conjugate of V
        Vt=conj(V);
    %Creates the Choi operator for the unitary U' and stock at a tensor for
    %the output
        Jout(:,:,i)=kron(eye(d),Vt) * max_ent_times_d * kron(eye(d),Vt');
end
tic;
%Find the maximal success probability of transforming 'k' isometries from
%d-dimensional space to D-dimensional space into its complex conjugate map
[maxp S F]=maxp_cvx(d, D, d, D, Jin, Jout, k, choose_protocol, isComplex);

maximal_success_probability=maxp
total_time_in_minutes=toc/60