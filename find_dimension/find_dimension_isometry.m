function dim=find_dimension_isometry(d,D,k)
%Function that evaluate the dimention of the linear space spanned by 'k'
%copies of isometry operations from 'd'-dimensional space to
%'D'-dimensional space
%Start generating 'n' random isometry operations
%initialise the 'n' parameter as 10^3
%This value is, in principle, arbitrary. We have chosen 10^3 because we have verified numerically that is relatively fast for the computer to analyse 10^3 unitaries of the form V^k 
    n=10^3;
 %Initialise a variable to stock a set of 'k' isometry operators
    list_V_k=zeros(D^k,d^k,n);
    for i=1:n 
        %Creates a random isometry from d-dimensional space to D-dimensional
        %space by truncating a random D-dimensional unitary sampled by the Haar measure
        V=RandomUnitary(D);
        for a=d+1:D
            V(:,d+1)=[];
        end
        %Add 'k' copies of the isometry operator into a list
        list_V_k(:,:,i)=Tensor(V,k);
    end
    %Evaluate the number of linear independent isometry operations
    nLI=isometry_linear_independent(list_V_k);
    %This while only stops when the number 'n' is larger nLI
    while n<=nLI
        %Generate more random isometry operations, the value 100 is arbitrary
        n=n+100;
        %Initialise a variable to stock a set of 'k' isometry operators
        list_V_k=zeros(D^k,d^k,n);
        for i=1:n 
            V=RandomUnitary(D);
            for a=d+1:D
                V(:,d+1)=[];
            end
            list_V_k(:,:,i)=Tensor(V,k);
        end
        %Evaluate the number of linear independent isometry operations
        nLI=isometry_linear_independent(list_V_k);
    end
%Show total time of this funcion
    dim=nLI;
end
