function nLI=isometry_linear_independent(list_V)
%This funcion calculate the number of linearly independent linear maps
%(in particular, isometry operations) in a list of operators 'list_V'.
%The variable 'list_V' should be in the form list_U(:,:,i), where i ranges
%from 1 to the number of linear operators

    N=size(list_V,3); %Evaluate the number of elements in the list
    D=size(list_V,1); %Evaluate the output dimension of the operators in the list
    d=size(list_V,2); %Evaluate the input dimension of the operators in the list
    A=nan(d^2*D^2,N); %Initialise a matrix that will stock the operators as vectors
    phiP=IsotropicState(d,1);  %Initialise a matrix that will stock the operators as vectors
    for i=1:N
        Vnow=list_V(:,:,i);
        %Create the Choi operator associated to the map induced by the
        %operator Vnow
        aux=kron(Vnow,eye(d))*phiP*kron(Vnow,eye(d))'; 
        %Stock the this choi operator in a matrix
        A(:,i)=aux(:);
    end
    %The number of linear independent maps is the rank of the matrix A
    nLI=rank(A);
end
