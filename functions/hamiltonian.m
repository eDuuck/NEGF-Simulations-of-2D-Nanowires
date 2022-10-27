function [H] = hamiltonian(sample, B)
    %HAMILTONIAN Returns the hamiltonian of a given sample.
    %   Returns the hamiltonian matrix for a sample, a magnetic field described
    %   by the matrix B can be applied to the sample. If the sample has a
    %   size of [M N] then B needs to be of size [M N 2] for the different
    %   coordinates of the vector field matching every point of the sample's
    %   units.

    hasMagField = exist('B','var') && ~isequal(B,0*B);

    q = 1.6022e-19;
    hbar = 1.0546e-34;
    qha = 1i*sample.a*q/hbar;
    wid = sample.width;
    len = sample.length;

    if lower(sample.arch) == "rectangular"

        t = sample.t;
        units = sample.getUnits();
        H = zeros(sample.M);
        if ~hasMagField
            %         H_elements_max = (wid+length(t)*(2*wid-length(t)-1))*len + ... %Main column elements.
            %             wid*length(t)*(2*len-length(t)-1);
            %
            %         H = spalloc(sample.M, sample.M,H_elements_max);
            H = zeros(sample.M);
            for j = 1:len
                H(((wid*(j-1))+1):(wid*j),((wid*(j-1))+1):(wid*j)) = ...
                    column(wid,units(:,j),t); %Fills column blocks.
            end

            for j = 1:length(t)
                H = H + sparse(t(j)*(diag(ones((len-j)*wid,1),wid*j)+diag(ones((len-j)*wid,1),-wid*j)));
            end

        else
            %This can be implemented a lot more elegantly with the use of
            %the built in matlab function kron(X,Y)
            if length(B) == 1           %Uniform magnetic field.
                B = ones(wid,len) * B;
            end
            A = anti_curl(B,sample.a);
            for x = 1:len
                alpha = diag(units(:,x));
                for j = 1:length(t)
                    alpha = alpha + diag(ones(wid-1,1)*t(j).*exp(qha * A(1:end-1,x,2)),1) ...
                        + diag(ones(wid-1,1)*t(j).*exp(-qha * A(2:end,x,2)),-1);
                end
                aI = [(x-1)*wid+1, x*wid]; %AlphaIndex.
                H(aI(1):aI(2),aI(1):aI(2)) = alpha;
                for j = 1:length(t)
                    beta = diag(ones(wid,1)*t(j).*exp(qha * A(:,x,1)));
                    if x + j - 1 < len  %%This checks that we're not out of matrix
                        bI = [aI(1), aI(2), aI(1)+j*wid,aI(2)+j*wid];
                        H(bI(1):bI(2),bI(3):bI(4)) = beta;
                    end
                    if x - (j - 1) > 1
                        bI = [aI(1), aI(2), aI(1)-j*wid,aI(2)-j*wid];
                        H(bI(1):bI(2),bI(3):bI(4)) = beta';
                    end
                end
            end
        end
    end

end

function A = column(y,eps,t)
    A = diag(eps);
    for j = 1:length(t)
        A = A + t(j)*(diag(ones(y-j,1),j) + diag(ones(y-j,1),-j));
    end
end