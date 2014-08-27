% @Author: Athul Vijayan
% @Date:   2014-08-23 14:52:29
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2014-08-23 14:52:37

clear all
numGrids = 5;                       % number of different grids over which solution is sought
phiVector = cell(numGrids,1);           % solution vectors
phi= cell(numGrids,1);           % solution vectors
xArray = cell(numGrids, 1);
yArray = cell(numGrids, 1);
for i = 2:numGrids
    % =========== First we produce the matrix A for the bottom rectangular part. ==============
    xNumFirstGridElems = 3;
    yNumFirstGridElems = 2; 
    xNumElem = xNumFirstGridElems*2^(i-1);      % x number of elements
    yNumElem = yNumFirstGridElems*2^(i-1);      % y number of elements
    

    deltaX = 1/yNumElem;                        % delta x = delta y

    % make B matrix., so that we can just pile this to form A matrix
    BdiagVec = -4*ones(xNumElem-1,1);           % exclude two boundaries
    BdiagVec(1) = -3+deltaX;                    %add derivative boundary condition, at all i = 1 for bottom part
    BdiagVec(xNumElem-1) = -3+deltaX;           %add derivative boundary condition, at all i = n-1 for bottom part
    if (xNumElem-1 > 1)
        B_bottom = gallery('tridiag',ones(xNumElem -2, 1),BdiagVec,ones(xNumElem -2, 1)); % matrix of coefficients
    else
        B_bottom = BdiagVec                     % case for $$\Delta x =1$$
    end

    A = B_bottom;                               %Initialize A matrix
    for k=1:yNumElem-2
        A = blkdiag(A, B_bottom);               %populate A matrix with B matrix as diagonal
    end
    % Now loop through the A matrix to add I matrix wherever necessary
    for k= 1: (xNumElem - 1)*(yNumElem-1)
        if (k + xNumElem -1 <= size(A, 2))      % leave the boundary, i=n and j=m
            A(k, k + xNumElem -1) = 1;
        end
        if (k > (xNumElem - 1))                 % leave the boundary, i=0 and j=0
            A(k, k -( xNumElem -1)) = 1;
        end
    end
    % Now the vector A has the coefficients of $$\phi$$s inside the bottom rectangle.
        % ====================================================================

    % Now we add to this matrix, the coefficients for the equations of top rectangular part
    % having different dimensions. So this will result in change of dimesion of B matrix and I matrix

    % =========== then we produce the matrix A for the bottom rectangular part. ==============
    xNumFirstGridElems = 2;
    yNumFirstGridElems = 2;
    yNumElemOld = yNumElem;
    xNumElemOld = xNumElem;
    xNumElem = xNumFirstGridElems*2^(i-1);      % number of x elements
    yNumElem = yNumFirstGridElems*2^(i-1);      % number of y elements

    BdiagVec = -4*ones(xNumElem-1,1);           % exclude two boundaries
    BdiagVec(1) = -3+deltaX;                    %add derivative boundary condition, at all i = 1 for bottom part
    if (xNumElem-1 > 1)
        B_top = gallery('tridiag',ones(xNumElem -2, 1),BdiagVec,ones(xNumElem -2, 1)); % matrix of coefficients
    else
        B_top = BdiagVec;
    end

    for k=1: yNumElem
        A = blkdiag(A, B_top);                  %populate A matrix with B matrix as diagonal
    end
    % Now loop through the A matrix to add I matrix wherever necessary
    for k= (yNumElemOld-1)*(xNumElemOld-1) + 1: (yNumElemOld-1)*(xNumElemOld-1) + (xNumElem - 1)*(yNumElem)
        if (k + xNumElem -1 <= size(A, 2))   
            A(k, k + xNumElem -1) = 1;
        end
        if (k > (yNumElemOld-1)*(xNumElemOld-1) + (xNumElem - 1))
            A(k, k -( xNumElem -1)) = 1;
        end
    end
        % ====================================================================

    % we still have not established the 'connection' between these matrices. 
    % i.e. in the sheet, we haven't considered the fact that there is no boundary 
    % between this top and bottom plate upto a point. right now the matrix assumes, zero boundary condition.
    % to fix this we can link by following code
    for k= (xNumElemOld - 1)*(yNumElemOld-2): (xNumElemOld - 1)*(yNumElemOld-2) + (xNumElem - 1)
        if k > 0
            A(k, k + xNumElemOld -1) = 1;
        end
    end

    for k= (xNumElemOld - 1)*(yNumElemOld-1) + 1: (xNumElemOld - 1)*(yNumElemOld-1) + xNumElem -1
        A(k, k - (xNumElemOld -1)) = 1;
    end
    % /*******************
    % Now the matrix is complete for the given plate with boundary condition as zero along boundaries.
                                                                    % ***************************/
    % We now add boundary conditions one by one

    % Now add derivative boundary condition at j=1
    y = 1;
    for index = y:y + (xNumElemOld-2)
        id = find(A(index, :) == -4);
        A(index, id) = -3 + deltaX;
    end
    % Now add derivative boundary condition at j=m-1
    y = (xNumElemOld - 1)*(yNumElemOld-1) + (xNumElem -1)*(yNumElem-1) + 1;
    for index = y :y + (xNumElem-2)
        id = find(A(index, :) == -4);
        A(index, id) = -3 + deltaX;
    end

    % 3 corner points, where two boundary terms will come
    A(1, 1) = A(1, 1) + (1+deltaX);
    A(xNumElemOld -1, xNumElemOld -1) = A(xNumElemOld -1, xNumElemOld -1) + (1+deltaX);
    A(size(A,1)-(xNumElem-2), size(A,1)-(xNumElem-2)) = A(size(A,1)-(xNumElem-2), size(A,1)-(xNumElem-2)) + (1+deltaX);
    clear('y', 'id', 'index');


    % initialize F
    F = zeros(size(A,1),1);
    % assign -100 to points corresponding to ones near inner boundary
    for k=(yNumElemOld-2)*(xNumElemOld-1)+ xNumElem:(yNumElemOld-2)*(xNumElemOld-1)+xNumElemOld-1
        F(k) = -100;
    end

    for y=1:yNumElem
        id = (yNumElemOld-1)*(xNumElemOld-1)+ y*(xNumElem-1);
        F(id) = -100;
    end

    % Now we have all the required matrices. we just need to invert A matrix to get temperature distribution
    phiVector{i} = inv(A)*F;
    phi{i} = NaN(yNumElemOld+yNumElem + 1, xNumElemOld +1);

    % Now we need to convert that row matrix into a 2D matrix for contour potting.
    % It is done by
    for y=1:size(phi{i}, 1)
        for x=1:size(phi{i}, 2)
            % update interior points
            if ((y < yNumElemOld + 1) && (y >1))
                phi{i}(y, 2: xNumElemOld) = phiVector{i}((y-1)*(xNumElemOld-1) +1: y*(xNumElemOld-1) );
            end

            % update internal boundary
            if (y > yNumElemOld)
                phi{i}(y, 2: xNumElem) = phiVector{i}((y-1)*(xNumElem-1) +1: y*(xNumElem-1) );
                if (x == xNumElem +1)
                    phi{i}(y, x) = -100;
                end
            end
            % update internal boundary
            if (y == yNumElemOld)
                phi{i}(y, xNumElem+1: end) = -100;
            end
        end
    end

    % update outer boundaries
    for y=1:size(phi{i}, 1)
        phi{i}(y, 1) = (1+deltaX)* phi{i}(y, 2);
        if y<yNumElemOld
            phi{i}(y, xNumElemOld + 1) = (1+deltaX)* phi{i}(y, xNumElemOld);
        end
    end

    % update outer boundaries
    for x=1:size(phi{i}, 2)
        phi{i}(1, x) = (1+deltaX)* phi{i}(2, x);
    end

    % Since we divided sheet to 4, remake that full matrix
    phi{i} = horzcat(phi{i}, flip(phi{i}(:, 2:end), 2));
    phi{i} = vertcat(phi{i}, flip(phi{i}(2:end, :)));

    % grid points
    xArray{i} = 0:2*deltaX:6;
    yArray{i} = 0:2*deltaX:8;

    % plot the contour
    figure;
    [X,Y] = meshgrid(xArray{i},yArray{i});
    contour(X, Y, phi{i}, 30); 
    title(['The contour plot of temperature for iteration ', num2str(i)]);
end

