function [T varargout] = myLanczos(mat,varargin)
%--------------------------------------------------------------------------
% Syntax:       T = myLanczos(mat);
%               T = myLanczos(mat,k);
%               T = myLanczos(mat,v);
%               T = myLanczos(mat,k,v);
%               [T Q] = myLanczos(mat);
%               [T Q] = myLanczos(mat,k);
%               [T Q] = myLanczos(mat,v);
%               [T Q] = myLanczos(mat,k,v);
%               
% Inputs:       mat is a symmetric or Hermitian N x N matrix
%
%               k is the number of Lanczos iterations to perform. The
%               default value is N
%
%               v is an N x 1 vector (internally forced to have unit norm)
%               specifying a particular starting vector to use. The default
%               vector is chosen uniformly at random from the unit sphere
%
% Outputs:      T is a k x k symmetric tridiagonal matrix. When k < N, T is
%               approximately similar to mat up to rank k. When k = N, T is
%               similar (in the linear algebra sense) to mat. That is, when
%               k = N, T has identical eigenvalues to mat
%
%               Q is the N x k similarity transformation such that
%               mat = Q * T * Q'. Note that equality holds only when k = N.
%               For k < N, mat ~ Q * T * Q'
%
% Description:  This function uses the Lanczos algorithm with full
%               reorthogonalization to compute k x k symmetric tridiagonal
%               matrix T that approximates mat up to rank k with respect to
%               transformation Q. That is, mat = Q * T * Q'. Note that
%               equality holds only for k = N. For k < N, mat ~ Q * T * Q'
%
%               NOTE: In particular, when k = N, T has identical
%               eigenvalues to mat
%
% Author:       Brian Moore
%               brimoor@umich.edu
%
% Date:         July 5, 2013
%--------------------------------------------------------------------------

% Knobs
symTol = 1e-8;

% Check input matrix size
[m,n] = size(mat);
if (m ~= n)
    error('Input matrix must be square');
end

% Make sure input matrix is symmetric
if max(max(abs(mat - mat'))) > symTol
    error('Input matrix is not symmetric to working precision');
end

% Parse user inputs
if (nargin == 2)
    if (length(varargin{1}) == n)
        k = n;
        v = varargin{1};
    else
        k = varargin{1};
        v = randn(n,1);
    end
elseif (nargin == 3)
        k = varargin{1};
        v = varargin{2};
else
    k = n;
    v = randn(n,1);
end

% Initialize variables
Q = nan(n,k);
q = v / norm(v);
Q(:,1) = q;
d = nan(k,1);
od = nan(k-1,1);

% Perform Lanczos iterations
for i = 1:k
    z = mat * q;
    d(i) = q' * z;
    
    %----------------------------------------------
    % Uncomment only one of the following 3 options
    %----------------------------------------------
    % Option 1: Full re-orthogonalization
    z = z - Q(:,1:i) * (Q(:,1:i)' * z);
    
    % Option 2: Full re-orthogonalization (x2)
%     z = z - Q(:,1:i) * (Q(:,1:i)' * z);
%     z = z - Q(:,1:i) * (Q(:,1:i)' * z);
    
    % Option 3: No re-orthogonalization
    %z = z - d(i) * q;
    %if (i > 1)
    %    z = z - od(i-1) * Q(:,i-1);
    %end
    %----------------------------------------------
    
    if (i ~= k)
        od(i) = norm(z);
        q = z / od(i);
        Q(:,i + 1) = q;
    end
end

% Construct T
T = diag(d) + diag(od,-1) + diag(od,1);

% Return user-requested information
if (nargout == 2)
    varargout{1} = Q;
end
