function [val, vec] = batcheig(mat)
% Compute the eigendecomposition of a batch of symmetric 2x2 or 3x3
% matrices (fast vectorized implementation).
% 
% mat - (N, 6)     ordered as [xx, yx, yy, zx, zy, zz] if D == 3, or
%       (N, 3)     ordered as [xx, yx, yy] if D == 2
% val - (N, D)     eigenvalues (sorted by increasing magnitude)
% vec - (N, D, D)  eigenvectors
%
% "Closed-form expressions of the eigen decomposition 
%  of 2x2 and 3x3 Hermitian matrices"
% Charles-Alban Deledalle, Loic Denis, Sonia Tabti, Florence Tupin
% https://hal.science/hal-01501221/document

% --------------
% Yael Balbastre
% 2024-01-23
% --------------

N = size(mat, 1);
DD = size(mat, 2);

if DD == 3

    % 2D
    D = 2;

    xx = mat(:,1);
    xy = mat(:,2);
    yy = mat(:,3);

    delta = sqrt(4 * (xy.*xy + (xx - yy).^2));
    val1 = (xx + yy - delta) * 0.5;
    val2 = (xx + yy + delta) * 0.5;
    val = cat(2, val1, val2);

    [~, idx] = sort(abs(val), 2);
    idx = sub2ind(size(val), repmat((1:N)', 1, D), idx);

    if nargout > 1
        vec = ones([N, D, D], class(mat));
        vec(:, 1, 1) = (val2 - yy) ./ xy;
        vec(:, 2, 1) = (val1 - yy) ./ xy;
        vec = vec ./ sqrt(sum(vec.*vec, 3));
        vec = reshape(vec, [], 2);
        vec = reshape(vec(idx, :), [], 2, 2);
    end
    val = reshape(val(idx), [], 3);


elseif DD == 6

    % 3D
    D = 3;

    xx = mat(:,1);
    xy = mat(:,2);
    yy = mat(:,3);
    xz = mat(:,4);
    yz = mat(:,5);
    zz = mat(:,6);

    x1 = xx.*xx + yy.*yy + zz.*zz ...
       - xx.*yy - xx.*zz - yy.*zz ...
       + (xy.*xy + xz.*xz + yz.*yz) * 3;

    x2 = -(2*xx - yy - zz).*(2*yy - xx - zz).*(2*zz - xx - yy) ...
       + 9 *(...
        (2*zz - xx - yy) .* (xy.*xy) + ...
        (2*yy - xx - zz) .* (xz.*xz) + ...
        (2*xx - yy - zz) .* (yz.*yz) ) ...
       - (xy.*xz.*yz) * 54;

    phi = atan2(sqrt(abs(4 * x1.^3 - x2.^2)), x2);

    val1 = (xx + yy + zz - 2 * sqrt(x1) .* cos(phi/3)) / 3;
    val2 = (xx + yy + zz + 2 * sqrt(x1) .* cos((phi-pi)/3)) / 3;
    val3 = (xx + yy + zz + 2 * sqrt(x1) .* cos((phi+pi)/3)) / 3;
    val = cat(2, val1, val2, val3);

    [~, idx] = sort(abs(val), 2);
    idx = sub2ind(size(val), repmat((1:N)', 1, D), idx);

    if nargout > 1
        m1 = (xy .*(zz - val1) - xz.* yz) ./ (xz.*(yy - val1) - xy.*yz);
        m2 = (xy .*(zz - val2) - xz.* yz) ./ (xz.*(yy - val2) - xy.*yz);
        m3 = (xy .*(zz - val3) - xz.* yz) ./ (xz.*(yy - val3) - xy.*yz);

        vec = ones([N, D, D], class(mat));
        vec(:, 1, 1) = (val1 - zz - yz.*m1) ./ xz;
        vec(:, 2, 1) = (val2 - zz - yz.*m2) ./ xz;
        vec(:, 3, 1) = (val3 - zz - yz.*m3) ./ xz;
        vec(:, 1, 2) = m1;
        vec(:, 2, 2) = m2;
        vec(:, 3, 2) = m3;
        vec = vec ./ sqrt(sum(vec.*vec, 3));
        vec = reshape(vec, [], 3);
        vec = reshape(vec(idx, :), [], 3, 3);
    end
    val = reshape(val(idx), [], 3);

else
    error("Matrix size should be 2 or 3")
end


end