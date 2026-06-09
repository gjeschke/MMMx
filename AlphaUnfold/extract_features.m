function features = extract_features(matrix, vector)

% fitrensemble, fitrsvm, fitrgp

    % matrix: n×n, vector: 1×n
    % returns a row vector of fixed length (e.g., 15 features)
    
    % matrix = (matrix + matrix')/2;
    % Matrix features
    matrix = (matrix + matrix')/2;
    mat = matrix(:);                 % all entries
    mat_mean = mean(mat);
    mat_std  = std(mat);
    mat_skew = skewness(mat);
    mat_kurt = kurtosis(mat);
    mat_norm = norm(matrix, 'fro');
    % Add eigenvalues (first 3 sorted descending, pad with zeros if fewer)
    e = sort(eig(matrix), 'descend');
    e_pad = [e(:)', zeros(1, max(0,3-length(e)))];
    e1 = e_pad(1); e2 = e_pad(2); e3 = e_pad(3);
    
    % Vector features
    vec = vector(:);
    v_mean = mean(vec);
    v_std  = std(vec);
    v_skew = skewness(vec);
    v_kurt = kurtosis(vec);
    v_min  = min(vec);
    v_max  = max(vec);
    
    % Combine
    features = [mat_mean, mat_std, mat_skew, mat_kurt, mat_norm,  ...
                e1, e2, e3, v_mean, v_std, v_skew, v_kurt, v_min, v_max];
end