% surface_r_20.m
%
% want to deal with NaN here, replace NaN with right value, failed
%
% zwy 20190126

% replace NaN with the value of above or below surface,i.e., 0.216/100 or
% 0, repectively.
% zwy 20190126
m_r = length(x2_test); m_c = length(x1_test);
m = reshape(m, m_r, m_c); m_grid = m;
nan_flags = isnan(m);
% How to determine which value a NaN position should be? 
for i = 1:m_r
    for j = 1:m_c
        if nan_flags(i,j) == 1
            for k1 = i:m_r
                if nan_flags(k1,j) == 0
                    % replace NaN with a unNaN value in the other row
                    m_grid(i,j) = m_grid(k1,j);
                end
            end
            if isnan(m_grid(i,j))
                for k2 = j:m_c
                    m_grid(i,j) = m_grid(i,k2);
                end
            end
            if isnan(m_grid(i,j))
                m_grid(i,j)=(rand>0.5)*0.216/100;
            end
        end
    end
end
