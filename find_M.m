function answer = find_M (M, m_E, sigma_E)
    answer = f((M - m_E) / sigma_E) + f((M + m_E) / sigma_E) - 1.8;
end