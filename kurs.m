function [k_mas, k_0_mas, k_1_mas, m_E_mas, sigma_E_mas, M_mas] = kurs
    T1 = 0.5; T2 = 3; g = 1; a = 5; d = 10; % параметры системы
    
    k_mas       = zeros(100);
    k_0_mas     = zeros(100);
    k_1_mas     = zeros(100);
    m_E_mas     = zeros(100);
    M_mas       = zeros(100);
    sigma_E_mas = zeros(100);
    
    m_E = 1; sigma_E = 1;       % начальные значения для вычислений
    i = 1; k = 0.5;
    
    while k <= 20
        m_E=1; sigma_E=1;
        k_0 = 2; k_1 = 1; eps = 1;
        while eps > 0.0001  % находим значения параметров для данного k
            m_E_tmp = g / (1 + k_0 * k); %m_E
            sigma_E_tmp = sqrt(k ^ 2 * k_1 ^ 2 * 2 * a * d * (T1 * T2 * ...
                (a * T1 * T2 + T1 + T2)) /  (2 * T1 * T2 * a * ...
                (1 + k_1 * k) * ((a * T1 * T2 + T1 + T2) * (a * T1 + ...
                a * T2 + 1 + k_1 * k) - T1 * T2 * (a + a * k_1 * k))));
            k_0_tmp = 2 / m_E * f(m_E / sigma_E);
            k_1_0 = 1 / sigma_E * sqrt((1 - 4 * f(m_E / sigma_E) ^ 2));
            k_1_1 = 2 / (sigma_E * sqrt(2 * pi)) * exp(-(m_E ^ 2 / ...
                ( 2 * sigma_E ^ 2)));
            k_1_tmp = (k_1_0 + k_1_1) / 2;
            eps = abs(m_E - m_E_tmp) + abs(sigma_E - sigma_E_tmp) + ...
                abs(k_0 - k_0_tmp) + abs(k_1 - k_1_tmp);
            m_E = m_E_tmp;
            sigma_E = sigma_E_tmp;
            k_0 = k_0_tmp;
            k_1 = k_1_tmp;
        end
        k_mas(i)    = k;
        k_0_mas(i)  = k_0;
        k_1_mas(i)  = k_1;
        
        m_E_mas(i)  = m_E;
        sigma_E_mas(i) = sigma_E;
        
        M = fzero(@(M) find_M(M, m_E, sigma_E), 1.0)
        M_mas(i) = M;
        
        i = i + 1;
        k = k + 0.1;
    end
end

