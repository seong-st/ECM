function params = parameters()
    % Parameters for the 2013 Nissan Leaf battery cell

    % Surface area of the battery cell [m²]
    a_surf = 0.067569;

    % Heat capacity of the battery cell [J/(kg K)]
    cp_cell = 1600;

    % Coulombic efficiency of the battery cell for charge and discharge [-]
    eta_chg = 0.98;
    eta_dis = 1.00;

    % Convective heat transfer coefficient [W/(m² K)]
    h_conv = 9.5;

    % Mass of a single battery cell [kg]
    m_cell = 0.799;

    % Total capacity of the battery cell [Ah]
    q_cell = 26;

    % Ambient temperature [K]
    tinf = 298.15;

    % Save parameters in a structure
    params.eta_chg = eta_chg;
    params.eta_dis = eta_dis;
    params.q_cell = q_cell;
    params.a_surf = a_surf;
    params.cp_cell = cp_cell;
    params.h_conv = h_conv;
    params.m_cell = m_cell;
    params.tinf = tinf;
end
