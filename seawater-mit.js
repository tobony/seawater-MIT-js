/**
 * SEAWATER THERMOPHYSICAL PROPERTIES LIBRARY
 * https://github.com/tobony/seawater-MIT-js
 */


/**
 * Boiling point elevation of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} BPE - Boiling point elevation [K]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_BPE(T, S) {
    if (T < 0 || T > 200) {
        throw new Error("Temperature is out of range for boiling point elevation function 0 < T < 200 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for boiling point elevation function 0 < S < 120 g/kg");
    }

    S = S / 1000;
    const a1 = -0.00045838530457;
    const a2 = 0.28230948284;
    const a3 = 17.945189194;
    const a4 = 0.00015361752708;
    const a5 = 0.052669058133;
    const a6 = 6.5604855793;

    const A = a1 * Math.pow(T, 2) + a2 * T + a3;
    const B = a4 * Math.pow(T, 2) + a5 * T + a6;

    return A * Math.pow(S, 2) + B * S;
}

/**
 * Chemical potential of salt in seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Chemical potential of salt [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_ChemPot_s(T, S, P) {
    if (T < 10 || T > 80) {
        throw new Error("Temperature is out of range for the chemical potential of salt function 10 < T < 80 C");
    }

    if (S < 0.1 || S > 120) {
        throw new Error("Salinity is out of range for the chemical potential of salt function 0.1 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for the chemical potential of salt function P_sat < P < 12 MPa");
    }

    let P0;
    if (T < 100) {
        P0 = 0.101325;
    } else {
        P0 = SW_Psat(T, S) / 1e6;
    }

    if (P > P0 && (S > 42 || T > 40)) {
        throw new Error("Salinity is out of range for the chemical potential of salt function 10 < T < 40 C; 0.1 < S < 42 g/kg; P_sat < P < 12 MPa");
    }

    const b1 = -2.4176e2;
    const b2 = -6.2462e-1;
    const b3 = 7.4761e-3;
    const b4 = 1.3836e-3;
    const b5 = -6.7157e-6;
    const b6 = 5.1993e-4;
    const b7 = 9.9176e-9;
    const b8 = 6.6448e1;
    const b9 = 2.0681e-1;

    const dg_ds_P0 = b1 + b2 * T + b3 * Math.pow(T, 2) + 
                     2 * b4 * S * T + 2 * b5 * S * Math.pow(T, 2) + 
                     3 * b6 * Math.pow(S, 2) + 3 * b7 * Math.pow(S, 2) * Math.pow(T, 2) + 
                     b8 * (Math.log(S) + 1) + b9 * T * (Math.log(S) + 1);

    const c5 = -7.2431e-1;
    const c6 = 1.5712e-3;
    const c7 = -1.8919e-5;
    const c8 = 2.5939e-8;

    const dg_ds_P = (P - P0) * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3));
    const dg_ds = dg_ds_P0 + dg_ds_P;
    const mu_s = SW_Gibbs(T, S, P) + (1000 - S) * dg_ds;

    return mu_s;
}

/**
 * Chemical potential of water in seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Chemical potential of water [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_ChemPot_w(T, S, P) {
    if (T < 10 || T > 80) {
        throw new Error("Temperature is out of range for the chemical potential of water function 10 < T < 80 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for the chemical potential of water function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for the chemical potential of water function P_sat < P < 12 MPa");
    }

    let P0;
    if (T < 100) {
        P0 = 0.101325;
    } else {
        P0 = SW_Psat(T, S) / 1e6;
    }

    if (P > P0 && (S > 42 || T > 40)) {
        throw new Error("Salinity is out of range for the chemical potential of water function 10 < T < 40 C; 0 < S < 42 g/kg; P_sat < P < 12 MPa");
    }

    const b1 = -2.4176e2;
    const b2 = -6.2462e-1;
    const b3 = 7.4761e-3;
    const b4 = 1.3836e-3;
    const b5 = -6.7157e-6;
    const b6 = 5.1993e-4;
    const b7 = 9.9176e-9;
    const b8 = 6.6448e1;
    const b9 = 2.0681e-1;

    let dg_ds, Sdg_dS;
    if (S > 0) {
        const dg_ds_P0 = b1 + b2 * T + b3 * Math.pow(T, 2) + 
                         2 * b4 * S * T + 2 * b5 * S * Math.pow(T, 2) + 
                         3 * b6 * Math.pow(S, 2) + 3 * b7 * Math.pow(S, 2) * Math.pow(T, 2) + 
                         b8 * (Math.log(S) + 1) + b9 * T * (Math.log(S) + 1);

        const c5 = -7.2431e-1;
        const c6 = 1.5712e-3;
        const c7 = -1.8919e-5;
        const c8 = 2.5939e-8;

        const dg_ds_P = (P - P0) * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3));
        dg_ds = dg_ds_P0 + dg_ds_P;
        Sdg_dS = S * dg_ds;
    } else {
        Sdg_dS = 0;
    }

    return SW_Gibbs(T, S, P) - Sdg_dS;
}

/**
 * Thermal conductivity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Thermal conductivity [W/m K]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Conductivity(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for thermal conductivity function 0 < T < 180 C");
    }

    if (S < 0 || S > 160) {
        throw new Error("Salinity is out of range for thermal conductivity function 0 < S < 160 g/kg");
    }

    const T68 = 1.00024 * T;    // convert from T_90 to T_68
    const SP = S / 1.00472;     // convert from S to S_P

    return 0.001 * Math.pow(10, (Math.log10(240 + 0.0002 * SP) + 0.434 * 
           (2.3 - (343.5 + 0.037 * SP) / (T68 + 273.15)) * 
           Math.pow((1 - (T68 + 273.15) / (647.3 + 0.03 * SP)), 1/3)));
}

/**
 * Pressure-dependent thermal conductivity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Thermal conductivity [W/m K]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_ConductivityP(T, S, P) {
    if (T < 10 || T > 90) {
        throw new Error("Temperature is out of range for pressure-dependent thermal conductivity function 10 < T < 90 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for pressure-dependent thermal conductivity function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for pressure-dependent thermal conductivity function P_sat < P < 12 MPa");
    }

    const T_star = (T + 273.15) / 300;
    const P_star = (P - 0.1) / 139.9;

    const k_fw0 = 0.797015135 * Math.pow(T_star, -0.193823894) - 
                  0.251242021 * Math.pow(T_star, -4.7166384) + 
                  0.0964365893 * Math.pow(T_star, -6.38463554) - 
                  0.0326956491 * Math.pow(T_star, -2.13362102);

    const A = 13.464 * Math.pow(T_star, 4) - 
              60.727 * Math.pow(T_star, 3) + 
              102.81 * Math.pow(T_star, 2) - 
              77.387 * T_star + 21.942;

    const k_fw = k_fw0 * (1 + A * P_star);
    const B = 0.00022;
    
    return k_fw / (B * S + 1);
}

/**
 * Density of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Density [kg/m^3]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_Density(T, S, P) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for density function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for density function 0 < S < 150 g/kg");
    }

    let P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for density function P_sat < P < 12 MPa");
    }

    const S_kgkg = S / 1000;

    // Pure water density coefficients
    const a1 = 9.9992293295e2;
    const a2 = 2.0341179217e-2;
    const a3 = -6.1624591598e-3;
    const a4 = 2.2614664708e-5;
    const a5 = -4.6570659168e-8;

    // Seawater density coefficients
    const b1 = 8.0200240891e2;
    const b2 = -2.0005183488;
    const b3 = 1.6771024982e-2;
    const b4 = -3.0600536746e-5;
    const b5 = -1.6132224742e-5;

    const rho_w = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 3) + a5 * Math.pow(T, 4);
    const D_rho = b1 * S_kgkg + b2 * S_kgkg * T + b3 * S_kgkg * Math.pow(T, 2) + 
                  b4 * S_kgkg * Math.pow(T, 3) + b5 * Math.pow(S_kgkg, 2) * Math.pow(T, 2);
    const rho_sw_sharq = rho_w + D_rho;

    // Pressure dependence coefficients
    const c1 = 5.0792e-4;
    const c2 = -3.4168e-6;
    const c3 = 5.6931e-8;
    const c4 = -3.7263e-10;
    const c5 = 1.4465e-12;
    const c6 = -1.7058e-15;
    const c7 = -1.3389e-6;
    const c8 = 4.8603e-9;
    const c9 = -6.8039e-13;
    const d1 = -1.1077e-6;
    const d2 = 5.5584e-9;
    const d3 = -4.2539e-11;
    const d4 = 8.3702e-9;

    const kT = c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
               c5 * Math.pow(T, 4) + c6 * Math.pow(T, 5) + 
               P * (c7 + c8 * T + c9 * Math.pow(T, 3)) + 
               S * (d1 + d2 * T + d3 * Math.pow(T, 2) + d4 * P);

    const F_P = Math.exp((P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + 
                c4 * Math.pow(T, 3) + c5 * Math.pow(T, 4) + c6 * Math.pow(T, 5) + 
                S * (d1 + d2 * T + d3 * Math.pow(T, 2))) + 
                0.5 * (Math.pow(P, 2) - Math.pow(P0, 2)) * 
                (c7 + c8 * T + c9 * Math.pow(T, 3) + d4 * S));

    return rho_sw_sharq * F_P;
}

/**
 * Thermal diffusivity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Thermal diffusivity [m^2/s]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Diffusivity(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for diffusivity function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for diffusivity function 0 < S < 150 g/kg");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const rho = SW_Density(T, S, P0);
    const cp = SW_SpcHeat(T, S, P0);
    const K = SW_Conductivity(T, S);

    return K / (rho * cp);
}

/**
 * Specific enthalpy of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific enthalpy [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_Enthalpy(T, S, P) {
    if (T < 10 || T > 120) {
        throw new Error("Temperature is out of range for enthalpy function 10 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for enthalpy function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for enthalpy function P_sat < P < 12 MPa");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const S_kgkg = S / 1000;

    // Pure water enthalpy
    const h_w = 141.355 + 4202.07 * T - 0.535 * Math.pow(T, 2) + 0.004 * Math.pow(T, 3);

    // Seawater coefficients
    const b1 = -2.34825e4;
    const b2 = 3.15183e5;
    const b3 = 2.80269e6;
    const b4 = -1.44606e7;
    const b5 = 7.82607e3;
    const b6 = -4.41733e1;
    const b7 = 2.1394e-1;
    const b8 = -1.99108e4;
    const b9 = 2.77846e4;
    const b10 = 9.72801e1;

    // Pressure dependence coefficients
    const c1 = 996.7767;
    const c2 = -3.2406;
    const c3 = 0.0127;
    const c4 = -4.7723e-5;
    const c5 = -1.1748;
    const c6 = 0.01169;
    const c7 = -2.6185e-5;
    const c8 = 7.0661e-8;

    const h_sw_P = (P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
                   S * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3)));

    return h_w - S_kgkg * (b1 + b2 * S_kgkg + b3 * Math.pow(S_kgkg, 2) + b4 * Math.pow(S_kgkg, 3) + 
           b5 * T + b6 * Math.pow(T, 2) + b7 * Math.pow(T, 3) + b8 * S_kgkg * T + 
           b9 * Math.pow(S_kgkg, 2) * T + b10 * S_kgkg * Math.pow(T, 2)) + h_sw_P;
}

/**
 * Specific entropy of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific entropy [J/kg-K]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_Entropy(T, S, P) {
    if (T < 10 || T > 120) {
        throw new Error("Temperature is out of range for entropy function 10 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for entropy function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for entropy function P_sat < P < 12 MPa");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const S_kgkg = S / 1000;

    // Pure water entropy coefficients
    const a1 = 1.543226508e-1;
    const a2 = 1.5382700241e1;
    const a3 = -2.9963211781e-2;
    const a4 = 8.1929151062e-5;
    const a5 = -1.3699640311e-7;

    const s_w = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 3) + a5 * Math.pow(T, 4);

    // Seawater coefficients
    const b1 = -4.2307343871e2;
    const b2 = 1.4630334922e4;
    const b3 = -9.8796297642e4;
    const b4 = 3.0946224962e5;
    const b5 = 2.5623880831e1;
    const b6 = -1.4432346624e-1;
    const b7 = 5.8790568541e-4;
    const b8 = -6.110676427e1;
    const b9 = 8.0408001971e1;
    const b10 = 3.0354282687e-1;

    // Pressure dependence coefficients
    const c1 = -4.4786e-3;
    const c2 = -1.1654e-2;
    const c3 = 6.1154e-5;
    const c4 = -2.0696e-7;
    const c5 = -1.5531e-3;
    const c6 = 4.0054e-5;
    const c7 = -1.4193e-7;
    const c8 = 3.3142e-10;

    const s_sw_P = (P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
                   S * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3)));

    return s_w - S_kgkg * (b1 + b2 * S_kgkg + b3 * Math.pow(S_kgkg, 2) + b4 * Math.pow(S_kgkg, 3) + 
           b5 * T + b6 * Math.pow(T, 2) + b7 * Math.pow(T, 3) + b8 * S_kgkg * T + 
           b9 * Math.pow(S_kgkg, 2) * T + b10 * S_kgkg * Math.pow(T, 2)) + s_sw_P;
}

/**
 * Specific flow exergy of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @param {number} T0 - Total dead state temperature [°C] (ITS-90)
 * @param {number} S0 - Total dead state salinity [g/kg] (reference-composition salinity)
 * @param {number} P0 - Total dead state pressure [MPa]
 * @returns {number} Specific flow exergy [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_FlowExergy(T, S, P, T0, S0, P0) {
    if (T < 10 || T > 80) {
        throw new Error("Temperature is out of range for flow exergy function 10 < T < 80 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for flow exergy function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for flow exergy function P_sat < P < 12 MPa");
    }

    // Setting default reference states
    T0 = T0 || 25;
    S0 = S0 || 35;
    P0 = P0 || 0.101325;

    if (S0 < 0.1) {
        throw new Error("Reference salinity is out of allowed range for flow exergy function 0.1 < S0 < 120");
    }

    const h_sw = SW_Enthalpy(T, S, P);
    const s_sw = SW_Entropy(T, S, P);

    // Restricted Dead State
    const h_sw_star = SW_Enthalpy(T0, S, P0);
    const s_sw_star = SW_Entropy(T0, S, P0);
    const mu_w_star = SW_ChemPot_w(T0, S, P0);
    const Smu_s_star = SW_SChemPot_s(T0, S, P0);

    // Total Dead State
    const mu_w_0 = SW_ChemPot_w(T0, S0, P0);
    const S0mu_s_0 = SW_SChemPot_s(T0, S0, P0);
    const Smu_s_0 = SW_SChemPot_s(T0, S0, P0) * (S / S0);

    return (h_sw - h_sw_star) - (T0 + 273.15) * (s_sw - s_sw_star) + 
           (1 - 0.001 * S) * (mu_w_star - mu_w_0) + 0.001 * (Smu_s_star - Smu_s_0);
}

/**
 * Specific Gibbs energy of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific Gibbs energy [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_Gibbs(T, S, P) {
    if (T < 10 || T > 120) {
        throw new Error("Temperature is out of range for Gibbs function 10 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for Gibbs function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for Gibbs function P_sat < P < 12 MPa");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;

    // Pure water Gibbs coefficients
    const a1 = 1.0677e2;
    const a2 = -1.4303;
    const a3 = -7.6139;
    const a4 = 8.3627e-3;
    const a5 = -7.8754e-6;

    const g_w = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 3) + a5 * Math.pow(T, 4);

    // Seawater coefficients
    const b1 = -2.4176e2;
    const b2 = -6.2462e-1;
    const b3 = 7.4761e-3;
    const b4 = 1.3836e-3;
    const b5 = -6.7157e-6;
    const b6 = 5.1993e-4;
    const b7 = 9.9176e-9;
    const b8 = 6.6448e1;
    const b9 = 2.0681e-1;

    let g_sw_P0 = 0;
    if (S > 0) {
        g_sw_P0 = b1 * S + b2 * S * T + b3 * S * Math.pow(T, 2) + 
                  b4 * Math.pow(S, 2) * T + b5 * Math.pow(S, 2) * Math.pow(T, 2) + 
                  b6 * Math.pow(S, 3) + b7 * Math.pow(S, 3) * Math.pow(T, 2) + 
                  b8 * S * Math.log(S) + b9 * S * T * Math.log(S);
    }

    // Pressure dependence coefficients
    const c1 = 996.1978;
    const c2 = 3.491e-2;
    const c3 = 4.7231e-3;
    const c4 = -6.9037e-6;
    const c5 = -7.2431e-1;
    const c6 = 1.5712e-3;
    const c7 = -1.8919e-5;
    const c8 = 2.5939e-8;

    const g_sw_P = (P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
                   S * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3)));

    return g_w + g_sw_P0 + g_sw_P;
}

/**
 * Specific internal energy of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific internal energy [J/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_IntEnergy(T, S, P) {
    if (T < 10 || T > 120) {
        throw new Error("Temperature is out of range for internal energy function 10 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for internal energy function 0 < S < 120 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for internal energy function P_sat < P < 12 MPa");
    }

    const rho = SW_Density(T, S, P);
    return SW_Enthalpy(T, S, P) - (P * 1e6 / rho);
}

/**
 * Isobaric thermal expansivity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Isobaric expansivity [1/K]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_IsobExp(T, S, P) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for isobaric expansivity function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for isobaric expansivity function 0 < S < 150 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for isobaric expansivity function P_sat < P < 12 MPa");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const S_kgkg = S / 1000;

    // Pure water density coefficients
    const a1 = 9.9992293295e2;
    const a2 = 2.0341179217e-2;
    const a3 = -6.1624591598e-3;
    const a4 = 2.2614664708e-5;
    const a5 = -4.6570659168e-8;

    // Seawater density coefficients
    const b1 = 8.0200240891e2;
    const b2 = -2.0005183488;
    const b3 = 1.6771024982e-2;
    const b4 = -3.0600536746e-5;
    const b5 = -1.6132224742e-5;

    const drho_wdT = a2 + 2 * a3 * T + 3 * a4 * Math.pow(T, 2) + 4 * a5 * Math.pow(T, 3);
    const dD_rhodT = b2 * S_kgkg + 2 * b3 * S_kgkg * T + 3 * b4 * S_kgkg * Math.pow(T, 2) + 
                     2 * b5 * Math.pow(S_kgkg, 2) * T;

    const rho_w = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 3) + a5 * Math.pow(T, 4);
    const D_rho = b1 * S_kgkg + b2 * S_kgkg * T + b3 * S_kgkg * Math.pow(T, 2) + 
                  b4 * S_kgkg * Math.pow(T, 3) + b5 * Math.pow(S_kgkg, 2) * Math.pow(T, 2);

    const rho_sw_sharq = rho_w + D_rho;
    const drho_sw_sharqdT = drho_wdT + dD_rhodT;

    // Pressure dependence coefficients
    const c1 = 5.0792e-4;
    const c2 = -3.4168e-6;
    const c3 = 5.6931e-8;
    const c4 = -3.7263e-10;
    const c5 = 1.4465e-12;
    const c6 = -1.7058e-15;
    const c7 = -1.3389e-6;
    const c8 = 4.8603e-9;
    const c9 = -6.8039e-13;
    const d1 = -1.1077e-6;
    const d2 = 5.5584e-9;
    const d3 = -4.2539e-11;
    const d4 = 8.3702e-9;

    const F_P = Math.exp((P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + 
                c4 * Math.pow(T, 3) + c5 * Math.pow(T, 4) + c6 * Math.pow(T, 5) + 
                S * (d1 + d2 * T + d3 * Math.pow(T, 2))) + 
                0.5 * (Math.pow(P, 2) - Math.pow(P0, 2)) * 
                (c7 + c8 * T + c9 * Math.pow(T, 3) + d4 * S));

    const dF_PdT = F_P * ((P - P0) * (c2 + 2 * c3 * T + 3 * c4 * Math.pow(T, 2) + 
                  4 * c5 * Math.pow(T, 3) + 5 * c6 * Math.pow(T, 4) + 
                  S * (d2 + 2 * d3 * T)) + 
                  0.5 * (Math.pow(P, 2) - Math.pow(P0, 2)) * 
                  (c8 + 3 * c9 * Math.pow(T, 2)));

    const rho = rho_sw_sharq * F_P;
    const drho_dT = drho_sw_sharqdT * F_P + rho_sw_sharq * dF_PdT;

    return -drho_dT / rho;
}

/**
 * Isothermal compressibility of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Isothermal compressibility [1/MPa]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_IsothComp(T, S, P) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for isothermal compressibility function 0 < T < 180 C");
    }

    if (S < 0 || S > 160) {
        throw new Error("Salinity is out of range for isothermal compressibility function 0 < S < 160 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for isothermal compressibility function P_sat < P < 12 MPa");
    }

    const c1 = 5.0792e-4;
    const c2 = -3.4168e-6;
    const c3 = 5.6931e-8;
    const c4 = -3.7263e-10;
    const c5 = 1.4465e-12;
    const c6 = -1.7058e-15;
    const c7 = -1.3389e-6;
    const c8 = 4.8603e-9;
    const c9 = -6.8039e-13;
    const d1 = -1.1077e-6;
    const d2 = 5.5584e-9;
    const d3 = -4.2539e-11;
    const d4 = 8.3702e-9;

    return c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
           c5 * Math.pow(T, 4) + c6 * Math.pow(T, 5) + 
           P * (c7 + c8 * T + c9 * Math.pow(T, 3)) + 
           S * (d1 + d2 * T + d3 * Math.pow(T, 2) + d4 * P);
}

/**
 * Kinematic viscosity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Kinematic viscosity [m^2/s]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Kviscosity(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for kinematic viscosity function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for kinematic viscosity function 0 < S < 150 g/kg");
    }

    const P0 = T >= 100 ? SW_Psat(T, S) / 1e6 : 0.101325;
    const rho = SW_Density(T, S, P0);
    const mu = SW_Viscosity(T, S);

    return mu / rho;
}

/**
 * Latent heat of vaporization of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Latent heat of vaporization [J/kg]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_LatentHeat(T, S) {
    if (T < 0 || T > 200) {
        throw new Error("Temperature is out of range for latent heat function 0 < T < 200 C");
    }

    if (S < 0 || S > 240) {
        throw new Error("Salinity is out of range for latent heat function 0 < S < 240 g/kg");
    }

    const a1 = 2500899.1412;
    const a2 = -2369.1806479;
    const a3 = 0.26776439436;
    const a4 = -0.0081027544602;
    const a5 = -0.000020799346624;

    const hfg_w = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 3) + a5 * Math.pow(T, 4);
    return hfg_w * (1 - 0.001 * S);
}

/**
 * Osmotic coefficient of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Osmotic coefficient [-]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_OsmCoeff(T, S) {
    if (T < 0 || T > 120) {
        throw new Error("Temperature is out of range for osmotic coefficient function 0 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for osmotic coefficient function 0 < S < 120 g/kg");
    }

    const a1 = 0.89453233003;
    const a2 = 0.00041560737424;
    const a3 = -0.0000046262121398;
    const a4 = 0.000000000022211195897;
    const a5 = -0.00011445456438;
    const a6 = -0.0000014783462366;
    const a7 = -0.000000000013526263499;
    const a8 = 0.0000070132355546;
    const a9 = 0.000000056960486681;
    const a10 = -0.00000000028624032584;

    const MW_s = 31.4038218;

    if (S <= 10) {
        const S_eq = 10;  // Correlation matches function at S_equivalent = 10

        const Phi_corr_eq = a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 4) + 
                           a5 * S_eq + a6 * T * S_eq + a7 * S_eq * Math.pow(T, 3) + 
                           a8 * Math.pow(S_eq, 2) + a9 * Math.pow(S_eq, 2) * T + 
                           a10 * Math.pow(S_eq, 2) * Math.pow(T, 2);

        const dPhi_corr_eq = a5 + a6 * T + a7 * Math.pow(T, 3) + 
                            2 * a8 * S_eq + 2 * a9 * S_eq * T + 
                            2 * a10 * S_eq * Math.pow(T, 2);

        const m_sum_eq = S_eq / (1000 - S_eq) * (1000 / MW_s);
        const dmds_eq = (1000 / MW_s) * (1 / (1000 - S_eq) + S_eq / Math.pow(1000 - S_eq, 2));

        // Pitzer-Bronsted equation constants
        const beta = -2 * (Math.pow(m_sum_eq, -0.5) * (Phi_corr_eq - 1) - 
                    dPhi_corr_eq * Math.pow(m_sum_eq, 0.5) / dmds_eq);
        const lambda = (Phi_corr_eq + beta * Math.pow(m_sum_eq, 0.5) - 1) / m_sum_eq;

        // Define molality as a function of salinity
        const m_sum = S / (1000 - S) * (1000 / MW_s);
        return 1 - beta * Math.pow(m_sum, 0.5) + lambda * m_sum;
    } else {
        return a1 + a2 * T + a3 * Math.pow(T, 2) + a4 * Math.pow(T, 4) + 
               a5 * S + a6 * T * S + a7 * S * Math.pow(T, 3) + 
               a8 * Math.pow(S, 2) + a9 * Math.pow(S, 2) * T + 
               a10 * Math.pow(S, 2) * Math.pow(T, 2);
    }
}

/**
 * Osmotic pressure of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Osmotic pressure [MPa]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_OsmPress(T, S) {
    if (T < 0 || T > 120) {
        throw new Error("Temperature is out of range for osmotic pressure function 0 < T < 120 C");
    }

    if (S < 0 || S > 120) {
        throw new Error("Salinity is out of range for osmotic pressure function 0 < S < 120 g/kg");
    }

    const Phi = SW_OsmCoeff(T, S);
    const T_K = T + 273.15;
    const Mw_sw = 31.4038218;  // Weighted mol. weight in g/mol
    const R = 8.3144598;       // J/mol-K

    const m_sum = 1000 * S / ((1000 - S) * Mw_sw);  // Define molality as a function of salinity

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, 0) / 1e6;
    const rho_w_kgm3 = SW_Density(T, 0, P0);

    const Pi = Phi * m_sum * R * T_K * rho_w_kgm3;
    return Pi / 1e6;
}

/**
 * Vapor pressure of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Vapor pressure [N/m^2]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Psat(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for vapor pressure function 0 < T < 180 C");
    }

    if (S < 0 || S > 160) {
        throw new Error("Salinity is out of range for vapor pressure function 0 < S < 160 g/kg");
    }

    const T_K = T + 273.15;

    // Pure water vapor pressure coefficients
    const a1 = -5800.2206;
    const a2 = 1.3914993;
    const a3 = -0.048640239;
    const a4 = 0.000041764768;
    const a5 = -0.000000014452093;
    const a6 = 6.5459673;

    const Pv_w = Math.exp((a1 / T_K) + a2 + a3 * T_K + a4 * Math.pow(T_K, 2) + 
                 a5 * Math.pow(T_K, 3) + a6 * Math.log(T_K));

    // Seawater coefficients
    const b1 = -4.5818e-4;
    const b2 = -2.0443e-6;

    return Pv_w * Math.exp(b1 * S + b2 * Math.pow(S, 2));
}

/**
 * Prandtl number of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Prandtl number [-]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Prandtl(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for Prandtl function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for Prandtl function 0 < S < 150 g/kg");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;
    const cp = SW_SpcHeat(T, S, P0);
    const mu = SW_Viscosity(T, S);
    const K = SW_Conductivity(T, S);

    return cp * mu / K;
}

/**
 * Specific heat capacity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific heat capacity [J/kg-K]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_SpcHeat(T, S, P) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for specific heat capacity function 0 < T < 180 C");
    }

    if (S < 0 || S > 180) {
        throw new Error("Salinity is out of range for specific heat capacity function 0 < S < 180 g/kg");
    }

    const P_sat = SW_Psat(T, S) / 1e6;

    if (P < P_sat || P > 12) {
        throw new Error("Pressure is out of range for specific heat capacity function P_sat < P < 12 MPa");
    }

    const P0 = T < 100 ? 0.101325 : SW_Psat(T, S) / 1e6;

    // Convert temperature from T_90 to T_68
    const T68 = 1.00024 * T;
    // Convert from S to S_P
    const SP = S / 1.00472;

    const A = 5.328 - 9.76e-2 * SP + 4.04e-4 * Math.pow(SP, 2);
    const B = -6.913e-3 + 7.351e-4 * SP - 3.15e-6 * Math.pow(SP, 2);
    const C = 9.6e-6 - 1.927e-6 * SP + 8.23e-9 * Math.pow(SP, 2);
    const D = 2.5e-9 + 1.666e-9 * SP - 7.125e-12 * Math.pow(SP, 2);

    const cp_sw_P0 = 1000 * (A + B * T68 + C * Math.pow(T68, 2) + D * Math.pow(T68, 3));

    // Pressure dependence coefficients
    const c1 = -3.1118;
    const c2 = 0.0157;
    const c3 = 5.1014e-5;
    const c4 = -1.0302e-6;
    const c5 = 0.0107;
    const c6 = -3.9716e-5;
    const c7 = 3.2088e-8;
    const c8 = 1.0119e-9;

    const cp_sw_P = (P - P0) * (c1 + c2 * T + c3 * Math.pow(T, 2) + c4 * Math.pow(T, 3) + 
                   S * (c5 + c6 * T + c7 * Math.pow(T, 2) + c8 * Math.pow(T, 3)));

    return cp_sw_P0 + cp_sw_P;
}

/**
 * Surface tension of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Surface tension [mN/m]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_SurfaceTension(T, S) {
    if (T < 0 || T > 100) {
        throw new Error("Temperature is out of range for surface tension function 0 < T < 100 C");
    }

    if (S < 0 || S > 131) {
        throw new Error("Salinity is out of range for surface tension function 0 < S < 131 g/kg");
    }

    const T_K = T + 273.15;
    const gamma_w = 235.8 * Math.pow(1 - (T_K / 647.096), 1.256) * 
                   (1 - 0.625 * (1 - (T_K / 647.096)));

    return gamma_w * (1 + 3.766e-4 * S + 2.347e-6 * S * T);
}

/**
 * Dynamic viscosity of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @returns {number} Dynamic viscosity [kg/m-s]
 * @throws {Error} If temperature or salinity is out of range
 */
function SW_Viscosity(T, S) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for viscosity function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for viscosity function 0 < S < 150 g/kg");
    }

    const S_kgkg = S / 1000;

    // Pure water viscosity coefficients
    const a1 = 0.15700386464;
    const a2 = 64.99262005;
    const a3 = -91.296496657;
    const a4 = 0.000042844324477;

    const mu_w = a4 + 1 / (a1 * Math.pow(T + a2, 2) + a3);

    // Seawater viscosity coefficients
    const a5 = 1.540913604;
    const a6 = 0.019981117208;
    const a7 = -0.000095203865864;
    const a8 = 7.9739318223;
    const a9 = -0.075614568881;
    const a10 = 0.00047237011074;

    const A = a5 + a6 * T + a7 * Math.pow(T, 2);
    const B = a8 + a9 * T + a10 * Math.pow(T, 2);

    return mu_w * (1 + A * S_kgkg + B * Math.pow(S_kgkg, 2));
}

/**
 * Specific volume of seawater
 * @param {number} T - Temperature [°C] (ITS-90)
 * @param {number} S - Salinity [g/kg] (reference-composition salinity)
 * @param {number} P - Pressure [MPa]
 * @returns {number} Specific volume [m^3/kg]
 * @throws {Error} If temperature, salinity or pressure is out of range
 */
function SW_Volume(T, S, P) {
    if (T < 0 || T > 180) {
        throw new Error("Temperature is out of range for specific volume function 0 < T < 180 C");
    }

    if (S < 0 || S > 150) {
        throw new Error("Salinity is out of range for specific volume function 0 < S < 150 g/kg");
    }

    const rho = SW_Density(T, S, P);
    return 1 / rho;
}

// Export functions
module.exports = {
    SW_BPE,
    SW_ChemPot_s,
    SW_ChemPot_w,
    SW_Conductivity,
    SW_ConductivityP,
    SW_Density,
    SW_Diffusivity,
    SW_Enthalpy,
    SW_Entropy,
    SW_FlowExergy,
    SW_Gibbs,
    SW_IntEnergy,
    SW_IsobExp,
    SW_IsothComp,
    SW_Kviscosity,
    SW_LatentHeat,
    SW_OsmCoeff,
    SW_OsmPress,
    SW_Psat,
    SW_Prandtl,
    SW_SpcHeat,
    SW_SurfaceTension,
    SW_Viscosity,
    SW_Volume
};


// Run the example if this file is executed directly
if (typeof require !== "undefined" && require.main === module) {
    // Code runs directly
    const T = 25; // Temperature in °C
    const S = 35; // Salinity in g/kg
    const P = 0.1; // Pressure in MPa

    console.log(`Density: ${SW_Density(T, S, P)} kg/m^3`);
    console.log(`Enthalpy: ${SW_Enthalpy(T, S, P)} J/kg`);
}
