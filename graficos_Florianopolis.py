import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sin, cos, tan, acos, pi, radians, log, exp

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    PARÂMETROS DE ENTRADA DO MODELO                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ─── DADOS DA SIMULAÇÃO ───
CIDADE_ALVO = "Florianópolis"
N_DIA = 318              # Dia do ano (1-365)
V_VENTO = 2.24           # Velocidade do vento [m/s]
USAR_TROCADOR = True     # True = com trocador de calor, False = sem

# ─── CONSTANTES FÍSICAS ───
G_SC = 1367.0            # Constante solar [W/m²]
SIGMA = 5.670374e-8      # Constante de Stefan-Boltzmann [W/(m²·K⁴)]
ALBEDO = 0.20            # Refletividade do solo

# ─── PARÂMETROS DO VIDRO (Duffie & Beckman - vidro padrão) ───
N_VIDRO = 1.526          # Índice de refração
KL = 4.0                 # Coeficiente de extinção [1/m]
DELTA_V_ESP = 0.004      # Espessura do vidro [m]
EPS_VIDRO = 0.88         # Emissividade do vidro
RHO_VIDRO = 2200.0       # Densidade [kg/m³]
CP_VIDRO = 670.0         # Calor específico [J/(kg·K)]

# ─── PARÂMETROS DO PV ───
ALPHA_PV = 0.94          # Absortividade
F_EPC = 0.804            # Fração de empacotamento
ETA_REF = 0.17           # Eficiência de referência @ 25°C
BETA_PV = -0.0045        # Coeficiente de temperatura [1/°C]
T_REF = 25.0             # Temperatura de referência [°C]
EPS_PV = 0.96            # Emissividade
RHO_PV = 2330.0          # Densidade [kg/m³]
CP_PV = 700.0            # Calor específico [J/(kg·K)]
DELTA_PV_ESP = 0.0002    # Espessura das células [m]
K_PV = 148.0             # Condutividade térmica [W/(m·K)]

# ─── PARÂMETROS DO ABSORVEDOR TÉRMICO ───
RHO_AT = 8920.0          # Densidade (cobre) [kg/m³]
CP_AT = 350.0            # Calor específico [J/(kg·K)]
DELTA_AT = 0.003         # Espessura [m]
K_AT = 380.0             # Condutividade térmica [W/(m·K)]
K_AD = 0.35              # Condutividade do adesivo [W/(m·K)]
DELTA_AD = 0.00046       # Espessura do adesivo [m]

# ─── PARÂMETROS DA TUBULAÇÃO ───
N_TUBOS = 10             # Número de tubos
D_EXT = 0.01             # Diâmetro externo [m]
ESP_TUBO = 0.0001        # Espessura da parede do tubo [m]
RHO_TUBO = 8920.0        # Densidade (cobre) [kg/m³]
CP_TUBO = 350.0          # Calor específico [J/(kg·K)]

# ─── PARÂMETROS DO ISOLAMENTO ───
RHO_ISO = 20.0           # Densidade [kg/m³]
CP_ISO = 670.0           # Calor específico [J/(kg·K)]
DELTA_ISO = 0.05         # Espessura [m]
K_ISO = 0.034            # Condutividade térmica [W/(m·K)]
DELTA_ISO_BORDA = 0.025  # Espessura do isolamento lateral [m]

# ─── PARÂMETROS DA ÁGUA ───
RHO_W = 999.0            # Densidade [kg/m³]
CP_W = 4180.0            # Calor específico [J/(kg·K)]
K_W = 0.61               # Condutividade térmica [W/(m·K)]
M_DOT_NOMINAL = 0.005 * 10  # Vazão mássica total [kg/s] (0.005 kg/s por tubo)

# ─── PARÂMETROS DO AR ───
K_AR = 0.02763           # Condutividade térmica [W/(m·K)]
RHO_AR = 1.2041          # Densidade [kg/m³]
CP_AR = 1005.0           # Calor específico [J/(kg·K)]
MU_AR = 1.85e-5          # Viscosidade dinâmica [kg/(m·s)]

# ─── GEOMETRIA DO COLETOR ───
L_COLETOR = 2.0          # Comprimento [m]
W_COLETOR = 1.0          # Largura [m]
ESPACAMENTO = W_COLETOR / N_TUBOS  # Espaçamento entre tubos [m]
DELTA_GAP = 0.020        # Distância vidro-PV [m]

# ─── PARÂMETROS DO TANQUE ───
D_INT_TK = 0.62          # Diâmetro interno [m]
L_TK = 1.0               # Altura [m]
ESP_ISO_TANQUE = 0.05    # Espessura isolamento tanque [m]
K_ISO_TK = 0.037         # Condutividade isolamento tanque [W/(m·K)]
N_TK = 12                # Número de camadas do tanque

# ─── PARÂMETROS DA SIMULAÇÃO ───
T_FIM = 22 * 3600        # Tempo final (22h) [s]
DT = 0.1                 # Passo de tempo [s]
SAVE_EVERY = int(60 / DT)  # Salvar dados a cada 60 segundos

# ─── DADOS CLIMÁTICOS DIÁRIOS ───
T_MAX = 32.0             # Temperatura máxima diária [°C]
T_MIN = 24.0             # Temperatura mínima diária [°C]
T_MED = (T_MAX + T_MIN) / 2.0

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                       VARIÁVEIS GLOBAIS AUXILIARES                        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

H_h_dia = None
phi = None
beta_rad = None
Ho = None
longitude = None
L_st = None

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    CÁLCULOS GEOMÉTRICOS DERIVADOS                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

D_INT = D_EXT - 2 * ESP_TUBO
A_V = L_COLETOR * W_COLETOR
A_PV_AT = A_V * ((W_COLETOR - D_EXT) / W_COLETOR)
X_P = ESPACAMENTO / 4.0
X_AT = (ESPACAMENTO - D_EXT) / 4.0
NU_AR = MU_AR / RHO_AR
ALPHA_AR = 25.164e-6

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║             CÁLCULO DAS PERDAS LATERAIS (KALOGIROU EQ. 3.48-3.49)        ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

ESP_TOTAL_COLETOR = (DELTA_V_ESP + DELTA_GAP + DELTA_PV_ESP + 
                     DELTA_AT + D_EXT + DELTA_ISO)

PERIMETRO_COLETOR = 2.0 * (L_COLETOR + W_COLETOR)

A_LATERAL = PERIMETRO_COLETOR * ESP_TOTAL_COLETOR

def calcular_U_borda(h_vento_atual):
    """Calcula o coeficiente de perda pelas bordas [W/m²K]."""
    R_iso_borda = DELTA_ISO_BORDA / K_ISO
    R_conv_borda = 1.0 / h_vento_atual
    U_borda = 1.0 / (R_iso_borda + R_conv_borda)
    return U_borda

print(f"\n📏 GEOMETRIA DAS PERDAS LATERAIS:")
print(f"   Espessura total do coletor: {ESP_TOTAL_COLETOR*1000:.1f} mm")
print(f"   Perímetro do coletor: {PERIMETRO_COLETOR:.2f} m")
print(f"   Área lateral total: {A_LATERAL:.4f} m²")
print(f"   Área superficial do coletor: {A_V:.2f} m²")
print(f"   Relação Área Lateral / Área Superficial: {A_LATERAL/A_V*100:.2f}%")
print(f"   Isolamento lateral: {DELTA_ISO_BORDA*1000:.1f} mm")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                         FUNÇÕES GLOBAIS DO MODELO                         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def calcular_LST(longitude):
    """
    Calcula o Meridiano Padrão (Standard Meridian) baseado na longitude.
    Fusos horários são centrados em múltiplos de 15°.
    """
    fuso_teorico = round(longitude / 15.0) * 15.0
    return fuso_teorico

def declinacao(n):
    """Calcula a declinação solar para o dia n do ano (D&B Eq. 1.6.1)."""
    return radians(23.45 * sin(2 * pi * (284 + n) / 365.0))

def equacao_tempo(n):
    """Calcula a Equação do Tempo (ET) em minutos (D&B Eq. 1.5.3)."""
    B = 360.0 * (n - 81) / 364.0
    ET = 9.87 * sin(radians(2 * B)) - 7.53 * cos(radians(B)) - 1.5 * sin(radians(B))
    return ET

def hora_solar(tsec, n, L_st, longitude):
    """Calcula o Tempo Solar Aparente (TSA) em horas (D&B Eq. 1.5.2)."""
    hora_padrao = tsec / 3600.0
    ET = equacao_tempo(n)
    TSA = hora_padrao + ET / 60.0 + 4.0 * (L_st - longitude) / 60.0
    return TSA

def omega_from_t_solar(tsec, n, L_st, longitude):
    """Ângulo horário baseado no tempo solar (D&B Eq. 1.5.2)."""
    t_solar = hora_solar(tsec, n, L_st, longitude)
    return (t_solar - 12.0) * pi / 12.0

def ws_sunset(phi, delta):
    """Ângulo horário do pôr do sol [rad] (D&B Eq. 1.6.10)."""
    c = np.clip(-tan(phi) * tan(delta), -1.0, 1.0)
    return acos(c)

def Ho_dia_horizontal(lat_deg, n):
    """Radiação extraterrestre diária em kWh/m² (D&B Eq. 1.10.3)."""
    lat = radians(lat_deg)
    delta = declinacao(n)
    dr = 1 + 0.033 * cos(2 * pi * n / 365.0)
    cosws = np.clip(-tan(lat) * tan(delta), -1.0, 1.0)
    ws = acos(cosws)
    Ho_val = (24 * 3600 / pi) * G_SC * dr * (cos(lat) * cos(delta) * sin(ws) + ws * sin(lat) * sin(delta))
    return Ho_val / 3600.0

def fracao_difusa_erbs(Kt_horario):
    """
    Calcula a fração difusa usando o modelo de Erbs et al. (1982) - Duffie & Beckman.
    Modelo recomendado para fração horária de radiação difusa.
    """
    if Kt_horario <= 0.22:
        Hd_H = 1.0 - 0.09 * Kt_horario
    elif Kt_horario <= 0.80:
        Hd_H = 0.9511 - 0.1604 * Kt_horario + 4.388 * Kt_horario**2 - 16.638 * Kt_horario**3 + 12.336 * Kt_horario**4
    else:
        Hd_H = 0.165
    
    return np.clip(Hd_H, 0.0, 1.0)

def fracao_difusa_diaria_erbs(Kt_diario):
    """
    Calcula a fração difusa diária usando o modelo de Erbs et al. (1982).
    """
    if Kt_diario <= 0.22:
        Hd_H = 1.0 - 0.09 * Kt_diario
    elif Kt_diario <= 0.80:
        Hd_H = 0.9511 - 0.1604 * Kt_diario + 4.388 * Kt_diario**2 - 16.638 * Kt_diario**3 + 12.336 * Kt_diario**4
    else:
        Hd_H = 0.165
    
    return np.clip(Hd_H, 0.0, 1.0)

def Rb_DB(omega, phi, beta_rad):
    """
    Fator de correção para radiação direta em superfície inclinada (D&B Eq. 1.8.2).
    ADAPTADO PARA HEMISFÉRIO SUL: Usa (phi - beta_rad) para superfície voltada ao Equador.
    """
    delta = declinacao(N_DIA)
    
    # Para o Hemisfério Sul, a superfície voltada ao Equador (Norte) tem inclinação beta
    # O ângulo de incidência corrigido usa (phi - beta)
    num = sin(delta) * sin(phi - beta_rad) + cos(delta) * cos(phi - beta_rad) * cos(omega)
    den = sin(delta) * sin(phi) + cos(delta) * cos(phi) * cos(omega)
    
    if den <= 0.001:
        return 0.0
    return max(0.0, num / den)

def transmissao_vidro_db(theta1, n):
    """
    Calcula a transmitância do vidro usando as equações de Fresnel (D&B Cap. 5).
    Retorna transmitância para polarização perpendicular e paralela.
    """
    if theta1 >= pi/2:
        return 0.0, 0.0
    
    # Lei de Snell (D&B Eq. 5.1.1)
    theta2 = np.arcsin(np.sin(theta1) / n)
    
    # Reflectâncias para polarização perpendicular e paralela (D&B Eq. 5.1.2)
    r_perp = (np.sin(theta2 - theta1) / np.sin(theta2 + theta1))**2
    r_par = (np.tan(theta2 - theta1) / np.tan(theta2 + theta1))**2
    
    # Transmitâncias
    tau_perp = 1 - r_perp
    tau_par = 1 - r_par
    
    return tau_perp, tau_par

def transmitancia_absorcao_db(KL, theta2):
    """
    Calcula a transmitância devido à absorção (D&B Eq. 5.3.1).
    """
    if abs(np.cos(theta2)) < 1e-6:
        return 0.0
    return np.exp(-KL / np.cos(theta2))

def calcular_propriedades_opticas_DB(tsec):
    """
    Calcula propriedades ópticas do vidro usando metodologia Duffie & Beckman.
    Considera múltiplas reflexões internas (D&B Eqs. 5.2.1 - 5.2.4).
    """
    global phi, beta_rad, L_st, longitude
    
    delta = declinacao(N_DIA)
    omega = omega_from_t_solar(tsec, N_DIA, L_st, longitude)
    
    # Ângulo de incidência para superfície inclinada (D&B Eq. 1.6.2)
    # Para Hemisfério Sul: superfície inclinada para Norte (azimute = 0°)
    cos_theta = sin(delta) * sin(phi - beta_rad) + cos(delta) * cos(phi - beta_rad) * cos(omega)
    cos_theta = np.clip(cos_theta, 0.0, 1.0)
    theta1 = np.arccos(cos_theta)
    
    # Se ângulo de incidência > 90°, sem radiação
    if theta1 >= pi/2 or cos_theta <= 0:
        return 0.0, 0.0, 1.0, cos_theta
    
    # Calcular transmitância pela lei de Bouguer (D&B Eq. 5.3.1)
    # Para vidro com índice de refração n e coeficiente de extinção KL
    theta2 = np.arcsin(np.sin(theta1) / N_VIDRO)
    tau_a = transmitancia_absorcao_db(KL * DELTA_V_ESP, theta2)
    
    # Calcular transmitância por reflexão usando método das múltiplas reflexões (D&B Eq. 5.2.1)
    tau_r_perp, tau_r_par = transmissao_vidro_db(theta1, N_VIDRO)
    
    # Para vidro sem revestimento, usar Eqs. 5.2.3 e 5.2.4
    # Reflectâncias para múltiplas reflexões
    r_perp = 1 - tau_r_perp
    r_par = 1 - tau_r_par
    
    # Transmitância considerando múltiplas reflexões (D&B Eq. 5.2.3)
    tau_r_perp_total = tau_r_perp * tau_a / (1 - (1 - tau_r_perp) * tau_a**2)
    tau_r_par_total = tau_r_par * tau_a / (1 - (1 - tau_r_par) * tau_a**2)
    
    # Transmitância total para radiação não polarizada (média)
    tau_v = (tau_r_perp_total + tau_r_par_total) / 2.0
    
    # Absortância do vidro (D&B Eq. 5.3.2)
    alpha_v = (1 - tau_a) * (1 - (tau_r_perp + tau_r_par)/2) / (1 - (tau_r_perp + tau_r_par)/2 * tau_a)
    
    # Refletância
    rho_v = 1.0 - alpha_v - tau_v
    rho_v = np.clip(rho_v, 0.0, 1.0)
    
    return alpha_v, tau_v, rho_v, cos_theta

def Gt_instantaneo_DB(tsec):
    """
    Calcula a irradiância solar incidente no plano do coletor [W/m²] usando Duffie & Beckman.
    """
    global H_h_dia, phi, beta_rad, Ho, longitude, L_st
    h = tsec / 3600.0
    
    delta = declinacao(N_DIA)
    ws = ws_sunset(phi, delta)
    hora_por_sol = 12.0 + ws * 12.0 / pi
    hora_nascer_sol = 12.0 - ws * 12.0 / pi
    
    if h < hora_nascer_sol or h > hora_por_sol:
        return 0.0
    
    omega = omega_from_t_solar(tsec, N_DIA, L_st, longitude)
    
    # Calcular Kt horário a partir do diário (método de Collares-Pereira & Rabl)
    ws_graus = ws * 180.0 / pi
    Kt_diario = H_h_dia / (Ho * 3600.0)
    
    # Parâmetros para distribuição horária (Collares-Pereira & Rabl 1979)
    a = 0.4090 + 0.5016 * np.sin(ws - 1.047)
    b = 0.6609 - 0.4767 * np.sin(ws - 1.047)
    
    # Radiação horária total
    r_t = (np.pi / 24.0) * (a + b * np.cos(omega)) * (np.cos(omega) - np.cos(ws)) / (np.sin(ws) - ws * np.cos(ws))
    r_t = np.clip(r_t, 0.0, 1.0)
    I_hor = H_h_dia * r_t
    
    # Fração difusa diária
    Hd_H_dia = fracao_difusa_diaria_erbs(Kt_diario)
    H_d_dia = Hd_H_dia * H_h_dia
    
    # Fração difusa horária (usando correlação de Erbs)
    Kt_horario = I_hor / G_SC  # Simplificado para radiação horária
    Hd_H_hor = fracao_difusa_erbs(Kt_horario)
    I_d_hor = Hd_H_hor * I_hor
    I_b_hor = np.clip(I_hor - I_d_hor, 0.0, I_hor)
    
    # Componentes no plano inclinado (D&B Eqs. 2.12.1 - 2.12.3)
    G_b_tilt = I_b_hor * Rb_DB(omega, phi, beta_rad)
    G_d_tilt = I_d_hor * (1 + np.cos(beta_rad)) / 2.0
    G_r_tilt = I_hor * ALBEDO * (1 - np.cos(beta_rad)) / 2.0
    
    return np.clip(G_b_tilt + G_d_tilt + G_r_tilt, 0.0, None)

def T_amb_fun(tsec):
    """Temperatura ambiente [°C] - Modelo de Ephrath et al. (1996)."""
    global L_st, longitude
    omega = omega_from_t_solar(tsec, N_DIA, L_st, longitude)
    T_amb = T_MED + (T_MAX - T_MIN) / 2.0 * np.cos(omega - pi / 4.0)
    return T_amb

def eta_el_T(Tcell_C):
    """Eficiência elétrica em função da temperatura da célula [%]."""
    return max(0.0, ETA_REF * (1.0 + BETA_PV * (Tcell_C - T_REF)))

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                  FUNÇÃO PRINCIPAL DE SIMULAÇÃO                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

def simular(use_hex=True):
    """
    Simulação do coletor PV/T com modelo Duffie & Beckman.
    use_hex: True = com trocador de calor, False = sem
    """
    global H_h_dia, phi, beta_rad, Ho, longitude, L_st
    
    # ── Fator de forma para condução tubo-PV ──
    FT_tubo = (DELTA_PV_ESP * L_COLETOR) / ((X_P/(2.0*K_PV)) + ((DELTA_AD*DELTA_PV_ESP)/(K_AD*D_EXT)))
    
    # ── Capacitâncias térmicas ──
    V_v = A_V * DELTA_V_ESP
    C_v = V_v * RHO_VIDRO * CP_VIDRO
    
    V_pv = A_V * DELTA_PV_ESP
    C_pv = V_pv * RHO_PV * CP_PV
    
    V_at = A_V * DELTA_AT
    C_at = V_at * RHO_AT * CP_AT
    
    V_iso_total = DELTA_ISO * L_COLETOR * W_COLETOR
    V_tubo_externo = N_TUBOS * ((D_EXT**2)/2.0) * (1.0 + pi/4.0) * L_COLETOR
    V_iso = V_iso_total - V_tubo_externo
    C_iso = V_iso * RHO_ISO * CP_ISO
    
    V_tubo = N_TUBOS * pi * ((D_EXT/2)**2 - (D_INT/2)**2) * L_COLETOR
    C_tubo = V_tubo * RHO_TUBO * CP_TUBO
    
    v_tubo_int = N_TUBOS * (pi/4.0) * D_INT**2 * L_COLETOR
    M_water = RHO_W * v_tubo_int
    C_agua = M_water * CP_W
    A_tubo_water = N_TUBOS * pi * D_INT * L_COLETOR
    Nua_water = 4.36
    h_water = Nua_water * K_W / D_INT
    A_tubo_iso = N_TUBOS * ((pi/2.0) + 1.0) * D_EXT * L_COLETOR
    
    def h_gap_vpv(Tg, Tpv, beta_rad):
        eps_eff = 1.0 / (1.0/EPS_VIDRO + 1.0/EPS_PV - 1.0)
        TgK, TpK = Tg + 273.15, Tpv + 273.15
        h_rad = eps_eff * SIGMA * (TgK**2 + TpK**2) * (TgK + TpK)
        Tm = (Tg + Tpv) / 2.0
        beta_exp = 2.0 / (TgK + TpK)
        dT = max(abs(Tg - Tpv), 1e-6)
        Ra = 9.81 * beta_exp * dT * (DELTA_GAP**3) / (ALPHA_AR * NU_AR)
        Ra = max(Ra, 1e-12)
        X = Ra * np.cos(beta_rad)
        X = max(X, 1e-12)
        def pos(z):
            return max(z, 0.0)
        term1 = pos(1.0 - 1708.0 / X)
        sin_term = abs(np.sin(1.8 * beta_rad))
        term2 = pos(1.0 - 1708.0 * (sin_term)**1.6 / X)
        term3 = pos((X / 5830.0)**(1.0/3.0) - 1.0)
        Nu = 1.0 + 1.44 * term1 * term2 * term3
        h_conv = Nu * K_AR / DELTA_GAP
        return h_conv + h_rad
    
    # ── Leitura dos dados da cidade ──
    df = pd.read_csv("irradiação_global_horizontal.csv", sep=";", encoding="utf-8")
    df["NAME"] = df["NAME"].str.strip()
    row = df[df["NAME"].str.lower() == CIDADE_ALVO.lower()].iloc[0]
    H_h_dia = float(row["NOV"])
    if H_h_dia < 100.0:
        H_h_dia *= 1000.0
    
    latitude = float(row["LAT"])
    longitude = float(row["LON"])
    phi = radians(latitude)
    
    L_st = calcular_LST(longitude)
    beta_rad = radians(abs(latitude))
    Ho = Ho_dia_horizontal(latitude, N_DIA)
    
    print(f"\n📍 DADOS DA CIDADE:")
    print(f"   Cidade: {CIDADE_ALVO}")
    print(f"   Latitude: {latitude:.2f}°")
    print(f"   Longitude: {longitude:.2f}° (do CSV)")
    print(f"   Meridiano Padrão (L_st) calculado: {L_st}°")
    
    # Calcular hora do nascer e pôr do sol (apenas para referência)
    delta_ref = declinacao(N_DIA)
    ws_ref = ws_sunset(phi, delta_ref)
    hora_nascer_sol = 12.0 - ws_ref * 12.0 / pi
    hora_por_sol = 12.0 + ws_ref * 12.0 / pi
    
    # Iniciar simulação no nascer do sol
    T_INI = hora_nascer_sol * 3600
    
    # Ajustar T_FIM para garantir que cubra pelo menos até o pôr do sol
    T_FIM = max(22 * 3600, (hora_por_sol + 2) * 3600)
    
    print(f"  Nascer do sol: {hora_nascer_sol:.2f}h")
    print(f"  Pôr do sol: {hora_por_sol:.2f}h")
    print(f"  Controle da bomba: desliga quando T_pv ≤ T_agua")
    
    # ── Inicialização do tanque ──
    T_tk = None
    if use_hex:
        D_ext_tk = D_INT_TK + 2 * ESP_ISO_TANQUE
        A_cs = 0.25 * np.pi * D_INT_TK**2
        dz = L_TK / N_TK
        V_i = A_cs * dz
        C_tk_i = RHO_W * V_i * CP_W
        A_lat_seg = np.pi * D_ext_tk * dz
        A_topo = 0.25 * np.pi * D_ext_tk**2
        A_fundo = 0.25 * np.pi * D_ext_tk**2
        
        def U_lateral(Ts, Ta):
            Tm = (Ts + Ta) / 2.0
            betaT = 1.0 / max(Tm + 273.15, 1.0)
            dT = max(abs(Ts - Ta), 1e-6)
            Ra = 9.81 * betaT * dT * (D_ext_tk**3) / (ALPHA_AR * NU_AR)
            Ra = max(Ra, 1e-12)
            Pr = 0.71
            Nu = 0.36 + (0.518 * Ra**0.25) / ((1 + (0.559/Pr)**(9/16))**(4/9))
            h_ext = Nu * K_AR / D_ext_tk
            R_iso = log(D_ext_tk / D_INT_TK) / (2.0 * np.pi * K_ISO_TK * L_TK)
            R_ext = 1.0 / (h_ext * np.pi * D_ext_tk * L_TK)
            U = 1.0 / (R_iso + R_ext)
            return max(U, 0.01)
        
        def U_topo_fundo(Ts, Ta):
            dT = max(abs(Ts - Ta), 1e-6)
            Lc = D_ext_tk / 4.0
            Tm = (Ts + Ta) / 2.0
            betaT = 1.0 / max(Tm + 273.15, 1.0)
            Ra = 9.81 * betaT * dT * (Lc**3) / (ALPHA_AR * NU_AR)
            Ra = max(Ra, 1e-12)
            if dT > 0:
                Nu = 0.54 * Ra**0.25 if Ra <= 1e7 else 0.15 * Ra**(1/3)
            else:
                Nu = 0.27 * Ra**0.25 if Ra <= 1e7 else 0.15 * Ra**(1/3)
            h_ext = Nu * K_AR / Lc
            R_iso = ESP_ISO_TANQUE / K_ISO_TK
            U = 1.0 / (R_iso + 1.0/max(h_ext, 1e-6))
            return max(U, 0.05)
        
        T_tk = np.ones(N_TK) * 25.0
    
    # ── Condições iniciais ──
    T_v = 26.0
    T_pv = 26.0
    T_at = 26.0
    T_tubo = 26.0
    T_iso = 26.0
    T_agua = 25.0
    
    # ── Histórico ──
    t_hist, Gt_hist, Ta_hist = [], [], []
    T_v_hist, T_pv_hist, T_at_hist, T_tubo_hist, T_iso_hist, T_agua_hist = [], [], [], [], [], []
    T_tk_hist = [[] for _ in range(N_TK)]
    m_dot_hist, eta_hist, P_el_hist = [], [], []
    Q_bordas_hist = []
    tau_v_hist = []
    alpha_v_hist = []
    theta_inc_hist = []
    
    m_dot_atual = M_DOT_NOMINAL if use_hex else 0.0
    t = T_INI
    step = 0
    hora_desliga_bomba = None
    
    print(f"  Simulando {'COM' if use_hex else 'SEM'} trocador (modelo Duffie & Beckman)...")
    
    while t <= T_FIM:
        hora_atual = t / 3600.0
        Gt = Gt_instantaneo_DB(t)
        
        # ============================================================
        # CONTROLE DA BOMBA:
        # Desliga quando a temperatura do módulo PV (T_pv) for 
        # menor ou igual à temperatura de saída da água (T_agua)
        # ============================================================
        if use_hex:
            if m_dot_atual > 0 and T_pv <= T_agua and hora_desliga_bomba is None:
                hora_desliga_bomba = hora_atual
                print(f"  ⏹️ Bomba desligada às {hora_desliga_bomba:.2f}h (T_pv = {T_pv:.1f}°C ≤ T_agua = {T_agua:.1f}°C)")
            
            if T_pv <= T_agua:
                m_dot_atual = 0.0
            else:
                m_dot_atual = M_DOT_NOMINAL
        
        # Usar modelo Duffie & Beckman para propriedades ópticas
        alpha_v, tau_v, rho_v_val, costh = calcular_propriedades_opticas_DB(t)
        Ta = T_amb_fun(t)
        
        # Armazenar propriedades ópticas para análise
        if step % SAVE_EVERY == 0:
            tau_v_hist.append(tau_v)
            alpha_v_hist.append(alpha_v)
            if costh > 0:
                theta_deg = np.arccos(costh) * 180.0 / np.pi
            else:
                theta_deg = 90.0
            theta_inc_hist.append(theta_deg)
        
        eta_now = eta_el_T(T_pv)
        P_el = eta_now * Gt * A_V
        
        # Produto óptico efetivo para PV (considerando múltiplas reflexões)
        # D&B Eq. 5.5.1
        if costh > 0.01:
            # Para células encapsuladas, considerar absorção no PV
            alpha_tau_pv = tau_v * ALPHA_PV / (1 - (1 - ALPHA_PV) * rho_v_val)
        else:
            alpha_tau_pv = 0.0
        
        h_gap = h_gap_vpv(T_v, T_pv, beta_rad)
        h_vento = 2.8 + 3.0 * V_VENTO
        
        # ╔══════════════════════════════════════════════════════════════════════╗
        # ║              PERDAS LATERAIS (KALOGIROU EQ. 3.48-3.49)            ║
        # ║              DISTRIBUIÇÃO PROPORCIONAL À TEMPERATURA              ║
        # ╚══════════════════════════════════════════════════════════════════════╝
        
        T_media_coletor = (T_v + T_pv + T_at + T_tubo + T_iso) / 5.0
        U_borda = calcular_U_borda(h_vento)
        Q_dot_bordas_total = U_borda * A_LATERAL * (T_media_coletor - Ta)
        
        # Distribuição proporcional às diferenças de temperatura
        delta_T_at = max(0, T_at - Ta)
        delta_T_iso = max(0, T_iso - Ta)
        delta_T_total = delta_T_at + delta_T_iso
        
        if delta_T_total > 0:
            frac_at = delta_T_at / delta_T_total
            frac_iso = delta_T_iso / delta_T_total
        else:
            # Fallback: valores padrão (caso de temperaturas iguais ou menores que ambiente)
            frac_at = 0.30
            frac_iso = 0.70
        
        Q_dot_bordas_at = Q_dot_bordas_total * frac_at
        Q_dot_bordas_iso = Q_dot_bordas_total * frac_iso
        
        # ── Balanços de energia ──
        # Vidro
        Q_dot_v = A_V * Gt * alpha_v
        Q_dot_conv_v_ar = h_vento * A_V * (Ta - T_v)
        Q_dot_rad_v_ar = EPS_VIDRO * A_V * SIGMA * ((Ta+273.15)**4 - (T_v+273.15)**4)
        Q_dot_conv_v_pv = h_gap * A_V * (T_pv - T_v)
        epsilon_pv_v = 1.0 / (1.0/EPS_VIDRO + 1.0/EPS_PV - 1.0)
        Q_dot_rad_v_pv = A_V * SIGMA * ((T_pv+273.15)**4 - (T_v+273.15)**4) / epsilon_pv_v
        
        # Células PV
        E_p = Gt * A_V * F_EPC * eta_now
        Q_dot_pv = A_V * Gt * alpha_tau_pv - E_p
        Q_dot_rad_pv_v = A_V * SIGMA * ((T_v+273.15)**4 - (T_pv+273.15)**4) / epsilon_pv_v
        Q_dot_conv_pv_v = h_gap * A_V * (T_v - T_pv)
        Q_dot_cond_pv_at = (K_AD/DELTA_AD) * A_PV_AT * (T_at - T_pv)
        Q_dot_cond_pv_tubo = FT_tubo * N_TUBOS * (T_tubo - T_pv)
        
        # Absorvedor térmico
        Q_dot_cond_at_pv = (K_AD/DELTA_AD) * A_PV_AT * (T_pv - T_at)
        Q_dot_cond_at_tubo = (2.0*K_AT/X_AT) * (DELTA_AT*L_COLETOR) * N_TUBOS * (T_tubo - T_at)
        Q_dot_cond_at_iso = (2.0*K_ISO/DELTA_ISO) * A_PV_AT * (T_iso - T_at)
        
        # Tubulação
        Q_dot_cond_tubo_at = (2.0*K_AT/X_AT) * (DELTA_AT*L_COLETOR) * N_TUBOS * (T_at - T_tubo)
        Q_dot_cond_tubo_pv = FT_tubo * N_TUBOS * (T_pv - T_tubo)
        Q_dot_cond_tubo_iso = (2.0*K_ISO/DELTA_ISO) * A_tubo_iso * (T_iso - T_tubo)
        Q_dot_conv_tubo_water = h_water * A_tubo_water * (T_agua - T_tubo)
        
        # Isolamento
        Q_dot_cond_iso_at = (2.0*K_ISO/DELTA_ISO) * A_PV_AT * (T_at - T_iso)
        Q_dot_cond_iso_tubo = (2.0*K_ISO/DELTA_ISO) * A_tubo_iso * (T_tubo - T_iso)
        h_iso_ar = 1.0 / ((DELTA_ISO/(2.0*K_ISO)) + (1.0/h_vento))
        Q_dot_conv_iso_ar = h_iso_ar * A_V * (Ta - T_iso)
        
        # Água
        Q_dot_conv_water_tubo = h_water * A_tubo_water * (T_tubo - T_agua)
        
        if m_dot_atual > 0 and T_tk is not None:
            T_in_HX = T_tk[-1]
            Q_dot_water = m_dot_atual * CP_W * (T_in_HX - T_agua)
        else:
            Q_dot_water = 0.0
        
        # ── Atualização das temperaturas ──
        T_v += (DT / C_v) * (Q_dot_v + Q_dot_conv_v_ar + Q_dot_rad_v_ar + Q_dot_conv_v_pv + Q_dot_rad_v_pv)
        T_pv += (DT / C_pv) * (Q_dot_pv + Q_dot_rad_pv_v + Q_dot_conv_pv_v + Q_dot_cond_pv_at + Q_dot_cond_pv_tubo)
        T_at += (DT / C_at) * (Q_dot_cond_at_pv + Q_dot_cond_at_tubo + Q_dot_cond_at_iso - Q_dot_bordas_at)
        T_tubo += (DT / C_tubo) * (Q_dot_cond_tubo_at + Q_dot_cond_tubo_pv + Q_dot_cond_tubo_iso + Q_dot_conv_tubo_water)
        T_iso += (DT / C_iso) * (Q_dot_cond_iso_at + Q_dot_cond_iso_tubo + Q_dot_conv_iso_ar - Q_dot_bordas_iso)
        T_agua += (DT / C_agua) * (Q_dot_conv_water_tubo + Q_dot_water)
        
        # ── Limites físicos ──
        T_v = np.clip(T_v, 10.0, 150)
        T_pv = np.clip(T_pv, 10.0, 150)
        T_at = np.clip(T_at, 10.0, 150)
        T_tubo = np.clip(T_tubo, 10.0, 150)
        T_iso = np.clip(T_iso, 10.0, 150)
        T_agua = np.clip(T_agua, 10.0, 100)
        
        # ╔══════════════════════════════════════════════════════════════════════╗
        # ║           ATUALIZAÇÃO DO TANQUE ESTRATIFICADO                      ║
        # ╚══════════════════════════════════════════════════════════════════════╝
        
        if use_hex and T_tk is not None and m_dot_atual > 0:
            # Bomba ligada: advecção + condução + perdas
            T_tk_n = T_tk.copy()
            v_tk = m_dot_atual / (RHO_W * A_cs)
            adv_fac = v_tk * DT / dz
            cond_fac = K_W * A_cs * DT / (RHO_W * CP_W * V_i * dz)
            
            for i in range(N_TK):
                if i == 0:
                    Qperda_topo = U_topo_fundo(T_tk_n[i], Ta) * A_topo * (Ta - T_tk_n[i])
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_topo + Qperda_lat
                elif i == N_TK - 1:
                    Qperda_fundo = U_topo_fundo(T_tk_n[i], Ta) * A_fundo * (Ta - T_tk_n[i])
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_fundo + Qperda_lat
                else:
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_lat
                
                T_up = T_agua if i == 0 else T_tk_n[i-1]
                adv_term = adv_fac * (T_up - T_tk_n[i])
                T_im1 = T_tk_n[i] if i == 0 else T_tk_n[i-1]
                T_ip1 = T_tk_n[i] if i == N_TK-1 else T_tk_n[i+1]
                cond_term = cond_fac * (T_im1 - 2*T_tk_n[i] + T_ip1)
                dTi = adv_term + cond_term + Qloss / C_tk_i
                T_tk[i] = T_tk_n[i] + dTi
            
            for _ in range(5):
                for i in range(N_TK - 1):
                    if T_tk[i] < T_tk[i + 1]:
                        T_media = (T_tk[i] + T_tk[i + 1]) / 2.0
                        T_tk[i] = T_media
                        T_tk[i + 1] = T_media
            
            T_tk = np.clip(T_tk, 10.0, 100)
        
        elif use_hex and T_tk is not None and m_dot_atual == 0:
            # Bomba desligada: apenas condução + perdas (sem advecção)
            T_tk_n = T_tk.copy()
            cond_fac = K_W * A_cs * DT / (RHO_W * CP_W * V_i * dz)
            
            for i in range(N_TK):
                if i == 0:
                    Qperda_topo = U_topo_fundo(T_tk_n[i], Ta) * A_topo * (Ta - T_tk_n[i])
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_topo + Qperda_lat
                elif i == N_TK - 1:
                    Qperda_fundo = U_topo_fundo(T_tk_n[i], Ta) * A_fundo * (Ta - T_tk_n[i])
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_fundo + Qperda_lat
                else:
                    Qperda_lat = U_lateral(T_tk_n[i], Ta) * A_lat_seg * (Ta - T_tk_n[i])
                    Qloss = Qperda_lat
                
                T_im1 = T_tk_n[i] if i == 0 else T_tk_n[i-1]
                T_ip1 = T_tk_n[i] if i == N_TK-1 else T_tk_n[i+1]
                cond_term = cond_fac * (T_im1 - 2*T_tk_n[i] + T_ip1)
                dTi = cond_term + Qloss / C_tk_i
                T_tk[i] = T_tk_n[i] + dTi
            
            for _ in range(5):
                for i in range(N_TK - 1):
                    if T_tk[i] < T_tk[i + 1]:
                        T_media = (T_tk[i] + T_tk[i + 1]) / 2.0
                        T_tk[i] = T_media
                        T_tk[i + 1] = T_media
            
            T_tk = np.clip(T_tk, 10.0, 100)
        
        # ── Salvamento de dados ──
        if step % SAVE_EVERY == 0:
            t_hist.append(t / 3600.0)
            Gt_hist.append(Gt)
            Ta_hist.append(Ta)
            T_v_hist.append(T_v)
            T_pv_hist.append(T_pv)
            T_at_hist.append(T_at)
            T_tubo_hist.append(T_tubo)
            T_iso_hist.append(T_iso)
            T_agua_hist.append(T_agua)
            m_dot_hist.append(m_dot_atual)
            eta_hist.append(eta_now * 100)
            P_el_hist.append(P_el)
            Q_bordas_hist.append(Q_dot_bordas_total)
            if use_hex and T_tk is not None:
                for i in range(N_TK):
                    T_tk_hist[i].append(T_tk[i])
        
        t += DT
        step += 1
    
    # ── Montagem do dicionário de resultados ──
    resultados = {
        't': np.array(t_hist), 'Gt': np.array(Gt_hist), 'Ta': np.array(Ta_hist),
        'T_v': np.array(T_v_hist), 'T_pv': np.array(T_pv_hist),
        'T_at': np.array(T_at_hist), 'T_tubo': np.array(T_tubo_hist),
        'T_iso': np.array(T_iso_hist), 'T_agua': np.array(T_agua_hist),
        'm_dot': np.array(m_dot_hist), 'eta': np.array(eta_hist),
        'P_el': np.array(P_el_hist), 'Q_bordas': np.array(Q_bordas_hist),
        'hora_desliga_bomba': hora_desliga_bomba,
        'tau_v': np.array(tau_v_hist) if tau_v_hist else None,
        'alpha_v': np.array(alpha_v_hist) if alpha_v_hist else None,
        'theta_inc': np.array(theta_inc_hist) if theta_inc_hist else None
    }
    
    if use_hex and T_tk is not None:
        for i in range(N_TK):
            resultados[f'T_tk_{i}'] = np.array(T_tk_hist[i])
        resultados['T_tk_top'] = np.array(T_tk_hist[0])
        resultados['T_tk_bot'] = np.array(T_tk_hist[-1])
    
    return resultados


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║                    EXECUÇÃO PRINCIPAL E GRÁFICOS                           ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print(f"MODELO PV/T - SIMULAÇÃO TÉRMICA DETALHADA")
    print(f"Modelo: Duffie & Beckman")
    print(f"Cidade: {CIDADE_ALVO}")
    print(f"Dia: {N_DIA}")
    print(f"Vento: {V_VENTO} m/s")
    print(f"Controle da bomba: desliga quando T_pv ≤ T_agua")
    print("=" * 70)
    
    print("\n[1/2] Simulando COM trocador...")
    res_com = simular(use_hex=True)
    
    print("\n[2/2] Simulando SEM trocador...")
    res_sem = simular(use_hex=False)
    
    from numpy import trapz
    
    energia_com_J = trapz(res_com['P_el'], res_com['t'] * 3600)
    energia_sem_J = trapz(res_sem['P_el'], res_sem['t'] * 3600)
    energia_com_kWh = energia_com_J / 3_600_000
    energia_sem_kWh = energia_sem_J / 3_600_000
    ganho = ((energia_com_kWh - energia_sem_kWh) / energia_sem_kWh * 100) if energia_sem_kWh > 0 else 0
    
    Q_bordas_total_com = trapz(res_com['Q_bordas'], res_com['t'] * 3600) / 1_000_000
    Q_bordas_total_sem = trapz(res_sem['Q_bordas'], res_sem['t'] * 3600) / 1_000_000
    
    print(f"\nEnergia COM trocador: {energia_com_kWh:.4f} kWh")
    print(f"Energia SEM trocador: {energia_sem_kWh:.4f} kWh")
    
    # ══════════════════════════════════════════════════════════════════════════
    #                    GRÁFICOS 2D COM LINHA DA BOMBA                         
    # ══════════════════════════════════════════════════════════════════════════
    
    # ── GRÁFICO 1: Radiação Solar + Temperatura Ambiente ──
    fig, ax1 = plt.subplots(figsize=(12, 6))
    color = 'orange'
    ax1.set_xlabel('Horário (h)', fontsize=12)
    ax1.set_ylabel('Radiação (W/m²)', fontsize=12, color=color)
    ax1.plot(res_com['t'], res_com['Gt'], color=color, linewidth=2, label='Gt - Radiação Incidente')
    ax1.fill_between(res_com['t'], 0, res_com['Gt'], alpha=0.2, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim(0, max(res_com['Gt']) * 1.2)
    ax1.set_xlim(res_com['t'][0], 22)
    ax1.grid(True, alpha=0.3)
    
    # Adicionar linha vertical do desligamento da bomba
    if res_com['hora_desliga_bomba'] is not None:
        HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
        ax1.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.0f}h)')
    
    ax2 = ax1.twinx()
    color2 = 'blue'
    ax2.set_ylabel('Temperatura Ambiente (°C)', fontsize=12, color=color2)
    ax2.plot(res_com['t'], res_com['Ta'], color=color2, linewidth=2, linestyle='--', label='T_amb')
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylim(20, 35)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=10)
    plt.title('Radiação Solar Incidente e Temperatura Ambiente', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'1_radiacao_temperatura_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 2: Temperatura do Módulo PV ──
    plt.figure(figsize=(12, 6))
    plt.plot(res_com['t'], res_com['T_pv'], 'b-', linewidth=2, label='PV - Com Trocador')
    plt.plot(res_sem['t'], res_sem['T_pv'], 'r--', linewidth=2, label='PV - Sem Trocador')
    plt.plot(res_com['t'], res_com['Ta'], 'g:', linewidth=1.5, label='Temperatura Ambiente')
    plt.axhline(y=25, color='gray', linestyle=':', alpha=0.5, label='T_inicial (25°C)')
    
    # Adicionar linha vertical do desligamento da bomba
    if res_com['hora_desliga_bomba'] is not None:
        HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
        plt.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.2f}h)')
    
    plt.ylabel('Temperatura (°C)', fontsize=12)
    plt.xlabel('Horário (h)', fontsize=12)
    plt.title('Temperatura do Módulo PV - Com e Sem Trocador', fontsize=14)
    plt.legend(loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(res_com['t'][0], 22)
    
    reducao = max(res_sem['T_pv']) - max(res_com['T_pv'])
    plt.figtext(0.5, 0.02, f'Redução máxima de temperatura PV com trocador: {reducao:.1f}°C', 
                ha='center', fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.tight_layout()
    plt.savefig(f'2_temp_pv_comparacao_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 3: Temperaturas das 6 Camadas (COM) ──
    plt.figure(figsize=(12, 7))
    plt.plot(res_com['t'], res_com['T_v'], 'g-', linewidth=1.5, label='Vidro')
    plt.plot(res_com['t'], res_com['T_pv'], 'r-', linewidth=2, label='Células PV')
    plt.plot(res_com['t'], res_com['T_at'], 'b-', linewidth=1.5, label='Absorvedor Térmico')
    plt.plot(res_com['t'], res_com['T_tubo'], 'c-', linewidth=1.5, label='Tubulação')
    plt.plot(res_com['t'], res_com['T_iso'], 'm-', linewidth=1.5, label='Isolamento')
    plt.plot(res_com['t'], res_com['T_agua'], 'darkorange', linewidth=2.5, label='Água (saída)')
    plt.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='T_inicial (25°C)')
    
    # Adicionar linha vertical do desligamento da bomba
    if res_com['hora_desliga_bomba'] is not None:
        HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
        plt.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.2f}h)')
    
    plt.ylabel('Temperatura (°C)', fontsize=12)
    plt.xlabel('Horário (h)', fontsize=12)
    plt.title('Temperaturas das 6 Camadas - Com Trocador', fontsize=14)
    plt.legend(loc='upper left', fontsize=9, ncol=2)
    plt.grid(True, alpha=0.3)
    plt.xlim(res_com['t'][0], 22)
    plt.tight_layout()
    plt.savefig(f'3_temperaturas_com_trocador_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 4: Temperaturas das Camadas (SEM) ──
    plt.figure(figsize=(12, 7))
    plt.plot(res_sem['t'], res_sem['T_v'], 'g-', linewidth=1.5, label='Vidro')
    plt.plot(res_sem['t'], res_sem['T_pv'], 'r-', linewidth=2, label='Células PV')
    plt.plot(res_sem['t'], res_sem['T_at'], 'b-', linewidth=1.5, label='Absorvedor Térmico')
    plt.plot(res_sem['t'], res_sem['T_tubo'], 'c-', linewidth=1.5, label='Tubulação')
    plt.plot(res_sem['t'], res_sem['T_iso'], 'm-', linewidth=1.5, label='Isolamento')
    plt.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='T_inicial (25°C)')
    plt.ylabel('Temperatura (°C)', fontsize=12)
    plt.xlabel('Horário (h)', fontsize=12)
    plt.title('Temperaturas das Camadas - Sem Trocador', fontsize=14)
    plt.legend(loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(res_sem['t'][0], 22)
    plt.tight_layout()
    plt.savefig(f'4_temperaturas_sem_trocador_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 5: Eficiência Elétrica ──
    plt.figure(figsize=(10, 6))
    plt.plot(res_com['t'], res_com['eta'], 'b-', linewidth=2, label='Com Trocador')
    plt.plot(res_sem['t'], res_sem['eta'], 'r--', linewidth=2, label='Sem Trocador')
    
    # Adicionar linha vertical do desligamento da bomba
    if res_com['hora_desliga_bomba'] is not None:
        HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
        plt.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.2f}h)')
    
    plt.ylabel('Eficiência (%)', fontsize=12)
    plt.xlabel('Horário (h)', fontsize=12)
    plt.title('Eficiência Elétrica do Sistema PV/T', fontsize=14)
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(res_com['t'][0], 22)
    plt.tight_layout()
    plt.savefig(f'5_eficiencia_eletrica_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 6: Potência Elétrica ──
    plt.figure(figsize=(10, 6))
    plt.fill_between(res_com['t'], 0, res_com['P_el'], alpha=0.3, color='blue', 
                     label=f'Com Trocador: {energia_com_kWh:.4f} kWh')
    plt.fill_between(res_sem['t'], 0, res_sem['P_el'], alpha=0.3, color='red', 
                     label=f'Sem Trocador: {energia_sem_kWh:.4f} kWh')
    plt.plot(res_com['t'], res_com['P_el'], 'b-', linewidth=2)
    plt.plot(res_sem['t'], res_sem['P_el'], 'r--', linewidth=2)
    
    # Adicionar linha vertical do desligamento da bomba
    if res_com['hora_desliga_bomba'] is not None:
        HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
        plt.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                   label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.0f}h)')
    
    plt.ylabel('Potência (W)', fontsize=12)
    plt.xlabel('Horário (h)', fontsize=12)
    plt.title(f'Potência Elétrica Gerada (Ganho: {ganho:.1f}%)', fontsize=14)
    plt.legend(loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(res_com['t'][0], 22)
    plt.tight_layout()
    plt.savefig(f'6_potencia_eletrica_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # ── GRÁFICO 7: Tanque Estratificado ──
    if 'T_tk_0' in res_com:
        plt.figure(figsize=(10, 6))
        plt.plot(res_com['t'], res_com['T_tk_0'], 'r-', linewidth=2, label='Topo (0%)')
        plt.plot(res_com['t'], res_com[f'T_tk_{N_TK//4}'], 'orange', linewidth=1.5, label=f'25% (camada {N_TK//4})')
        plt.plot(res_com['t'], res_com[f'T_tk_{N_TK//2}'], 'g-', linewidth=1.5, label=f'50% (camada {N_TK//2})')
        plt.plot(res_com['t'], res_com[f'T_tk_{3*N_TK//4}'], 'b-', linewidth=1.5, label=f'75% (camada {3*N_TK//4})')
        plt.plot(res_com['t'], res_com[f'T_tk_{N_TK-1}'], 'purple', linewidth=2, label=f'Fundo (camada {N_TK-1})')
        plt.axhline(y=25, color='gray', linestyle='--', alpha=0.5, label='T_inicial (25°C)')
        
        # Adicionar linha vertical do desligamento da bomba
        if res_com['hora_desliga_bomba'] is not None:
            HORA_DESLIGA_BOMBA = res_com['hora_desliga_bomba']
            plt.axvline(x=HORA_DESLIGA_BOMBA, color='red', linestyle='--', linewidth=2, alpha=0.7, 
                       label=f'Desliga bomba ({HORA_DESLIGA_BOMBA:.2f}h)')
        
        plt.ylabel('Temperatura (°C)', fontsize=12)
        plt.xlabel('Horário (h)', fontsize=12)
        plt.title('Tanque de Armazenamento Estratificado', fontsize=14)
        plt.legend(loc='upper left', fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.xlim(res_com['t'][0], 22)
        plt.tight_layout()
        plt.savefig(f'7_tanque_estratificado_DB_{CIDADE_ALVO.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # ╔══════════════════════════════════════════════════════════════════════════╗
    # ║                         RELATÓRIO FINAL                                 ║
    # ╚══════════════════════════════════════════════════════════════════════════╝
    
    print(f"\n{'='*60}")
    print(f"RESULTADOS FINAIS - MODELO DUFFIE & BECKMAN")
    print(f"{'='*60}")
    
    print(f"\nCOM TROCADOR:")
    if res_com['hora_desliga_bomba'] is not None:
        print(f"  ⏱️ Bomba desligada às: {res_com['hora_desliga_bomba']:.2f}h")
    print(f"  T_v máx: {max(res_com['T_v']):.1f}°C")
    print(f"  T_pv máx: {max(res_com['T_pv']):.1f}°C")
    print(f"  T_at máx: {max(res_com['T_at']):.1f}°C")
    print(f"  T_tubo máx: {max(res_com['T_tubo']):.1f}°C")
    print(f"  T_iso máx: {max(res_com['T_iso']):.1f}°C")
    print(f"  T_agua máx: {max(res_com['T_agua']):.1f}°C")
    print(f"  Eficiência máx: {max(res_com['eta']):.2f}%")
    print(f"  Potência máx: {max(res_com['P_el']):.1f} W")
    print(f"  Energia: {energia_com_kWh:.4f} kWh")
    print(f"  Perdas laterais totais: {Q_bordas_total_com:.4f} MJ")
    print(f"  Perdas laterais máx: {max(res_com['Q_bordas']):.1f} W")
    
    print(f"\nSEM TROCADOR:")
    print(f"  T_v máx: {max(res_sem['T_v']):.1f}°C")
    print(f"  T_pv máx: {max(res_sem['T_pv']):.1f}°C")
    print(f"  T_at máx: {max(res_sem['T_at']):.1f}°C")
    print(f"  T_tubo máx: {max(res_sem['T_tubo']):.1f}°C")
    print(f"  T_iso máx: {max(res_sem['T_iso']):.1f}°C")
    print(f"  Eficiência máx: {max(res_sem['eta']):.2f}%")
    print(f"  Potência máx: {max(res_sem['P_el']):.1f} W")
    print(f"  Energia: {energia_sem_kWh:.4f} kWh")
    print(f"  Perdas laterais totais: {Q_bordas_total_sem:.4f} MJ")
    print(f"  Perdas laterais máx: {max(res_sem['Q_bordas']):.1f} W")
    
    print(f"\nANÁLISE COMPARATIVA:")
    print(f"  Ganho de energia: {ganho:.2f}%")
    print(f"  Redução T_pv: {max(res_sem['T_pv']) - max(res_com['T_pv']):.1f}°C")
    
    if 'T_tk_0' in res_com:
        print(f"\nTANQUE:")
        print(f"  T_topo final: {res_com['T_tk_0'][-1]:.1f}°C")
        print(f"  T_fundo final: {res_com['T_tk_bot'][-1]:.1f}°C")
    
    print(f"\nMODELO DUFFIE & BECKMAN - PRINCIPAIS CARACTERÍSTICAS:")
    print(f"  • Coeficiente de extinção KL: {KL} m⁻¹ (valor padrão D&B)")
    print(f"  • Correlação de radiação difusa: Erbs et al. (1982)")
    print(f"  • Método de múltiplas reflexões internas")
    print(f"  • Ângulo de incidência corrigido para Hemisfério Sul")
    
    print(f"\n{'='*60}")
    print(f"GRÁFICOS GERADOS:")
    print(f"  1_radiacao_temperatura_DB_*.png")
    print(f"  2_temp_pv_comparacao_DB_*.png")
    print(f"  3_temperaturas_com_trocador_DB_*.png")
    print(f"  4_temperaturas_sem_trocador_DB_*.png")
    print(f"  5_eficiencia_eletrica_DB_*.png")
    print(f"  6_potencia_eletrica_DB_*.png")
    print(f"  7_tanque_estratificado_DB_*.png")
    print(f"  8_propriedades_opticas_DB_*.png")
    print(f"{'='*60}")
    print(f"\n✅ SIMULAÇÃO COM MODELO DUFFIE & BECKMAN CONCLUÍDA COM SUCESSO!")