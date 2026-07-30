# PV/T Simulation - Brazilian Cities

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)]()

> Modelagem e simulação de um sistema fotovoltaico-térmico (PV/T) acoplado a um tanque de armazenamento estratificado para diferentes regiões do Brasil.

---

## 📋 Descrição

Este projeto implementa a modelagem matemática completa de um sistema **Fotovoltaico-Térmico (PV/T)** com:

| Recurso | Descrição |
|---------|-----------|
| ✅ Coletor PV/T | Multicamada com 6 nós térmicos |
| ✅ Trocador de calor | Escoamento laminar (Nu = 4,36) |
| ✅ Tanque | Estratificado com 12 camadas |
| ✅ Radiação solar | Modelos Duffie & Beckman, Liu & Jordan, Erbs |
| ✅ Controle | Termostato diferencial |
| ✅ Gráficos | Geração automática em alta resolução |
| ✅ Cidades | Suporte para 5 regiões brasileiras |

---

## 🚀 Instalação

```bash
# 1. Clone o repositório
git clone https://github.com/almeidamobius/PVT-Simulation-Brazil
cd PVT-Simulation-Brazil

# 2. Crie um ambiente virtual
# Windows
python -m venv venv
venv\Scripts\activate

# Linux/macOS
python3 -m venv venv
source venv/bin/activate

# 3. Instale as dependências
pip install -r requirements.txt
```
##  📊 Como Usar
Configuração
Edite os parâmetros no início do arquivo main.py:
```
# ─── DADOS DA SIMULAÇÃO ───
CIDADE_ALVO = "Recife"           # Santana, Recife, Brasília, São Paulo, Florianópolis
N_DIA = 321                      # Dia do ano (1-365) - 321 = 17/nov
V_VENTO = 2.24                   # Velocidade do vento [m/s]
USAR_TROCADOR = True             # True = com trocador, False = sem

# ─── DADOS CLIMÁTICOS DIÁRIOS ───
T_MAX = 32.0                     # Temperatura máxima diária [°C]
T_MIN = 24.0                     # Temperatura mínima diária [°C]
```
##  Execução
```
python main.py
```
##  📁 Estrutura do Projeto

```
PVT-Simulation-Brazil/
│
├── 📄 main.py                      # Código principal
├── 📄 requirements.txt             # Dependências
├── 📄 README.md                    # Este arquivo
├── 📄 LICENSE                      # Licença MIT
│
├── 📁 data/
│   └── 📄 irradiacao_global_horizontal.csv  # Dados INPE (2017)
│
├── 📁 resultados/                  # Gráficos gerados
│   └── 🖼️ *.png                    # 7 gráficos por cidade
│
└── 📁 docs/                        # Documentação adicional
```
### Gráficos Gerados (7 por cidade)

| # | Arquivo | Descrição |
|---|---------|-----------|
| 1 | `1_radiacao_temperatura_*.png` | Radiação solar e temperatura ambiente |
| 2 | `2_temp_pv_comparacao_*.png` | Temperatura PV (com/sem trocador) |
| 3 | `3_temperaturas_com_trocador_*.png` | 6 camadas (com trocador) |
| 4 | `4_temperaturas_sem_trocador_*.png` | 6 camadas (sem trocador) |
| 5 | `5_eficiencia_eletrica_*.png` | Eficiência elétrica |
| 6 | `6_potencia_eletrica_*.png` | Potência elétrica gerada |
| 7 | `7_tanque_estratificado_*.png` | Perfil do tanque |

### Relatório no Console

| Item | Descrição |
|------|-----------|
| 📊 Temperaturas máximas | Por camada do coletor (Vidro, PV, Absorvedor, Tubulação, Isolamento, Água) |
| ⚡ Energia gerada | Em kWh/dia para cada cidade |
| 📈 Ganho percentual | Com resfriamento ativo em relação ao cenário sem circulação |
| 🔥 Perdas laterais | Totais e máximas ao longo do dia |
| ⏱️ Horário de desligamento | Momento em que a bomba é desligada por termostato diferencial |

---

## 🛠️ Dependências

| Biblioteca | Versão | Uso |
|------------|--------|-----|
| numpy | ≥ 1.21.0 | Cálculos matemáticos e vetoriais |
| pandas | ≥ 1.3.0 | Leitura e processamento de dados |
| matplotlib | ≥ 3.5.0 | Geração de gráficos |

---

## 📚 Referências

| # | Referência |
|---|------------|
| 1 | Duffie, J. A.; Beckman, W. A. *Solar Engineering of Thermal Processes*. 4th ed. Wiley, 2013. |
| 2 | Pereira, E. B. et al. *Atlas Brasileiro de Energia Solar*. 2nd ed. INPE, 2017. |
| 3 | Kalogirou, S. A. *Solar Energy Engineering*. 3rd ed. Academic Press, 2024. |
| 4 | Rubio Ospina, L. M. *Modelagem e análise de um coletor fotovoltaico térmico*. UFPE, 2017. |
| 5 | Liu, B. Y. H.; Jordan, R. C. *Solar Energy*, 1960. |
| 6 | Erbs, D. G. et al. *Solar Energy*, 1982. |

---

## 📝 Licença

| Licença | Descrição |
|---------|-----------|
| MIT License | Veja o arquivo [LICENSE](LICENSE) para detalhes |

---

## 👥 Autores

| Nome | Função |
|------|--------|
| **Matheus Almeida de Menezes** | Desenvolvimento e modelagem |
| **Prof. Alisson Cocci de Souza** | Orientador |

---

## 🙏 Agradecimentos

| Instituição | Contribuição |
|-------------|--------------|
| **Universidade Federal Rural de Pernambuco (UFRPE)** | Apoio acadêmico e institucional |
| **Instituto Nacional de Pesquisas Espaciais (INPE)** | Dados meteorológicos |

