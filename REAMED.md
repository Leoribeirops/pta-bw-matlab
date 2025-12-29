# Pan–Tompkins QRS Detection (Bandpass BW – MATLAB)

Este repositório contém uma implementação em **MATLAB** do algoritmo de detecção de complexos **QRS** baseado no método de **Pan–Tompkins**, utilizando um **filtro Bandpass (BW)** e processamento clássico em etapas:

1. Bandpass filtering  
2. Differenciação  
3. Elevação ao quadrado  
4. Moving Window Integration (MWI)  
5. Thresholding adaptativo  

O projeto utiliza sinais do **MIT-BIH Arrhythmia Database** para validação e comparação com a implementação original do Pan–Tompkins.

---

## 📂 Estrutura do Projeto

```text
pta-bw-matlab/
│
├── data/
│   └── ECG_MIT_01.mat         # Exemplo de sinal ECG (MIT-BIH)
│
├── Thresholding_BP_MWI.m      # Thresholding adaptativo (entrada: BW + MWI)
├── plot_signal_window.m       # Função genérica de plotagem
├── pan_tompkins_original.m   # Implementação original (referência)
│
├── main.m                     # Script principal de execução
└── README.md
