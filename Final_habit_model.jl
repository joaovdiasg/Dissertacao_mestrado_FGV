#######################################################################
########### Final version: Verdelhan's and Wachter's models ###########
#######################################################################

# A Model with external habits formation with time varying risk free rate

using Pkg, LinearAlgebra, Distributions, FastGaussQuadrature, BasicInterpolators, Printf, Expectations, Random, Plots, GLM, XLSX, Tables, DataFrames

cd("G:/Meu Drive/JOAO/Mestrado/Segundo ano/Dissertacao");

include("Final_functions11.jl")

setprecision(15)

# Escolha do método de solução de P/C (Ponto Fixo ou Séries)
method = 0  # 0 método por ponto fixo;
            # 1 método por séries de Wachter (2005)
            # 2 comparar os dois métodos

# Escolha a calibração
calib = 7   # 1 - Paper Campbell & Cochrane (1999)
            # 2 - Paper Wachter (2006)
            # 3 - Verdelhan (2010)
            # 4 - Verdelhan sem λ constante
            # 5 - João com B > 0
            # 6 - João com B < 0
            # 7 - João com B < 0 e λ constante 


grid  = 1   # 0 - grid scheindorfer 'padrão'
            # 1 - grid João
            # 2 - grid Wachter (2006)
            # 3 - grid Jivago

max_error = 1e-8
max_iter  = 1e+5

# Parâmetros Campbell & Cochrane (1999)
if calib == 1
    tsc = 12
    g = 0.0189 / tsc
    sig = 0.015 / sqrt(tsc)
    rf0 = 0.0094 / tsc
    phi = 0.87^(1/tsc)
    gamma = 2
    B = 0
    verd = 0
    ann = 0
    
# Parâmetros Wachter (2005)
elseif calib == 2
    tsc = 4
    g = 0.022 / tsc
    sig = 0.0086 / sqrt(tsc)
    rf0 = 0.0147 / tsc
    phi = (0.89)^(1/tsc)
    gamma = 2
    B = 0.011
    verd = 0
    ann = 0
    
# Parâmetros Verdelhan (WP)
elseif calib == 3
    tsc = 4
    g = 0.02115 / tsc
    sig = 0.0102 / sqrt(tsc)
    rf0 = 0.0136 / tsc
    phi = 0.97^(1/tsc)
    gamma = 2
    B = -0.01
    verd = 1
    ann = 1

# Parâmetros Verdelhan (λ variável)
elseif calib == 4
    tsc = 4
    g = 0.02115 / tsc
    sig = 0.0102 / sqrt(tsc)
    rf0 = 0.0136 / tsc
    phi = 0.97^(1/tsc)
    gamma = 2
    B = -0.000025
    verd = 0
    ann = 1

# Calibração João 1 (B>0)
elseif calib == 5
    tsc = 4
    g = 0.01903962 / tsc
    sig = 0.01797305 / sqrt(tsc)
    rf0 = 0.007010546 / tsc
    phi = (0.8421438)^(1/tsc)
    gamma = 2
    B = 0.011
    verd = 0
    ann = 0

# Calibração João 2 (B<0)
elseif calib == 6
    tsc = 4
    g = 0.01903962 / tsc
    sig = 0.01797305 / sqrt(tsc)
    rf0 = 0.007010546 / tsc
    phi = (0.8421438)^(1/tsc)
    gamma = 2
    B = -0.0001
    verd = 0
    ann = 0

# Calibração João 2 (B<0)
elseif calib == 7
    tsc = 4
    g = 0.01903962 / tsc
    sig = 0.01797305 / sqrt(tsc)
    rf0 = 0.007010546 / tsc
    phi = (0.8421438)^(1/tsc)
    gamma = 2
    B = -0.011
    verd = 1
    ann = 0
end    

## Parâmetros derivados do modelo
rho   = 0.15#-1:.1:1 # correlação países simulados
S_bar = sig*sqrt(gamma/(1-phi-B/gamma))
s_bar = log(S_bar)
s_max = s_bar + (1-S_bar^2)/2
S_max = exp(s_max)
delta  = exp(-rf0 + gamma*g -.5*(gamma*(1-phi)-B)) # Equação (12) do paper C&C1999. # pagina 15 da wachter
#(1-S_bar)/S_bar

if grid == 0 # scheindorfer
    S_grid      = LinRange(0,S_max,14)';
    S_grid      = S_grid[2:end];
    S_grid      = sort([S_grid;collect(0.09:0.001:0.093)]', dims =2)';
    logS_grid   = log.(S_grid);

elseif grid == 1 # Joao
    logS_grid   = mkgrids(215, 0);
    S_grid      = exp.(logS_grid);

elseif grid == 2 # Wachter
    logS_grid   = mkgrids(15, 1);
    S_grid      = exp.(logS_grid);

elseif grid == 3 # Jivago
    logS_grid   = mkgrids(15, 0);
    S_grid      = exp.(logS_grid);
end    


szgrid = length(S_grid) # Tamanho do grid
ncalc = 100000 # Número de simulações
bondsel = [1, 2, 3, 5, 10, 20] # Maturidades dos bonds a serem calculados
#bondsel = [1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 20]
maxcb = maximum(bondsel)
seedval = 123 # Número do gerador randômico

chk = 1  # Se 1 irá simular dados pós-guerra.
         # Se 0, irá simular dados a partir de 1926.1
flag1 = 0  # Simulará dados fictícios.
flag2 = 1  # 1 simulará dados anuais. 0 seria trimestral.
# ann = 0  # 1 anualiza dados de consumo e risk-free
           # da maneira do Verdelhan.
con = 0  # Usa na simulação dos países o rf vindo de simulação.
         # Se con = 0 usa rf interpolado na curva calibrada.

# Integração numérica
#a = abs(sig)*(-8) # limites
#b = abs(sig)*(8)
quad_points = 40 # quantos pontos serão usados

dist = Normal(0, sig);
E = expectation(dist, n = quad_points);

# Encontrando o valor de P/C e das perpetuidades pelo método de ponto fixo
if method == 0 || method == 2
    PC_ratio = findlpc(logS_grid)
    plot_PC = plot(S_grid, PC_ratio / tsc)  # plotando o P/C ratio anualizado (dividido por tsc)
    lnpca_pf = log.(PC_ratio)
    lnpca = lnpca_pf
    
# Encontrando o valor de P/C e das perpetuidades pelo método de séries
elseif method == 1 || method == 2
    PC_ratio_s = findlpcSeries(logS_grid) 
    plot_PC_s = plot(S_grid, PC_ratio_s / tsc)
    lnpca_s = log.(PC_ratio_s)
    lnpca = lnpca_s
    
# Comparando gráficos dos dois métodos
elseif method == 2
    plot_PC_2 = plot(S, PC_ratio_s / tsc, "r", S, PC_ratio / tsc, "g")
    legend(["Series method", "Fixed-point method"], loc=2)
    xlabel("Consumption surplus ratio (S{t})")
    ylabel("Price-consumption ratio (P{t}/C{t})")
    comp = maximum(abs.((W_PC_ratio - PC_ratio) ./ PC_ratio))
end

if method == 1
    PC_ratio = PC_ratio_s
end

#PC_verd_1 = PC_ratio_s
#PC_verd_0 = PC_ratio_s
#PC_joao_1 = PC_ratio_s
#PC_joao_0 = PC_ratio_s

#plot_PC_s

#savefig(plot_PC_s, "plot_PC_s.png")

verd = 0

# Pelo método de ponto fixo
if method == 0 || method == 2
    er_pf, elnr_pf, sdr_pf, sdlnr_pf, lnrf_pf, lnrf1_pf, lnycb_pf, elnrcb_pf, 
    sdlnrcb_pf, slpmv_pf = finders(logS_grid);

    # Pelo método de Séries
elseif method == 1 || method == 2
    er_s, elnr_s, sdr_s, sdlnr_s, lnrf_s, lnrf1_s, lnycb_s, elnrcb_s, 
    sdlnrcb_s, slpmv_s = finders(logS_grid);
end

#lny_verd_1 = lnycb_s
#lny_verd_0 = lnycb_s
#lny_joao_1 = lnycb_s
#lny_joao_0 = lnycb_s

# rates for each grid value

rate_plot = plot(S_grid, lnycb_s[:,1]*tsc, ylimits=(-0.05,0.35), label="3 months real")
            plot!(S_grid, lnycb_s[:,60]*tsc, ylimits=(-0.05,0.35), label="5 years real")
            
#savefig(rate_plot, "plot_rates_s.png")            

####### Simular dados do modelo

## Simulando séries temporais e anualizando-as
# (Caso não se simule através de dados reais, isto é, flag1 == 0)
# verd=0;

dc = 0 # Simulating consumption data

# Pelo método de ponto fixo
if method == 0 || method == 2
    Random.seed!(seedval); # Ajustando o seed randômico
    alndctsim_pf, astsim_pf, alnpctsim_pf, alnrtsim_pf, alnrfsim_pf, asdlnrtsim_pf, alnchpsim_pf, alnysim_pf,
    aelnrcbsim_pf, asdlnrcbsim_pf, atesterfsim_pf = annvars(dc, lnpca_pf,
    er_pf, elnr_pf, sdr_pf, sdlnr_pf, elnrcb_pf, sdlnrcb_pf, lnycb_pf, lnrf1_pf);

# Pelo método de Séries
elseif method == 1 || method == 2
    Random.seed!(seedval); 
    alndctsim_s, astsim_s, alnpctsim_s, alnrtsim_s, alnrfsim_s, asdlnrtsim_s, alnchpsim_s, alnysim_s,
    aelnrcbsim_s, asdlnrcbsim_s, atesterfsim_s = annvars(dc, lnpca_s,
    er_s, elnr_s, sdr_s, sdlnr_s, elnrcb_s, sdlnrcb_s, lnycb_s, lnrf1_s);
end

term_structure_plot = plot(mean(alnysim_s, dims=1)'.*100, label="Real yield")

#savefig(term_structure_plot, "plot_term_structure_s.png") 


if method == 0 || method == 2
   if ann == 1
    Edc_pf = tsc*mean(alndctsim_pf);
    Stdc_pf = sqrt(tsc)*std(alndctsim_pf);
   else
    Edc_pf = mean(alndctsim_pf); # Esperança da evolução do log do consumo
    Stdc_pf = std(alndctsim_pf); # Desvio da evolução do log do consumo
   end
   Erf_pf = mean(alnrfsim_pf); # Esperança do log da taxa livre de risco
   Stdrf_pf = std(alnrfsim_pf); # Desvio do log da taxa livre de risco
   Erfinterp_pf = mean(atesterfsim_pf);
   Stdrfinterp_pf = std(atesterfsim_pf);
   exrett_pf = alnrtsim_pf - alnrfsim_pf; # Excesso de retonos
   exrettinterp_pf = alnrtsim_pf - atesterfsim_pf;
   Shpr_pf = mean(exrett_pf)/std(exrett_pf); # Razão de Sharpe para log dos retornos
   ShpR_pf = mean(exp.(alnrtsim_pf)-exp.(alnrfsim_pf))/std(exp.(alnrtsim_pf)- exp.(alnrfsim_pf));
   Shprinterp_pf = mean(exrettinterp_pf)/std(exrettinterp_pf);
   ShpRinterp_pf = mean(exp.(alnrtsim_pf)- exp.(atesterfsim_pf))/std(exp.(alnrtsim_pf)-exp.(atesterfsim_pf));
   Eexrett_pf = mean(exrett_pf); # Média do excesso de retornos (em log)
   Stdexrett_pf = std(exrett_pf); # Desvio do excesso de retornos (em log)
   Eexrettinterp_pf = mean(exrettinterp_pf);
   Stdexrettinterp_pf = std(exrettinterp_pf);
   Ep_d_pf = mean(alnpctsim_pf); # Média do log da razão preço-consumo simulada
   Stdp_d_pf = std(alnpctsim_pf); # Desvio do log da razão preço-consumo simulada
   table = zeros(13,1);
   table[1, 1] = Edc_pf
   table[2, 1] = Stdc_pf
   table[3, 1] = Erf_pf
   table[4, 1] = Stdrf_pf
   table[5, 1] = Shpr_pf
   table[6, 1] = ShpR_pf
   table[7, 1] = Eexrett_pf
   table[8, 1] = Stdexrett_pf
   table[9, 1] = Ep_d_pf
   table[10, 1] = Stdp_d_pf
   table[11, 1] = S_max
   table[12, 1] = S_bar
   table[13, 1] = delta^tsc
elseif method == 1 || method == 2
   if ann == 1
    Edc_s = tsc*mean(alndctsim_s);
    Stdc_s = sqrt(tsc)*std(alndctsim_s);
   else
    Edc_s = mean(alndctsim_s); # Esperança da evolução do log do consumo
    Stdc_s = std(alndctsim_s); # Desvio da evolução do log do consumo
   end
   Erf_s = mean(alnrfsim_s); # Esperança do log da taxa livre de risco
   Stdrf_s = std(alnrfsim_s); # Desvio do log da taxa livre de risco
   Erfinterp_s = mean(atesterfsim_s);
   Stdrfinterp_s = std(atesterfsim_s);
   exrett_s = alnrtsim_s - alnrfsim_s; # Excesso de retonos
   exrettinterp_s = alnrtsim_s - atesterfsim_s;
   Shpr_s = mean(exrett_s)/std(exrett_s); # Razão de Sharpe para log dos retornos
   ShpR_s = mean(exp.(alnrtsim_s)-exp.(alnrfsim_s))/std(exp.(alnrtsim_s)-exp.(alnrfsim_s));
   Shprinterp_s = mean(exrettinterp_s)/std(exrettinterp_s);
   ShpRinterp_s = mean(exp.(alnrtsim_s)-exp.(atesterfsim_s))/std(exp.(alnrtsim_s)-exp.(atesterfsim_s));
   Eexrett_s = mean(exrett_s); # Média do excesso de retornos (em log)
   Stdexrett_s = std(exrett_s); # Desvio do excesso de retornos (em log)
   Eexrettinterp_s = mean(exrettinterp_s);
   Stdexrettinterp_s = std(exrettinterp_s);
   Ep_d_s = mean(alnpctsim_s); #     simulada
   Stdp_d_s = std(alnpctsim_s); # Desvio do log da razão preço-consumo simulada
   table = zeros(13,1);
   table[1,1] = Edc_s
   table[2,1] = Stdc_s
   table[3,1] = Erf_s
   table[4,1] = Stdrf_s
   table[5,1] = Shpr_s
   table[6,1] = ShpR_s
   table[7,1] = Eexrett_s
   table[8,1] = Stdexrett_s
   table[9,1] = Ep_d_s
   table[10,1] = Stdp_d_s
   table[11,1] = S_max
   table[12,1] = S_bar
   table[13,1] = delta^tsc
elseif method == 2
   table = zeros(13,2);
   table[1,1] = Edc_pf; table[1,2] = Edc_s
   table[2,1] = Stdc_pf; table[2,2] = Stdc_s
   table[3,1] = Erf_pf; table[3,2] = Erf_s
   table[4,1] = Stdrf_pf; table[4,2] = Stdrf_s
   table[5,1] = Shpr_pf; table[5,2] = Shpr_s
   table[6,1] = ShpR_pf; table[6,2] = ShpR_s
   table[7,1] = Eexrett_pf; table[7,2] = Eexrett_s
   table[8,1] = Stdexrett_pf; table[8,2] = Stdexrett_s
   table[9,1] = Ep_d_pf; table[9,2] = Ep_d_s
   table[10,1] = Stdp_d_pf; table[10,2] = Stdp_d_s
   table[11,1] = S_max; table[11,2] = S_max
   table[12,1] = S_bar; table[12,2] = S_bar
   table[13,1] = delta^tsc; table[13,2] = delta^tsc
end

bondsel2 = [1,3,5,10]
rates = zeros(length(bondsel2),2);
for (index, i) in enumerate(bondsel2*4)
    rates[index,1] = mean(alnysim_pf, dims=1)[i]
    rates[index,2] = std(alnysim_pf, dims=1)[i]
end

rates.*100

#XLSX.writetable("table112.xlsx", Tables.table(table))

# Simulação dos SDF's
# Simulando os dados novamente sob a mesma semente aleatória
if method == 0 || method == 2
    Random.seed!(seedval) # Ajustando o seed randômico
    
    stsim, vtsim, lndctsim, lnpctsim, lnrtsim, lnrfsim, ertsim, elnrtsim, sdrtsim, sdlnrtsim, elnrcbsim, sdlnrcbsim,
    lnysim, testerfsim = simvars(dc, lnpca_pf, er_pf, elnr_pf, sdr_pf, sdlnr_pf, elnrcb_pf,
    sdlnrcb_pf, lnycb_pf, lnrf1_pf);
    
    # Fator estocástico de desconto dos USA
    SDFus = delta .* exp.(-g .* gamma) .* exp.(-gamma .* vtsim) .* exp.(-gamma .* (stsim[2:end] .- stsim[1:end-1]))
    
    # Initialize SDFx matrix with zeros
    SDFx = zeros(length(stsim)-1, length(rho))
    stsimx = zeros(length(stsim), length(rho))
    vtsimx = zeros(length(stsim)-1, length(rho))
    lndctsimx = zeros(length(stsim)-1, length(rho))
    lnrfsimx = zeros(length(stsim), length(rho))
    
    for k in 1:length(rho)
        #stx, vtx, lndctx, lnrfx = simulacorr(rho[k], lnrf1_pf)
        stx, vtx, lndctx, lnrfx = simulacorr(rho[k])
        
        # delta.*exp(-gamma(Δc_t+1 + Δs_t+1))
        SDFx[:, k] = delta .* exp.(-g .* gamma) .* exp.(-gamma .* vtx) .* exp.(-gamma .* (stx[2:end] .- stx[1:end-1]))
        
        stsimx[:, k] = stx
        vtsimx[:, k] = vtx
        lndctsimx[:, k] = lndctx
        lnrfsimx[:, k] = lnrfx
    end

    # Encontrando as taxas de câmbio real
    # log(M^i) - log(M)
    deltaq = log.(SDFx) .- repeat(log.(SDFus), 1, size(SDFx, 2))

    # Regressões da UIP
    Betas_pf = zeros(length(rho), 3)
    Interv_pf = zeros(2, 1)  

    for k in 1:length(rho)
        dif = lnrfsim - lnrfsimx[:, k]
        #regressor = hcat(ones(size(diff, 1)), diff)
            
        data = DataFrame(Y = deltaq[:, k], X1 = dif[1:end-1])
        model = lm(@formula(Y ~ 1 + X1), data)
        Interv_pf = hcat(Interv_pf, confint(model))

        Betas_pf[k, 1] = rho[k]
        Betas_pf[k, 2] = GLM.coef(model)[1]
        Betas_pf[k, 3] = GLM.coef(model)[2]
    end 
    Interv_pf = Interv_pf[:, 2:end]

elseif method == 1 || method == 2
    Random.seed!(seedval) # Ajustando o seed randômico
    
    stsim, vtsim, lndctsim, lnpctsim, lnrtsim, lnrfsim, ertsim, elnrtsim, sdrtsim, sdlnrtsim, elnrcbsim, sdlnrcbsim,
    lnysim, testerfsim = simvars(dc, lnpca_s, er_s, elnr_s, sdr_s, sdlnr_s, elnrcb_s,
    sdlnrcb_s, lnycb_s, lnrf1_s);
    
    # Fator estocástico de desconto dos USA
    SDFus = delta .* exp.(-g .* gamma) .* exp.(-gamma .* vtsim) .* exp.(-gamma .* (stsim[2:end] .- stsim[1:end-1]))
    
    # Initialize SDFx matrix with zeros
    SDFx = zeros(length(stsim)-1, length(rho))
    stsimx = zeros(length(stsim), length(rho))
    vtsimx = zeros(length(stsim)-1, length(rho))
    lndctsimx = zeros(length(stsim)-1, length(rho))
    lnrfsimx = zeros(length(stsim), length(rho))
    lnysimx = zeros(length(stsim), length(rho))
    
    for k in 1:length(rho)
        #stx, vtx, lndctx, lnrfx = simulacorr(rho[k], lnrf1_pf)
        stx, vtx, lndctx, lnrfx = simulacorr(rho[k])
        
        # delta.*exp(-gamma(Δc_t+1 + Δs_t+1))
        SDFx[:, k] = delta .* exp.(-g .* gamma) .* exp.(-gamma .* vtx) .* exp.(-gamma .* (stx[2:end] .- stx[1:end-1]))
        
        stsimx[:, k] = stx
        vtsimx[:, k] = vtx
        lndctsimx[:, k] = lndctx
        lnrfsimx[:, k] = lnrfx
    end

    # Encontrando as taxas de câmbio real
    # log(M^i) - log(M)
    deltaq = log.(SDFx) .- repeat(log.(SDFus), 1, size(SDFx, 2))

    # Regressões da UIP
    Betas_s = zeros(length(rho), 3)
    Interv_s = zeros(2, 1)  

    for k in 1:length(rho)
        dif = lnrfsim - lnrfsimx[:, k]
        #regressor = hcat(ones(size(diff, 1)), diff)
            
        data = DataFrame(Y = deltaq[:, k], X1 = dif[2:end])
        model = lm(@formula(Y ~ 1 + X1), data)
        Interv_s = hcat(Interv_s, confint(model))

        Betas_s[k, 1] = rho[k]
        Betas_s[k, 2] = GLM.coef(model)[1]
        Betas_s[k, 3] = GLM.coef(model)[2]
    end
    Interv_s = Interv_s[:, 2:end]

elseif method == 2
    Betas = hcat(Betas_pf, fill(NaN, size(Betas_pf, 1), 1))
    Betas = hcat(Betas, Betas_s)
    Intervalo = hcat(Interv_pf, fill(NaN, size(Interv_pf, 1), 1))
    Intervalo = hcat(Intervalo, Interv_s)
end

Betas_pf'
##############################################################################################################################################
##############################################################################################################################################

if verd == 1
    lambda_st = (1-S_bar)/S_bar

    # Define the expressions
    expected_M_t_plus_1 = delta .* exp.(gamma*(-g + 0.5*gamma*sig^2*(lambda_st + 1)^2 .+ (phi - 1).*(s_bar .- logS_grid)))
    var_M_t_plus_1 = delta^2 * (exp(gamma^2*sig^2*(lambda_st + 1)^2) - 1) * exp.(gamma*(-2*g + gamma*sig^2*(lambda_st + 1)^2 .+ 2*(phi - 1).*(s_bar .- logS_grid)))
    std_dev_M_t_plus_1 = sqrt.(var_M_t_plus_1)

    # Calculate the ratio
    ratio_std_dev_to_exp_value = std_dev_M_t_plus_1 ./ expected_M_t_plus_1
    mean(ratio_std_dev_to_exp_value)
    std(ratio_std_dev_to_exp_value)

    verd = 0

    # Letting λ(s_t) varies 
    expected_M_t_plus_2 = delta .* exp.(gamma*(-g .+ 0.5*gamma*sig^2 .*(lambda.(logS_grid) .+ 1).^2 .+ (phi - 1).*(s_bar .- logS_grid)))
    var_M_t_plus_2 = delta^2 .* (exp.(gamma^2*sig^2 .*(lambda.(logS_grid) .+ 1).^2) .- 1) .* exp.(gamma*(-2*g .+ gamma*sig^2 .*(lambda.(logS_grid) .+ 1).^2 .+ 2*(phi - 1).*(s_bar .- logS_grid)))
    std_dev_M_t_plus_2 = sqrt.(var_M_t_plus_1)

    ratio_std_dev_to_exp_value_2 = std_dev_M_t_plus_2 ./ expected_M_t_plus_2
    mean(ratio_std_dev_to_exp_value_2)
    std(ratio_std_dev_to_exp_value_2)
end

Real_tips =[2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20;
0.443954722	0.637316103	0.791997997	0.922977301	1.036641344	1.136490939	1.224616279	1.302659078	1.371530524	1.432622381	1.485182762	1.533577189	1.57478153	1.610108187	1.641383529	1.66704969	1.68859113	1.706136026	1.719593705]

plot(Real_tips[1,:], Real_tips[2,:], title="American real yield curve",label = "TIPS yield", lw=2, lc=:black)
xlabel!("years")
ylabel!("%")
#savefig("TIPS.png")


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# Teste da hipótese das expectativas
bn = zeros(length(bondsel)-1, 2)
bnconf = zeros(length(bondsel)-1, 4)

for i = 2:length(bondsel)
    chgyield = lnysim[2:end, i-1] .- lnysim[1:end-1, i]
    spreadyield = (lnysim[1:end-1, i] .- lnysim[1:end-1, 1]) ./ (bondsel[i] - 1)

    # You need to replace the regress function with the appropriate Julia code
    # For regression in Julia, you can use the GLM.jl package or similar
    # Example using GLM.jl for a linear regression:
    # model = lm(@formula(YieldChange ~ Spread), data)
    # aux1 = coef(model)
    # aux2 = confint(model)

    bn[i-1, :] = aux1'
    bnconf[i-1, :] = reshape(aux2, 1, 4)
end

# Conditional check for method == 2
if method == 2
    Betas = hcat(Betas_pf, fill(NaN, size(Betas_pf, 1), 1))
    Betas = hcat(Betas, Betas_s)
    Intervalo = hcat(Interv_pf, fill(NaN, size(Interv_pf, 1), 1))
    Intervalo = hcat(Intervalo, Interv_s)
end

# Teste da hipótese das expectativas
bn = zeros(length(bondsel)-1, 2)
bnconf = zeros(length(bondsel)-1, 4)

for i = 2:length(bondsel)
    chgyield = lnysim[2:end, i-1] .- lnysim[1:end-1, i]
    spreadyield = (lnysim[1:end-1, i] .- lnysim[1:end-1, 1]) ./ (bondsel[i] - 1)

    # You need to replace the regress function with the appropriate Julia code
    # For regression in Julia, you can use the GLM.jl package or similar
    # Example using GLM.jl for a linear regression:
    # model = lm(@formula(YieldChange ~ Spread), data)
    # aux1 = coef(model)
    # aux2 = confint(model)

    bn[i-1, :] = aux1'
    bnconf[i-1, :] = reshape(aux2, 1, 4)
end

#means = Float64[]  # Create an empty vector to store the means

#for i in bondsel
    # Calculate the mean of the i-th column and multiply by 400, then append to the vector
#    push!(means, mean(lnysim[:, i]) * 400)
#end
