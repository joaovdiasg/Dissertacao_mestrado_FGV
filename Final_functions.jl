#######################################################################
########### Final version: Verdelhan's and Wachter's models ###########
#######################################################################

### functions file

function mkgrids(szgrid,flag)
    # Construirá o grid de s de uma maneira eficiente. 
    # global s_max S_max s_bar

    if flag == 0
        sg = range(start = 0, stop = S_max, length = szgrid)
        aux = [sg[szgrid]-.01, sg[szgrid]-.02, sg[szgrid]-.03, sg[szgrid]-.04]
        sg = [sg; aux]
        sg = sort(sg)
        sg = log.(sg[2:length(sg)])
        if maximum(sg .== s_bar) == 0 # Certificando-se que s_bar vai estar no grid.
            sg = [sg; s_bar]
            sg = sort(sg)
        end
        if maximum(sg .== s_max) == 0
            sg = [sg; s_max]
            sg = sort(sg)
        end
        # Colocar mais densidade no início do grid para melhorar a 
        # iteração durante o procedimento de ponto fixo. 
        idens=log.([.0005, 0.0015, .0025, .0035, .0045])
        sg=[sg; idens]
        sg = sort(sg) # Ordenando os valores presentes no grid
    end

    # Grid 3 da Wachter (2005)
    if flag == 1
        sg = range(start = 0, stop = S_max, length = szgrid)
        sg = log.(sg[2:length(sg)])
        if maximum(sg .== s_bar) == 0 # Certificando-se que s_bar vai estar no grid.
            sg = [sg; s_bar]
            sg = sort(sg)
        end
        if maximum(sg .== s_max) == 0
            sg = [sg; s_max]
            sg = sort(sg)
        end
        u = minimum(sg)
        aux = range(start = -300, stop = u, length = 200)
        sg = [sg; aux]
        sg = sort(sg)
    end
    return sg
end 

function strans(s,v)
    # Esta função retornará o valor de s_{t+1} = log(S_{t+1}) 
    # s_{t+1} = (1-phi)*s_bar + phi*s_t + lambda(s_t)*v_{t+1}; 

    news = (1-phi)*s_bar .+ phi.*s .+ lambda.(s).*v';
    return news
end

function lambda(s)
    # Retorna o valor da função sensibilidade λ(s_t)

    if verd == 0
        if s <= s_max
            y = (1 / S_bar)*sqrt(maximum([0; 1-2*(s-s_bar)]))-1;
        else
            y = 0;
        end
    elseif verd == 1
        y=(1-S_bar)/S_bar;
    end
    return y
end

function mrsinsd(s, v)
    # Retorna a taxa marginal de substituição intertemporal no modelo (M_{t+1}).
    # delta*[(S_{t+1}/S_t)*(C_{t+1}/C_t)]^(-gamma)
    # delta*exp{-gamma[g-(1-phi)(s_t-s_bar) + (1+ λ(s_t))*v_{t+1}]}

    out = delta*exp.(-gamma*(g.-(1-phi).*(s.-s_bar) .+ (1 .+ lambda.(s)).*v'))
    return out
end

function recursionPC(v)
    # função usada para no processo iterativo do ponto fixo de P/C
    s1          = strans(logS_grid, v);
    M           = mrsinsd(logS_grid, v); # [nS,quad_points] M_{t+1}
    delta_c     = exp.(g .+ v); # E_t[C_{t+1}/C_t] = exp(g + (σ^2)/2)
    Fit         = LinearInterpolator(vec(logS_grid), vec(log.(PC)), NoBoundaries());
    PC_new      = (M.*delta_c'.*(1 .+exp.(Fit.(s1))))#*pdf(Normal(0,sig), v);
    return PC_new
end

function findlpc(logS_grid)
    # Este é o procedimento que irá calcular o ponto fixo P/C

   global PC, delta_c

    PC       = ones(length(logS_grid));
    d_PC     = 10000 # só pra inicializar o loop

    iter        = 0;
    for iter in 1:max_iter
        #global Fit
        #Fit    = linear_interpolation(vec(logS_grid), vec(log.(PC)), extrapolation_bc = Line());
        PC_new = E(recursionPC);
        d_PC   = PC_new-PC;
        PC     = PC_new;
        iter   = iter + 1;
        println(maximum(abs.(d_PC)))
        
        if maximum(abs.(d_PC)) < max_error
            break
        end
    end
    return PC
end

function seriesPC(v)
    # função usada para no processo séries da Wachter para encontrar P/C

    s1          = strans(logS_grid, v);
    M           = mrsinsd(logS_grid, v); # [nS,quad_points] M_{t+1}
    delta_c     = exp.(g .+ v); # E_t[C_{t+1}/C_t] = exp(g + (σ^2)/2)
    Fit         = LinearInterpolator(vec(logS_grid), vec(log.(F)), NoBoundaries());
    Fit_        = s -> max.(0,exp.(Fit.(s)));
    F_new       = (M.*delta_c'.*Fit_(s1))#*pdf(Normal(0,sig), v);
    return F_new
end

function findlpcSeries(logS_grid)
    # Este é o procedimento que irá calcular o ponto fixo P/C através do método de séries da Wachter.
    
    global F

    PC_s     = ones(szgrid,1)#Interpolations.deduplicate_knots!(ones(nS,1));
    F        = ones(szgrid,1);

    iter        = 0;    
    for iter=1:max_iter
        F    = E(seriesPC);
        PC_s = PC_s + F;
        iter   = iter + 1;
        println(maximum(abs.(F./PC_s)))

        if maximum(abs.(F./PC_s)) < max_error
            break
        end
    
    end
    return PC_s
end

function bondsP(v)
    s1          = strans(logS_grid, v);
    M           = mrsinsd(logS_grid, v);
    Fit_        = s -> max.(0,Fit.(s));
    P_new       = (M.*Fit_.(s1))
    return P_new
end

function finders(sg)
    # Procedimento que calcula os retornos esperados dos ativos de consumo. ele
    # nos fornecerá E(R), SD(R), lnrf, a curvatura da fronteira média variância
    # dada pela variável (slpmv) e integrará a razão de Sharpe do vetor de ativos de consumo. 

    global lnpcb
    
    # Inclinação da Fronteira de Média-Variância 
    slpmv = (exp.((gamma*sig)^2 .*(1 .+lambda.(sg)).^2).-1).^(0.5);
    
    ## Estrutura termo da taxa de juros dada por: 
    # lnrf = -ln(delta) + gamma*g - 
    # gamma*(1-phi)*(s{t}-s_bar)-.5((gamma*sig)^2)*(1+lambda(s{t}))^2
    # - B*(sg - s_bar) --> Variável no estado
    
    # eq (6) in Wachter 2006
    lnrf = -log(delta) + gamma*g .+ gamma*(1-phi).*(s_bar.-sg) .- 0.5.*(gamma*sig.*(1 .+ lambda.(sg))).^2;
    
    lnpcb = ones(szgrid, maxcb*tsc+1); # price of real bonds

    for j = 2:(maxcb*tsc+1)
        global Fit
        # real bonds
        Fit        = LinearInterpolator(vec(logS_grid), vec(lnpcb[:,j-1]), NoBoundaries());
        lnpcb[:,j] = E(bondsP);
    end  

    lnycb = zeros(szgrid, maxcb*tsc+1); # yields of real bonds
    for j = 1:(maxcb*tsc+1)
        lnycb[:,j] = -(1/(j-1))*log.(lnpcb[:,j])#*tsc # multiplying by tsc makes it annual 
    end  
    
    lnycb = lnycb[:, 2:end]

    #lnycb[:,1] = lnrf;

    ## Retornos e Desvios Padrão Esperados 
    
    # Comsumption claim:
    lnrf1 = - log.(E(intemrs)) # E_t[M_{t+1}]
    er = E(erinsd) # E_t[R_{t+1}]
    elnr = E(x -> log.(erinsd(x))) # E_t[r_{t+1}]
    sdr = E(x -> erinsd(x).^2) 
    sdr = (sdr .- er.^2).^(.5) # var_t[R_{t+1}]
    sdlnr = E(x -> log.(erinsd(x)).^2)
    sdlnr = (sdlnr .- elnr.^2).^(.5) # var_t[r_{t+1}]
    
    # Zero-coupon bonds:
    elnrcb = zeros(length(sg),maxcb*tsc);  
    sdlnrcb = zeros(length(sg),maxcb*tsc);

    elnrcb[:,1] = lnrf
        
    for k = 2:maxcb*tsc
        p_aux = (lnpcb[:, (k-1)]);
        p_auxt = (lnpcb[:, k]);
        global p_aux, p_auxt
        elnrcb[:,k] = E(x -> log.(ercbin(x)));
        sdlnrcb[:,k] = E(x -> log.(ercbin(x)).^2);
        for i = 1:length(logS_grid)
            if sdlnrcb[i,k] .- elnrcb[i,k].^2 > 0 
                sdlnrcb[i,k] = (sdlnrcb[i,k] .- elnrcb[i,k].^2).^(.5);
            else
                sdlnrcb[i,k] = 0
            end
        end
    end
    return er, elnr, sdr, sdlnr, lnrf, lnrf1, lnycb, elnrcb, sdlnrcb, slpmv
end


function intemrs(v)
    # Retornará a taxa marginal de substituição de tal forma que possa ser 
    # utilizada na integração do ponto fixo

    #out = pdf.(Normal(0,sig), v).*mrsinsd(s, v);
    out = mrsinsd(logS_grid, v);
    return out 
end

function erinsd(v)
    # Procedimento para calcular os retornos do consumption claim

    Fit = LinearInterpolator(vec(logS_grid), vec(log.(PC_ratio)), NoBoundaries());
    er = zeros(length(logS_grid), length(v));
    for i = 1:length(v)
        er[:, i] = ((1 .+ exp.(Fit.(strans(logS_grid,v[i]))))./exp.(Fit.(logS_grid))).*exp.(g .+ v);
    #er = ((1 .+exp.(Fit.(strans(logS_grid,v))))./exp.(Fit.(logS_grid))).*exp.(g.+v);
    end
    return er
end

function ercbin(v)
    # Retornos esperados da estrutura a termo 
    #global s sg lnpcb matur

    Fit = LinearInterpolator(vec(logS_grid), vec(log.(p_aux)), NoBoundaries());
    #Fit2 = linear_interpolation(vec(logS_grid), vec((p_aux[:,end])), extrapolation_bc = Line());
    out = zeros(length(logS_grid), length(v));
    for i = 1:length(v)
        #out[:, i] = exp.(Fit.(strans(logS_grid, v[i])))./exp.(p_auxt);
        out[:, i] = exp.(Fit.(strans(logS_grid, v)))./p_auxt;
    end
    return out
end
# com exp() ou sem exp()?
# acho q sem

function annvars(dc, lnpca, er, elnr, sdr, sdlnr, elnrcb, sdlnrcb, lny, lnrf1)
    # This program works with artificial time series data generated outside real data
    # It annualizes expected returns and still needs to annualize simulated bond returns.
#    global tsc bondsel ann
    
    # Running artificial time series
    stsim, vtsim, lndctsim, lnpctsim, lnrtsim, lnrfsim, ertsim, elnrtsim, sdrtsim, sdlnrtsim, elnrcbsim, sdlnrcbsim, 
    lnysim, testerf = simvars(dc, lnpca, er, elnr, sdr, sdlnr, elnrcb, sdlnrcb, lny, lnrf1)
    
    T = length(stsim)
    
    # Consumption
    if ann == 1
        alndctsim = lndctsim
    else
        alndctsim = chgfreq(lndctsim, true, tsc, 1)
    end
    
    # State variable of the model
    if T > 1
        astsim = chgfreq(stsim, false, tsc, 1)
        #astsim = astsim[2:end]
    end
    
    # P/C ratio
    if length(lnpctsim) > 1
        alnpctsim = chgfreq(lnpctsim, false, tsc, 1) .- log(tsc)
        #alnpctsim = alnpctsim[2:end]
    end
    
    # Annualize the returns
    if length(lnrtsim) > 1
        alnrtsim = chgfreq(lnrtsim, true, tsc, 1)
        #alnrtsim = alnrtsim[2:end]
    end
    
    # Simulated risk-free return
    if length(lnrfsim) > 1
        alnrfsim = chgfreq(lnrfsim, true, tsc, 1)
        #alnrfsim = alnrfsim[2:end]
    end


    # Interpolated risk-free return
    if length(testerf) > 1
        atesterf = chgfreq(testerf, true, tsc, 1)
        #atesterf = atesterf[2:end]
    end
    
    # Conditional return deviations
    if length(sdlnrtsim) > 1
        asdlnrtsim = chgfreq(sdlnrtsim, true, tsc, 1)
        #asdlnrtsim = asdlnrtsim[2:end]
    end
    
    # Price evolution
    if length(lnpctsim) > 1
        lnchpsim = lnpctsim[2:T] - lnpctsim[1:T-1] + lndctsim
        alnchpsim = chgfreq(lnchpsim, true, tsc, 1)
        #alnchpsim = alnchpsim[2:end]
    end
    
    # Yields
    if length(lnysim) > 1
        #alnysim = similar(lnysim)
        alnysim = zeros(Int(round(T/tsc)), maxcb*tsc)
        for i = 1:maxcb*tsc
            alnysim[:, i] = chgfreq(lnysim[:, i], true, tsc, 1)
        end
        #alnysim = alnysim[2:end, :]
    end
        
    # Bonds - Average Return
    if length(elnrcbsim) > 1
        #aelnrcbsim = similar(elnrcbsim)
        aelnrcbsim = zeros(Int(round(T/tsc)), maxcb*tsc)
        for i = 1:maxcb*tsc
            aelnrcbsim[:, i] = chgfreq(elnrcbsim[:, i], true, tsc, 1)
        end
        #aelnrcbsim = aelnrcbsim[2:end, :]
    end

    # Bonds - Standard Deviation
    if length(sdlnrcbsim) > 1
        #asdlnrcbsim = similar(sdlnrcbsim)
        asdlnrcbsim = zeros(Int(round(T/tsc)), maxcb*tsc)
        for i = 1:maxcb*tsc
            asdlnrcbsim[:, i] = chgfreq(sdlnrcbsim[:, i], true, tsc, 1)
        end
        #asdlnrcbsim = asdlnrcbsim[2:end, :]
    end

    return alndctsim, astsim, alnpctsim, alnrtsim, alnrfsim, asdlnrtsim, alnchpsim, 
           alnysim, aelnrcbsim, asdlnrcbsim, atesterf
end


function simvars(dc, lnpca, er, elnr, sdr, sdlnr, elnrcb, sdlnrcb, lny, lnrf1)
    
    # Esta rotina simulará séries temporais das principais grandezas deste 
    # modelo e terá a opção de utilizar dados reais ou fictícios. 

    # Simulará: 
    # - processo de s=log(S); 
    # - P/C; 
    # - R{t+1}; 
    # - E{t}[R{t+1}]; 
    # - SD{t}[R{t+1}]; 
    # - USA Rf{t+1} calibrado; 
    # - Rf{t+1} com consumo correlacionado com os USA 
    # - Bonds 
    
    # Todas as séries simuladas serão sufixadas por 'tsim' e terão a mesma 
    # nomeação herdada dos programas anteriores (findlpc e finders). 

    # Se o input de simvars é nulo, que dizer que este programa irá simular o 
    # processo de consumo do modelo. Caso o usuário queira inserir os dados 
    # reais da economia americana, as outras variáveis serão simuladas através
    # destes dados de consumo. 
    
    if dc == 0  # Simulating consumption data
        T = ncalc
        vtsim = sig * randn(T)
        lndctsim = g .+ vtsim
    else  # Using real consumption data
        if minimum(dc) <= 0
            println("simvars: You entered the log of consumption growth.")
            println("I need you to enter the data for gross consumption growth,")
            println("i.e., neither log nor net growth.")
        end
        T = length(dc)
        lndctsim = log.(dc)
        vtsim = lndctsim .- g
    end
    
    stsim = zeros(T + 1)
    stsim[1] = s_bar
    
    # Simultating s_{t+n} for n = 1, 2, ... T+1
    for i = 2:T + 1
        if (1-phi)*s_bar .+ phi.*stsim[i - 1] .+ lambda.(stsim[i - 1]).*vtsim[i-1] <= s_max
            #strans(stsim[i - 1], vtsim[i - 1]) <= s_max
            #stsim[i] = strans(stsim[i - 1], vtsim[i - 1])
            stsim[i] = (1-phi)*s_bar .+ phi.*stsim[i - 1] .+ lambda.(stsim[i - 1]).*vtsim[i-1]
        else
            stsim[i] = (1 - phi) * s_bar + phi * stsim[i - 1]
        end
    end
    
    # Simulating P/C
    Fit = LinearInterpolator(vec(logS_grid), vec(lnpca), NoBoundaries());
    lnpctsim = Fit.(stsim)
    
    # Retornos ex-post
    # R = (C'/C){(1+(P/C)')/(P/C)} %
    lnrtsim = lndctsim .+ log.(1 .+ exp.(lnpctsim[2:T+1])) .- lnpctsim[1:T]
    
    # Log de Rf variante no tempo
    # lnrfsim = -log(delta) + gamma*g - gamma*(1-phi)*(stsim-s_bar)...
    # - 0.5*(gamma*sig*(1+lambda(stsim))).^2;
    lnrfsim = -log(delta) .+ gamma .* g .- gamma .* (1 - phi) .* (stsim .- s_bar) .- 0.5 .* (gamma^2 .* sig^2) .* (1 .+ lambda.(stsim)).^2
    
    # Retornos esperados e desvio condicional
    Fit1 = LinearInterpolator(vec(logS_grid), vec(lnrf1), NoBoundaries());
    testerf = Fit1.(stsim)

    Fit2 = LinearInterpolator(vec(logS_grid), vec(er), NoBoundaries());
    ertsim = Fit2.(stsim)

    Fit3 = LinearInterpolator(vec(logS_grid), vec(elnr), NoBoundaries());
    elnrtsim = Fit3.(stsim)

    Fit4 = LinearInterpolator(vec(logS_grid), vec(sdr), NoBoundaries());
    sdrtsim = Fit4.(stsim)

    Fit5 = LinearInterpolator(vec(logS_grid), vec(sdlnr), NoBoundaries());
    sdlnrtsim = Fit5.(stsim)
    
    # Retornos esperados das T-Bill 90
    Fit6 = LinearInterpolator(vec(logS_grid), vec(elnrcb[:, 1]), NoBoundaries());
    elnrcbsim = Fit6.(stsim) # Retornos esperados do bonds reais
    
    Fit7 = LinearInterpolator(vec(logS_grid), vec(sdlnrcb[:, 1]), NoBoundaries());
    sdlnrcbsim = Fit7.(stsim)

    Fit10 = LinearInterpolator(vec(logS_grid), vec(lny[:, 1]), NoBoundaries());
    lnysim = Fit10.(stsim) #  Real bond yields
    lny2sim = zeros(length(lnysim));

    #for i in bondsel*tsc
    for i = 2:maxcb*tsc

        Fit12 = LinearInterpolator(vec(logS_grid), vec(lny[:, i]), NoBoundaries());
        lnysim = hcat(lnysim, Fit12.(stsim))
        Fit13 = LinearInterpolator(vec(logS_grid), vec(lny[:, i-1]), NoBoundaries());
        lny2sim = hcat(lny2sim, Fit13.(stsim))
        Fit14 = LinearInterpolator(vec(logS_grid), vec(elnrcb[:, i]), NoBoundaries());
        elnrcbsim = hcat(elnrcbsim, Fit14.(stsim))
        Fit15 = LinearInterpolator(vec(logS_grid), vec(sdlnrcb[:, i]), NoBoundaries());
        sdlnrcbsim = hcat(sdlnrcbsim, Fit15.(stsim))
    end
        
    # Retornos das maturidades ajustadas em bondsel = [1 2 4 8 12 16 20]
    #lnrcbsim = vcat(0,lnysim[1:T-1,1]);
    #lnrcbNsim = vcat(0,lnyNsim[1:T-1,1]);
    #for i = 2:length(bondsel)
    #    lnrcbsim = hcat(lnrcbsim, vcat(0,(-lny2sim[2:T,i]*((bondsel[i]-1/tsc)) + lnysim[1:T-1,i]*bondsel[i])/tsc));
    #    lnrcbNsim = hcat(lnrcbNsim, vcat(0,(-lny2Nsim[2:T,i]*((bondsel[i]-1/tsc)) + lnyNsim[1:T-1,i]*bondsel[i])/tsc));
    #end
    return stsim, vtsim, lndctsim, lnpctsim, lnrtsim, lnrfsim, ertsim, elnrtsim, sdrtsim, sdlnrtsim, 
           elnrcbsim, sdlnrcbsim, lnysim, testerf #, lnrcbsim, lnrcbNsim
end

function chgfreq(series::Vector{Float64}, sum_logs::Bool, original_freq::Int, target_freq::Int)
    # Check if the length of the series is appropriate
    if length(series) % original_freq != 0
        
    num_to_remove = length(series) % tsc

    # Trim the vector
    series = series[1+num_to_remove:end]
    end

    # Reshape the series to have `original_freq` periods per row
    reshaped_series = reshape(series, (original_freq, :))

    # Process the series based on whether the logs can be summed
    if sum_logs
        # If logs can be summed, simply sum over the periods
        converted_series = sum.(eachcol(reshaped_series))
    else
        # If logs cannot be summed, exponentiate, average, and then log the mean
        converted_series = log.(mean.(eachcol(exp.(reshaped_series))))
    end

    return converted_series
end


function simulacorr(rho)
    T = ncalc
    Random.seed!(seedval)
    x = sig * randn(T)
    y = sig * randn(T)
    vtsim = rho .* x .+ sqrt(1 - rho^2) .* y  # Controlling the correlation of the consumption processes
    
    lndctsim = g .+ vtsim

    # Simulate the log state variable S
    stsim = zeros(T + 1)
    stsim[1] = s_bar  # The economy starts at its stationary state.

    for i = 2:T + 1
        stsim[i] = strans(stsim[i - 1], vtsim[i - 1])
    end

    # Log of time-varying Rf
    lnrfsim = -log(delta) .+ gamma .* g .- 0.5.*(gamma * (1 - phi) - B) .- B .* (stsim .- s_bar)

    return stsim, vtsim, lndctsim, lnrfsim
end



#e, w = gausslegendre(quad_points)  ## using 40 quadrature points as watcher
    
#fx = w' * sig_c *[pdf.(dist, v)'; pdf.(dist, v)'] * f.((b-a)/2 * e .+ (a+b)/2)*(b-a)/2



############################
##### não estou usando #####
############################

function EE(f)
    
    rawNodes, rawWeights = gausslegendre(quad_points)
    # Transform nodes to proper interval.
    nodes = map(x -> (0.5(b - a)) * x + (a + b) / 2, rawNodes)
    # Add pdf to weights.
    compoundWeights = [rawWeights[i] * sig * pdf(dist, nodes[i]) for i in 1:quad_points]
    # Add scale factor to weights.
    weights = (b - a) / 2 * compoundWeights

    E = f(nodes[1]) * weights[1]
    @inbounds for i in 2:length(nodes)
        n, w = nodes[i], weights[i]
        E += f(n) * w
    end
    return E
end

function intgl1(f, a, b)
    # Rotina que desenvolve o método de integração de Gauss-Legendre.
    # Nesta rotina, o vetor 'e' compõe-se de valores entre 0 e 1 (e ~ (0,1)) e o 
    # vetor 'w' das raízes do polinômio de Legendre cujos coeficientes vem de 'e'
    
    e, w = gausslegendre(quad_points)  ## using 40 quadrature points as watcher
    
    fx = w' * f.((b-a)/2 * e .+ (a+b)/2)*(b-a)/2
    # função do luan
    
    return fx
end


function inter(v)
    # Fornece o valor do retorno esperado de acordo com uma distribuição normal

    #out = pdf.(Normal(0,sig), v).*erinsd(s, v);
    out = erinsd(logS_grid, v);
    return out
end

function intelnr(v)
    # Fornece o valor do retorno esperado de acordo com uma distribuição normal

    #out = pdf.(Normal(0,sig), v).*log.(erinsd(s, v));
    aux = map(s -> erinsd(s, v)[1], logS_grid)
    out = log.(aux);
    return out
end

function inter2(v)
    # Fornece o valor do retorno esperado de acordo com uma distribuição normal

    #out = pdf.(Normal(0,sig), v).*erinsd(s, v).^2;
    aux = map(s -> erinsd(s, v)[1], logS_grid);
    out = aux.^2;
    return out
end

function intelnr2(v)
    # Fornece o valor do retorno esperado de acordo com uma distribuição normal

    #out = pdf.(Normal(0,sig), v).*log.(erinsd(s, v)).^2;
    aux = map(s -> erinsd(s, v)[1], logS_grid);
    out = log.(aux).^2;
    return out
end

function intercb(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*ercbin(v);
    #out = map(s -> ercbin(s, v)[1], logS_grid);
    out = ercbin(logS_grid, v);
    return out
end


function intercb2(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*ercbin(v).^2;
    #aux = map(s -> ercbin(s, v)[1], logS_grid);
    #out = aux.^2;
    out = ercbin(logS_grid, v).^2;
    return out
end

function intelnrcb(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*log.(ercbin(v));
    #aux = map(s -> ercbin(s, v)[1], logS_grid);
    #out = log.(aux);
    out = log.(ercbin(v));
    return out
end

function intelnrcb2(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*log.(ercbin(v)).^2;
    #aux = map(s -> ercbin(s, v)[1], logS_grid);
    #out = log.(aux).^2;
    out = log.(ercbin(logS_grid, v)).^2
    return out
end

function ercbinN(v)
    # Retornos esperados da estrutura a termo nominal 
    #global s sg lnpcb matur

    # não sei se esses exp e log estao certos, pq aqui usei o preço e nao a taxa

    Fit = linear_interpolation(vec(logS_grid), vec((p_auxN)), extrapolation_bc = Line());
    #Fit2 = linear_interpolation(vec(logS_grid), vec((p_aux[:,end])), extrapolation_bc = Line());
    out = zeros(length(logS_grid), length(v));
    for i = 1:length(v)
        out[:, i] = exp.(Fit.(strans(logS_grid, v[i])))./exp.(p_auxtN);
    end
    return out
end

function intercbN(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*ercbinN(v);
    #out = map(s -> ercbinN(s, v)[1], logS_grid);;
    out = ercbinN(logS_grid, v)
    return out
end

function intercb2N(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*ercbinN(v).^2;
    #aux = map(s -> ercbinN(s, v)[1], logS_grid);
    #out = aux.^2;
    out = ercbinN(logS_grid, v).^2
    return out
end

function intelnrcbN(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*log.(ercbinN(v));
    #aux = map(s -> ercbinN(s, v)[1], logS_grid);
    #out = log.(aux);
    out = log.(ercbinN(logS_grid, v))
    return out
end

function intelnrcb2N(v)
    # Integrando dos retornos esperados dos títulos públicos 
    # ou, da estrutura a termo. 

    #out = pdf.(Normal(0,sig), v).*log.(ercbinN(v)).^2;
    #aux = map(s -> ercbinN(s, v)[1], logS_grid);
    #out = log.(aux).^2;
    out = log.(ercbin(logS_grid, v)).^2
    return out
end
 
# nao estou usando
function intpcb(v)
    # Função que fornece o preço dos bonds para cada maturidade %
    # global s sg lnpcb sig

    d = Normal(0, sig);
    #out = pdf(d,v) .*mrsinsd(v).* exp.(interp(strans(s,v),sg,lnpcb[length(lnpcb)]))';
    Fit         = linear_interpolation(vec(logS_grid), vec(lnpcb[:, end]), extrapolation_bc = Line());
    Fit_        = s -> max.(0, Fit.(s))
    out = pdf(d,v) .* mrsinsd(s, v) .* exp.(Fit_.(strans.(s,v)));
    return out
end

# nao estou usando
function intpcbN(v)
    # Função que fornece o preço dos bonds para cada maturidade %
    # global s sg lnpcb sig

    d = Normal(0, sig);
    #out = pdf(d,v) .*mrsinsd(v).* exp.(interp(strans(s,v),sg,lnpcb[length(lnpcb)]))';
    Fit         = linear_interpolation(vec(logS_grid), vec(pcbN[:, end]), extrapolation_bc = Line());
    mrsinsdN = delta*exp.(-gamma*(g.-(1-phi).*(s.-s_bar)) .+ (xiX - gamma*(1 .+ lambda.(s))).*v)
    out = pdf(d,v) .* mrsinsdN .* exp.(Fit.(strans.(s,v)));
    return out
end

