%Estructura matrices W (HFB)
clc;clear;
format long
%Definición base de estados monoparticulares
level_sp=[1 1 2 2 3 3];
level   =[0 0 0 0 0 0];    % nivel de energía 
m=6; % Dimension matrices, grados de libertad 
spin=[-1 1 -1 1 -1 1];   % 1=up , -1=down 
N0=2; %Nº de particulas

%Matriz compleja aleatoria
%rng(1); %Semilla
%Z0 = (randn(m) + 1i*randn(m));
Z0=randn(m);
%Imponer antisimetria
Z0 = Z0 - Z0.';

% construir U
M = eye(m) + Z0'*Z0;
U0 = inv(sqrtm(M));
V0 = Z0*U0;

% V0=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0;0 0 0 0 0 0]; %Particulas en el nivel 2
% U0=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1];


% construir W con estructura HFB
 W = [U0 conj(V0); V0 conj(U0)];
% 
% % comprobar unitariedad
 error = norm(W'*W - eye(2*m));
% 
% %Comprobar condiciones U y V
 cond1 = norm(U0'*U0 + V0'*V0)-1;
 cond2 = norm(U0.'*V0 + V0.'*U0);
% 
% %Construcción de las matriz de densidad (rho) y el tensor de paridad (kappa)
 

 rho=conj(V0)*V0.';
 kappa=conj(V0)*U0.';

% %Comprobación de sus propiedades

 prop1=rho'-rho; %Hermiticidad de la matriz de densidad (prop1==zeros(N))
 prop2=kappa.'+kappa; %Antisimetria del tensor de paridad

% Definición de Delta y Gamma
%Definición tensor v(abcd)
G = 1;

v_barra = zeros(m,m,m,m);

for a=1:m
for b=1:m
for c=1:m
for d=1:m

    if level_sp(a)==level_sp(b) && level_sp(c)==level_sp(d) %Deltas de Kronecker
      if spin(a)==1 && spin(b)==-1  && spin(c)==1 && spin(d)==-1
        v_barra(a,b,c,d) = -G;
      elseif spin(a)==-1 && spin(b)==1 && spin(c)==1 && spin(d)==-1 
         v_barra(a,b,c,d) = G;
      elseif spin(a)==1 && spin(b)==-1 && spin(c)==-1 && spin(d)==1 
          v_barra(a,b,c,d) = G;
      elseif spin(a)==-1 && spin(b)==1 && spin(c)==-1 && spin(d)==1
          v_barra(a,b,c,d) = -G;
      end

    end

end
end
end
end



%Construccion gamma 

Gamma = zeros(m,m);
for a=1:m
for c=1:m
    
    sum_val = 0;
    
    for d=1:m
    for b=1:m
        
        sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);
        
    end
    end
    
    Gamma(a,c) = sum_val;
    
end
end

%Construccion Delta
Delta = zeros(m,m);
for a=1:m
for b=1:m
    
    sum_val = 0;
    
    for c=1:m
    for d=1:m
       
        sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);
        
    end
    end
    
    Delta(a,b) = 0.5*sum_val;
    
end
end


% Energías monoparticulares 
eps = level;  

t = diag(level);
h = t + Gamma;

%Def de H20

H20 = U0' * h * conj(V0) ...
    + U0' * Delta * conj(U0) ...
    - V0' * h.' * conj(U0) ...
    - V0' * conj(Delta) * conj(V0);


%Def de N20

N20 = U0' * conj(V0) - V0' * conj(U0);

%Comprobacion antisimetria
check_H20 = norm(H20.' + H20); %¿=0?
check_N20 = norm(N20.' + N20); %¿=0?

%Cálculo de E_HFB

E1 = trace(t * rho);
E2 = 0.5 * trace(Gamma * rho);
E3 = -0.5 * trace(Delta * conj(kappa));
E_HFB = E1 + E2 + E3;
N=trace(rho); %Hay que ajustarlo al valor deseado
%ES REAL?
check_EReal =norm(imag(E_HFB)); %¿=0?


%% Ajuste puro del constraint (Número de partículas)
% --- Parámetros de control ---
max_iter = 10;    % Máximo de iteraciones
tol = 1e-6;        % Tolerancia de convergencia

fprintf('Iteración | Num. Partículas | Error N \n');
fprintf('------------------------------------------\n');

U = U0;
V = V0;

for k = 1:max_iter
    % 1. Recalcular las densidades actuales
    rho = conj(V) * V.';   
    N_calc=real(trace(rho));
    N20 = U' * conj(V) - V' * conj(U);

    % 2. Evaluar el error
    error_N = N0 - N_calc;

    % Criterio de parada temprano
    if abs(error_N) < tol
        fprintf('Convergencia alcanzada en la iteración %d con N=%d\n', k,N_calc);
        break;
    end
    % 3. Cálculo del paso (eta_N)  (Ec. (B.8))
    eta_N = error_N / real(trace(N20 * N20'));

    % 4. Dirección de actualización Z (Solo usamos N20)
    Z = eta_N * N20;
    Z = 0.5 * (Z - Z.'); % Asegurar que la matriz Z sea estrictamente antisimétrica 

    % 5. Reconstrucción de U y V (Preservando unitariedad mediante Cholesky)
    A0 = eye(m) + Z.' * conj(Z);
    L0 = chol(A0, 'lower');

    Uold = U;
    Vold = V;

    % Actualizacion U,V (ec. 25)
    U = (Uold + conj(Vold) * conj(Z)) * (inv(L0))';
    V = (Vold + conj(Uold) * conj(Z)) * (inv(L0))';

    % 6. Verificación de unitariedad Bogoliubov
    tol_unitarity = 1e-7;
    condic1 = norm(U'*U + V'*V - eye(m));
    condic2 = norm(U.'*V + V.'*U);

    if condic1 > tol_unitarity || condic2 > tol_unitarity
        warning('Iteración %d: ¡Pérdida de unitariedad detectada!', k);
        fprintf('  Error Condición 1: %.2e\n', condic1);
        fprintf('  Error Condición 2: %.2e\n', condic2);
    end
    % Mostrar progreso
    fprintf('%9d | %15.4f | %10.2e \n', k, N_calc, error_N);
end





%% Minimización de Energía por Gradiente con Paso Adaptativo y Constraint
% --- Parámetros de Control ---
max_iter_E = 1000;     % Iteraciones máximas del bucle de energía (externo)
max_iter_N = 10;      % Iteraciones máximas del bucle de partículas (interno)
eta_E = 0.02;         % ¡Paso inicial agresivo! Se auto-ajustará solo.
tol_E = 1e-6;         % Tolerancia para la convergencia en energía
tol_N = 1e-6;         % Tolerancia para el número de partículas


fprintf('Iter Ext | Energía HFB | Num. Partículas | Delta E  | Iter Int\n');
fprintf('---------------------------------------------------------------------------------------\n');

% --- Inicialización para la gráfica ---
historial_E = [];    % Almacenará los valores de energía
historial_iter = []; % Almacenará el número de iteración exitosa
E_actual=E_HFB;

for iter_E = 1:max_iter_E
    % =========================================================================
    % 1. CÁLCULO DE GRADIENTES EN EL PUNTO ACTUAL
    % =========================================================================
    rho = conj(V) * V.';
    kappa = conj(V)*U.';
    % Gamma actual

    Gamma = zeros(m,m);

    for a=1:m
        for c=1:m

            sum_val = 0;

            for d=1:m
                for b=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);

                end
            end

            Gamma(a,c) = sum_val;

        end
    end

    % Delta actual
    Delta = zeros(m,m);

    for a=1:m

        for b=1:m

            sum_val = 0;

            for c=1:m
                for d=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);

                end
            end

            Delta(a,b) = 0.5*sum_val;

        end
    end
    

    N_final = real(trace(rho));
    h = t + Gamma;
    H20 = U' * h * conj(V) + U' * Delta * conj(U) - V' * h.' * conj(U) - V' * conj(Delta) * conj(V);
    N20 = U' * conj(V) - V' * conj(U); 
    E_old=E_actual;
    % =========================================================================
    % 2. PROYECCIÓN DEL MULTIPLICADOR DE LAGRANGE (LAMBDA)
    % =========================================================================
    denom_lambda = trace(N20 * N20.');
    if denom_lambda > 1e-12
        lambda = trace(H20 * N20.') / denom_lambda;
    else
        lambda = 0;
    end

    % =========================================================================
    % 3. PASO DE GRADIENTE PARA LA ENERGÍA
    % =========================================================================
    Z_E = -eta_E * (H20 - lambda * N20); %Direccion de gradiente
    Z_E = 0.5 * (Z_E - Z_E.'); %Unitariedad Z (forzada)

    A0 = eye(m) + Z_E.' * conj(Z_E);
    L0 = chol(A0, 'lower');

    U_old=U;
    V_old=V;

    U = (U_old + conj(V_old) * conj(Z_E)) * (inv(L0))';
    V = (V_old + conj(U_old) * conj(Z_E)) * (inv(L0))';

  
    % =========================================================================
    % 4. BUCLE INTERNO: RESTAURACIÓN ESTRICTA DEL NÚMERO DE PARTÍCULAS
    % =========================================================================
    sub_iter = 0; 
    for k_N = 1:max_iter_N
        rho_sub = conj(V) * V.';
        N_sub = real(trace(rho_sub));
        error_N = N0 - N_sub;

        if abs(error_N) < tol_N
            sub_iter = k_N - 1;
            break;
        end

        N20_sub = U' * conj(V) - V' * conj(U);
        norm_N20 = real(trace(N20_sub * N20_sub'));

        if norm_N20 < 1e-14, break; end

        eta_N = error_N / norm_N20;
        Z_N = eta_N * N20_sub;
        Z_N = 0.5 * (Z_N - Z_N.');

        A_sub = eye(m) + Z_N.' * conj(Z_N);
        L_sub = chol(A_sub, 'lower');

        Uold_sub = U; Vold_sub = V;
        U = (Uold_sub + conj(Vold_sub) * conj(Z_N)) * (inv(L_sub))';
        V = (Vold_sub + conj(Uold_sub) * conj(Z_N)) * (inv(L_sub))';
    end

    % =========================================================================
    % 5. EVALUACIÓN DE LA ENERGÍA DEL NUEVO ESTADO CANDIDATO
    % =========================================================================
    rho = conj(V) * V.';
    kappa =conj(V) * U.';
    N_final = real(trace(rho));

    % Recalcular Gamma y Delta en el nuevo punto
    Gamma = zeros(m,m);
    for a=1:m
        for c=1:m

            sum_val = 0;

            for d=1:m
                for b=1:m

                    sum_val = sum_val + v_barra(a,b,c,d)*rho(d,b);

                end
            end

            Gamma(a,c) = sum_val;

        end
    end

    % Delta actual
    Delta = zeros(m,m);
    for a=1:m
    for b=1:m
        sum_val = 0;

            for c=1:m
            for d=1:m

                sum_val = sum_val + v_barra(a,b,c,d)*kappa(c,d);

            end
            end

            Delta(a,b) = 0.5*sum_val;

    end
    end
    

    E_actual = trace(t * rho) + 0.5 * trace(Gamma * rho) - 0.5 * trace(Delta * conj(kappa)); 
   % Monitorear el progreso en consola
    fprintf('%8d | %12.6f | %15.4f | %10.2e | %d\n', ...
            iter_E, E_actual, N_final, E_old-E_actual, sub_iter);

    %Verificación de unitariedad Bogoliubov
    tol_unitarity = 1e-7;
    condic1 = norm(U'*U + V'*V - eye(m));
    condic2 = norm(U.'*V + V.'*U);

    if condic1 > tol_unitarity || condic2 > tol_unitarity
        warning('Iteración %d: ¡Pérdida de unitariedad detectada!', k);
        fprintf('  Error Condición 1: %.2e\n', condic1);
        fprintf('  Error Condición 2: %.2e\n', condic2);
    end       

    % Criterio de parada absoluto (convergencia real)
    if abs(E_old-E_actual) < tol_E && abs(N_final - N0) < tol_N
        fprintf('-----------------------------------------------------------------------------------\n');
        fprintf('¡Convergencia HFB absoluta alcanzada en la iteración %d con E = %.6f!\n', iter_E, E_actual);
        break;

    end
    historial_E = [historial_E, real(E_actual)];
    historial_iter = [historial_iter, iter_E];
end
%E3=- 0.5 * trace(Delta * conj(kappa));
% =========================================================================
% 7. REPRESENTACIÓN GRÁFICA DE LA CONVERGENCIA
% =========================================================================
% if ~isempty(historial_E)
%     close all
%     figure('Name', 'Convergencia HFB', 'Color', 'w');
% 
%     % Gráfica principal con marcadores
%     plot(historial_iter, historial_E, '-o', ...
%         'LineWidth', 1.5, ...
%         'MarkerSize', 6, ...
%         'MarkerEdgeColor', [0 0.4470 0.7410], ...
%         'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
%         'Color', [0 0.4470 0.7410]);
% 
%     % Estética y etiquetas
%     title('Minimización de Energía HFB por Gradiente', 'FontSize', 14, 'FontWeight', 'bold');
%     xlabel('Iteraciones Externas (Aceptadas)', 'FontSize', 12);
%     ylabel('Energía Total (E)', 'FontSize', 12);
%     grid on;
% 
%     % Ajustar márgenes 
%     ax = gca;
%     ax.GridAlpha = 0.3;
%     ax.FontSize = 11;
% else
%     disp('No se registraron iteraciones exitosas para graficar.');
% end