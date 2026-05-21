clc; clear; close all;

%% Q1

disp('========== Q1 ) System Idendification ==========');

% Part a
Ts = 0.01;
T_final = 5;
t = 0:Ts:T_final;
u = zeros(size(t));
u(1) = 1/Ts;
sim_u.time = t';
sim_u.signals.values = u';
sim_u.signals.dimensions = 1;

% Part b
simOut = sim('blackboxsystem', 'StartTime', '0', 'StopTime', num2str(T_final), 'FixedStep', num2str(Ts), 'Solver', 'ode4');
y = simOut.sim_y.Data;
if size(y,1) > size(y,2)
    y = y';
end
y = y * Ts * 2; 

% Part c
L = 250;
y_dyn = y(2:end); 
H = zeros(L, L);
Hp = zeros(L, L);
for i = 1:L
    for j = 1:L
        k = i + j - 1;
        H(i,j) = y_dyn(k);
        k_shift = i + j;
        if k_shift <= length(y_dyn)
            Hp(i,j) = y_dyn(k_shift);
        end
    end
end

[U, S, V] = svd(H);
sigma = diag(S);

disp('Top 10 Singular Values:');
disp(sigma(1:10));

figure(1);
semilogy(sigma, 'o-', 'LineWidth', 1.5);
title('Q1) Singular Values');
xlabel('Index');
ylabel('Magnitude (Log Scale)');
grid on;

% Part d
n = 4;

% Part e
Un = U(:, 1:n);
Sn = S(1:n, 1:n);
Vn = V(:, 1:n);
Obs = Un * sqrt(Sn);
Con = sqrt(Sn) * Vn';
Cd = Obs(1, :);
Bd = Con(:, 1);
Sn_sqrt_inv = inv(sqrt(Sn));
Ad = Sn_sqrt_inv * Un' * Hp * Vn * Sn_sqrt_inv;

disp('Ad (Discrete):');
disp(Ad);
disp('Bd (Discrete):');
disp(Bd);
disp('Cd (Discrete):');
disp(Cd);

% Part f
Dd = y(1);
sys_d = ss(Ad, Bd, Cd, Dd, Ts);
sys_c = d2c(sys_d, 'zoh');
[Ac, Bc, Cc, Dc] = ssdata(sys_c);
disp('Ac (Continuous):');
disp(Ac);
disp('Bc (Continuous):');
disp(Bc);
disp('Cc (Continuous):');
disp(Cc);
disp('Dc (Continuous):');
disp(Dc);

%% Q2

disp('========== Q2 ) Validation ==========');

% Part a & b
u_step = ones(size(t));
sim_u.time = t';
sim_u.signals.values = u_step';
sim_u.signals.dimensions = 1;
simOut_step = sim('blackboxsystem', 'StartTime', '0', 'StopTime', num2str(T_final), 'FixedStep', num2str(Ts), 'Solver', 'ode4');
y_real = simOut_step.sim_y.Data;
if size(y_real, 1) > size(y_real, 2)
    y_real = y_real';
end
[y_est, t_est] = step(sys_c, t);
y_est = y_est';

% Part c
figure(2);
plot(t, y_real, 'b', 'LineWidth', 2); hold on;
plot(t, y_est, 'r--', 'LineWidth', 1.5);
title('Q2) Step Response Comparison');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Black Box (Real)', 'Identified Model');
grid on;

% Part d
error_val = norm(y_real - y_est) / norm(y_real);
disp('Relative Error:');
disp(error_val);

%% Q3

disp('========== Q3 ) State Space Matrices ==========');
format short
disp('Matrix A (Continuous):');
disp(Ac);
disp('Matrix B (Continuous):');
disp(Bc);
disp('Matrix C (Continuous):');
disp(Cc);
disp('Matrix D (Continuous):');
disp(Dc);

%% Q4
disp('========== Q4 ) Stability Analysis ==========');

% Part a
Poles = eig(Ac);
disp('System Poles (Eigenvalues of Ac):');
disp(Poles);

% Part b
real_parts = real(Poles);
if all(real_parts < 0)
    disp('Stability Result: The system is STABLE (All poles are in LHP).');
else
    disp('Stability Result: The system is UNSTABLE.');
end

% Part c
figure(3);
plot(real(Poles), imag(Poles), 'bx', 'MarkerSize', 12, 'LineWidth', 2);
hold on;
yline(0, 'k--', 'LineWidth', 1); 
xline(0, 'k--', 'LineWidth', 1); 
title('Q4) Pole Locations');
xlabel('Real Axis');
ylabel('Imaginary Axis');
grid on;
hold off;

%% Q5
disp('========== Q5 ) Controllability & Observability ==========');

% Part a
Cnt_manual = [Bc, Ac*Bc, Ac^2*Bc, Ac^3*Bc]

% Part b
Obs_manual = [Cc; Cc*Ac; Cc*Ac^2; Cc*Ac^3]

% Part c
n = 4;
rank_C = rank(Cnt_manual);
rank_O = rank(Obs_manual);

disp('Rank of Controllability Matrix:');
disp(rank_C);

if rank_C == n
    disp('Result: The system is FULLY CONTROLLABLE.');
else
    disp('Result: The system is UNCONTROLLABLE.');
end

disp('Rank of Observability Matrix:');
disp(rank_O);

if rank_O == n
    disp('Result: The system is FULLY OBSERVABLE.');
else
    disp('Result: The system is UNOBSERVABLE.');
end

% Part d
Cnt_builtin = ctrb(Ac, Bc);
Obs_builtin = obsv(Ac, Cc);

disp('Built-in Controllability Matrix (ctrb):');
disp(Cnt_builtin);
disp('Rank (ctrb):');
disp(rank(Cnt_builtin));

disp('Built-in Observability Matrix (obsv):');
disp(Obs_builtin);
disp('Rank (obsv):');
disp(rank(Obs_builtin));

%% Q6

disp('========== Q6 ) Minimality Check ==========');

sys_min = minreal(sys_c);
sys_min
n_original = order(sys_c);
n_minimal = order(sys_min);

disp(['Order of Original System: ', num2str(n_original)]);
disp(['Order of Minimal System:  ', num2str(n_minimal)]);

%% Q7
disp('========== Q7 ) State Transition Matrix & Transfer Function ==========');

% Part a
disp('State Transition Matrix Phi(t) evaluated at t = 1s:');
Phi_at_1 = expm(Ac * 1);
disp(Phi_at_1);

% Part b
disp('Transfer Function G(s):');
sys_c_clean = sys_c;
sys_c_clean.D = 0;
G = minreal(tf(sys_c_clean), 1e-4);
G

% Part d
disp('Output of "tf" command:');
disp(tf(sys_c_clean));

disp('Output of "pole" command:');
disp(pole(sys_c_clean));

disp('Output of "zero" command:');
disp(zero(sys_c_clean));

%% Q8
disp('========== Q8 ) Pole Placement Design ==========');

P_slow = [-1, -1.6, -2.2, -2.8]; 
K_slow = place(Ac, Bc, P_slow);

N_slow = inv( -(Cc * inv(Ac - Bc*K_slow) * Bc) );

disp('--- Design 1 (Slow) ---');
disp('Desired Poles:'); disp(P_slow);
disp('Gain Matrix K_slow:'); disp(K_slow);
disp('Pre-filter Gain N_slow:'); disp(N_slow);

P_fast = [-5, -6, -7, -8];
K_fast = place(Ac, Bc, P_fast);

N_fast = inv( -(Cc * inv(Ac - Bc*K_fast) * Bc) );

disp('--- Design 2 (Fast) ---');
disp('Desired Poles:'); disp(P_fast);
disp('Gain Matrix K_fast:'); disp(K_fast);
disp('Pre-filter Gain N_fast:'); disp(N_fast);

%% Q10
disp('========== Q10 ) Integral Control Design ==========');

% Part a

n = 4;
Aa = [Ac, zeros(n, 1); -Cc, 0];
Ba = [Bc; 0];

disp('Part a) Augmented System Matrices Created.');

% Part b
Ra_aug = ctrb(Aa, Ba);
rank_aug = rank(Ra_aug);

disp('Rank of Augmented Controllability Matrix:');
disp(rank_aug);

if rank_aug == n + 1
    disp('Result: The Augmented System is FULLY CONTROLLABLE.');
else
    disp('Result: The Augmented System is NOT Controllable.');
end

% Part c
P_aug = [-2, -2.5, -3, -3.5, -4]; 

disp('Desired Poles for Augmented System:');
disp(P_aug);

K_aug = place(Aa, Ba, P_aug);

K_x = K_aug(1:n);
K_I = K_aug(n+1);

disp('Calculated Gains:');
disp('State Feedback Gain (K_x):'); disp(K_x);
disp('Integral Gain (K_I):');       disp(K_I);

%% Q11
disp('========== Q11 ) Observer Design ==========');

% Part b
P_obs = [-25, -30, -35, -40]; 

disp('Desired Observer Poles:');
disp(P_obs);

Lt = place(Ac', Cc', P_obs);
L = Lt';

disp('Observer Gain Matrix L:');
disp(L);

disp('Eigenvalues of (A - LC):');
disp(eig(Ac - L*Cc));

K_val = K_fast; 
disp('Using Controller Gain K (from Q8 Fast design):');
disp(K_val);

x0_sys = [0.1; -0.1; 0.05; -0.05]; 

x0_obs = [0; 0; 0; 0];

%% Q12
disp('========== Q12 ) Reduced-Order Observer ==========');

F_obs = diag([-25, -30, -35]);

L_red = [1; 1; 1];

T = sylvester(-F_obs, Ac, L_red * Cc);

H_obs = T * Bc;

disp('F_obs:'); disp(F_obs);
disp('L_red:'); disp(L_red);
disp('T matrix:'); disp(T);
disp('H_obs (Input Gain):'); disp(H_obs);

M = [Cc; T];
if rank(M) < 4
    error('Matrix [C; T] is singular.');
end

InvM = inv(M);
N1 = InvM(:, 1);
N2 = InvM(:, 2:end);

disp('Reconstruction Gain N1 (for y):'); disp(N1);
disp('Reconstruction Gain N2 (for z):'); disp(N2);

%% Q14
disp('========== Q14 ) Lyapunov Analysis ==========');

% Part a
Q1 = eye(4);

A_T = Ac';
I = eye(4);
M = kron(I, A_T) + kron(A_T, I);
vec_Q = Q1(:);

vec_P = -inv(M) * vec_Q;
P1_manual = reshape(vec_P, 4, 4);

disp('Part a) Matrix P Calculated:');
disp(P1_manual);

% Part b
eig_P1 = eig(P1_manual);
disp('Part b) Eigenvalues of P1:');
disp(eig_P1);

if all(eig_P1 > 0)
    disp('Result: P is Positive Definite.');
else
    disp('Result: P is NOT Positive Definite.');
end

% Part d
P1_command = lyap(Ac', Q1);

disp('Part d) Matrix P from "lyap" command:');
disp(P1_command);

% Part e
Q2 = diag([10, 5, 20, 2]); 
P2 = lyap(Ac', Q2);
eig_P2 = eig(P2);

disp('Part e) Analysis with new Q2:');
disp('New Q2 matrix:'); disp(Q2);
disp('New P2 matrix:'); disp(P2);
disp('Eigenvalues of New P2:'); disp(eig_P2);

%% Q15
disp('========== Q15 ) Lyapunov Energy Function ==========');

P_lyap = P1_manual; 
Q_lyap = eye(4);

x_test = [0.1; -0.1; 0.05; -0.05];

% Part a & b
V_value = x_test' * P_lyap * x_test;

disp('Value of V(x) = x''Px:');
disp(V_value);

if V_value > 0
    disp('Check: V(x) is POSITIVE (Correct).');
else
    disp('Check: V(x) is NOT positive (Error).');
end

% Part c
V_dot_value = -x_test' * Q_lyap * x_test;

disp('Value of V_dot(x) = -x''Qx:');
disp(V_dot_value);

if V_dot_value < 0
    disp('Check: V_dot(x) is NEGATIVE (Energy is decreasing).');
else
    disp('Check: V_dot(x) is NOT negative (Error).');
end

%% Q16
disp('========== Q16 ) Closed-Loop Stability ==========');

% Part a
A_cl = Ac - Bc * K_fast;

disp('Part a) Closed-Loop System Matrix (A_cl):');
disp(A_cl);

disp('Eigenvalues of A_cl (Verification):');
disp(eig(A_cl));

% Part b
Q_cl = eye(4);

P_cl = lyap(A_cl', Q_cl);

disp('Part b) Closed-Loop Lyapunov Matrix (P_cl):');
disp(P_cl);

% Part c
eig_P_cl = eig(P_cl);
disp('Eigenvalues of P_cl:');
disp(eig_P_cl);

if all(eig_P_cl > 0)
    disp('Result: P_cl is Positive Definite. Closed-Loop System is Asymptotically Stable.');
else
    disp('Result: P_cl is NOT Positive Definite.');
end

%% Q17
disp('========== Q17 ) Controllable Canonical Form ==========');

% Part a
p = poly(Ac);
a1 = p(2);
a2 = p(3);
a3 = p(4);
a4 = p(5);

disp('a) Characteristic Polynomial Coefficients:');
disp(['a1 = ', num2str(a1)]);
disp(['a2 = ', num2str(a2)]);
disp(['a3 = ', num2str(a3)]);
disp(['a4 = ', num2str(a4)]);

% Part b
M_c = [Bc, Ac*Bc, Ac^2*Bc, Ac^3*Bc];

W = [a3, a2, a1, 1;
     a2, a1, 1, 0;
     a1, 1, 0, 0;
     1,  0, 0, 0];

Tc = M_c * W;

disp('b) Transformation Matrix Tc:');
disp(Tc);

% Part c
A_bar = inv(Tc) * Ac * Tc;
B_bar = inv(Tc) * Bc;
C_bar = Cc * Tc;

disp('c) Converted Matrix A_bar:');
disp(A_bar);
disp('Converted Vector B_bar:');
disp(B_bar);
disp('Converted Vector C_bar:');
disp(C_bar);

% Part e
sys_original = ss(Ac, Bc, Cc, 0);
[sys_canon, T_matlab] = canon(sys_original, 'companion');

disp('e) Result of "canon" command (A matrix):');
disp(sys_canon.A);

%% Q18
disp('========== Q18 ) Observability Canonical Form ==========');

% Part a
Ob = obsv(Ac, Cc);
To = inv(W * Ob);

disp('a) Transformation Matrix To:');
disp(To);

% Part b
A_o = inv(To) * Ac * To;
B_o = inv(To) * Bc;
C_o = Cc * To;

disp('b) Converted Matrix A_o:');
disp(A_o);
disp('Converted Vector C_o:');
disp(C_o);
disp('Converted Vector B_o:');
disp(B_o);

% Part d
expected_Ao = A_bar';

error_duality = norm(A_o - expected_Ao);
disp('d) Difference between A_o and A_c transpose:');
disp(error_duality);

%% Q19
disp('========== Q19 ) Jordan Form & Modal Analysis ==========');

% Part a
[T_J, J] = eig(Ac);

disp('a) Jordan Matrix J (Diagonal of Eigenvalues):');
disp(J);
disp('Transformation Matrix Tj (Eigenvectors):');
disp(T_J);

% Part b
A_J = inv(T_J) * Ac * T_J;
B_J = inv(T_J) * Bc;
C_J = Cc * T_J;

disp('b) Input Matrix in Jordan Form (B_J):');
disp(B_J);
disp('Output Matrix in Jordan Form (C_J):');
disp(C_J);

% Part c
n = length(Ac);
eigenvalues = diag(J);
tol = 1e-5;

disp('c) Modal Controllability and Observability Tests (PBH):');
disp('-----------------------------------------------------');

for i = 1:n
    lambda = eigenvalues(i);
    
    M_con = [Ac - lambda*eye(n), Bc];
    rank_con = rank(M_con, tol);
    
    M_obs = [Ac - lambda*eye(n); Cc];
    rank_obs = rank(M_obs, tol);
    
    fprintf('Mode %d (lambda = %.4f + %.4fi):\n', i, real(lambda), imag(lambda));
    
    if rank_con == n
        fprintf('  -> Controllable (Rank = %d)\n', rank_con);
    else
        fprintf('  -> NOT Controllable (Rank = %d)\n', rank_con);
    end
    
    if rank_obs == n
        fprintf('  -> Observable (Rank = %d)\n', rank_obs);
    else
        fprintf('  -> NOT Observable (Rank = %d)\n', rank_obs);
    end
    disp('-----------------------------------------------------');
end