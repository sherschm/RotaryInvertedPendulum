import sympy as sp
import numpy as np
from sympy import Matrix, Function, symbols, simplify, diff
from parameters import rc,m,Ip,g,damping
from jax import jacfwd
import jax.numpy as jnp

# Define time variable and functions
t = symbols('t')
theta1 = Function('theta1')(t)
theta2 = Function('theta2')(t)
theta1d = Function('theta1d')(t)
theta2d = Function('theta2d')(t)
theta1dd = Function('theta1dd')(t)
theta2dd = Function('theta2dd')(t)
u = Function('u')(t)

# Pack variables
vars = (t, theta1, theta2, theta1d, theta2d, theta1dd, theta2dd, u)
x = Matrix([theta1, theta2, theta1d, theta2d])

# Rotation matrix function: assuming yaw-pitch rotation (Z-Y axes)
def calc_R20(theta1, theta2):
    Rz = Matrix([
        [sp.cos(theta1), -sp.sin(theta1), 0],
        [sp.sin(theta1),  sp.cos(theta1), 0],
        [0, 0, 1]
    ])
    Ry = Matrix([
        [sp.cos(theta2), 0, sp.sin(theta2)],
        [0, 1, 0],
        [-sp.sin(theta2), 0, sp.cos(theta2)]
    ])
    return Ry * Rz 

# Rotation matrix and derivative
R20 = calc_R20(theta1, theta2)

subs_dict = {
    sp.Derivative(theta1, t): theta1d,
    sp.Derivative(theta2, t): theta2d,
    sp.Derivative(theta1d, t): theta1dd,
    sp.Derivative(theta2d, t): theta2dd,
}

# Time derivative of R20
R20_dot = R20.diff(t).doit().subs(subs_dict)

# Angular velocity matrix
Omega = R20_dot * R20.T
omega = Matrix([
    Omega[2, 1],
    Omega[0, 2],
    Omega[1, 0]
])

# Rotational kinetic energy
# Ip is expected to be a 3x3 inertia matrix
if isinstance(Ip, sp.Symbol):
    Ip = sp.eye(3) * Ip  # default inertia matrix if not specified
T_rot = (1/2) * (omega.T * Ip * omega)[0]

# Center of mass velocity
v_com = R20_dot.T * Matrix(rc)

# Linear kinetic energy
T_lin = (1/2) * m * (v_com.T * v_com)[0]

# Total kinetic energy
T = simplify(T_rot + T_lin)

# Potential energy
V = simplify(m * g * (R20.T * Matrix(rc))[2])

# Display expressions
print("Kinetic Energy T:")
sp.pprint(T, use_unicode=True)

print("\nPotential Energy V:")
sp.pprint(V, use_unicode=True)

# Optional: create lambda functions for numerical use
T_func = sp.lambdify((theta1, theta2, theta1d, theta2d), T, modules='numpy')
V_func = sp.lambdify((theta1, theta2, theta1d, theta2d), V, modules='numpy')

# Function to compute total energy
def Total_energy(x_vals):
    return T_func(*x_vals) + V_func(*x_vals)

def Lagrangian_dynamics(T, V):
    """
    Compute symbolic Lagrangian dynamics: Mass matrix M and nonlinear vector N.
    Parameters:
        T: Kinetic energy expression (SymPy)
        V: Potential energy expression (SymPy)
        vars: Tuple of symbolic variables (t, θ1(t), θ2(t), θ1d(t), θ2d(t), θ1dd(t), θ2dd(t), u(t))
    Returns:
        M: 2x2 SymPy Matrix (mass matrix)
        N: 2x1 SymPy Matrix (nonlinear terms)
    """
    L = T - V  # Lagrangian
    
    # Partial derivatives of L
    dL_dθ = [sp.diff(L, theta1), sp.diff(L, theta2)]
    dL_dθd = [sp.diff(L, theta1d), sp.diff(L, theta2d)]

    # Total derivatives w.r.t. time
    dt_dL_dθd = sp.Matrix(dL_dθd).diff(t).doit()

    # Euler-Lagrange equations
    Eq = [dt_dL_dθd[i] - dL_dθ[i] for i in range(2)]

    # Substitute time derivatives
    subs = {
        sp.Derivative(theta1, t): theta1d,
        sp.Derivative(theta2, t): theta2d,
        sp.Derivative(theta1d, t): theta1dd,
        sp.Derivative(theta2d, t): theta2dd
    }
    Eq = [sp.simplify(eq.subs(subs)) for eq in Eq]
    # Mass matrix extraction
    thetadd_vec = sp.Matrix([theta1dd, theta2dd])
    M = sp.Matrix(Eq).jacobian(thetadd_vec)

    #M = sp.Matrix([[M11, M12], [M21, M22]])

    # Compute nonlinear terms: N = Eq - M*[θ1dd; θ2dd]
    accel_vec = sp.Matrix([theta1dd, theta2dd])
    N = sp.simplify(sp.expand(sp.Matrix(Eq) - M * accel_vec))

    return M, N

# Define T, V before calling
M, N = Lagrangian_dynamics(T, V)

x_syms = [theta1, theta2, theta1d, theta2d]
M_f = sp.lambdify(x, M.tolist(), modules="jax")
N_f = sp.lambdify(x, N.tolist(), modules="jax")

def dynamics_acc_ctrl_terms(x):

    """
    Compute modified dynamics with acceleration as control input.

    Args:
        M_f: function returning mass matrix M given state x (length-4 vector)
        N_f: function returning nonlinear vector N given state x (length-4 vector)
        x: state vector [theta1, theta2, theta1d, theta2d]
        Damping: scalar damping coefficient

    Returns:
        M_acc: Mass matrix (2×2)
        N_acc: Modified nonlinear vector (2×1)
        B_acc: Control input mapping matrix (2×1)
    """
    A = jnp.array([[1.0, 0.0]])  # constraint matrix
    M = jnp.array(M_f(x[0], x[1], x[2], x[3]))
    N = jnp.array(N_f(x[0], x[1], x[2], x[3])).flatten()

    D_mat = jnp.array([[0.0, 0.0],
                      [0.0, -damping]])
    Damping_force = (D_mat @ jnp.asarray(x)[2:]).flatten()

    M_inv = jnp.linalg.inv(M)
    AMinvAT_inv = jnp.linalg.inv(A @ M_inv @ A.T)  # scalar

    N_bar = N + Damping_force
    proj = A.T @ AMinvAT_inv @ A @ M_inv
    N_acc = N_bar - proj @ N_bar
    B_acc = A.T @ AMinvAT_inv
    M_acc = M

    return M_acc, N_acc, B_acc

# -----------------------------------------------------------------------------------
# 2. Returns dx/dt for first-order state-space dynamics
# -----------------------------------------------------------------------------------
def rot_pend_dynamics_num(x, u):
    """
    Evaluate the first-order ODE dynamics of the rotating pendulum system.

    Args:
        x: state vector [theta1, theta2, theta1d, theta2d] (length 4)
        u: control input (desired angular acceleration θ̈₁)
        M_f: function returning M(x)
        N_f: function returning N(x)
        Damping: scalar damping coefficient

    Returns:
        dxdt: derivative of state vector (length 4)
    """
    # Get modified mass matrix, nonlinear vector, input mapping
    M_a, N_a, B_a = dynamics_acc_ctrl_terms(x)

    # Compute acceleration from projected dynamics
    acc = jnp.linalg.inv(M_a) @ (B_a.flatten() * u - N_a.T)  # shape (2,)
    #print(acc)
    dxdt = [x[2], x[3], acc[0], acc[1]]
    return dxdt

def f_wrapped(xu):
    x = xu[:4]
    u = xu[4]
    return rot_pend_dynamics_num(x, u)

Jacobian = jacfwd(f_wrapped)(jnp.array([0.0, np.pi, 0.0, 0.0, 0.0]))
A_matrix = jnp.stack(Jacobian, axis=0)[:4,:4]
B_matrix = jnp.stack(Jacobian, axis=0)[:,5]

# M, N built previously from your Lagrangian_dynamics(...)
dxdt = rot_pend_dynamics_num([0.0, np.pi, 0.0, 0.0], 0.0)

#A_sym = rot_pend_dynamics_num.jacobian(x)

print("dx/dt =", dxdt)

print(A_matrix)
print(B_matrix)

