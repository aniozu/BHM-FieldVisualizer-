import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from functions import *

st.set_page_config(layout="wide")

# パラメータ設定（Streamlit用）
st.sidebar.title("パラメータ設定")
a = st.sidebar.slider("半径 (a)", min_value=0.5, max_value=2.0, value=1.5, step=0.1)
L = st.sidebar.slider("長さ (L)", min_value=1.0, max_value=5.0, value=2.0, step=0.1)
N = st.sidebar.slider("本数 (N)", min_value=10, max_value=20, value=14, step=1)

# グリッド設定
Ngrid = 200
rho_vals = np.linspace(-5, 5, Ngrid)
z_vals = np.linspace(-5, 5, Ngrid)
RHO, Z = np.meshgrid(rho_vals, z_vals)

# プロット
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# 長方形の座標
rect_x = [-a, a, a, -a, -a]
rect_y = [-L/2, -L/2, L/2, L/2, -L/2]

for i in range(3):
    axs[i].plot(rect_x, rect_y, 'k-', linewidth=2)
    axs[i].fill([-a, a, a, -a], [0, 0, L/2, L/2], color='red', alpha=0.2)
    axs[i].fill([-a, a, a, -a], [-L/2, -L/2, 0, 0], color='green', alpha=0.2)
    axs[i].set_xlim(-5, 5)
    axs[i].set_ylim(-5, 5)
    axs[i].set_aspect(1)
    axs[i].set_xticks([])
    axs[i].set_yticks([])

# H_rho計算
H_rho_vals = H_rho(RHO, Z, a, L)
H_z_vals = H_z(RHO, Z, a, L)

# Hの流線 (streamplot)
start_points_H = getStartPointsH(a, L, N)
stream_H = axs[1].streamplot(RHO, Z, H_rho_vals, H_z_vals, color='k', density=20, arrowsize=2, start_points=start_points_H)
axs[1].set_title(r"$\mu _0H$")

# B計算
B_rho_vals = B_rho(RHO, Z, a, L)
B_z_vals = B_z(RHO, Z, a, L)

# Bの流線 (streamplot)
start_points_B = getStartPointsB(a, L, N, stream_H)
axs[0].streamplot(RHO, Z, B_rho_vals, B_z_vals, color='k', density=20, arrowsize=2, start_points=start_points_B)
axs[0].set_title(r"$B$")

# M計算
M_rho_vals = M_rho(RHO, Z, a, L)
M_z_vals = M_z(RHO, Z, a, L)

# Mの流線 (streamplot)
start_points_M = getStartPointsM(a, L, N)
axs[2].streamplot(RHO, Z, M_rho_vals, M_z_vals, color='k', density=100, arrowsize=2, start_points=start_points_M)
axs[2].set_title(r"$\mu _0M$")

st.pyplot(fig)