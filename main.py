import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
from functions import *

st.set_page_config(layout="wide")

# パラメータ設定（Streamlit用）
st.sidebar.title("パラメータ設定")
st.title(r'$B$と$H$に関するチュートリアル')
st.write('このページでは問の答えをシミュレーションを使って確認することができます。')
option = st.selectbox("問を選んでください", ["問2-3", "問2-5", "問4-1"])
a = st.sidebar.slider("半径 (a)", min_value=0.5, max_value=2.0, value=1.5, step=0.1)
L = st.sidebar.slider("長さ (L)", min_value=1.5, max_value=5.0, value=4.0, step=0.1)
N = st.sidebar.slider("本数 (N)", min_value=6, max_value=20, value=8, step=1)

# グリッド設定
Ngrid = 200
rho_vals = np.linspace(-5, 5, Ngrid)
z_vals = np.linspace(-5, 5, Ngrid)
Z, RHO = np.meshgrid(z_vals, rho_vals)


if option == "問2-3":
    # プロット
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # 長方形の座標
    rect_x = [-a, a, a, -a, -a]
    rect_y = [-L/2, -L/2, L/2, L/2, -L/2]
    
    axs[0].plot([L/2, L/2], [-a, a], 'r-', linewidth=4) 
    axs[0].plot([-L/2, -L/2], [-a, a], 'b-', linewidth=4) 
    
    axs[1].plot(rect_y, rect_x, 'k-', linewidth=2)
    axs[1].fill([0, 0, L/2, L/2], [-a, a, a, -a], color='red', alpha=0.2)
    axs[1].fill([-L/2, -L/2, 0, 0], [-a, a, a, -a], color='green', alpha=0.2)

    for i in range(1 + 1):
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
    
    stream_H = axs[0].streamplot(Z, RHO, H_z_vals, H_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_H)
    stream_H = axs[1].streamplot(Z, RHO, H_z_vals, H_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_H)
    
    buf = BytesIO()
    fig.savefig(buf, format="png")
    st.image(buf)


if option == "問2-5":
    # H_rho計算
    H_rho_vals = H_rho(RHO, Z, a, L)
    H_z_vals = H_z(RHO, Z, a, L)
    
    # Hの流線 (streamplot)
    start_points_H = getStartPointsH(a, L, N)
    plt.figure()
    stream_H = plt.gca().streamplot(Z, RHO, H_z_vals, H_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_H)
    plt.close('all')
    
    # プロット
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))

    # 長方形の座標
    rect_x = [-a, a, a, -a, -a]
    rect_y = [-L/2, -L/2, L/2, L/2, -L/2]

    axes.plot(rect_y, rect_x, 'k-', linewidth=2)
    axes.fill([0, 0, L/2, L/2], [-a, a, a, -a], color='red', alpha=0.2)
    axes.fill([-L/2, -L/2, 0, 0], [-a, a, a, -a], color='green', alpha=0.2)
    axes.set_xlim(-5, 5)
    axes.set_ylim(-5, 5)
    axes.set_aspect(1)
    axes.set_xticks([])
    axes.set_yticks([])
    
    # B計算
    B_rho_vals = B_rho(RHO, Z, a, L)
    B_z_vals = B_z(RHO, Z, a, L)

    # Bの流線 (streamplot)
    start_points_B = getStartPointsB(a, L, N, stream_H)
    axes.streamplot(Z, RHO, B_z_vals, B_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_B)
    
    buf = BytesIO()
    fig.savefig(buf, format="png")
    st.image(buf)


if option == "問4-1":
    # プロット
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # 長方形の座標
    rect_x = [-a, a, a, -a, -a]
    rect_y = [-L/2, -L/2, L/2, L/2, -L/2]

    for i in range(3):
        axs[i].plot(rect_y, rect_x, 'k-', linewidth=2)
        axs[i].fill([0, 0, L/2, L/2], [-a, a, a, -a], color='red', alpha=0.2)
        axs[i].fill([-L/2, -L/2, 0, 0], [-a, a, a, -a], color='green', alpha=0.2)
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
    
    stream_H = axs[0].streamplot(Z, RHO, H_z_vals, H_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_H)
    axs[0].set_title(r"$\mu _0H$", fontsize=20)
    
    # B計算
    B_rho_vals = B_rho(RHO, Z, a, L)
    B_z_vals = B_z(RHO, Z, a, L)

    # Bの流線 (streamplot)
    start_points_B = getStartPointsB(a, L, N, stream_H)
    axs[2].streamplot(Z, RHO, B_z_vals, B_rho_vals, color='k', density=20, arrowsize=2, start_points=start_points_B)
    axs[2].set_title(r"$B$", fontsize=20)

    # M計算
    M_rho_vals = M_rho(RHO, Z, a, L)
    M_z_vals = M_z(RHO, Z, a, L)

    # Mの流線 (streamplot)
    start_points_M = getStartPointsM(a, L, N)
    axs[1].streamplot(Z, RHO, M_z_vals, M_rho_vals, color='k', density=100, arrowsize=2, start_points=start_points_M)
    axs[1].set_title(r"$\mu _0M$", fontsize=20)
    
    buf = BytesIO()
    fig.savefig(buf, format="png")
    st.image(buf)