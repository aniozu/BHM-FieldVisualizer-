import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.special import elliprj

# Define constants
mu0 = 4 * np.pi * 1e-7  # 真空の透磁率
I = 1.0  # 電流

def ellippi(n, m):
  return sp.elliprf(                                                                                                                                                                                        
        0., 1. - m, 1.) + (n / 3.) * sp.elliprj(0., 1. - m, 1., 1. - n)

def f_rho(rho, zeta, a, L):
    k2 = 4 * a * rho / ((a + rho)**2 + zeta**2)
    k = np.sqrt(k2)
    Kk = sp.ellipk(k2)
    Ek = sp.ellipe(k2)
    return (k2 - 2) / k * Kk + 2 / k * Ek
    
# 磁場B成分の定義
def B_rho(rho, z, a, L):
    rho_sign = np.sign(rho)
    rho = np.abs(rho)
    zeta_p = z + L / 2
    zeta_m = z - L / 2
    factor = (mu0 * I) / (2 * np.pi * L) * np.sqrt(a / rho)
    return factor * (f_rho(rho, zeta_p, a, L) - f_rho(rho, zeta_m, a, L)) * rho_sign

def f_z(rho, zeta, a, L):
    k2 = 4 * a * rho / ((a + rho)**2 + zeta**2)
    k = np.sqrt(k2)
    Kk = sp.ellipk(k2)
    h2 = 4 * a * rho / (a + rho)**2
    h = np.sqrt(h2)
    #Pik = ellippi_array(h2, k2)
    Pik = ellippi(h2, k2)
    return zeta * k * (Kk + (a - rho) / (a + rho) * Pik)

def B_z(rho, z, a, L):
    rho = np.abs(rho)
    zeta_p = z + L / 2
    zeta_m = z - L / 2
    factor = (mu0 * I) / (2 * np.pi * L) * 1 / (2 * np.sqrt(a * rho))
    return factor * (f_z(rho, zeta_p, a, L) - f_z(rho, zeta_m, a, L))

# 磁場mu0H動径成分の定義
def H_rho(rho, z, a, L):
    rho_sign = np.sign(rho)
    rho = np.abs(rho)
    zeta_p = z + L / 2
    zeta_m = z - L / 2
    factor = (mu0 * I) / (2 * np.pi * L) * np.sqrt(a / rho)
    return factor * (f_rho(rho, zeta_p, a, L) - f_rho(rho, zeta_m, a, L)) * rho_sign

def H_z(rho, z, a, L):
    rho = np.abs(rho)
    zeta_p = z + L / 2
    zeta_m = z - L / 2
    factor = (mu0 * I) / (2 * np.pi * L)
    return factor * (1 / (2 * np.sqrt(a * rho)) * (f_z(rho, zeta_p, a, L) - f_z(rho, zeta_m, a, L)) - 2 * np.pi * np.where(np.logical_and(rho < a, np.abs(z) < L/2), 1, 0)) 

def M_z(rho, z, a, L):
    rho = np.abs(rho)
    factor = (mu0 * I) / (2 * np.pi * L)
    return factor * 2 * np.pi * np.where(np.logical_and(np.abs(rho) < a, np.abs(z) < L/2), 1, 0)

def M_rho(rho, z, a, L):
    return np.zeros_like(rho)
 
def calcStartRhos(a, L, N, showFig = False):
    rho_arr = np.linspace(-a, a, 1000 + 2)[1:-1]# rho の範囲（-a, aを踏まない）
    
    # Hz_out, Hz_in の計算
    Hz_out = H_z(rho_arr, np.full_like(rho_arr, L/2 + 0.000001), a, L)
    Hz_in = H_z(rho_arr, np.full_like(rho_arr, L/2 - 0.000001), a, L)

    if showFig:
        plt.figure()
        plt.plot(rho_arr, Hz_out, label="Hz_out")
        plt.plot(rho_arr, Hz_in, label="Hz_in")
        plt.plot(rho_arr, Hz_out + np.abs(Hz_in), label="Hz_out + |Hz_in|")
        plt.xlabel("rho")
        plt.ylabel("Hz")
        plt.axvline(-a, linestyle='--', color='gray')
        plt.axvline(a, linestyle='--', color='gray')
        plt.legend()
        plt.show()
    
    # Hz_out, Hz_in の絶対値総和
    sum_abs_Hz_out = np.nansum(np.abs(Hz_out))
    sum_abs_Hz_in = np.nansum(np.abs(Hz_in))
    
    print(sum_abs_Hz_out)
    print(sum_abs_Hz_in)
    
    # 合計
    total_sum = sum_abs_Hz_out + sum_abs_Hz_in
    
    # N を Hz_out と Hz_in の絶対値総和の比に応じて分割
    Nout = int(round(N * (sum_abs_Hz_out / total_sum))) if total_sum > 0 else N // 2
    Nin = N - Nout
    
    print(f"Nout: {Nout}, Nin: {Nin}")
    
    # 逆数を計算（ゼロと NaN を防ぐ）
    epsilon = 0#1e-10  # ゼロ除算防止用の微小値
    Hz_out_inv = np.where(Hz_out != 0, 1 / (np.abs(Hz_out) + epsilon), 0)
    Hz_in_inv = np.where(Hz_in != 0, 1 / (np.abs(Hz_in) + epsilon), 0)
    
    # 逆数の積分値（累積分布を作成）
    cum_weights_out = np.cumsum(Hz_out_inv)
    cum_weights_in = np.cumsum(Hz_in_inv)
    
    # 正規化（ゼロ除算を防ぐ）
    if cum_weights_out[-1] > 0:
        cum_weights_out /= cum_weights_out[-1]
    else:
        print("警告: cum_weights_out[-1] がゼロのため、正規化できません。")
        cum_weights_out[:] = np.linspace(0, 1, len(cum_weights_out))  # 線形補間で代用
    
    if cum_weights_in[-1] > 0:
        cum_weights_in /= cum_weights_in[-1]
    else:
        print("警告: cum_weights_in[-1] がゼロのため、正規化できません。")
        cum_weights_in[:] = np.linspace(0, 1, len(cum_weights_in))  # 線形補間で代用
    
    # 区間エッジの計算
    edges_out = np.interp(np.linspace(0, 1, Nout + 1), cum_weights_out, rho_arr)
    edges_in = np.interp(np.linspace(0, 1, Nin + 1), cum_weights_in, rho_arr)
    
    # 結果の出力
    print("範囲のエッジ (Hz_out):", edges_out)
    print("範囲のエッジ (Hz_in):", edges_in)
    
    # 各区間の中心を求める
    midpoints_out = (edges_out[:-1] + edges_out[1:]) / 2
    midpoints_in = (edges_in[:-1] + edges_in[1:]) / 2
    return midpoints_out, midpoints_in

def getStartPointsH(a, L, N):
    midpoints_out, midpoints_in = calcStartRhos(a, L, N)
    
    # ストリームラインの開始点（in/out, positive/negative の4種類）
    start_rho_out = midpoints_out
    start_rho_in = midpoints_in
    
    delta = 0.03
    start_z_pos_out = np.full(len(midpoints_out), L / 2 + delta)
    start_z_neg_out = np.full(len(midpoints_out), -L / 2 - delta)
    start_z_pos_in = np.full(len(midpoints_in), L / 2 - delta)
    # start_z_neg_in = np.full(len(midpoints_in), -L / 2 + delta)
    
    start_points_out_pos = np.vstack((start_z_pos_out, start_rho_out)).T
    start_points_out_neg = np.vstack((start_z_neg_out, start_rho_out)).T
    start_points_in_pos = np.vstack((start_z_pos_in, start_rho_in)).T
    # start_points_in_neg = np.vstack((start_rho_in, start_z_neg_in)).T
    
    # 全ての開始点を統合
    start_points = np.vstack((start_points_out_pos, start_points_out_neg, start_points_in_pos))
    return start_points

def find_nearest_zero_crossing_points(a, L, N, stream_H):
    """
    rho = ±a を横切る流線のうち、-L/2 < z < L/2 の範囲で横切るものを対象に、
    それぞれの流線上で z=0 に最も近い点を求める。

    Parameters:
    - stream_H: plt.streamplot の戻り値（StreamplotSet）
    - a (float): 横切るべき rho の絶対値
    - L (float): z の範囲制限 (-L/2 < z < L/2)

    Returns:
    - list of tuples: 各流線の (rho, z) で z=0 に最も近い点のリスト
    """
    crossing_points = []

    # すべての流線について処理
    for path in stream_H.lines.get_paths():
        verts = path.vertices  # (N, 2) の座標リスト
        rho_vals, z_vals = verts[:, 1], verts[:, 0]

        # rho = ±a を横切るか判定
        crosses_a = False
        for i in range(len(rho_vals) - 1):
            rho1, rho2 = rho_vals[i], rho_vals[i + 1]
            z1, z2 = z_vals[i], z_vals[i + 1]

            # rho = ±a を横切る & -L/2 < z < L/2 の範囲内
            if (-L/2 < z1 < L/2 and -L/2 < z2 < L/2) and \
               ((rho1 - a) * (rho2 - a) < 0 or (rho1 + a) * (rho2 + a) < 0):
                crosses_a = True
                break  # 1つでも横切ればその流線を対象とする

        if crosses_a:
            # z=0 に最も近い点を探す
            min_z_idx = np.argmin(np.abs(z_vals))  # z の絶対値が最小のインデックス
            nearest_point = (rho_vals[min_z_idx], z_vals[min_z_idx])
            crossing_points.append(nearest_point)
    
    try:
        return np.array(crossing_points)[:,[1,0]]
    except:
        return None
        
def getStartPointsB(a, L, N, stream_H):
    midpoints_out, midpoints_in = calcStartRhos(a, L, N)
    
    # ストリームラインの開始点
    start_rho = midpoints_out
    start_z_pos = np.full(len(midpoints_out), L / 2)
    start_points = np.vstack((start_z_pos, start_rho)).T
    
    addpoints = find_nearest_zero_crossing_points(a, L, N, stream_H)
    
    # 既存の start_points に追加
    if addpoints is not None:
        start_points = np.vstack((start_points, addpoints))
    return start_points

def getStartPointsM(a, L, N):
    # ストリームラインの開始点
    edges = np.linspace(-a, a, N+1)
    start_rho = (edges[:-1] + edges[1:]) / 2
    # start_z_pos = np.full(N, L / 2)
    start_z_neg = np.full(N, -L / 2)
    start_points = np.vstack((start_z_neg, start_rho)).T
    return start_points