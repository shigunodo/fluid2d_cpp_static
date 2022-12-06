#2次元での流体計算
#使用モジュール：
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gc


class Fluid2d:
    #------------初期化メソッド---------------------
    #引数：
    #Ni=第一方向のグリッド数
    #Nj=第二方向のグリッド数
    #Nb=境界条件（外挿または固定）用の列数（Ni,Nj に含まれる、デフォでは境界の外側にもう一列ダミーを用意）
    #KH不安定のための等分割グリッドを生成 & 初期条件を構成
    #-----------------------------------------------
    def __init__(self, Ni, Nj, Nb=4, gamma=1.4):
        #
        #計算空間でのインデックス
        self.Nb = Nb
        self.Ni = Ni
        self.Nj = Nj
        self.gam = gamma
        i = np.tile(np.arange(Ni, dtype=int), (Nj,1)).T
        j = np.tile(np.arange(Nj, dtype=int), (Ni,1))
        self.index = np.stack([i,j],-1)
        #
        #物理空間座標（計算空間と同じに初期化）
        i = np.tile(np.arange(Ni, dtype=float), (Nj,1)).T
        j = np.tile(np.arange(Nj, dtype=float), (Ni,1))
        self.x = np.stack([i,j],-1)
        # 基本変数
        self.q = np.zeros((Ni,Nj,4))
        # パラメータ
        rhoU = 0.5 # rhoU < rhoD
        rhoD = 1.
        uU   = 1.
        uD   = -1.
        eD   = 1.
        Ma   = 0.3
        p    = 1./self.gam/Ma/Ma
        eU   = p/(self.gam-1.)+0.5*rhoU*uU*uU
        eD   = p/(self.gam-1.)+0.5*rhoD*uD*uD
        # 格子生成と初期化
        dx = 1./(Ni-2*Nb)
        dy = 1./(Ni-2*Nb)
        for i in range(Ni):
            for j in range(Nj):
                self.x[i,j,0] = 0.+(i-Nb)*dx
                self.x[i,j,1] = -0.5 + (j-Nb)*dy
        N0 = 0
        for j in range(Nj):
            if self.x[0,j,1]<0.:
                N0 = N0+1
        self.q[:,:N0,0] = rhoD
        self.q[:,:N0,1] = uD
        self.q[:,:N0,2] = 0.
        self.q[:,:N0,3] = eD
        self.q[:,N0:,0] = rhoU
        self.q[:,N0:,1] = uU
        self.q[:,N0:,2] = 0.
        self.q[:,N0:,3] = eU
    #
    ############## Fortran I/O ############################
    def input_basic_fort(self, name):
        data = np.loadtxt(name)
        data = data.reshape((4,self.Nj,self.Ni))
        self.q[:,:,0] = data[0,:,:].T
        self.q[:,:,1] = data[1,:,:].T
        self.q[:,:,2] = data[2,:,:].T
        self.q[:,:,3] = data[3,:,:].T
    def output_basic_fort(self, name):
        rho = self.q[:,:,0]
        u   = self.q[:,:,1]
        v   = self.q[:,:,2]
        e   = self.q[:,:,3]
        np.savetxt(name, np.array([rho.T,u.T,v.T,e.T]).reshape(4*self.Ni*self.Nj))
    def output_coordinate_fort(self, name):
        x = self.x[:,:,0]
        y = self.x[:,:,1]
        np.savetxt(name, np.array([x.T,y.T]).reshape(2*self.Ni*self.Nj))
    def input_coordinate_fort(self, name):
        data = np.loadtxt(name)
        data = data.reshape((2,self.Nj,self.Ni))
        self.x[:,:,0] = data[0,:,:].T
        self.x[:,:,1] = data[1,:,:].T
        #ヤコビアンの逆とメトリクス
        self.S = self.calc_S()
        self.met = self.calc_metrix()
    def adjust_basic_fort(self, fname_of_old_coordinate, Ni_old,Nj_old, fname_of_basic):
        c_old = np.loadtxt(fname_of_old_coordinate).reshape((2,Nj_old*Ni_old))
        x_old = c_old[0,:]
        y_old = c_old[1,:]
        x_new = self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,0]
        y_new = self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,1]
        q_old = np.loadtxt(fname_of_basic).reshape(4,Nj_old*Ni_old)
        points = (x_old,y_old)
        for i in range(4):
            self.q[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,i] = griddata(points, q_old[i,:], xi=(x_new,y_new), method='cubic')
    ############## C++ I/O ############################
    def input_basic_cpp(self, name):
        data = np.loadtxt(name)
        data = data.reshape((4,self.Ni,self.Nj))
        self.q[:,:,0] = data[0,:,:]
        self.q[:,:,1] = data[1,:,:]
        self.q[:,:,2] = data[2,:,:]
        self.q[:,:,3] = data[3,:,:]
    def output_basic_cpp(self, name):
        rho = self.q[:,:,0]
        u   = self.q[:,:,1]
        v   = self.q[:,:,2]
        e   = self.q[:,:,3]
        np.savetxt(name, np.array([rho,u,v,e]).reshape(4*self.Ni*self.Nj))
    def output_coordinate_cpp(self, name):
        x = self.x[:,:,0]
        y = self.x[:,:,1]
        np.savetxt(name, np.array([x,y]).reshape(2*self.Ni*self.Nj))
    def input_coordinate_cpp(self, name):
        data = np.loadtxt(name)
        data = data.reshape((2,self.Ni,self.Nj))
        self.x[:,:,0] = data[0,:,:]
        self.x[:,:,1] = data[1,:,:]
    def adjust_basic_cpp(self, fname_of_old_coordinate, Ni_old,Nj_old, fname_of_basic):
        c_old = np.loadtxt(fname_of_old_coordinate).reshape((2,Nj_old*Ni_old))
        x_old = c_old[0,:]
        y_old = c_old[1,:]
        x_new = self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,0]
        y_new = self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,1]
        q_old = np.loadtxt(fname_of_basic).reshape(4,Nj_old*Ni_old)
        points = (x_old,y_old)
        for i in range(4):
            self.q[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,i] = griddata(points, q_old[i,:], xi=(x_new,y_new), method='cubic')
    #######################PLOT##########################
    def plot_rho(self, fname, time):
        # 解の図
        fig, ax1 = plt.subplots(figsize=(6,6),dpi=150, facecolor="w")

        #タイトル
        fig.subplots_adjust(hspace=0.2,left=0.05, right=0.9)
        title = "Kelvin-Helmholtz Instability \n Euler eq./Roe's FDS/MP5/SSPRK3 \n t = {:3.1f}".format(time)
        fig.suptitle(title)

        map1 = ax1.pcolormesh(self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,0], self.x[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,1], self.q[self.Nb:self.Ni-self.Nb,self.Nb:self.Nj-self.Nb,0], cmap="inferno",vmax=1.,vmin=0.5, zorder=1)
        #カラーマップ
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="2%", pad=0.1)
        plt.colorbar(map1, cax=cax1,label='density')


        #図の設定
        ax1.spines['left'].set(lw=1.5)
        ax1.spines['bottom'].set(lw=1.5)
        ax1.spines['right'].set(lw=1.5)
        ax1.spines['top'].set(lw=1.5)
        ax1.set_xlim(0,1)
        ax1.set_aspect('equal')

        #図の保存
        plt.savefig(fname, bbox_inches = 'tight',dpi=150, transparent=False)
        plt.clf()
        plt.close()
        gc.collect()