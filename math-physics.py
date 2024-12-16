import sys
import math
import numpy as np
from PyQt5.QtWidgets import (QApplication, QWidget, QHBoxLayout, QFormLayout, 
                             QDoubleSpinBox, QPushButton, QLabel, QSlider)
from PyQt5.QtCore import QTimer, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

class InvertedBowlSimulation(QWidget):
    def __init__(self):
        super().__init__()
        
        # Параметры по умолчанию
        self.m = 1.0        # масса шара
        self.r_ball = 0.05  # радиус шара
        self.R = 1.0        # радиус полусферы
        self.g = 9.81
        self.theta0_deg = 89.0   # почти на краю
        self.phi_dot0 = 0.1      # небольшой стартовый импульс по φ
        
        self.dt = 0.01      # шаг интегрирования
        self.speed_factor = 1.0  # множитель скорости симуляции
        
        # Коэффициенты трения (демпфирования)
        self.damping_theta = 0.2  # демпфирование по θ
        self.damping_phi = 0.2    # демпфирование по φ
        
        self.reset_sim_params()
        
        # Для отслеживания максимальной скорости
        self.max_speed = 0.0
        self.time_at_max_speed = 0.0
        
        # Состояние симуляции
        self.running = False
        
        self.init_ui()
        
    def init_ui(self):
        layout = QHBoxLayout(self)
        
        # Левая панель
        param_layout = QFormLayout()
        
        self.mass_spin = QDoubleSpinBox()
        self.mass_spin.setRange(0.1, 100.0)
        self.mass_spin.setValue(self.m)
        
        self.r_ball_spin = QDoubleSpinBox()
        self.r_ball_spin.setRange(0.001, 1.0)
        self.r_ball_spin.setValue(self.r_ball)
        
        self.R_spin = QDoubleSpinBox()
        self.R_spin.setRange(0.1, 10.0)
        self.R_spin.setValue(self.R)
        
        self.g_spin = QDoubleSpinBox()
        self.g_spin.setRange(1.0, 20.0)
        self.g_spin.setValue(self.g)
        
        self.theta0_spin = QDoubleSpinBox()
        self.theta0_spin.setRange(0.0,90.0)
        self.theta0_spin.setValue(self.theta0_deg)
        
        self.phi_dot0_spin = QDoubleSpinBox()
        self.phi_dot0_spin.setRange(-5.0,5.0)
        self.phi_dot0_spin.setValue(self.phi_dot0)
        
        # Слайдер для скорости симуляции
        self.speed_slider = QSlider(Qt.Horizontal)
        self.speed_slider.setRange(1, 100)  # от 1 до 100
        self.speed_slider.setValue(int(self.speed_factor*10)) 
        self.speed_slider.valueChanged.connect(self.on_speed_change)
        
        param_layout.addRow("Масса (кг):", self.mass_spin)
        param_layout.addRow("Радиус шара (м):", self.r_ball_spin)
        param_layout.addRow("Радиус чаши (м):", self.R_spin)
        param_layout.addRow("g (м/с²):", self.g_spin)
        param_layout.addRow("Начальный угол θ0 (°):", self.theta0_spin)
        param_layout.addRow("Нач. φ_dot (рад/с):", self.phi_dot0_spin)
        param_layout.addRow("Скорость симуляции:", self.speed_slider)
        
        self.start_btn = QPushButton("Старт")
        self.start_btn.clicked.connect(self.start_sim)
        self.stop_btn = QPushButton("Пауза")
        self.stop_btn.clicked.connect(self.stop_sim)
        self.reset_btn = QPushButton("Сброс")
        self.reset_btn.clicked.connect(self.reset_sim)
        
        param_layout.addRow(self.start_btn, self.stop_btn)
        param_layout.addRow(self.reset_btn)
        
        self.info_label = QLabel("Время: 0.00 c\nСкорость: 0.00 м/с\nМакс. скорость: 0.00 м/с при 0.00 c")
        param_layout.addRow(self.info_label)
        
        left_panel = QWidget()
        left_panel.setLayout(param_layout)
        
        # Правая панель
        self.fig = Figure(figsize=(5,5))
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        self.ax.set_xlabel("X (м)")
        self.ax.set_ylabel("Y (м)")
        self.ax.set_zlabel("Z (м)")
        
        layout.addWidget(left_panel)
        layout.addWidget(self.canvas)
        
        self.setLayout(layout)
        
        # Таймер
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_sim)
        
        self.draw_sphere()
        self.update_plot()
        
    def on_speed_change(self, value):
        self.speed_factor = value / 10.0
        if self.running:
            self.timer.start(int(self.dt*1000/self.speed_factor))
        
    def draw_sphere(self):
        self.ax.clear()
        self.ax.set_xlabel("X (м)")
        self.ax.set_ylabel("Y (м)")
        self.ax.set_zlabel("Z (м)")
        
        u = np.linspace(0, 2*np.pi, 30)
        v = np.linspace(0, np.pi/2, 30)
        x_sphere = self.R * np.outer(np.sin(v), np.cos(u))
        y_sphere = self.R * np.outer(np.sin(v), np.sin(u))
        z_sphere = -self.R * np.outer(np.cos(v), np.ones_like(u))
        self.ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='b', alpha=0.3)
        
        marker_size = 50*self.r_ball
        self.ball_plot, = self.ax.plot([0],[0],[0], 'ro', markersize=marker_size)
        
        self.ax.set_xlim(-self.R-0.2, self.R+0.2)
        self.ax.set_ylim(-self.R-0.2, self.R+0.2)
        self.ax.set_zlim(-self.R-0.2, 0.2)
    
    def start_sim(self):
        self.m = self.mass_spin.value()
        self.r_ball = self.r_ball_spin.value()
        self.R = self.R_spin.value()
        self.g = self.g_spin.value()
        self.theta0_deg = self.theta0_spin.value()
        self.phi_dot0 = self.phi_dot0_spin.value()
        
        self.reset_sim_params()
        self.draw_sphere()
        
        self.running = True
        self.timer.start(int(self.dt*1000/self.speed_factor))
        
    def stop_sim(self):
        self.running = False
        self.timer.stop()
        
    def reset_sim(self):
        self.stop_sim()
        self.reset_sim_params()
        self.draw_sphere()
        self.update_plot()
        
    def reset_sim_params(self):
        self.t = 0.0
        self.theta = math.radians(self.theta0_deg)
        self.theta_dot = 0.0
        self.phi = 0.0
        self.phi_dot = self.phi_dot0
        
        self.R_eff = self.R - self.r_ball
        
        # Изначально у нас был Lphi, но при трении угловой момент не сохраняется.
        # Поэтому просто работаем с phi_dot напрямую и будем демпфировать её.
        
        self.max_speed = 0.0
        self.time_at_max_speed = 0.0
        
    def update_sim(self):
        if not self.running:
            return
        
        sin_th = math.sin(self.theta)
        cos_th = math.cos(self.theta)
        
        # θ_ddot = sinθ cosθ φ_dot² - (5/7)*(g/R_eff)*sinθ
        # Добавим демпфирование по θ: θ_ddot -= damping_theta*θ_dot
        if sin_th != 0:
            # phi_dot не константен. Добавим демпфирование по φ:
            self.phi_dot = self.phi_dot * (1 - self.damping_phi*self.dt)
        else:
            # Если sin_th=0, шар в самом дне, просто затухаем φ_dot
            self.phi_dot *= (1 - self.damping_phi*self.dt)
        
        theta_ddot = sin_th*cos_th*(self.phi_dot**2) - (5.0/7.0)*(self.g/self.R_eff)*sin_th
        theta_ddot -= self.damping_theta*self.theta_dot  # демпфирование θ
        
        # Интеграция θ по РК2
        def f_theta(th, th_dot, phi_d):
            s = math.sin(th)
            c = math.cos(th)
            th_dd = s*c*(phi_d**2) - (5.0/7.0)*(self.g/self.R_eff)*s - self.damping_theta*th_dot
            return (th_dot, th_dd)
        
        k1 = f_theta(self.theta, self.theta_dot, self.phi_dot)
        th_mid = self.theta + k1[0]*self.dt*0.5
        th_dot_mid = self.theta_dot + k1[1]*self.dt*0.5
        k2 = f_theta(th_mid, th_dot_mid, self.phi_dot)
        
        self.theta += k2[0]*self.dt
        self.theta_dot += k2[1]*self.dt
        
        # Интегрируем φ (просто φ += φ_dot*dt)
        self.phi += self.phi_dot * self.dt
        self.t += self.dt
        
        # Скорость шара
        v = self.R_eff * math.sqrt(self.theta_dot**2 + (sin_th**2)*(self.phi_dot**2))
        
        if v > self.max_speed:
            self.max_speed = v
            self.time_at_max_speed = self.t
        
        self.info_label.setText(
            f"Время: {self.t:.3f} c\n"
            f"Скорость: {v:.3f} м/с\n"
            f"Макс. скорость: {self.max_speed:.3f} м/с при {self.time_at_max_speed:.3f} c"
        )
        
        # Условие остановки: если θ около 0 и скорости малы
        if self.theta < 0.01 and abs(self.theta_dot) < 0.001 and abs(self.phi_dot) < 0.001:
            self.stop_sim()
        
        self.update_plot()
        
    def update_plot(self):
        x = self.R_eff * math.sin(self.theta)*math.cos(self.phi)
        y = self.R_eff * math.sin(self.theta)*math.sin(self.phi)
        z = -self.R_eff * math.cos(self.theta)
        
        self.ball_plot.set_data([x],[y])
        self.ball_plot.set_3d_properties([z])
        
        self.ax.set_xlim(-self.R-0.2, self.R+0.2)
        self.ax.set_ylim(-self.R-0.2, self.R+0.2)
        self.ax.set_zlim(-self.R-0.2, 0.2)
        
        self.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = InvertedBowlSimulation()
    w.setWindowTitle("3D Симуляция с трением")
    w.show()
    sys.exit(app.exec_())