#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from datetime import datetime
import argparse
import simulate


pi = 3.14
temp_s = 0
R = 52
alpha_r = 5.2
alpha_m = 233
a = 0.2
epsilon = 0.9
kr = 0.7
kp = 0.38
ks = 0.093
ri = 0.0175
ro = 0.022
cw = 4184
Le = 2267000
Lm = 3333500
ro_i = 917
ro_w = 1000
vw = 10
Qp = 300



class control():

	# snowfall sensor
	def judge_1(self, F):
		if(float(F)>0):	switch = 'on'
		else:			switch = 'off'
		return switch



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='control mode')
	parser.add_argument('mode')
	args = parser.parse_args()


	simulate = simulate.simulate()

	num = 0
	weather = simulate.get_weather(num)
	print(weather)


	datetime = datetime.strptime(weather[0], "%m-%d %H:%M")
	F = float(weather[1])
	temp_o = float(weather[2])
	RH = float(weather[3])
	I = float(weather[4])
	V = float(weather[5])
	P = float(weather[6].split('\n')[0])

	on_t = 0
	off_t = 0


	Xo = simulate.calc_Xo(RH, temp_o, P)

	Xr = simulate.calc_Xr(temp_o, P)
	Xs = simulate.calc_Xs(temp_o, P)
	ro_f = simulate.calc_ro_f(temp_o)
	alpha_c = simulate.calc_alpha_c(V)
	alpha_x = simulate.calc_alpha_x(V)
	alpha_i = simulate.calc_alpha_i()
	rs = simulate.calc_rs(temp_o)
	temp_w = simulate.calc_temp_w(on_t)

	# 相対外気温
	temp_e = simulate.calc_temp_e(temp_o, I, alpha_c)

	# 融雪量
	M = 0

	# 雪氷量
	S = 0

	# 水分の路面からの浸透高さ
	dw = 0

	# 雪面の路面からの高さ
	ds = 0

	# 路面からの蒸発量
	Er = 0

	# 路面から外界への潜熱伝熱量
	Qer = 0

	# 雪から外界への顕熱伝熱量
	Qas = 0

	# 雪面からの蒸発量
	Es = 0

	# 雪から外界への潜熱伝熱量
	Qes = 0

	# 融雪熱量
	Qm = 0

	# 路面から相変化中の雪への伝熱量
	Qs = 0

	# 路面および路盤内部温度
	temp_r = simulate.calc_temp_r(Qs, alpha_m)

	# 路面から外界への顕熱伝熱量
	Qar = 0

	# 路面のうち雪氷あるいは水分の相変化に寄与する面積の割合
	f = 0
	E = 0

	# 路面放熱量
	Qr = simulate.calc_Qr(f, Qs, Qar, Qer)

	# 水分量
	W = 0


	if(args.mode == str(1)):
		switch = control().judge_1(F)
	else:
		switch = 'off'
	print(switch)


	b_time = datetime
	S_1 = S
	W_1 = W

	simulate.result_output(datetime, M, S, ds, Qm, temp_r, W)


#	roopnum = len(simulate.weathers)
	roopnum = 3
	for n in range(roopnum-1):

		num += 1
		weather = simulate.get_weather(num)
		print(weather)


		datetime = datetime.strptime(weather[0], "%m-%d %H:%M")
		F = float(weather[1])
		temp_o = float(weather[2])
		RH = float(weather[3])
		I = float(weather[4])
		V = float(weather[5])
		P = float(weather[6].split('\n')[0])

		elapsed_time = datetime - b_time
		elapsed_t = float(elapsed_time.seconds/60)

		if(switch=='on'):
			on_t += elapsed_t
			off_t = 0
		else:
			on_t = 0
			off_t += elapsed_t


		Xo = simulate.calc_Xo(RH, temp_o, P)

		Xr = simulate.calc_Xr(temp_o, P)
		Xs = simulate.calc_Xs(temp_o, P)
		ro_f = simulate.calc_ro_f(temp_o)
		alpha_c = simulate.calc_alpha_c(V)
		alpha_x = simulate.calc_alpha_x(V)
		alpha_i = simulate.calc_alpha_i()
		rs = simulate.calc_rs(temp_o)
		temp_w = simulate.calc_temp_w(on_t)

		# 相対外気温
		temp_e = simulate.calc_temp_e(temp_o, I, alpha_c)

		# 融雪量
		M = simulate.calc_M(S_1, elapsed_t, rs, F)

		# 雪氷量
		S = simulate.calc_S(S_1, rs, F, M, elapsed_t)

		if(S>0.0):
			# 積雪の密度
			ro_s = simulate.calc_ro_s(temp_o, on_t, S_1, F, M, elapsed_t, ro_s_1, ro_f)

		# 水分の路面からの浸透高さ
		dw = simulate.calc_dw()

		if(S>0.0):
			# 雪面の路面からの高さ
			ds = simulate.calc_ds(S, ro_s)

		# 路面からの蒸発量
		Er = simulate.calc_Er(alpha_x, Xr, Xo, S)

		# 路面から外界への潜熱伝熱量
		Qer = simulate.calc_Qer(Er)

		if(S>0.0):
			# 雪から外界への顕熱伝熱量
			Qas = simulate.calc_Qas(alpha_c, ds, dw, temp_e)

			# 雪面からの蒸発量
			Es = simulate.calc_Es(alpha_x, Xs, Xo)

			# 雪から外界への潜熱伝熱量
			Qes = simulate.calc_Qes(Es)

		# 融雪熱量
		Qm = simulate.calc_Qm(M)

		if(S>0.0):
			# 路面から相変化中の雪への伝熱量
			Qs = simulate.calc_Qs(Qas, Qes, Qm)

		# 路面および路盤内部温度
		temp_r = simulate.calc_temp_r(Qs, alpha_m)

		# 路面から外界への顕熱伝熱量
		Qar = simulate.calc_Qar(temp_r, temp_e, alpha_c, ds)

		# 路面のうち雪氷あるいは水分の相変化に寄与する面積の割合
		f, M = simulate.calc_f(M, W_1, elapsed_t, rs, F, temp_r, temp_e, Er, Es)

		# 路面放熱量
		Qr = simulate.calc_Qr(f, Qs, Qar, Qer)

		# 水分量
		W = simulate.calc_W(W_1, rs, F, M, E, elapsed_t)


		if(args.mode == '1'):
			switch = control().judge_1(F)
		else:
			switch = 'off'
		print(switch)


		b_time = datetime
		S_1 = S
		W_1 = W

		simulate.result_output(datetime, M, S, ds, Qm, temp_r, W)
