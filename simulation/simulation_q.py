#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
import re
import math
import time
import datetime
import argparse
import os
import os.path
import linecache
import sys

import q_control


interval = 3			# [min] calculatioin interval


CK      = 1.0		# 緩和係数(0.7~1.5程度)
maxT    = 70		# [℃ ] maximum temperature of heat source
maxPene = 0.05		# [m] maximum penetration height of water
abrate0 = 0.2		# solar radiatioin absorption rate of snow(0.2)
windC   = 0.5		# wind speed correction coefficient(~1)
BT      = 8.0		# temperature?
C       = 0.052629437
area    = 0.02783305
Hfusion = 80.0		# heat of fusion

fn1 = open('smap.run', 'r')
smap = fn1.readlines()[4].split('\t')
Qs   = float(smap[0])		# [kcal/h] Q of heat source
fn1.close()

fn2 = open('file.net', 'r')
net = re.split(" +", fn2.readlines()[4])
npn = int(net[1])
HR = []
for j in range(npn):
	target_line = linecache.getline('file.net', 6+j)
	netj = re.split( " +", target_line)
	HR.append(float(netj[3]))
fn2.close()
linecache.clearcache()


class sim:

	def __init__(self):
		self.net  = open('file.net', 'r')
		self.logf = open('q_logs_'+str(MODE)+'.txt', 'a')

		self.control = q_control.control()


	# absolute humidity
	def abshumid(self, Temp):
		P = 51.7063 * 10**(6.21147-2886.37/(1.8*Temp+491.69)\
							- 337269.46/(1.8*Temp+491.69)**2)
		funx = 0.622*P / (760-P)
		return funx


	# 対流熱伝達率convection heat transfer coefficient
	def funa(self, Wspeed):
		funa = 5 + 3.4*Wspeed
		return funa


	# absolute humidity, nighttime radiation
	def absoluteHumid(self, vaporP, cloud, temp_o, nightR):
		# when there is water vapor in the air
		if(vaporP >= 0):
			absH = 0.622*vaporP/(1013.25-vaporP)	# [kg/kg]
				
			# nighttime radiation when there is a cloud
			if(cloud >= 0):
				RH = math.sqrt(vaporP*760/1013.25)
				nightR = 0.0000000488 * (temp_o+273.16)**4 \
						* (1-0.62*cloud/10) * (0.49-0.076*RH)	# [W/m^2]

		# when no vapor in the air
		else:
			absH = 0.7 * sim().abshumid(temp_o)		# [kg/kg]

		return absH, nightR


	# density of snowfall (max:50) [kg/m^3]
	def snowfall_density(self, temp_o):
		sfdens = 1000 / (0.091*temp_o**2 - 1.81*temp_o + 9.47)
		if(sfdens < 50):	sfdens = 50
		return sfdens


	# amount of snow and rain in precipitation
	def calc_plus(self, temp_o, pre):
		if(temp_o < 0):
			snow_plus = pre						# [kg/m^2/min]
			rain_plus = 0						# [kg/m^2/min]
		elif(temp_o>=0 and temp_o<2):
			snow_plus = pre * (2-temp_o)/2		# [kg/m^2/min]
			rain_plus = pre * temp_o / 2		# [kg/m^2/min]
		elif(temp_o>=2):
			snow_plus = 0						# [kg/m^2/min]
			rain_plus = pre						# [kg/m^2/min]
		return snow_plus, rain_plus


	# penetration height
	def penetration_height(self, cover, snow, Water):
		peneH = 0
		# when there is snow cover
		if(cover > 0):
			dens  = snow / cover		# [kg/m^2] / [m] -> [kg/m^3]
			peneH = Water / (1000-dens)	# [m]
		# when no snow cover
		else:
			peneH = Water / 1000		# [m]

		if(peneH > maxPene):
			peneH = maxPene				# [m]

		return peneH



class QL:
	
	def __init__(self):
		alpha = 0.1
		gamma = 0.99
		self.Qlearn = q_control.Qlearning(MODE, alpha, gamma, interval)


	def next_Tlevel(self, all_data, shift, data_cnt):
		if(shift):
			nextTemp = float(all_data[data_cnt+1].split(', ')[5])
		else:
			nextTemp = float(all_data[data_cnt].split(', ')[5])

		nextTLv = self.Qlearn.Tlevel(nextTemp)
		return nextTLv


	def next_Plevel(self, all_data, shift, data_cnt):
		if(shift):
			nextPre = float(all_data[data_cnt+1].split(', ')[9]) / 10
		else:
			nextPre = float(all_data[data_cnt].split(', ')[9]) / 10

		nextPLv = self.Qlearn.Plevel(nextPre)
		return nextPLv


	def next_Slevel(self, act, cover, temp_o, nightR, Wspeed, Water, \
					rain_plus, snow, snow_plus, sfdens, tset, ilp, Qr):
		ict = 0
		E   = 0

		abrate = 0.8 - 30*cover
		if(abrate < abrate0):	abrate = abrate0

		sat = temp_o + (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4)

		# total water
		wat = Water + rain_plus

		mlt = 0
		TS = BT
		# 路面に雪があるとき
		if( (snow+snow_plus)>0 ):
			if( TS>0 or sat>0 ):
				mlt = 1
		# 路面に水があるとき
		if(wat > 0):
			if(TS < 0):
				mlt = 1

		# penetration height [m]
		peneH = sim().penetration_height(cover, snow, Water)

		if(act==1):
			Q = Qs
		else:
			Q = 0

		melt = 0
		evaporate = 0
		BF = 0
		if(mlt==1):
			DH = cover - peneH
			if(DH<=0 or sat>0):
				DH  = 0
				tsv = 0
				# evaporation amount
				evaporate = 4 * sim().funa(Wspeed) * (sim().abshumid(tsv)-absH)
				# total water [kg/m^2]
				wat = Water + rain_plus - evaporate*(interval/60.0)
				# abnormal value correction
				if(wat < 0):
					evaporate = Water/(interval/60.0) + rain_plus
					wat = 0
			htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)

			# calc amount of snow melting
			melt = (200*TS + htrm*sat - 590*evaporate) \
						* (interval/60.0) / (Hfusion*180)

			BF = 1

			if(melt < -1*wat):
				melt = -1 * wat
			elif( melt > (snow+snow_plus) ):
				melt = snow + snow_plus
			if(melt < 0):
				melt = 0

		else:
			tsv = BT
			evaporate = 4 * sim().funa(Wspeed) * (sim().abshumid(tsv) - absH)
			if( evaporate > (Water/(interval/60.0) + rain_plus) ):
				evaporate = Water/(interval/60.0) + rain_plus

		htr  = 1 / (1/(sim().funa(Wspeed)+4) + cover/0.08)
		Qe   = -590*evaporate*(interval/60.0) * area * (1-BF)
		Snow = snow - melt + snow_plus

		# snow accumulation [m]
		if(Snow > 0):
			if(melt > 0):
				Scover = cover + snow_plus/sfdens - melt/916
			else:
				Scover = cover + snow_plus/sfdens
		else:
			Scover = 0

		T = BT

		while(True):
			if(ict==1):
				T = tset
				E = 0
			S1 = (Q+Qe+E)*(interval/60.0) + BT*C
			S2 = C

			for j in range(npn):
				S1 += HR[j] * T * (interval/60.0)
				S2 += HR[j] * (interval/60.0)

			tmp = (1-BF) * htr * sat
			S1 += tmp * (interval/60.0) * area
			S2 += (BF*200 + (1-BF)*htr) * (interval/60.0) * area
			if(ict==1):
				E = (S2*tset - S1) / (interval/60.0)
			else:
				ER = S1/S2 - T
				T = CK*ER + T

			if(ilp < 2):
				if( T>maxT and Q>0 ):
					ilp += 1
					tset = maxT
					ict = 1
					Q = 0
					continue
				if(E < 0):
					ilp += 1
					ict = 0
					E = 0
					continue
				else:	break
			else:	break

		# abnormal valu correction
		if(Snow < 0):
			Snow = 0

		if( Scover>0 and Snow>0 ):
			gm = Snow / Scover
			ee = 16 * math.exp(0.021*gm)
			gm *= math.exp( Snow/2/ee*(interval/60.0)/24 )
			if(gm > 916):	gm = 916
			Scover = Snow / gm


		if( pre>0.0 and T>0 ):
			TS = T
			Qr_plus = BF*200*TS + (1-BF)*htr*(TS-sat)
			Qr += Qr_plus * area

		snow = Snow

		nextSLv = self.Qlearn.Slevel(snow)
		return nextSLv



if __name__ == '__main__':

	start = time.time()
	print('start time :', time.ctime())

	parser = argparse.ArgumentParser(description='MODE')
	parser.add_argument('MODE', help='S, P, ST, SP' )
	parser.add_argument('loop_num', default=1000)
	args = parser.parse_args()

	global MODE, Qloop
	Qloop = int(args.loop_num)
	MODE  = args.MODE
	# MODE S : snow accumulation
	# MODE T : temperature
	# MODE P : precipitation
	if(MODE!='S' and MODE!='P' and MODE!='ST' and MODE!='SP'):
		print('invalid MODE ERROR')
		sys.exit()
	elif( Qloop < 1 ):
		print('NO LOOP')
		sys.exit()

	snow_minusT = 0
	wet_minusT  = 0

	print('interval :', interval, '[min]')
	print('  MODE   :', MODE)

	# initialize
	temp_o = 0		# temperature outside
	Qsup   = 0		# supplied Q
	onT    = 0		# operating time
	heater = 0		# on(1) / off(0)
	noPreT = 0		# no precipitation time
	lpmx   = 0
	ntime  = np.zeros((2, 3))
	Qr     = 0		# Q from surface(snowing, more than 0℃ )
	Qe     = 0
	BF     = 0
	Water  = 0		# [kg/m^2]
	Snow   = 0
	ww     = 0		# [kg/m^2]
	htr    = 0
	sat    = 0		# 相当外気温 Sol-Air Temperature

	snow   = 0.0	# [kg/m^2] snow mass
	Scover = 0.0	# [m] snow volume
	cover  = 0.0	# [m] snow volume

	tset = 0

	data     = open('sapporo2017.csv', 'r')
	all_data = data.readlines()
	data_num = len(all_data) - 1
	data.close()

	day_cnt  = 0
	day_1    = 0
	oldmin   = -1

	# initialize Q table
	Qtable, comp = QL().Qlearn.initializeQ()


	try:
		""" Q-learn loop """
		for qnum in range(Qloop):
			percent = float(qnum)/Qloop * 100
			data_cnt = 0

#			if(os.path.exists('q_logs_'+str(MODE)+'.txt')):
#				os.remove('q_logs_'+str(MODE)+'.txt')

			data1 = all_data[1].split(', ')
			month  = int(data1[1])
			day    = data1[2]
			Hour   = data1[3]
			minute = int(data1[4])
			date = '2017-'+str(month)+'-'+day+' '+Hour+':'+str(minute)
			date = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M')

			### interval loop ###
			while( data_cnt < data_num-1 ):
				if(day != day_1):
					day_cnt += 1
				day_1 = day
#				sim().logf.write('\n' + str(date))

				shift = False

				# next 10 minutes
				if( int(date.minute)//10 != oldmin ):
					shift = True
					data_cnt += 1
					t0 = temp_o
					percent += 100.0/Qloop /data_num

					data1 = all_data[data_cnt].split(', ')
					month   = int(data1[1])
					day     = data1[2]
					Hour    = data1[3]
					minute  = int(data1[4])
					temp_o  = float(data1[5])		# temperature
					vaporP  = float(data1[6])		# [hPa] vapor pressure
					Wspeed0 = float(data1[7])		# [m/s] wind speed
					sun     = float(data1[8])		# [MJ/m^2] solar radiatioin
					pre     = float(data1[9])/10	# [mm/min] precipitation
					cloud   = float(data1[10])		# cloud cover
					nightR  = 45					# [W/m^2] nighttime radiation

					sun = sun/4.186*1000		# [kcal/m^2]

					# absolute humidity, night R
					absH, nightR = sim().absoluteHumid(vaporP,cloud,temp_o,nightR)

					oldmin = minute//10


				elif( (int(date.minute)+interval)//10 != oldmin ):
					shift = True


				Wspeed = Wspeed0 * windC		# [m/sec] after ccorrection

				# amount of snow and rain in precipitation (mass)
				snow_plus, rain_plus = sim().calc_plus(temp_o, pre)
				snow_plus = snow_plus * interval		# [kg/m^2]
				rain_plus = rain_plus * interval		# [kg/m^2]

				# density of snowfall [kg/m^3]
				sfdens = sim().snowfall_density(temp_o)


				water = 0		# [kg/m^2] water amount
				SE    = 0		# [kg/m^2] heat source calorific value
				erot  = 0
				ilp = 0
				erx = 0

				# no precipitation time [min]
				noPreT = noPreT + interval
				if(pre > 0):	noPreT = 0


				# Novenber -> April
				if( month<=4 or month>=11 ):
					# when snowing
					if(pre > 0):
						level = 0
					# when not snowing
					else:
						# wet
						if( (snow+Water) > 0 ):
							level = 1
						# dry
						else:
							level = 2


				""" Q learning """
				# levels
				Slevel = QL().Qlearn.Slevel(snow)
				if( MODE.find('S') >= 0 ):
					onSLv  = QL().next_Slevel(1, cover, temp_o, nightR, \
												Wspeed, Water, rain_plus, snow,\
												snow_plus, sfdens, tset, ilp, Qr)
					offSLv = QL().next_Slevel(0, cover, temp_o, nightR, \
												Wspeed, Water, rain_plus, snow,\
												snow_plus, sfdens, tset, ilp, Qr)
				if( MODE.find('T') >= 0 ):
					Tlevel = QL().Qlearn.Tlevel(temp_o)
					nextTLv = QL().next_Tlevel(all_data, shift, data_cnt)
				if( MODE.find('P') >= 0 ):
					Plevel = QL().Qlearn.Plevel(pre)
					nextPLv = QL().next_Plevel(all_data, shift, data_cnt)

				# decide action
				if( MODE == 'S' ):
					heater = QL().Qlearn.select_act_S(Qtable,Slevel,onSLv,offSLv)
				if( MODE == 'P' ):
					heater = QL().Qlearn.select_act_P(Qtable,Slevel,Plevel,nextPLv)
				elif( MODE=='ST' or MODE=='TS' ):
					heater = QL().Qlearn.select_act_ST(Qtable, Slevel, Tlevel,\
														onSLv, offSLv, nextTLv)
				elif( MODE=='SP' or MODE=='PS' ):
					heater = QL().Qlearn.select_act_SP(Qtable, Slevel, Plevel,\
														onSLv, offSLv, nextPLv)

				# update Q table
				if( MODE == 'S' ):
					Qtable, comp = QL().Qlearn.updateQ_S(Qtable, comp, heater,\
														Slevel, onSLv, offSLv)
				elif( MODE=='ST' or MODE=='TS' ):
					Qtable, comp = QL().Qlearn.updateQ_ST(Qtable, comp, heater,\
														Slevel, Tlevel, onSLv,\
														offSLv, nextTLv)
				elif( MODE == 'P' ):
					Qtable, comp = QL().Qlearn.updateQ_P(Qtable, comp, heater,\
														Slevel, Plevel, nextPLv)
				elif( MODE=='SP' or mode=='PS' ):
					Qtable, comp = QL().Qlearn.updateQ_SP(Qtable, comp, heater,\
														Slevel, Plevel, onSLv,\
														offSLv, nextPLv)


				ntime[heater][level] += 1

	
				E   = 0
				NC  = 0
				mlt = 0
				ict = 0

				# when heater on
				if(heater==1):
					# Q from heat source [kcal/h]
					Q = Qs
					# sum of Q
					erot = erot + Q*interval/60.0		# [kcal]
				# when heater off
				else:
					Q = 0

	
				# solar radiatioin absorption rate of snow
				abrate = 0.8 - 30*cover
				if(abrate < abrate0):	abrate = abrate0


				# 相当外気温
				sat = temp_o + (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4)


				TS  = BT			# temperature(?)
				# 路面に雪があるとき
				if( (snow+snow_plus)>0 ):
					if( TS>0 or sat>0 ):
						mlt = 1
				# 路面に水があるとき
				if( (Water+rain_plus)>0 ):
					if( TS<0 ):
						mlt = 1


				# penetration height [m]
				peneH = sim().penetration_height(cover, snow, Water)


				# total water
				wat = Water + rain_plus		# [kg/m^2]


				melt  = 0
				evaporate = 0
				BF = 0			# (?)
				if(mlt==1):
					DH = cover - peneH
					if(DH<=0 or sat>0):
						DH  = 0
						tsv = 0			# (?)
						# evaporation amount
						evaporate = 4 * sim().funa(Wspeed) \
									* (sim().abshumid(tsv)-absH)	# [kg/m^2/h]
						# total water [kg/m^2]
						wat = Water + rain_plus - evaporate*(interval/60.0)
						# abnormal value correction
						if(wat < 0):
							evaporate  = Water/(interval/60.0) + rain_plus
							wat = 0
					htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)

					# calc amount of snow melting
					melt = (200*TS+htrm*sat-590*evaporate) \
								* (interval/60.0)/(Hfusion*180)
					BF   = 1			# (?)

					# 算出された融雪量が水分量より多い時
					if(melt < -1*wat):
						melt = -1 * wat
					# 算出された融雪量が存在していた雪より多い時
					elif( melt > (snow+snow_plus) ):
						melt = snow + snow_plus		# [kg/m^2]
					if(melt < 0):
						melt = 0

				else:
					tsv = BT			# (?)
					evaporate = 4 * sim().funa(Wspeed)\
								* (sim().abshumid(tsv)- absH)
					if( evaporate > (Water/(interval/60.0) + rain_plus) ):
						evaporate = Water/(interval/60.0) + rain_plus

				htr  = 1 / (1/(sim().funa(Wspeed)+4) + cover/0.08)
				Qe   = -590*evaporate*(interval/60.0) * area * (1-BF)
				Snow = snow - melt + snow_plus						# [kg/m^2]
				water_plus = rain_plus - evaporate*(interval/60.0)	# [kg/m^2]
				ww   = Water + melt + water_plus					# [kg/m^2]

				# snow accumulation (volume) [m]
				if(Snow > 0):
					if(melt > 0):
						Scover = cover + snow_plus/sfdens - melt/916
					else:
						Scover = cover + snow_plus/sfdens
				else:
					Scover = 0


				T = BT


				while(True):
					if(ict==1):
						T = tset
						E = 0
					S1 = (Q+Qe+E)*interval/60.0 + BT*C
					S2 = C

					for j in range(npn):
						tmp = T
						S1 = S1 + HR[j]*(tmp)*interval/60.0
						S2 = S2 + HR[j]*interval/60.0

					tmp = (1-BF)*htr*sat
					S1 = S1 + ( tmp )*interval/60.0*area
					S2 = S2 + (BF*200+(1-BF)*htr)*interval/60.0*area
					if(ict==1):
						E = (S2*tset - S1) / interval/60.0
					else:
						TT = S1 / S2
						ER = TT - T
						aer = abs(ER)
						if(aer > erx):
							erx = aer
						T = CK*ER + T

					if(erx > 0.0001):
						erx = 0

					if(ilp < 2):
						if( T>maxT and Q>0 ):
							ilp += 1
							tset = maxT
							ict = 1
							Q = 0
							erx = 0
							continue
						if(E < 0):
							ilp += 1
							ict = 0
							E = 0
							erx = 0
							continue
						else:	break
					else:	break


				if( level==0 and T<0 ):
					snow_minusT += 1
				elif( level==1 and T<0 ):
					wet_minusT += 1

				# abnormal value correction
				if(Snow < 0):
					Snow = 0
				if(ww < 0):
					ww = 0		# [kg/m^2]


				water = water + ww		# [kg/m^2]

				if( Scover>0 and Snow>0 ):
					gm = Snow / Scover
					ee = 16 * math.exp(0.021*gm)
					gm = gm * math.exp( Snow/2/ee*interval/60.0/24 )
					if(gm > 916):	gm = 916
					Scover = Snow / gm		# [m]


				SE = E + Q + erot

				if(SE > 0):
					onT += 1
					Qsup += SE

				if( pre>0.0 and T>0 ):
					TS = T
					Qr_plus = BF*200*TS + (1-BF)*htr*(TS-sat)
					Qr = Qr + (Qr_plus)*area

				if( (Snow+pre)==0.0 and ww>0.0 ):
					ww = 0.0		# [kg/m^2]
				Water = ww			# [kg/m^2]
				snow  = Snow		# [kg/m^2]
#				sim().logf.write('{:>5}℃\t' .format(temp_o))
#				sim().logf.write('+{:.2f}[kg/m^2]  ' .format(snow_plus))
#				sim().logf.write('-> {:.3f}[kg/m^2]' .format(snow))
				cover = Scover		# [m]

				BT = T

				Tsnow_off = ntime[0][0] * interval
				Tsnow_on  = ntime[1][0] * interval
				Twet_off  = ntime[0][1] * interval
				Twet_on   = ntime[1][1] * interval
				Tdry_off  = ntime[0][2] * interval
				Tdry_on   = ntime[1][2] * interval

				str_snow = '{:.3f}[kg/m^2]' .format(snow)
				sys.stderr.write('\r{} {:>6}℃   {:.3f}[mm/min] {:>15}  {}  {:.2f}%'
						.format(date, temp_o, pre, str_snow, heater, percent))

				date = date + datetime.timedelta(minutes=interval)

#				time.sleep(0.5)

			sim().logf.write('\n'+str(Qtable)+'\n')

	finally:
		onT         = onT * interval			# [min]
		snow_minusT = snow_minusT * interval	# [min]
		wet_minusT  = wet_minusT * interval		# [min]
		Qsup        = Qsup * interval/60.0
		Qr          = Qr * interval/60.0

		Tsnow_off = ntime[0][0] * interval
		Tsnow_on  = ntime[1][0] * interval
		Twet_off  = ntime[0][1] * interval
		Twet_on   = ntime[1][1] * interval
		Tdry_off  = ntime[0][2] * interval
		Tdry_on   = ntime[1][2] * interval

		sim().logf.write('\n'+str(Qtable))

		time = time.time() - start
		print('\ntime{0:4d}:{0:02d}' .format(int(time)//60, int(time)%60))
