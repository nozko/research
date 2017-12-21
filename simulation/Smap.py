#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
import re
import math
import time
import argparse


class sim:

	def __init__(self):
		fn1 = 'file.net'
		self.net = open(fn1, 'r')
		fn3 = 'smap.run'
		self.run = open(fn3, 'r')


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



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='weather data file')
	parser.add_argument('weather')
	args = parser.parse_args()

	snow_minusT = 0
	wet_minusT  = 0


	# initialize
	temp_o = 0		# temperature outside
	Qsup   = 0		# supplied Q
	onT    = 0		# operating time
	height = 0		# height of snow
	slevel = 0		# level of snow accumulation
	heater = 0		# on(1) / off(0)
	irt    = 0
	noPreT = 0		# no precipitation time
	lpmx   = 0
	ntime  = np.zeros((2, 3))
	Qr  = 0		# Q from surface(snowing, plus temperature)
	Qr2 = 0		# Q from surface(snowing, plus temperature) @moist sensor node


	interval = 0.1		# calculatioin interval
	diffM    = 1			# difference method(1:後退差分)
	CK       = 1.0		# 緩和係数(0.7~1.5程度)
	maxloop  = 100		# maximum number of loops(200)
	print('interval :', interval)

	idx = int( 1/(interval+0.4) )

	smap4 = sim().run.readlines()[4].split('\t')
	Qs      = float(smap4[0])	# Q of heat source
	maxT    = 70				# maximum temperature of heat source
	timer   = float(smap4[5])	# delay time
	circuit = int(smap4[7])		# number of rotatioin circuit

	smap5 = sim().run.readlines()[5].split('\t')
	snowT  = float(smap5[0])	# operating temperature during snowfall
	B1     = float(smap5[1])	# 外気温にかかる係数1
	stopT1 = float(smap5[2])	# stop temperature(difference from operating one)
	wetT   = float(smap5[3])	# operating temperature with wet surface
	B2     = float(smap5[4])	# 外気温にかかる係数2
	stopT2 = float(smap5[5])	# stop temperature
	dryT   = float(smap5[6])	# operating temperature with dry surface
	B3     = float(smap5[7])	# 外気温にかかる係数3
	stopT3 = float(smap5[8])	# stop temperature
	offT   = float(smap5[9])	# not driving temperature

	remainW = 0.0		# remained water after melting
	maxPene = 0.05		# maximum penetration height of water
	abrate0 = 0.2		# solar radiatioin absorption rate of snow(0.2)
	windC   = 0.5		# wind speed correction coefficient(~1)
	Scover0 = 0.0		# initial snow cover
	dens0   = 100		# initial snow density
	snow0 = Scover0 * dens0

	sim().run.close()


	Qe    = 0
	BF    = 0
	Water = 0
	Snow  = 0
	ww    = 0
	Evapo = 0		# eveporate
	htr   = 0
	sat   = 0		# 相当外気温 Sol-Air Temperature

	BT = 8.0
	C  = 0.052629437

	area = 0.02783305
	QQ  = 0.0

	snow   = snow0
	Scover = Scover0
	cover  = Scover0

	net = re.split(" +", sim().net.readlines()[4])
	npn = int(net[1])
	HR = []
	for j in range(int(net[1])):
		netj = re.split(" +", sim().net.readlines()[5+j])
		HR.append(float(netj[3]))

	sim().net.close()


	data     = open(args.weather, 'r')
	all_data = data.readlines()
	day_num  = int( math.floor( len(all_data)/24 ) )


	### day loop ###
	print('\nperiod :', day_num, 'days')
	for d in range( day_num ):
		print('\n----- day', d+1, '------')
		time.sleep(1)


		### time loop ###
		for hour in range(24):
			t0 = temp_o
			data1 = all_data[ d*24+hour+1 ].split(', ')
			month   = int(data1[1])
			day     = int(data1[2])
			Hour    = data1[3]
			temp_o  = float(data1[4])	# temperature
			vaporP  = float(data1[5])	# vapor pressure
			Wspeed0 = float(data1[6])	# wind speed
			sun     = float(data1[7])	# solar radiatioin
			pre     = float(data1[8])	# precipitation
			cloud   = float(data1[9])	# (maybe) cloud cover
			nightR  = 45				# nighttime radiation
			print('\n'+str(month)+'/'+str(day)+'  '+Hour+':00')
			print('temperature :',temp_o,'℃ ', '\tvapor pressure      :',vaporP)
			print('wind speed  :',Wspeed0,     '\tsolar radiation     :',sun)
			print('cloud cover :',cloud,       '\tnighttime radiation :',nightR)
			data.close()

			sun = sun/4.186*1000


			# absolute humidity when there is vapor in the air
			if(vaporP >= 0):
				absH = 0.622*vaporP/(1013.25-vaporP)
				
				# nighttime radiation when there is a cloud
				if(cloud >= 0):
					RH = vaporP*760/1013.25
					nightR = 0.0000000488 * (temp_o+273.16)**4 \
							* (1-0.62*cloud/10) * (0.49-0.076*math.sqrt(RH))

			# absolute humidity when no vapor in the air
			else:
				absH = 0.7 * sim().abshumid(tempA)


			for idt in range(idx):
				#temp_o = t0 + (tempA-t0)*interval*idt
				Wspeed = Wspeed0 * windC

				# amount of snow and rain in precipitation
				if(temp_o < 0):
					snow_plus = pre
					rain_plus = 0
				elif(temp_o>=0 and temp_o<2):
					snow_plus = pre * (2-temp_o)/2
					rain_plus = pre * temp_o / 2
				elif(temp_o>=2):
					snow_plus = 0
					rain_plus = pre

				# density of snowfall (max:50)
				sfdens = 1000 / (0.091*temp_o**2 - 1.81*temp_o + 9.47)
				if(sfdens < 50):	sfdens = 50


				sca   = 0
				FD    = 0
				water = 0		# water amount
				SE    = 0		# heat source calorific value
				erot  = 0
				evaporate = 0	# evaporation amount
				ilp = 0
				erx = 0
				L1 = 0
				L2 = 0

				noPreT = noPreT + interval
				if(pre > 0):	noPreT = 0


				# Novenber -> April
				if( month<=4 or month>=11 ):
					# when snowing
					if(pre > 0):
						level = 0
						TN = snowT + B1*temp_o
						TF = TN + stopT1
					# when not snowing
					else:
						# wet
						if( (snow+Water) > 0 ):
							level = 1
							TN = wetT + B2*temp_o
							TF = TN + stopT2
						# dry
						else:
							level = 2
							TN = dryT + B3*temp_o
							TF = TN + stopT3


					# decide on/off
					if(BT < TN):				heater = 1
					elif(BT > TF):				heater = 0
					if(noPreT > timer+0.0001):	heater = 0
					if(temp_o > offT):			heater = 0

					ntime[heater][level] += 1


				E = 0
				Q = 0
				NC = 0
				mlt = 0
				ict = 0

				# when heater on
				if(heater==1):
					Q = Qs
					erot = erot + Q


				abrate = 0.8 - 30*cover
				if(abrate < abrate0):	abrate = abrate0
				sat = temp_o + (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4)
				Evapo  = 0
				TS  = BT
				lps = 0
				dps = cover
				if( (snow+snow_plus)>0 ):
					if(TS<0 or sat>0):
						mlt = 1
				if( (Water+rain_plus)>0 ):
					if( TS<0 or (Snow+rain_plus)>0 ):
						mlt = 1


				# penetration height
				peneH = 0
				if(cover > 0):
					peneH = Water / (1000 - snow/cover)
					# snow/cover : density of snow cover
				else:
					peneH = Water / 1000
				if(peneH > maxPene):
					peneH = maxPene

				wat = Water + rain_plus*interval


				while(True):
					melt  = 0
					Evapo = 0
					BF = 0
					if(mlt==1):
						DH = dps - peneH
						if(DH<=0 or sat>0):
							DH  = 0
							tsv = 0
							Evapo = 4 * sim().funa(Wspeed) \
									* (sim().abshumid(tsv)-absH)
							wat = Water + (rain_plus-Evapo)*interval
							if(wat < 0):
								Evapo  = Water/interval + rain_plus
								wat = 0
						htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)

						# calc amount of snow melting
						melt  = (200*TS + htrm*sat - 590*Evapo) * interval/80
						BF   = 1

						# 算出された融雪量が水分量より多い時
						if(melt < -1*wat):
							melt = -1 * wat
						# 算出された融雪量が存在する雪より多い時
						elif( melt > (snow+snow_plus*interval) ):
							melt = snow + snow_plus*interval
						print('melt :', melt)

					elif(mlt==0):
						tsv = BT
						Evapo = 4 * sim().funa(Wspeed)\
								* (sim().abshumid(tsv)- absH)
						if( Evapo > (Water/interval + rain_plus) ):
							Evapo = Water/interval + rain_plus

					htr  = 1 / (1/(sim().funa(Wspeed)+4) + dps/0.08)
					Qe   = -590*Evapo * area * (1-BF)
					Snow = snow - melt + snow_plus*interval
					water_plus = rain_plus-Evapo
					ww   = Water + melt + water_plus*interval

					if(Snow > 0):
						if(melt > 0):
							G = (snow + snow_plus*interval) \
								/ (cover + snow_plus*interval/sfdens)
							Scover = Snow / G
						else:
							Scover = cover + snow_plus*interval/sfdens - melt/916
					else:
						Scover = 0

					ddps = (Scover+cover)/2 - dps
					if( abs(ddps)>0.001 and lps<10 ):
						dps = 0.7*ddps + dps
						lps += 1
					else:
						break


				T = BT

				while(True):
					lps = 0

					while(True):
						if(ict==1):
							T = tset
							E = 0
						S1 = (Q+Qe+E)*interval + BT*C
						S2 = C

						for j in range(npn):
							tmp = diffM * T
							S1 = S1 + HR[j]*(tmp)*interval
							S2 = S2 + diffM*HR[j]*interval

						tmp = (1-BF)*htr*sat - BF*200*(1-diffM)*BT \
								- (1-BF)*htr*(1-diffM)*BT
						S1 = S1 + ( tmp )*interval*area
						S2 = S2 + (BF*200+(1-BF)*htr)*interval*area*diffM
						if(ict==1):
							E = (S2*tset - S1) / interval
						else:
							if(diffM<=0.01):
								T = S1 / S2
								erx = 0
							else:
								TT = S1 / S2
								ER = TT - T
								aer = abs(ER)
								if(aer > erx):
									erx = aer
									ier = 0
								T = CK*ER + T

						lps += 1
						if(lps > maxloop):
							break
						else:
							if(erx > 0.0001):
								erx = 0
							else:
								if(lps > lpmx):	lpmx = lps
								break

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
					L1 = 1
				elif( level==1 and T<0 ):
					L2 = 1

				# abnormal value correction
				if(Snow < 0):
					Snow = 0
				if(ww < 0):
					ww = 0

				if(Scover > FD):
					FD = Scover

				sca = sca + Snow
				water = water + ww
				evaporate = evaporate + Evapo*interval

				if( Scover>0 and Snow>0 ):
					gm = Snow / Scover
					ee = 16 * math.exp(0.021*gm)
					gm = gm * math.exp( Snow/2/ee*interval/24 )
					if(gm > 916):
						gm = 916
					Scover = Snow / gm


				if(FD > 0):
					ifd = int(FD*100)
					if(ifd > height):
						height = ifd
						print('snow height :', height)

				if(sca > 0):
					isca = int(sca)
					nsca = np.zeros(10)
					for k in range(isca+1):
						nsca[k] += 1
					if(isca > slevel):
						slevel = isca

				snow_minusT += L1
				wet_minusT  += L2
				SE = E + Q + erot

				if(SE > 0):
					onT += 1
					Qsup += SE
					rse = Qsup / onT

				if( pre>0 and T>0 ):
					TS = diffM*T + (1-diffM)*BT
					Qr2 = Qr2 + BF*200*TS + (1-BF)*htr*(TS-sat)

				if( (sca+pre)==0 and ww>remainW ):
					ww = remainW
				Water = ww
				snow = Snow
				print('snow  :', snow)
				cover = Scover
				print('cover :', cover)
				if( pre>0 and T>0 ):
					TS = diffM*T + (1-diffM)*BT
					Qr_plus = BF*200*TS + (1-BF)*htr*(TS-sat)
					Qr = Qr + (Qr_plus)*area

				BT = T

			Tsnow_off = ntime[0][0] * interval
			Tsnow_on  = ntime[1][0] * interval
			Twet_off  = ntime[0][1] * interval
			Twet_on   = ntime[1][1] * interval
			Tdry_off  = ntime[0][2] * interval
			Tdry_on   = ntime[1][2] * interval
			print('snow.on time :', Tsnow_on, '\tsnow.off time :', Tsnow_off, \
				'\nwet.on time  :', Twet_on,  '\twet.off time  :', Twet_off, \
				'\ndry.on time  :', Tdry_on,  '\tdry.off time  :', Tdry_off)
			if(heater==0):	print('heater : off')
			else:			print('heater : on')

			time.sleep(2)


		sim().run.close()

		onT         = onT * interval
		snow_minusT = snow_minusT * interval
		wet_minusT  = wet_minusT * interval
		Qsup        = Qsup * interval
		Qr          = Qr * interval
		Qr2         = Qr2 * interval

		print(height)
		for i in range(height+1):
			Theight = nsca[i] * interval
		for i in range(slevel+1):
			Tcover = nsca[i] * interval
	
		Tsnow_off = ntime[0][0] * interval
		Tsnow_on  = ntime[1][0] * interval
		Twet_off  = ntime[0][1] * interval
		Twet_on   = ntime[1][1] * interval
		Tdry_off  = ntime[0][2] * interval
		Tdry_on   = ntime[1][2] * interval
		print('snow.on :', Tsnow_on, '\tsnow.off :', Tsnow_off, \
			'\nwet.on  :', Twet_on,  '\twet.off  :', Twet_off, \
			'\ndry.on  :', Tdry_on,  '\tdry.off  :', Tdry_off)
