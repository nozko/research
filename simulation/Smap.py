#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
import re
import math
import time

RD   = 3.1416 / 180
pido = 43.04
pked = 141.33
pcth = 0
pctl = 0.2
tvm  = (pcth + pctl) / 2
cth  = 1.67


class sim:

	def __init__(self):
		fn1 = 'file.net'
		self.net = open(fn1, 'r')
		fn3 = 'smap.run'
		self.run = open(fn3, 'r')
		fn5 = 'asahikawa.csv'
		self.data = open(fn5, 'r')


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


	def sunpo(self, TM, SH, CH, SA, CA):
		w = 2 * 3.1416 * Lday / 366
		Del = 0.362213 - 23.2476*math.cos(w+0.153231) \
			- 0.336891*math.cos(2*w+0.207099) - 0.185265*math.cos(3*w+0.620129)
		AE = -0.000278641 + 0.122772*math.cos(w+1.49831) \
			- 0.165458*math.cos(2*w-1.26155) - 0.00535383*math.cos(3*w-1.1571)
		sjik = 15 * (TM-12+AE) + pked - 135
		SH = math.sin(pido*RD)*math.sin(Del*RD) \
			+ math.cos(pido*RD)*math.cos(Del*RD)*math.cos(sjik*RD)
		if(SH<=0):
			SH = 0
			CH = 0
			SA = 0
			CA = 0
		else:
			CH = math.sqrt(1-SH**2)
			SA = math.cos(Del*RD)*math.sin(sjik*RD)/CH
			CA = (SH*math.sin(pido*RD)-math.sin(Del*RD)) / CH / math.cos(pido*RD)
		return SH, CH, SA, CA


	def qds(SH, CH, SA, CA, rdn, rsh, rtv):
		if(SH > 0):
			rdv = rdn * SH
		else:
			rdv = 0
		rfv = rsh
		rtv = rdv + rfv
		return rdv, rfv, rtv



if __name__ == '__main__':

	nsca = np.zeros(100)
	HR = []

	snow_minusT = 0
	wet_minusT  = 0

	area = 0	#?


	# initialize
	temp_o = 0		# temperature
	Qsup   = 0		# supply calorie
	onT     = 0		# operating time
	height = 0		# height of snow
	slevel = 0		# level of snow accumulation
	heater = 0		# on(1) / off(0)
	irt    = 0
	melt   = 0		# melted snow
	noPreT = 0		# no precipitation time
	lpmx   = 0
	lpr    = 1
	ntime  = np.zeros((2, 3))
	Qr  = 0		# Q from surface(snowing, plus temperature)
	Qr2 = 0		# Q from surface(snowing, plus temperature) @moist sensor node


	interval = 0.1		# calculatioin interval
	diffM    = 1			# difference method(1:後退差分)
	CK       = 1.0		# 緩和係数(0.7~1.5程度)
	maxloop  = 100		# maximum number of loops(200)
	print('interval :', interval)

	idx = int( 1/(interval+0.4) )

	Dstart = 1		# day to start calculation(days from top of the file)
	lyear  = 1		# number of years to calculate(1~)
	Dend   = 181	# last day of calculation(days from top of the file)
	print('years :', lyear)
	Lyear  = 0
	Lday   = Dstart - 1

	smap4 = sim().run.readlines()[4].split('\t')
	Qs      = float(smap4[0])	# Q of heat source
	flowR   = 40				# flow rate of antifreeze
	maxT    = 70				# maximum temperature of heat source
	timer   = float(smap4[5])	# delay time
	step    = int(smap4[6])		# rotation step
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


	line_num = 0

	# array initialize
	QE    = 0
	BF    = 0
	W     = 0		# water
	SS    = 0
	WW    = 0
	EV    = 0
	htr   = 0
	sat   = 0		# 相当外気温 Sol-Air Temperature
	smf   = 0
	cover = 0		# snow cover

	line_num += 1
	net1 = re.split(" +", sim().net.readlines()[line_num])
	BT = float(net1[2])
	C  = float(net1[4])
	CL = float(net1[5])

	line_num += 1
	net2 = re.split(" +", sim().net.readlines()[line_num])
	nsd = float(net2[1])
	AS  = float(net2[2])
	QQ  = float(net2[4])

	Scover = 0
	snow   = 0

	line_num += 1

	if(int(net2[1])==1):
		area   = area + float(net2[2])
		snow   = snow0
		Scover = Scover0
		cover  = Scover0

	line_num += 1
	net4 = re.split(" +", sim().net.readlines()[line_num])
	npn = int(net4[1])
	for j in range(int(net4[1])):
		line_num += 1
		net4j = re.split(" +", sim().net.readlines()[line_num])
		HR.append(float(net4j[3]))

	flr = 0

	sim().net.close()


	### day loop ###
	for day in range(Dend-Dstart):
		Lday += 1
		if(Lday > 365):	Lday = 1
		if(Lday==Dstart):	Lyear += 1
		print('\n----- year', Lyear, ', day', Lday, '------\n')
		print('start day :', Dstart, ', end day :', Dend)


		### time loop ###

		for hour in range(24):
			t0 = temp_o
			data1 = sim().data.readlines()[ (Lday-1)*24+hour ].split(',')
			month   = int(data1[1])
			day     = int(data1[2])
			Hour    = data1[3]
			tempA   = float(data1[4])	# temperature
			vaporP  = float(data1[5])	# vapor pressure
			Wspeed0 = float(data1[6])	# wind speed
			sun     = float(data1[7])	# solar radiatioin
			pre     = float(data1[8])	# precipitation降水量
			cloud   = float(data1[9])	# (maybe) cloud cover
			nightR  = 45				# nighttime radiation
			print('\n'+str(month)+'/'+str(day)+'  '+Hour+':00')
			print('temperature :', tempA,   '\tvapor pressure      :', vaporP)
			print('wind speed  :', Wspeed0, '\tsolar radiation     :', sun)
			print('cloud cover :', cloud,   '\tnighttime radiation :', nightR)

			sun = sun/4.186*1000


			# there is vapor in the air
			if(vaporP >= 0):
				absH = 0.622*vaporP/(1013.25-vaporP)	# absolute humidity
				# nighttime radiation when there is a cloud
				if(cloud >= 0):
					RH = vaporP*760/1013.25
					nightR = 0.0000000488 * (tempA+273.16)**4 \
							* (1-0.62*cloud/10) * (0.49-0.076*math.sqrt(RH))
			# no vapor in the air
			else:
				absH = 0.7 * sim().abshumid(tempA)


			for idt in range(idx):
				TM = hour - 1 + interval*idt

				temp_o = t0 + (tempA-t0)*interval*idt
				print('temp_o :', temp_o, '℃ ')
				Wspeed = Wspeed0 * windC

				# amount of snow and rain in precipitation
				if(temp_o < 0):
					Fs = pre
					Fr = 0
				elif(temp_o>=0 and temp_o<2):
					Fs = pre * (2-temp_o)/2
					Fr = pre * temp_o / 2
				elif(temp_o>=2):
					Fs = 0
					Fr = pre

				# density of snowfall (max:50)
				sfdens = 1000 / (0.091*temp_o**2 - 1.81*temp_o + 9.47)
				if(sfdens < 50):	sfdens = 50


				# set
				sca = 0
				ssmf = 0
				FD = 0
				water = 0		# water amount
				SE = 0			# heat source calorific value
				erot = 0
				evaporate = 0	# evaporation amount
				ilp = 0
				erx = 0
				L1 = 0
				L2 = 0

				noPreT = noPreT + interval
				if(pre > 0):	noPreT = 0


				# Novenber -> April
				if(month<=4 or month>=11):
					# snowing
					if(pre > 0):
						level = 0
						TN = snowT + B1*temp_o
						TF = TN + stopT1
					# not snowing
					else:
						# wet
						if( (snow+W) > 0 ):
							level = 1
							TN = wetT + B2*temp_o
							TF = TN + stopT2
						# dry
						else:
							level = 2
							TN = dryT + B3*temp_o
							TF = TN + stopT3


					# decide on/off
					if(BT < TN):		heater = 1
					elif(BT > TF):		heater = 0
					if(noPreT > timer+0.0001):	heater = 0
					if(temp_o > offT):			heater = 0

					# when heater off
					if(heater==0):
						irt = -1
					# when heater on
					elif(heater==1):
						irt += 1
						IA = irt % (circuit*step)	# circuit*step : rotation time
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

				CR = C
				if(CL > 0):
					if(BT > tvm):
						CR = C + CL
					else:
						CR = C + CL*0.5
					Q = flr
					if(IL==1):
						tic = BT
						pic = (pcth-tic) / (pcth-pctl)
						CR = C + CL*(1-pic) \
								+ CL*0.5*pic + CL*80/abs(pcth-pctl)

				if(nsd==1):
					abrate = 0.8 - 30*cover
					if(abrate < abrate0):	abrate = abrate0
					sat = ( temp_o \
							+ (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4) )
					EV = 0
					TS = BT
					lps = 0
					dps = cover
					if( (snow+Fs)>0 ):
						if(TS<0 or sat>0):
							mlt = 1
					if( (W+Fr)>0 ):
						if(TS < 0):
							mlt = 1
					if( (snow+Fs)>0 and (W+Fr)>0 ):
						mlt = 1
					peneH = 0			# penetration height
					if(cover > 0):
						peneH = W / (1000 - snow/cover)
						# snow/cover : density of snow cover
					else:
						peneH = W / 1000
					if(peneH > maxPene):
						peneH = maxPene
					wat = W + Fr*interval


					while(True):
						SM = 0		# amount of melted snow(?)
						EV = 0
						BF = 0
						if(mlt==1):
							DH = dps - peneH
							if(DH<=0 or sat>0):
								DH = 0
								tsv = 0
								EV = 4 * sim().funa(Wspeed) \
										* (sim().abshumid(tsv)-absH)
								wat = W + (Fr-EV)*interval
								if(wat < 0):
									EV = W/interval + Fr
									wat = 0
							htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)
							SM0 = (200*TS + htrm*sat - 590*EV) \
									* interval/80
							SM = SM0
							BF = 1
							if(SM0 < -1*wat):
								SM = -1 * wat
								BF = SM / SM0
							elif( SM0 > (snow+Fs*interval) ):
								SM = snow + Fs*interval
								BF = SM / SM0
						elif(mlt==0):
							tsv = BT
							EV = 4 * sim().funa(Wspeed)\
									* (sim().abshumid(tsv)- absH)
							if( EV > (W/interval + Fr) ):
								EV = W/interval + Fr

						htr = 1 / (1/(sim().funa(Wspeed)+4) + dps/0.08)
						QE  = -590*EV * AS * (1-BF)
						SS  = snow - SM + Fs*interval
						WW  = W + SM + (Fr-EV)*interval
						smf = SM

						if(SS > 0):
							if(SM > 0):
								G = (snow + Fs*interval) \
									/ (cover + Fs*interval/sfdens)
								Scover = SS / G
							else:
								Scover = cover + Fs*interval/sfdens \
												- SM/916
						else:
							Scover = 0

						ddps = (Scover+cover)/2 - dps
						if(abs(ddps)>0.001 and lps<10):
							dps = 0.7*ddps + dps
							lps += 1
						else:
							break


				# matrix
				T = BT

				while(True):
					lps = 0

					while(True):
						if(ict==1):
							T = tset
							E = 0
						S1 = (Q+QE+E)*interval + BT*CR
						S2 = CR

						for j in range(npn):
							tmp = diffM * T
							S1 = S1 + HR[j]*(tmp)*interval
							S2 = S2 + diffM*HR[j]*interval

						if(nsd==1):
							F = BF
							tmp = (1-F)*htr*sat \
								- F*200*(1-diffM)*BT \
								- (1-F)*htr*(1-diffM)*BT
							S1 = S1 + ( tmp )*interval*AS
							S2 = S2 \
									+ (F*200+(1-F)*htr)\
										*interval*AS*diffM
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

				if(CL > 0):
					IL = 0
					flr = 0
					if( T<pcth and T>pctl ):
						IL = 1
					elif( BT>pcth and T<pcth ):
						flr = (C+CL) * (T-pcth)
						T = pcth
						IL = 1
					elif( BT>pctl and T<pctl ):
						flr = (C + CL*80/abs(pcth-pctl)) \
									* (T-pctl)
						T = pctl
					elif( BT<pctl and T>pctl ):
						flr = (C + CL*0.5) * (T-pctl)
						T = pctl
						IL = 1
					elif( BT<pcth and T>pcth ):
						flr = (C + CL*80/abs(pcth-pctl)) \
									* (T-pcth)
						T = pcth

				if(nsd==1):
					if( level==0 and T<0 ):
						L1 = 1
					elif( level==1 and T<0 ):
						L2 = 1
					if(SS < 0):
						SS = 0
					if(WW < 0):
						WW = 0
					if(Scover > FD):
						FD = Scover

					sca = sca + SS*AS/area
					water = water + WW*AS/area
					evaporate = evaporate + EV*interval*AS/area
					ssmf = ssmf + smf*AS/area

					if( Scover>0 and SS>0 ):
						GM = SS / Scover
						EE = 16 * math.exp(0.021*GM)
						GM = GM * math.exp( SS/2/EE*interval/24 )
						if(GM > 916):
							GM = 916
						Scover = SS / GM
				

				if(FD > 0):
					#ifd = int(FD*100)
					ifd = int(FD)
					if(ifd > height):
						height = ifd
						print('height :', height)

				if(sca > 0):
					isca = int(sca)
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

				melt += ssmf

				lpr += 1

				if( pre>0 and T>0 ):
					TS = diffM*T + (1-diffM)*BT
					Qr2 = Qr2 + BF*200*TS \
							+ (1-BF)*htr*(TS-sat)

				if(nsd==1):
					if( (sca+pre)==0 and WW>remainW ):
						WW = remainW	# remained water
					W = WW
					snow = SS
					print('snow  :', snow)
					cover = Scover
					print('cover :', cover)
					if( pre>0 and T>0 ):
						TS = diffM*T + (1-diffM)*BT
						Qr_plus = BF*200*TS \
									+ (1-BF)*htr*(TS-sat)
						Qr = Qr + (Qr_plus)*AS

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

			time.sleep(3)

		if( Lday!=Dend or Lyear!=lyear ):	continue

		sim().net.close()
		sim().run.close()
		sim().data.close()

		onT = onT * interval
		snow_minusT = snow_minusT * interval
		wet_minusT  = wet_minusT * interval
		Qsup = Qsup * interval
		Qr = Qr * interval
		Qr2 = Qr2 * interval

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
