#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import print_function
import math
import numpy as np
import re
import math

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
		SH = math.sin(pido*RD)*math.sin(Del*RH) \
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

	BT = []
	Q = []
	E = []
	tset = []
	ID = []
	QQ = []
	npn = []
	nsd = []
	nqd = []
	C = []
	CR = []
	AS = []
	iopt = []
	iope = []
	iopq = []
	CL = []
	ict = []
	nsca = np.zeros(100)
	NP = []
	nkr = []
	HR = []
	S1 = []
	S2 = []
	TT = []
	ER = []
	flr = []
	IL = []
	snow = []			# amount of snow and ice
	Scover = []
	mlt = []
	NC = []
	Theight = []
	Tcover = []

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
	tmr    = 0
	lpmx   = 0
	lpr    = 1
	ntime  = np.zeros((2, 3))
	Qr  = 0		# Q from surface(snowing, plus temperature)
	Qr2 = 0		# Q from surface(snowing, plus temperature) @moist sensor node


	interval   = 0.1		# calculatioin interval
	diffM      = 1			# difference method(1:後退差分)
	CK         = 1.0		# 緩和係数(0.7~1.5程度)
	maxloop    = 100		# maximum number of loops(200)
	print('interval :', interval)

	idx    = int( 1/(interval+0.4) )

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
	Tsensor = int(smap4[3])		# number of node(temperature sensor)
	Wsensor = 77				# number of node(moist sensor)
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
	ipx = int(sim().net.readlines()[line_num])		# data num

	# array initialize
	QE    = np.zeros(ipx)
	BF    = np.zeros(ipx)
	W     = np.zeros(ipx)		# water
	SS    = np.zeros(ipx)
	WW    = np.zeros(ipx)
	EV    = np.zeros(ipx)
	htr   = np.zeros(ipx)
	sat   = np.zeros(ipx)		# 相当外気温 Sol-Air Temperature
	smf   = np.zeros(ipx)
	cover = np.zeros(ipx)		# snow cover

	### read all data in 'file.net'
	for ip in range(ipx):
		line_num += 1
		net1 = re.split(" +", sim().net.readlines()[line_num])
		BT.append(float(net1[2]))
		ID.append(float(net1[3]))
		C.append(float(net1[4]))
		CL.append(float(net1[5]))
		print(ID)

		line_num += 1
		net2 = re.split(" +", sim().net.readlines()[line_num])
		nsd.append(float(net2[1]))
		AS.append(float(net2[2]))
		nqd.append(float(net2[3]))
		QQ.append(float(net2[4]))

		Scover.append(0)
		snow.append(0)

		line_num += 1
		net3 = re.split(" +", sim().net.readlines()[line_num])
		iopt.append(int(net3[1]))
		iopq.append(float(net3[2]))
		iope.append(float(net3[3]))
		if(int(net2[1])==1):
			area = area + float(net2[2])
			snow.append(snow0)
			Scover.append(Scover0)
			cover[ip] = Scover0

		line_num += 1
		net4 = re.split(" +", sim().net.readlines()[line_num])
		npn.append(int(net4[1]))
		_NP = []
		_nkr = []
		_HR = []
		for j in range(int(net4[1])):
			line_num += 1
			net4j = re.split(" +", sim().net.readlines()[line_num])
			_NP.append(int(net4j[1]))
			_nkr.append(float(net4j[2]))
			_HR.append(float(net4j[3]))
		NP.append(_NP)
		nkr.append(_nkr)
		HR.append(_HR)

		flr.append(0)

	sim().net.close()


	ih = ipx - 1


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

			if(vaporP >= 0):
				absH = 0.622*vaporP/(1013.25-vaporP)	# absolute humidity
				RH = vaporP*760/1013.25
				if(cloud >= 0):
					nightR = 0.0000000488 * (tempA+273.16)**4 \
							* (1-0.62*cloud/10) * (0.49-0.076*math.sqrt(RH))
			else:
				absH = 0.7 * sim().abshumid(tempA)

			for idt in range(idx):
				TM = hour - 1 + interval*idt

				temp_o = t0 + (tempA-t0)*interval*idt
				print('temp_o :', temp_o, '℃ ')
				Wspeed = Wspeed0 * windC

				if(temp_o < 0):
					Fs = pre		# amount of snow in precipitation
					Fr = 0			# smount of rain in precipitation
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
				if(heater==0):
					pmx = 0
					ptl = 0
					ptm = 0
				if(pre > 0):
					ptl = ptl + pre*interval
					ptm = ptm + interval
				if(pre > pmx):
					pmx = pre
				if(snowT==99):
					wetT = pmx * 2
					if(wetT < 2):	wetT = 2
					B1 = wetT
					B2 = wetT
				if(snowT==999):
					wetT = pmx + 2
					dryT = wetT
					B1 = wetT
					B2 = wetT
					B3 = wetT
				tmr = tmr + interval

				if(pre > 0):	tmr = 0
				if(month<=4 or month>=11):
					if(pre > 0):
						level = 0
					else:
						if( (snow[Wsensor]+W[Wsensor]) > 0 ):
							level = 1
						else:
							level = 2
					if(level==0):
						TN = snowT + B1*temp_o
						TF = TN + stopT1
					elif(level==1):
						TN = wetT + B2*temp_o
						TF = TN + stopT2
					elif(level==2):
						TN = dryT + B3*temp_o
						TF = TN + stopT3

					# decide on/off
					if(BT[Tsensor] < TN):		heater = 1
					elif(BT[Tsensor] > TF):		heater = 0
					if(tmr > timer+0.0001):		heater = 0
					if(temp_o > offT):			heater = 0

					if(heater==0):
						irt = -1
						imo = 0
					elif(heater==1):
						irt += 1
						IA = irt % (circuit*step)	# circuit*step : rotation time
						if(IA < step*3):	imo = 1
					ntime[heater][level] += 1

				for ip in range(ipx):
					E.append(0)
					Q.append(0)
					NC.append(0)
					mlt.append(0)
					ict.append(0)

					if(ID[ip]==2):
						tset.append(temp_o)
						ict[ip] = 1
					elif(ID[ip]==9):
						tset.append(BT[ip])
						ict[ip] = 1

					if(nqd[ip]==1):
						Q[ip] = QQ[ip]
					elif(nqd[ip]==11):
						if(heater==1):
							ih = ip - 1
							Q[ip] = Qs
							ict[ip] = 0

					if(heater==1):
						if(ID[ip]==10):		# normal heater
							Q[ip] = Qs
							erot = erot + Q[ip]
						elif(ID[ip]>=11 and ID[ip]<=13):
							Q[ip] = Qs * 0.5
							if(imo==1):	Q[ip] = Qs*2
							erot = erot + Q[ip]

					for j in range(npn[ip]):
						if(nkr[ip][j]==11):
							HR[ip][j] = 0
							if(heater==1):
								HR[ip][j] = flowR
							flow = HR[ip][j]

					CR.append(C[ip])
					if(CL[ip] > 0):
						if(BT[ip] > tvm):
							CR[ip] = C[ip] + CL[ip]
						else:
							CR[ip] = C[ip] + CL[ip]*0.5
						Q[ip] = flr[ip]
						if(IL[ip]==1):
							tic = BT[ip]
							pic = (pcth-tic) / (pcth-pctl)
							CR[ip] = C[ip] + CL[ip]*(1-pic) \
									+ CL[ip]*0.5*pic + CL[ip]*80/abs(pcth-pctl)

					if(nsd[ip]==1):
						abrate = 0.8 - 30*cover[ip]
						if(abrate < abrate0):	abrate = abrate0
						sat[ip] = ( temp_o \
								+ (abrate*sun-0.9*nightR)/(sim().funa(Wspeed)+4) )
						EV[ip] = 0
						TS = BT[ip]
						lps = 0
						dps = cover[ip]
						if( (snow[ip]+Fs)>0 ):
							if(TS<0 or sat[ip]>0):
								mlt[ip] = 1
						if( (W[ip]+Fr)>0 ):
							if(TS < 0):
								mlt[ip] = 1
						if( (snow[ip]+Fs)>0 and (W[ip]+Fr)>0 ):
							mlt[ip] = 1
						peneH = 0			# penetration height
						if(cover[ip] > 0):
							peneH = W[ip] / (1000 - snow[ip]/cover[ip])
							# snow[ip]/cover[ip] : density of snow cover
						else:
							peneH = W[ip] / 1000
						if(peneH > maxPene):
							peneH = maxPene
						wat = W[ip] + Fr*interval


						while(True):
							SM = 0		# amount of melted snow(?)
							EV[ip] = 0
							BF[ip] = 0
							if(mlt[ip]==1):
								DH = dps - peneH
								if(DH<=0 or sat[ip]>0):
									DH = 0
									tsv = 0
									EV[ip] = 4 * sim().funa(Wspeed) \
											* (sim().abshumid(tsv)-absH)
									wat = W[ip] + (Fr-EV[ip])*interval
									if(wat < 0):
										EV[ip] = W[ip]/interval + Fr
										wat = 0
								htrm = 1 / (1/(sim().funa(Wspeed)+4) + DH/0.08)
								SM0 = (200*TS + htrm*sat[ip] - 590*EV[ip]) \
										* interval/80
								SM = SM0
								BF[ip] = 1
								if(SM0 < -1*wat):
									SM = -1 * wat
									BF[ip] = SM / SM0
								elif( SM0 > (snow[ip]+Fs*interval) ):
									SM = snow[ip] + Fs*interval
									BF[ip] = SM / SM0
							elif(mlt[ip]==0):
								tsv = BT[ip]
								EV[ip] = 4 * sim().funa(Wspeed)\
										* (sim().abshumid(tsv)- absH)
								if( EV[ip] > (W[ip]/interval + Fr) ):
									EV[ip] = W[ip]/interval + Fr

							htr[ip] = 1 / (1/(sim().funa(Wspeed)+4) + dps/0.08)
							QE[ip]  = -590*EV[ip] * AS[ip] * (1-BF[ip])
							SS[ip]  = snow[ip] - SM + Fs*interval
							WW[ip]  = W[ip] + SM + (Fr-EV[ip])*interval
							smf[ip] = SM

							if(SS[ip] > 0):
								if(SM > 0):
									G = (snow[ip] + Fs*interval) \
										/ (cover[ip] + Fs*interval/sfdens)
									Scover[ip] = SS[ip] / G
								else:
									Scover[ip] = cover[ip] + Fs*interval/sfdens \
													- SM/916
							else:
								Scover[ip] = 0

							ddps = (Scover[ip]+cover[ip])/2 - dps
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
						for ip in range(ipx):
							if(ict[ip]==1):
								T[ip] = tset[ip]
								E[ip] = 0
							S1.append( (Q[ip]+QE[ip]+E[ip])*interval \
											+ BT[ip]*CR[ip] )
							S2.append( CR[ip] )

							for j in range(npn[ip]):
								np_n = NP[ip][j]
								tmp = diffM*T[np_n-1] \
									+ (1-diffM)*(BT[np_n-1]-BT[ip])
								S1[ip] = S1[ip] + HR[ip][j]*(tmp)*interval
								S2[ip] = S2[ip] + diffM*HR[ip][j]*interval

							if(nsd[ip]==1):
								F = BF[ip]
								tmp = (1-F)*htr[ip]*sat[ip] \
									- F*200*(1-diffM)*BT[ip] \
									- (1-F)*htr[ip]*(1-diffM)*BT[ip]
								S1[ip] = S1[ip] + ( tmp )*interval*AS[ip]
								S2[ip] = S2[ip] \
										+ (F*200+(1-F)*htr[ip])\
											*interval*AS[ip]*diffM
							if(ict[ip]==1):
								E[ip] = (S2[ip]*tset[ip] - S1[ip]) / interval
							else:
								if(diffM<=0.01):
									T[ip] = S1[ip] / S2[ip]
									erx = 0
								else:
									TT.append( S1[ip] / S2[ip] )
									ER.append( TT[ip] - T[ip] )
									aer = abs(ER[ip])
									if(aer > erx):
										erx = aer
										ier = ip
									T[ip] = CK*ER[ip] + T[ip]

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
						if( T[ih]>maxT and Q[ih]>0 ):
							ilp += 1
							tset[ih] = maxT
							ict[ih] = 1
							Q[ih] = 0
							erx = 0
							continue
						if(E[ih] < 0):
							ilp += 1
							ict[ih] = 0
							E[ih] = 0
							erx = 0
							continue
						else:	break
					else:	break

				for ip in range(ipx):
					if(CL[ip] > 0):
						IL.append(0)
						flr[ip] = 0
						if( T[ip]<pcth and T[ip]>pctl ):
							IL[ip] = 1
						elif( BT[ip]>pcth and T[ip]<pcth ):
							flr[ip] = (C[ip]+CL[ip]) * (T[ip]-pcth)
							T[ip] = pcth
							IL[ip] = 1
						elif( BT[ip]>pctl and T[ip]<pctl ):
							flr[ip] = (C[ip] + CL[ip]*80/abs(pcth-pctl)) \
										* (T[ip]-pctl)
							T[ip] = pctl
						elif( BT[ip]<pctl and T[ip]>pctl ):
							flr[ip] = (C[ip] + CL[ip]*0.5) * (T[ip]-pctl)
							T[ip] = pctl
							IL[ip] = 1
						elif( BT[ip]<pcth and T[ip]>pcth ):
							flr[ip] = (C[ip] + CL[ip]*80/abs(pcth-pctl)) \
										* (T[ip]-pcth)
							T[ip] = pcth

					if(nsd[ip]==1):
						if( level==0 and T[ip]<0 ):
							L1 = 1
						elif( level==1 and T[ip]<0 ):
							L2 = 1
						if(SS[ip] < 0):
							SS[ip] = 0
						if(WW[ip] < 0):
							WW[ip] = 0
						if(Scover[ip] > FD):
							FD = Scover[ip]

						sca = sca + SS[ip]*AS[ip]/area
						water = water + WW[ip]*AS[ip]/area
						evaporate = evaporate + EV[ip]*interval*AS[ip]/area
						ssmf = ssmf + smf[ip]*AS[ip]/area

						if( Scover[ip]>0 and SS[ip]>0 ):
							GM = SS[ip] / Scover[ip]
							EE = 16 * math.exp(0.021*GM)
							GM = GM * math.exp( SS[ip]/2/EE*interval/24 )
							if(GM > 916):
								GM = 916
							Scover[ip] = SS[ip] / GM
				

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
				SE = E[ih] + Q[ih] + erot

				if(SE > 0):
					onT += 1
					Qsup += se
					rse = Qsup / onT

				melt += ssmf

				lpr += 1

				if( pre>0 and T[ip]>0 ):
					TS = diffM*T[Wsensor] + (1-diffM)*BT[Wsensor]
					Qr2 = Qr2 + BF[Wsensor]*200*TS \
							+ (1-BF[Wsensor])*htr[Wsensor]*(TS-sat[Wsensor])

				for ip in range(ipx):
					if(nsd[ip]==1):
						if( (sca+pre)==0 and WW[ip]>remainW ):
							WW[ip] = remainW	# remained water
						W[ip] = WW[ip]
						snow[ip] = SS[ip]
						cover[ip] = Scover[ip]
						if( pre>0 and T[ip]>0 ):
							TS = diffM*T[ip] + (1-diffM)*BT[ip]
							Qr_plus = BF[ip]*200*TS \
										+ (1-BF[ip])*htr[ip]*(TS-sat[ip])
							Qr = Qr + (Qr_plus)*AS[ip]

					BT[ip] = T[ip]

			Tsnow_off = ntime[0][0] * interval
			Tsnow_on  = ntime[1][0] * interval
			Twet_off  = ntime[0][1] * interval
			Twet_on   = ntime[1][1] * interval
			Tdry_off  = ntime[0][2] * interval
			Tdry_on   = ntime[1][2] * interval
			print('snow.on time :', Tsnow_on, '\tsnow.off time :', Tsnow_off, \
				'\nwet.on time :', Twet_on, '\twet.off time :', Twet_off, \
				'\ndry.on time :', Tdry_on, '\tdry.off time :', Tdry_off)
			if(heater==0):	print('heater : off')
			else:			print('heater : on')

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
			Theight.append( nsca[i] * interval )
		for i in range(slevel+1):
			Tcover.append( nsca[i] * interval )
	
		Tsnow_off = ntime[0][0] * interval
		Tsnow_on  = ntime[1][0] * interval
		Twet_off  = ntime[0][1] * interval
		Twet_on   = ntime[1][1] * interval
		Tdry_off  = ntime[0][2] * interval
		Tdry_on   = ntime[1][2] * interval
		print('TIME\t', 'snow.on :', Tsnow_on, 'snow.off :', Tsnow_off, \
						', wet.on :', Twet_on, 'wet.off :', Twet_off, \
						', dry.on :', Tdry_on, 'dry.off :', Tdry_off)
